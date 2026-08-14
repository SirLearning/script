"""Full-lineage syntenic sites: cached site lists, overlap join, proportional sampling."""

from __future__ import annotations

import gzip
import os
import random
import re
import subprocess
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from infra.utils.graph import plot_categorical_bar_chart
from infra.utils.io import save_df_to_tsv, save_thresholds
from infra.wheat.ref_v1 import get_chromosome, get_pos_on_chromosome

LINEAGE_SPECIES: dict[str, tuple[str, ...]] = {
    "A": ("traesA", "trboe", "trmon", "trura", "trdicA", "trdurA"),
    "B": ("traesB", "trdicB", "trdurB", "aespe"),
    "D": ("traesD", "aetau", "aecom", "aeuni", "aelon", "aebic"),
}

CHR_LABELS: tuple[str, ...] = tuple(f"{i}{sg}" for i in range(1, 8) for sg in "ABD")
PLINK_CHR_FILE_IDS: tuple[str, ...] = tuple(f"chr{i:03d}" for i in range(1, 43))

_MAF_RE = re.compile(r"(?:^|;)MAF=([^;]+)")
_SCRIPT_DIR = Path(__file__).resolve().parents[5] / "workflow/Genetics/subworkflows/tmp/scripts/collinear"


@dataclass(frozen=True)
class CollinearCandidate:
    chr_label: str
    subgenome: str
    ref_pos_1based: int
    chr_file_id: str
    variant_id: str
    maf: float


def _parse_maf(info_field: str) -> float | None:
    match = _MAF_RE.search(info_field)
    if not match:
        return None
    return float(match.group(1))


def _plink_ids_for_chr_label(chr_label: str) -> list[int]:
    subgenome = chr_label[-1]
    genome_num = int(chr_label[:-1])
    if subgenome == "A":
        return [(genome_num - 1) * 6 + 1, (genome_num - 1) * 6 + 2]
    if subgenome == "B":
        return [(genome_num - 1) * 6 + 3, (genome_num - 1) * 6 + 4]
    if subgenome == "D":
        return [(genome_num - 1) * 6 + 5, (genome_num - 1) * 6 + 6]
    raise ValueError(f"Unsupported chr_label: {chr_label}")


def collinear_layout(root: str | Path) -> dict[str, Path]:
    """Standard layout under ``process/test_collinear/collinear/``."""
    root = Path(root)
    layout = {
        "root": root,
        "full_lineage_sites": root / "full_lineage_sites",
        "common_variant_sites": root / "common_variant_sites",
        "extract": root / "extract",
        "info": root / "info",
        "plots": root / "plots",
    }
    for path in layout.values():
        if isinstance(path, Path):
            path.mkdir(parents=True, exist_ok=True)
    return layout


def plot_full_lineage_subgenome_fraction(
    summary_tsv: str | Path,
    output_png: str | Path,
) -> pd.DataFrame:
    """Bar plot of full-lineage fraction by subgenome (A/B/D) from summary TSV."""
    df = pd.read_csv(summary_tsv, sep="\t")
    sub = df[df["scope"].isin(["A", "B", "D"])].copy()
    if sub.empty:
        raise ValueError(f"No A/B/D rows in {summary_tsv}")

    categories = sub["scope"].tolist()
    fractions = sub["fraction_full_lineage"].astype(float).tolist()
    colors = {"A": "#4c72b0", "B": "#c44e52", "D": "#55a868"}
    plot_categorical_bar_chart(
        categories=categories,
        values=fractions,
        title="Full-lineage syntenic site fraction by subgenome",
        ylabel="Fraction full-lineage (k≥2 sites)",
        filename=str(output_png),
        colors=[colors.get(c, "steelblue") for c in categories],
        ylim=(0.0, max(0.4, max(fractions) * 1.15)),
        value_formatter=lambda v: f"{v:.1%}",
    )
    return sub


def extract_common_variant_sites_chr(
    chr_file_id: str,
    main_raw_dir: str | Path,
    output_gz: str | Path,
    *,
    maf_threshold: float = 0.05,
) -> int:
    """Write MAF-filtered variant coordinates for one PLINK segment file."""
    plink_chr = int(chr_file_id.replace("chr", ""))
    chr_label = get_chromosome(plink_chr)
    if chr_label is None:
        return 0

    pvar_path = Path(main_raw_dir) / f"{chr_file_id}.plink2.pvar"
    output_gz = Path(output_gz)
    output_gz.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    with open(pvar_path, encoding="utf-8") as src, gzip.open(output_gz, "wt", encoding="utf-8") as out:
        out.write("chr_label\tref_pos_1based\tvariant_id\tmaf\n")
        for line in src:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            local_pos = int(fields[1])
            maf = _parse_maf(fields[5])
            if maf is None or not (maf > maf_threshold):
                continue
            ref_pos = int(get_pos_on_chromosome(plink_chr, local_pos))
            if get_chromosome(plink_chr) != chr_label:
                continue
            out.write(f"{chr_label}\t{ref_pos}\t{fields[2]}\t{maf}\n")
            n_written += 1
    return n_written


def build_common_variant_sites(
    main_raw_dir: str | Path,
    out_dir: str | Path,
    *,
    maf_threshold: float = 0.05,
    workers: int = 7,
    force: bool = False,
) -> dict[str, int]:
    """Parallel extract of common-variant site lists (cached under ``common_variant_sites/``)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tasks = []
    for chr_file_id in PLINK_CHR_FILE_IDS:
        out_path = out_dir / f"{chr_file_id}.tsv.gz"
        if not force and out_path.stat().st_size > 0 if out_path.exists() else False:
            continue
        pvar_path = Path(main_raw_dir) / f"{chr_file_id}.plink2.pvar"
        if not pvar_path.exists():
            continue
        tasks.append((chr_file_id, main_raw_dir, out_path, maf_threshold))

    counts: dict[str, int] = {}
    if not tasks:
        for path in sorted(out_dir.glob("chr*.tsv.gz")):
            counts[path.stem] = _count_gzip_data_lines(path)
        return counts

    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(extract_common_variant_sites_chr, cid, mdir, op, maf_threshold=mt): cid
            for cid, mdir, op, mt in tasks
        }
        for fut in as_completed(futures):
            cid = futures[fut]
            n = fut.result()
            counts[cid] = n
            print(f"common_variant_sites {cid}: {n:,} variants (MAF>{maf_threshold})")

    for path in sorted(out_dir.glob("chr*.tsv.gz")):
        counts.setdefault(path.stem, _count_gzip_data_lines(path))
    return counts


def build_full_lineage_sites(
    syntenic_matrix_dir: str | Path,
    out_dir: str | Path,
    *,
    workers: int = 7,
    force: bool = False,
) -> None:
    """Parallel mawk extract of full-lineage coordinates (cached under ``full_lineage_sites/``)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not force and all((out_dir / f"{label}.tsv.gz").exists() for label in CHR_LABELS):
        print(f"Reuse existing full_lineage_sites in {out_dir}")
        return

    script = _SCRIPT_DIR / "run_full_lineage_extract_all.sh"
    if not script.exists():
        raise FileNotFoundError(script)
    env = {
        **os.environ,
        "MATRIX_DIR": str(syntenic_matrix_dir),
        "OUT_DIR": str(out_dir),
        "WORKERS": str(workers),
    }
    subprocess.run(["bash", str(script)], check=True, env=env)


def _count_gzip_data_lines(path: Path) -> int:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return sum(1 for _ in handle) - 1


def _load_common_variant_dict(
    common_dir: Path,
    chr_label: str,
) -> dict[int, tuple[str, str, float]]:
    index: dict[int, tuple[str, str, float]] = {}
    for plink_chr in _plink_ids_for_chr_label(chr_label):
        chr_file_id = f"chr{plink_chr:03d}"
        path = common_dir / f"{chr_file_id}.tsv.gz"
        if not path.exists():
            continue
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            header = handle.readline()
            if not header:
                continue
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4 or parts[0] != chr_label:
                    continue
                ref_pos = int(parts[1])
                index[ref_pos] = (chr_file_id, parts[2], float(parts[3]))
    return index


def _count_overlaps(
    full_lineage_dir: Path,
    common_dir: Path,
) -> tuple[Counter, list[dict[str, int | str]]]:
    totals: Counter = Counter()
    per_chr: list[dict[str, int | str]] = []

    for chr_label in CHR_LABELS:
        fl_path = full_lineage_dir / f"{chr_label}.tsv.gz"
        if not fl_path.exists():
            continue
        common = _load_common_variant_dict(common_dir, chr_label)
        n_fl = 0
        n_ov = 0
        with gzip.open(fl_path, "rt", encoding="utf-8") as handle:
            handle.readline()
            for line in handle:
                n_fl += 1
                ref_pos = int(line.rstrip("\n").split("\t")[1])
                if ref_pos in common:
                    n_ov += 1
        totals["n_full_lineage"] += n_fl
        totals["n_overlap_maf"] += n_ov
        per_chr.append(
            {
                "chr_label": chr_label,
                "subgenome": chr_label[-1],
                "n_full_lineage": n_fl,
                "n_overlap_variant": n_ov,
                "n_overlap_maf": n_ov,
                "n_maf_index": len(common),
            }
        )
        print(f"{chr_label}: full_lineage={n_fl:,} overlap_maf>{0.05}={n_ov:,}")
    totals["n_overlap_variant"] = totals["n_overlap_maf"]
    return totals, per_chr


def _reservoir_sample(
    full_lineage_dir: Path,
    common_dir: Path,
    target_n: int,
    seed: int,
    *,
    totals: Counter | None = None,
    per_chr_rows: list[dict[str, int | str]] | None = None,
) -> list[CollinearCandidate]:
    if totals is None or per_chr_rows is None:
        totals, per_chr_rows = _count_overlaps(full_lineage_dir, common_dir)
    n_candidates = int(totals["n_overlap_maf"])
    if n_candidates == 0:
        return []
    if n_candidates <= target_n:
        return _collect_all_candidates(full_lineage_dir, common_dir)

    by_sg_count = Counter()
    for row in per_chr_rows:
        by_sg_count[str(row["subgenome"])] += int(row["n_overlap_maf"])

    alloc = _allocate_proportional(by_sg_count, target_n)
    rng = random.Random(seed)
    reservoirs: dict[str, list[CollinearCandidate]] = {sg: [] for sg in ("A", "B", "D")}
    seen: Counter = Counter()

    for chr_label in CHR_LABELS:
        subgenome = chr_label[-1]
        k = alloc.get(subgenome, 0)
        if k == 0:
            continue
        fl_path = full_lineage_dir / f"{chr_label}.tsv.gz"
        if not fl_path.exists():
            continue
        common = _load_common_variant_dict(common_dir, chr_label)
        with gzip.open(fl_path, "rt", encoding="utf-8") as handle:
            handle.readline()
            for line in handle:
                ref_pos = int(line.rstrip("\n").split("\t")[1])
                hit = common.get(ref_pos)
                if hit is None:
                    continue
                seen[subgenome] += 1
                chr_file_id, variant_id, maf = hit
                cand = CollinearCandidate(
                    chr_label=chr_label,
                    subgenome=subgenome,
                    ref_pos_1based=ref_pos,
                    chr_file_id=chr_file_id,
                    variant_id=variant_id,
                    maf=maf,
                )
                pool = reservoirs[subgenome]
                if len(pool) < k:
                    pool.append(cand)
                else:
                    j = rng.randint(0, seen[subgenome] - 1)
                    if j < k:
                        pool[j] = cand

    selected: list[CollinearCandidate] = []
    for sg in ("A", "B", "D"):
        selected.extend(reservoirs[sg])
    return selected


def _collect_all_candidates(
    full_lineage_dir: Path,
    common_dir: Path,
) -> list[CollinearCandidate]:
    out: list[CollinearCandidate] = []
    for chr_label in CHR_LABELS:
        fl_path = full_lineage_dir / f"{chr_label}.tsv.gz"
        if not fl_path.exists():
            continue
        common = _load_common_variant_dict(common_dir, chr_label)
        with gzip.open(fl_path, "rt", encoding="utf-8") as handle:
            handle.readline()
            for line in handle:
                ref_pos = int(line.rstrip("\n").split("\t")[1])
                hit = common.get(ref_pos)
                if hit is None:
                    continue
                chr_file_id, variant_id, maf = hit
                out.append(
                    CollinearCandidate(
                        chr_label=chr_label,
                        subgenome=chr_label[-1],
                        ref_pos_1based=ref_pos,
                        chr_file_id=chr_file_id,
                        variant_id=variant_id,
                        maf=maf,
                    )
                )
    return out


def _allocate_proportional(by_sg_count: Counter, target_n: int) -> dict[str, int]:
    total = sum(by_sg_count.values())
    if total == 0:
        return {"A": 0, "B": 0, "D": 0}
    alloc: dict[str, int] = {}
    remainders: list[tuple[float, str]] = []
    assigned = 0
    for sg in ("A", "B", "D"):
        n = by_sg_count.get(sg, 0)
        if not n:
            alloc[sg] = 0
            continue
        exact = target_n * n / total
        base = int(exact)
        alloc[sg] = base
        assigned += base
        remainders.append((exact - base, sg))
    remainder = target_n - assigned
    remainders.sort(reverse=True)
    for _, sg in remainders:
        if remainder <= 0:
            break
        if by_sg_count.get(sg, 0):
            alloc[sg] += 1
            remainder -= 1
    return alloc


def proportional_sample_candidates(
    candidates: list[CollinearCandidate],
    target_n: int,
    *,
    seed: int = 42,
) -> list[CollinearCandidate]:
    """Sample ``target_n`` candidates with subgenome proportions preserved."""
    if not candidates:
        return []
    if len(candidates) <= target_n:
        return list(candidates)

    by_sg: dict[str, list[CollinearCandidate]] = {"A": [], "B": [], "D": []}
    for item in candidates:
        by_sg[item.subgenome].append(item)

    by_sg_count = Counter({sg: len(by_sg[sg]) for sg in by_sg})
    alloc = _allocate_proportional(by_sg_count, target_n)
    rng = random.Random(seed)
    selected: list[CollinearCandidate] = []
    for sg in ("A", "B", "D"):
        pool = by_sg[sg]
        n = min(alloc.get(sg, 0), len(pool))
        if n:
            selected.extend(rng.sample(pool, n))
    return selected


def write_extract_lists(
    selected: list[CollinearCandidate],
    extract_dir: str | Path,
    *,
    all_chr_file_ids: tuple[str, ...] | None = None,
) -> dict[str, int]:
    """Write per-chr plink2 --extract ID lists; return counts per chr_file_id."""
    extract_dir = Path(extract_dir)
    extract_dir.mkdir(parents=True, exist_ok=True)
    by_chr: dict[str, list[str]] = {}
    for item in selected:
        by_chr.setdefault(item.chr_file_id, []).append(item.variant_id)

    if all_chr_file_ids is None:
        all_chr_file_ids = PLINK_CHR_FILE_IDS

    counts: dict[str, int] = {}
    for chr_file_id in all_chr_file_ids:
        ids = by_chr.get(chr_file_id, [])
        out_path = extract_dir / f"{chr_file_id}.collinear.txt"
        with open(out_path, "w", encoding="utf-8") as handle:
            for vid in ids:
                handle.write(f"{vid}\n")
        counts[chr_file_id] = len(ids)
    return counts


def prepare_collinear_test_dataset(
    *,
    syntenic_matrix_dir: str | Path,
    main_raw_dir: str | Path,
    full_lineage_summary_tsv: str | Path,
    output_dir: str | Path,
    target_n: int = 2_000_000,
    maf_threshold: float = 0.05,
    seed: int = 42,
    workers: int = 7,
    force: bool = False,
) -> dict[str, int | float | str]:
    """Build cached site lists, join, sample, and write plink extract lists."""
    layout = collinear_layout(output_dir)
    for path in layout.values():
        if isinstance(path, Path):
            path.mkdir(parents=True, exist_ok=True)

    plot_full_lineage_subgenome_fraction(
        full_lineage_summary_tsv,
        layout["plots"] / "full_lineage_subgenome_fraction.png",
    )

    print("Step 1/3: full_lineage_sites (parallel mawk)...")
    build_full_lineage_sites(
        syntenic_matrix_dir,
        layout["full_lineage_sites"],
        workers=workers,
        force=force,
    )

    print("Step 2/3: common_variant_sites (parallel pvar scan)...")
    build_common_variant_sites(
        main_raw_dir,
        layout["common_variant_sites"],
        maf_threshold=maf_threshold,
        workers=workers,
        force=force,
    )

    print("Step 3/3: join + proportional sample...")
    totals, per_chr_rows = _count_overlaps(layout["full_lineage_sites"], layout["common_variant_sites"])
    save_df_to_tsv(pd.DataFrame(per_chr_rows), str(layout["info"] / "collinear_overlap_by_chr.tsv"))

    summary = {
        "n_full_lineage": int(totals["n_full_lineage"]),
        "n_overlap_variant": int(totals["n_overlap_variant"]),
        "n_overlap_maf": int(totals["n_overlap_maf"]),
        "n_candidates": int(totals["n_overlap_maf"]),
        "target_n": int(target_n),
        "maf_threshold": maf_threshold,
    }
    save_df_to_tsv(pd.DataFrame([summary]), str(layout["info"] / "collinear_overlap_summary.tsv"))

    selected = _reservoir_sample(
        layout["full_lineage_sites"],
        layout["common_variant_sites"],
        target_n,
        seed,
        totals=totals,
        per_chr_rows=per_chr_rows,
    )
    sg_counts = Counter(item.subgenome for item in selected)
    extract_counts = write_extract_lists(selected, layout["extract"])

    sample_summary = {
        **summary,
        "n_selected": len(selected),
        "selected_A": sg_counts.get("A", 0),
        "selected_B": sg_counts.get("B", 0),
        "selected_D": sg_counts.get("D", 0),
        "n_extract_chr": sum(1 for n in extract_counts.values() if n > 0),
    }
    save_thresholds(sample_summary, str(layout["info"] / "collinear_sample.th.tsv"))

    selected_df = pd.DataFrame(
        [
            {
                "chr_label": c.chr_label,
                "subgenome": c.subgenome,
                "ref_pos": c.ref_pos_1based,
                "chr_file_id": c.chr_file_id,
                "variant_id": c.variant_id,
                "maf": c.maf,
            }
            for c in selected
        ]
    )
    save_df_to_tsv(selected_df, str(layout["info"] / "collinear_selected_variants.tsv"))
    print(f"Selected {len(selected):,} variants for test_collinear (target {target_n:,})")
    return sample_summary
