"""Helpers for SINGER ARG sampling from a prepared VCF window."""

from __future__ import annotations

import csv
import random
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

from infra.wheat.ref_v1 import (
    get_centromere_end,
    get_centromere_start,
    get_chromosomes,
    get_ref_v1_chr_name,
    pvar_chrom_to_plink_chr_id,
)


def nucleotide_diversity_in_region(
    vcf_path: str | Path,
    *,
    start: int,
    end: int,
) -> float:
    """Pairwise diversity in ``[start, end]``, matching SINGER ``singer_master``."""
    if end <= start:
        raise ValueError("end must be greater than start")
    total_diversity = 0.0
    path = Path(vcf_path)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            pos = int(fields[1])
            if pos < start:
                continue
            if pos > end:
                break
            alleles = [fields[3]] + fields[4].split(",")
            allele_counts = [0] * len(alleles)
            n_samples = 0
            for sample_field in fields[9:]:
                gt = sample_field.split(":", 1)[0]
                if "." in gt:
                    continue
                gt_alleles = gt.replace("|", "/").split("/")
                n_samples += 1
                for allele in gt_alleles:
                    allele_counts[int(allele)] += 1
            if n_samples == 0:
                continue
            allele_freqs = [count / (2 * n_samples) for count in allele_counts]
            total_diversity += sum(
                freq * (1.0 - freq) for freq in allele_freqs if 0.0 < freq < 1.0
            )
    return total_diversity / (end - start)


def estimate_ne_from_vcf(
    vcf_path: str | Path,
    *,
    start: int,
    end: int,
    mut_rate: float,
) -> float:
    """Diploid Ne from π = 4 Ne μ, as in SINGER ``singer_master``."""
    if mut_rate <= 0:
        raise ValueError("mut_rate must be positive")
    diversity = nucleotide_diversity_in_region(vcf_path, start=start, end=end)
    return diversity / 4.0 / mut_rate


def write_ne_file(path: str | Path, ne: float) -> None:
    Path(path).write_text(f"{ne}\n", encoding="utf-8")


def window_overlaps_centromere(ref_chrom: str, start_bp: int, end_bp: int) -> bool:
    """True when ``[start_bp, end_bp)`` overlaps the IWGSC v1 centromere of ``ref_chrom``."""
    token = str(ref_chrom).strip()
    chrom = token[3:] if token.startswith("chr") else token
    c_start = get_centromere_start(chrom)
    c_end = get_centromere_end(chrom)
    if c_start is None or c_end is None:
        return False
    return start_bp < c_end and end_bp > c_start


def count_pvar_window_bins(
    pvar_path: str | Path,
    *,
    subgenome: str,
    window_bp: int = 2_000_000,
) -> list[dict]:
    """Count SNPs in half-open ``window_bp`` bins of a pvar."""
    if window_bp <= 0:
        raise ValueError("window_bp must be positive")
    counts: dict[tuple[str, int], int] = defaultdict(int)
    with Path(pvar_path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 2)
            if len(fields) < 2:
                continue
            chrom = fields[0].strip()
            pos = int(fields[1])
            start = (pos // window_bp) * window_bp
            counts[(chrom, start)] += 1
    rows = []
    for (chrom, start), n_snps in counts.items():
        plink_chr = pvar_chrom_to_plink_chr_id(chrom) or chrom
        try:
            ref_chrom = get_ref_v1_chr_name(int(plink_chr))
        except (TypeError, ValueError):
            ref_chrom = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
        rows.append(
            {
                "subgenome": subgenome,
                "plink_chr": str(plink_chr),
                "ref_chrom": ref_chrom,
                "start_bp": int(start),
                "end_bp": int(start + window_bp),
                "n_snps": int(n_snps),
            }
        )
    return rows


def select_spread_windows(
    bins: list[dict],
    *,
    n_windows: int,
    min_snps: int = 50,
    skip_leading_telomere: bool = True,
) -> list[dict]:
    """Pick ``n_windows`` 2 Mb bins spread across chromosomes, densest first per chrom."""
    if n_windows <= 0:
        raise ValueError("n_windows must be positive")
    eligible: list[dict] = []
    for row in bins:
        start = int(row["start_bp"])
        end = int(row["end_bp"])
        if int(row["n_snps"]) < min_snps:
            continue
        if skip_leading_telomere and start == 0:
            continue
        if window_overlaps_centromere(str(row["ref_chrom"]), start, end):
            continue
        eligible.append(row)
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    for row in eligible:
        by_chrom[str(row["ref_chrom"])].append(row)
    for chrom, rows in by_chrom.items():
        rows.sort(key=lambda item: (-int(item["n_snps"]), int(item["start_bp"])))
        by_chrom[chrom] = rows

    chrom_order = [f"chr{name}" for name in get_chromosomes() if f"chr{name}" in by_chrom]
    extra = [chrom for chrom in sorted(by_chrom) if chrom not in chrom_order]
    chrom_order.extend(extra)

    picked: list[dict] = []
    cursors = {chrom: 0 for chrom in chrom_order}
    while len(picked) < n_windows:
        progressed = False
        for chrom in chrom_order:
            idx = cursors[chrom]
            rows = by_chrom[chrom]
            if idx >= len(rows):
                continue
            picked.append(rows[idx])
            cursors[chrom] = idx + 1
            progressed = True
            if len(picked) >= n_windows:
                break
        if not progressed:
            break
    if len(picked) < n_windows:
        raise ValueError(
            f"only {len(picked)} windows passed filters; need {n_windows}"
        )
    return picked


def write_windows_tsv(windows: list[dict], path: str | Path) -> None:
    columns = ["subgenome", "plink_chr", "ref_chrom", "start_bp", "end_bp", "n_snps"]
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in windows:
            handle.write("\t".join(str(row[col]) for col in columns) + "\n")


def prepare_wap_china_singer_inputs(
    *,
    germplasm_tsv: str | Path,
    psam_path: str | Path,
    pvar_by_subgenome: dict[str, str | Path],
    keep_path: str | Path,
    windows_path: str | Path,
    country: str = "China",
    n_windows: int = 64,
    window_bp: int = 2_000_000,
    min_snps: int = 50,
) -> tuple[int, int]:
    """Write a China-WAP keep file and a spread 2 Mb window table."""
    from genetics.germplasm.sample.wap_china import (
        china_wap_taxa,
        load_psam_iids,
        write_plink2_keep,
    )

    taxa = china_wap_taxa(
        germplasm_tsv,
        load_psam_iids(psam_path),
        country=country,
    )
    write_plink2_keep(taxa, keep_path)
    bins: list[dict] = []
    for subgenome, pvar in pvar_by_subgenome.items():
        bins.extend(count_pvar_window_bins(pvar, subgenome=subgenome, window_bp=window_bp))
    windows = select_spread_windows(bins, n_windows=n_windows, min_snps=min_snps)
    write_windows_tsv(windows, windows_path)
    return len(taxa), len(windows)


def singer_window_id(row: Mapping[str, object]) -> str:
    """Stable window id used by extract / Beagle / SINGER publish names."""
    return (
        f"{row['subgenome']}.{row['ref_chrom']}."
        f"{int(row['start_bp'])}_{int(row['end_bp'])}"
    )


def arg_nodes_path(singer_dir: str | Path, window_id: str) -> Path:
    return Path(singer_dir) / f"singer_{window_id}_nodes_0.txt"


def load_windows_tsv(path: str | Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def list_windows_missing_arg(
    windows: Sequence[Mapping[str, object]] | str | Path,
    singer_dir: str | Path,
) -> list[dict]:
    """Windows whose published ARG is missing ``singer_{id}_nodes_0.txt``."""
    rows: Sequence[Mapping[str, object]]
    if isinstance(windows, (str, Path)):
        rows = load_windows_tsv(windows)
    else:
        rows = windows
    missing: list[dict] = []
    for row in rows:
        window_id = singer_window_id(row)
        if not arg_nodes_path(singer_dir, window_id).is_file():
            missing.append({**dict(row), "window_id": window_id})
    return missing


def random_singer_seeds(
    n: int,
    *,
    exclude: Iterable[int] = (42,),
    rng: random.Random | None = None,
) -> list[int]:
    """Our own SINGER ``-seed`` values. Never reuse 42; not ``singer_master -debug``."""
    if n <= 0:
        raise ValueError("n must be positive")
    banned = set(exclude)
    generator = rng if rng is not None else random.SystemRandom()
    seeds: list[int] = []
    seen: set[int] = set()
    while len(seeds) < n:
        seed = int(generator.randint(1, 2_147_483_646))
        if seed in banned or seed in seen:
            continue
        seen.add(seed)
        seeds.append(seed)
    return seeds


def parse_estimated_ne(run_log: str | Path) -> float | None:
    """Read ``estimated_Ne=`` from a SINGER process ``*.run.log``."""
    path = Path(run_log)
    if not path.is_file():
        return None
    for line in path.open(encoding="utf-8", errors="replace"):
        if line.startswith("estimated_Ne="):
            return float(line.split("=", 1)[1].strip())
    return None


def _count_nonempty_lines(path: Path) -> int:
    if not path.is_file():
        return 0
    return sum(1 for line in path.open(encoding="utf-8", errors="replace") if line.strip())


def _node_times(path: Path) -> list[float]:
    times: list[float] = []
    for line in path.open(encoding="utf-8", errors="replace"):
        token = line.strip().split()
        if not token:
            continue
        times.append(float(token[0]))
    return times


def _summarize_nodes(path: Path) -> tuple[int, int, float | None]:
    times = _node_times(path)
    if not times:
        return 0, 0, None
    return len(times), sum(time == 0.0 for time in times), max(times)


def summarize_singer_windows(
    windows: Sequence[Mapping[str, object]] | str | Path,
    singer_dir: str | Path,
) -> list[dict]:
    """Per-window ARG preview: success flag, diversity Ne, and sample-0 ARG size."""
    rows: Sequence[Mapping[str, object]]
    if isinstance(windows, (str, Path)):
        rows = load_windows_tsv(windows)
    else:
        rows = windows
    out_dir = Path(singer_dir)
    summaries: list[dict] = []
    for row in rows:
        window_id = singer_window_id(row)
        start_bp = int(row["start_bp"])
        end_bp = int(row["end_bp"])
        nodes0 = arg_nodes_path(out_dir, window_id)
        rec: dict = {
            "window_id": window_id,
            "subgenome": str(row["subgenome"]),
            "ref_chrom": str(row["ref_chrom"]),
            "start_bp": start_bp,
            "end_bp": end_bp,
            "start_mb": start_bp / 1_000_000.0,
            "n_snps": int(row["n_snps"]),
            "arg_ok": nodes0.is_file(),
            "estimated_ne": parse_estimated_ne(out_dir / f"singer_{window_id}.run.log"),
            "n_arg_samples": 0,
            "n_nodes": None,
            "n_leaves": None,
            "t_max": None,
            "n_recomb": None,
            "n_muts": None,
            "mean_t_max": None,
            "mean_n_recomb": None,
        }
        if not rec["arg_ok"]:
            summaries.append(rec)
            continue
        n_nodes, n_leaves, t_max = _summarize_nodes(nodes0)
        rec["n_nodes"] = n_nodes
        rec["n_leaves"] = n_leaves
        rec["t_max"] = t_max
        rec["n_recomb"] = _count_nonempty_lines(out_dir / f"singer_{window_id}_recombs_0.txt")
        rec["n_muts"] = _count_nonempty_lines(out_dir / f"singer_{window_id}_muts_0.txt")
        t_maxes: list[float] = []
        n_recombs: list[int] = []
        sample = 0
        while (out_dir / f"singer_{window_id}_nodes_{sample}.txt").is_file():
            _, _, sample_t = _summarize_nodes(out_dir / f"singer_{window_id}_nodes_{sample}.txt")
            if sample_t is not None:
                t_maxes.append(sample_t)
            n_recombs.append(
                _count_nonempty_lines(out_dir / f"singer_{window_id}_recombs_{sample}.txt")
            )
            sample += 1
        rec["n_arg_samples"] = sample
        rec["mean_t_max"] = sum(t_maxes) / len(t_maxes) if t_maxes else None
        rec["mean_n_recomb"] = sum(n_recombs) / len(n_recombs) if n_recombs else None
        summaries.append(rec)
    return summaries


def plot_singer_arg_preview(
    windows: Sequence[Mapping[str, object]] | str | Path,
    singer_dir: str | Path,
    output_prefix: str | Path,
) -> dict:
    """Write window-status tables and preview plots from published SINGER texts."""
    import pandas as pd

    from infra.utils.graph import (
        plot_distribution_with_stats,
        plot_scatter_with_thresholds,
        plot_stacked_bar_chart,
    )
    from infra.utils.io import save_df_to_tsv, save_thresholds

    rows = summarize_singer_windows(windows, singer_dir)
    frame = pd.DataFrame(rows)
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    save_df_to_tsv(frame, f"{prefix}.window.info.tsv", float_format="%.6g")

    n_ok = int(frame["arg_ok"].sum()) if not frame.empty else 0
    n_fail = int(len(frame) - n_ok)
    ok = frame.loc[frame["arg_ok"]].copy() if n_ok else frame.iloc[0:0].copy()
    stats = {
        "n_windows": len(frame),
        "n_ok": n_ok,
        "n_fail": n_fail,
        "median_estimated_ne": float(ok["estimated_ne"].median()) if n_ok else "",
        "median_n_recomb": float(ok["n_recomb"].median()) if n_ok else "",
        "median_t_max": float(ok["t_max"].median()) if n_ok else "",
    }
    save_thresholds(stats, f"{prefix}.window.th.tsv")

    status = frame.copy()
    status["status"] = status["arg_ok"].map({True: "ok", False: "fail"})
    counts = pd.crosstab(status["subgenome"], status["status"])
    for col in ("ok", "fail"):
        if col not in counts.columns:
            counts[col] = 0
    plot_stacked_bar_chart(
        counts[["ok", "fail"]],
        title="SINGER ARG windows by subgenome",
        filename=f"{prefix}.status.subgenome.png",
        x_label="Subgenome",
        y_label="Windows",
        stack_col_title="ARG",
        stack_colors={"ok": "steelblue", "fail": "crimson"},
        rotate_xlabels=0,
    )
    plot_scatter_with_thresholds(
        status,
        "start_mb",
        "n_snps",
        title="Window SNP count vs start (ok vs fail)",
        filename=f"{prefix}.nsnps.start.png",
        xlabel="Window start (Mb)",
        ylabel="SNPs in 2 Mb window",
        group_col="status",
    )
    if n_ok:
        plot_distribution_with_stats(
            ok,
            "estimated_ne",
            title="Diversity Ne in windows with an ARG",
            filename=f"{prefix}.ne.dist.png",
            x_label="Estimated Ne (pi / 4 mu)",
            y_label="Windows",
            bins=min(20, max(n_ok, 1)),
            mean_val=float(ok["estimated_ne"].mean()),
            median_val=float(ok["estimated_ne"].median()),
        )
        plot_distribution_with_stats(
            ok,
            "t_max",
            title="Max node time in ARG sample 0",
            filename=f"{prefix}.tmax.dist.png",
            x_label="Max node time (SINGER units)",
            y_label="Windows",
            bins=min(20, max(n_ok, 1)),
            mean_val=float(ok["t_max"].mean()),
            median_val=float(ok["t_max"].median()),
        )
        plot_scatter_with_thresholds(
            ok,
            "n_snps",
            "n_recomb",
            title="ARG recombinations vs SNPs (sample 0)",
            filename=f"{prefix}.nrecomb.nsnps.png",
            xlabel="SNPs in 2 Mb window",
            ylabel="Recombination records (sample 0)",
            group_col="subgenome",
        )
    return {"n_ok": n_ok, "n_fail": n_fail}

