"""Random genome site panel for test_site mod: popdep grid + variant overlap."""

from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.variant.popdep import (
    _tabix_popdep_lookup,
    ensure_popdep_ref_tabix,
)
from infra.utils.graph import plot_distribution_with_stats, plot_scatter_two_categories
from infra.utils.io import save_df_to_tsv
from infra.wheat.ref_v1 import get_chr_id_length

PLINK_CHR_FILE_IDS: tuple[str, ...] = tuple(f"chr{i:03d}" for i in range(45))
TABIX_BATCH_SIZE = 50_000


@dataclass(frozen=True)
class TestSiteLayout:
    root: Path
    sites: Path
    extract: Path
    variant: Path
    info: Path

    @classmethod
    def under(cls, process_dir: str | Path) -> TestSiteLayout:
        root = Path(process_dir)
        layout = cls(
            root=root,
            sites=root / "sites",
            extract=root / "sites" / "extract",
            variant=root / "variant",
            info=root / "info",
        )
        for path in (layout.sites, layout.extract, layout.variant, layout.info):
            path.mkdir(parents=True, exist_ok=True)
        return layout


def plink_chr_int(chr_file_id: str) -> int:
    return int(chr_file_id.replace("chr", ""))


def list_plink_chr_file_ids(*, chr_exclude: set[int] | None = None) -> list[str]:
    exclude = chr_exclude or set()
    return [cid for cid in PLINK_CHR_FILE_IDS if plink_chr_int(cid) not in exclude]


def _allocate_proportional(counts: dict[str, int], target_n: int) -> dict[str, int]:
    total = sum(counts.values())
    if total <= 0 or target_n <= 0:
        return {k: 0 for k in counts}
    if target_n >= total:
        return dict(counts)

    raw = {k: target_n * counts[k] / total for k in counts}
    alloc = {k: int(raw[k]) for k in counts}
    remainder = target_n - sum(alloc.values())
    order = sorted(counts.keys(), key=lambda k: (raw[k] - alloc[k]), reverse=True)
    for key in order:
        if remainder <= 0:
            break
        if alloc[key] < counts[key]:
            alloc[key] += 1
            remainder -= 1
    return alloc


def sample_random_genome_sites(
    target_bp: int = 10_000_000,
    *,
    seed: int = 42,
    chr_file_ids: list[str] | None = None,
    chr_exclude: set[int] | None = None,
) -> pd.DataFrame:
    """Sample ``target_bp`` random 1 bp positions proportional to chromosome length."""
    chr_ids = chr_file_ids or list_plink_chr_file_ids(chr_exclude=chr_exclude)
    lengths: dict[str, int] = {}
    for cid in chr_ids:
        chrom = plink_chr_int(cid)
        length = int(get_chr_id_length(chrom))
        if length > 0:
            lengths[cid] = length
    if not lengths:
        raise ValueError("No chromosomes available for site sampling")

    alloc = _allocate_proportional(lengths, target_bp)
    rng = random.Random(seed)
    rows: list[dict[str, int | str]] = []
    for cid in chr_ids:
        n = alloc.get(cid, 0)
        if n <= 0:
            continue
        chrom = plink_chr_int(cid)
        length = lengths[cid]
        if n >= length:
            positions = list(range(1, length + 1))
        else:
            positions = rng.sample(range(1, length + 1), n)
        for pos in positions:
            rows.append({"Chr_id": cid, "CHROM": chrom, "Position": pos})

    df = pd.DataFrame(rows)
    if len(df) < target_bp:
        deficit = target_bp - len(df)
        existing = {(r["Chr_id"], r["Position"]) for r in rows}
        pool: list[tuple[str, int, int]] = []
        for cid, length in lengths.items():
            chrom = plink_chr_int(cid)
            taken = alloc.get(cid, 0)
            if taken >= length:
                continue
            for pos in range(1, length + 1):
                if (cid, pos) not in existing:
                    pool.append((cid, chrom, pos))
        rng.shuffle(pool)
        for cid, chrom, pos in pool[:deficit]:
            rows.append({"Chr_id": cid, "CHROM": chrom, "Position": pos})
        df = pd.DataFrame(rows)

    return df.sort_values(["CHROM", "Position"], ignore_index=True)


def scan_variant_positions(
    main_raw_dir: str | Path,
    chr_file_id: str,
    wanted_positions: set[int],
) -> dict[int, str]:
    """Return ``pos -> variant_id`` for positions present in ``{chr}.plink2.pvar``."""
    if not wanted_positions:
        return {}
    pvar_path = Path(main_raw_dir) / f"{chr_file_id}.plink2.pvar"
    if not pvar_path.exists():
        return {}
    max_pos = max(wanted_positions)
    found: dict[int, str] = {}
    with open(pvar_path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 4)
            pos = int(parts[1])
            if pos > max_pos:
                break
            if pos in wanted_positions:
                found[pos] = parts[2]
    return found


def _tabix_lookup_chunked(chrom: int, positions: list[int], popdep_dir: str) -> dict[int, tuple]:
    bgz = ensure_popdep_ref_tabix(f"chr{chrom:03d}", popdep_dir)
    result: dict[int, tuple] = {}
    for start in range(0, len(positions), TABIX_BATCH_SIZE):
        chunk = positions[start : start + TABIX_BATCH_SIZE]
        result.update(_tabix_popdep_lookup(chrom, chunk, bgz))
    return result


def extract_popdep_rows_for_chr(
    chr_file_id: str,
    sites_df: pd.DataFrame,
    popdep_dir: str,
) -> pd.DataFrame:
    """Attach popdep metrics to sampled sites on one PLINK segment."""
    chrom = plink_chr_int(chr_file_id)
    sub = sites_df[sites_df["Chr_id"] == chr_file_id].copy()
    if sub.empty:
        return sub
    positions = sub["Position"].astype(int).tolist()
    try:
        popdep_map = _tabix_lookup_chunked(chrom, positions, popdep_dir)
    except (FileNotFoundError, RuntimeError, OSError) as exc:
        print(f"[Warning] Popdep lookup failed for {chr_file_id}: {exc}")
        popdep_map = {}

    depth_means = []
    depth_sds = []
    rel_means = []
    rel_sds = []
    depth_cvs = []
    rel_cvs = []
    for pos in positions:
        depth = popdep_map.get(pos)
        if depth is None:
            depth_means.append(np.nan)
            depth_sds.append(np.nan)
            rel_means.append(np.nan)
            rel_sds.append(np.nan)
            depth_cvs.append(np.nan)
            rel_cvs.append(np.nan)
            continue
        mean, sd, rel_mean, rel_sd = depth
        depth_means.append(mean)
        depth_sds.append(sd)
        rel_means.append(rel_mean)
        rel_sds.append(rel_sd)
        depth_cvs.append(sd / mean if mean > 0 else np.nan)
        rel_cvs.append(rel_sd / rel_mean if rel_mean > 0 else np.nan)

    sub["Depth_Mean"] = depth_means
    sub["Depth_SD"] = depth_sds
    sub["Depth_CV"] = depth_cvs
    sub["RelativeDepth_Mean"] = rel_means
    sub["RelativeDepth_SD"] = rel_sds
    sub["RelativeDepth_CV"] = rel_cvs
    return sub


def write_variant_extract_list(variant_ids: list[str], output_path: str | Path) -> int:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        for vid in variant_ids:
            handle.write(f"{vid}\n")
    return len(variant_ids)


def prepare_test_site_panel(
    *,
    popdep_dir: str | Path,
    main_raw_dir: str | Path,
    output_dir: str | Path,
    target_bp: int = 10_000_000,
    seed: int = 42,
    chr_exclude: set[int] | None = None,
    force: bool = False,
) -> dict[str, int | float | str]:
    """Sample sites, extract popdep, mark variants, write per-chr artefacts."""
    layout = TestSiteLayout.under(output_dir)
    panel_path = layout.sites / "panel.sites.tsv"
    summary_path = layout.info / "panel.summary.tsv"

    if panel_path.exists() and not force:
        summary = pd.read_csv(summary_path, sep="\t")
        if not summary.empty:
            return summary.iloc[0].to_dict()

    chr_ids = list_plink_chr_file_ids(chr_exclude=chr_exclude)
    sites = sample_random_genome_sites(
        target_bp,
        seed=seed,
        chr_file_ids=chr_ids,
        chr_exclude=chr_exclude,
    )

    all_frames: list[pd.DataFrame] = []
    for cid in chr_ids:
        sub = sites[sites["Chr_id"] == cid]
        if sub.empty:
            write_variant_extract_list([], layout.extract / f"{cid}.variant.txt")
            continue

        pos_set = set(sub["Position"].astype(int).tolist())
        variant_map = scan_variant_positions(main_raw_dir, cid, pos_set)
        write_variant_extract_list(
            list(variant_map.values()),
            layout.extract / f"{cid}.variant.txt",
        )

        popdep_sub = extract_popdep_rows_for_chr(cid, sites, str(popdep_dir))
        popdep_sub["Is_Variant"] = popdep_sub["Position"].astype(int).isin(variant_map)
        popdep_sub["Variant_ID"] = popdep_sub["Position"].astype(int).map(variant_map).fillna(".")

        save_df_to_tsv(popdep_sub, str(layout.sites / f"{cid}.popdep.sites.tsv"))
        all_frames.append(popdep_sub)

    panel = pd.concat(all_frames, ignore_index=True) if all_frames else sites.copy()
    save_df_to_tsv(panel, str(panel_path))

    n_variant = int(panel["Is_Variant"].sum()) if "Is_Variant" in panel.columns else 0
    summary = {
        "target_bp": int(target_bp),
        "n_sites": int(len(panel)),
        "n_variant_sites": n_variant,
        "n_non_variant_sites": int(len(panel) - n_variant),
        "seed": int(seed),
        "popdep_dir": str(popdep_dir),
        "main_raw_dir": str(main_raw_dir),
        "variant_fraction": float(n_variant / len(panel)) if len(panel) else 0.0,
    }
    save_df_to_tsv(pd.DataFrame([summary]), str(summary_path))
    return summary


def load_test_site_panel_sites(process_dir: str | Path) -> pd.DataFrame:
    panel_path = Path(process_dir) / "sites" / "panel.sites.tsv"
    if not panel_path.exists():
        raise FileNotFoundError(f"Missing test_site panel: {panel_path}")
    df = pd.read_csv(panel_path, sep="\t")
    if "Is_Variant" in df.columns:
        df["Is_Variant"] = df["Is_Variant"].astype(str).str.lower().isin(
            {"true", "1", "yes", "t"}
        )
    return df


def _metric_specs() -> tuple[dict, ...]:
    return (
        {
            "mean": "Depth_Mean",
            "cv": "Depth_CV",
            "tag": "depth",
            "label_mean": "Mean depth",
            "label_cv": "Depth CV",
        },
        {
            "mean": "RelativeDepth_Mean",
            "cv": "RelativeDepth_CV",
            "tag": "reldepth",
            "label_mean": "Relative mean depth",
            "label_cv": "Relative depth CV",
        },
    )


def ana_test_site_popdep_plots(
    process_dir: str | Path,
    output_prefix: str,
    *,
    max_along_chr_sites: int = 200_000,
) -> None:
    """Popdep-only QC for test_site: all sampled sites, variants highlighted in plots."""
    df = load_test_site_panel_sites(process_dir)
    req = {"Position", "Depth_Mean", "RelativeDepth_Mean", "Is_Variant"}
    if not req.issubset(df.columns):
        missing = req - set(df.columns)
        raise ValueError(f"panel.sites.tsv missing columns: {sorted(missing)}")

    df = df.replace([np.inf, -np.inf], np.nan)
    df_valid = df.dropna(subset=["Depth_Mean", "RelativeDepth_Mean"], how="any")
    if df_valid.empty:
        print("[Error] No valid popdep rows for test_site plots.")
        return

    non_var = df_valid[~df_valid["Is_Variant"]]
    variants = df_valid[df_valid["Is_Variant"]]

    for spec in _metric_specs():
        mean_col = spec["mean"]
        if mean_col not in df_valid.columns:
            continue
        for subset, suffix in ((non_var, "nonvariant"), (variants, "variant"), (df_valid, "all")):
            if subset.empty:
                continue
            plot_distribution_with_stats(
                subset,
                col=mean_col,
                title=f"{spec['label_mean']} — {suffix} (test_site panel)",
                x_label=spec["label_mean"],
                filename=f"{output_prefix}.{spec['tag']}.{suffix}.mean_dist.png",
            )

    plot_scatter_two_categories(
        df_valid,
        x_col="Depth_Mean",
        y_col="RelativeDepth_Mean",
        category_col="Is_Variant",
        category_labels={False: "Non-variant site", True: "Variant site"},
        category_colors={False: "#bdbdbd", True: "#c44e52"},
        title="Population depth — all sampled sites (test_site)",
        filename=f"{output_prefix}.depth_vs_reldepth_by_variant.png",
        xlabel="Mean depth",
        ylabel="Relative mean depth",
        alpha=0.35,
        s=4,
    )

    for spec in _metric_specs():
        x_col = spec["mean"]
        y_col = spec["cv"]
        if y_col not in df_valid.columns:
            continue
        plot_scatter_two_categories(
            df_valid.dropna(subset=[x_col, y_col]),
            x_col=x_col,
            y_col=y_col,
            category_col="Is_Variant",
            category_labels={False: "Non-variant site", True: "Variant site"},
            category_colors={False: "#bdbdbd", True: "#c44e52"},
            title=f"{spec['label_mean']} vs {spec['label_cv']} (test_site)",
            filename=f"{output_prefix}.{spec['tag']}.mean_vs_cv_by_variant.png",
            xlabel=spec["label_mean"],
            ylabel=spec["label_cv"],
            alpha=0.35,
            s=4,
        )

    for chrom, sub in df_valid.groupby("CHROM", sort=True):
        if len(sub) > max_along_chr_sites:
            sub = sub.sample(n=max_along_chr_sites, random_state=42)
        plot_scatter_two_categories(
            sub.sort_values("Position"),
            x_col="Position",
            y_col="Depth_Mean",
            category_col="Is_Variant",
            category_labels={False: "Non-variant site", True: "Variant site"},
            category_colors={False: "#bdbdbd", True: "#c44e52"},
            title=f"Mean depth along chr{int(chrom):03d} (test_site panel)",
            filename=f"{output_prefix}.chr{int(chrom):03d}.depth_along_coord.png",
            xlabel="Position",
            ylabel="Mean depth",
            alpha=0.4,
            s=3,
        )

    save_df_to_tsv(
        pd.DataFrame(
            [{
                "n_sites_plotted": int(len(df_valid)),
                "n_variant_sites": int(variants.shape[0]),
                "n_non_variant_sites": int(non_var.shape[0]),
            }]
        ),
        f"{output_prefix}.plot_summary.tsv",
    )
