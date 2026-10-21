"""Random genome site panel for test_site mod: popdep grid + variant overlap."""

from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.variant.popdep import (
    _stream_popdep_lookup,
    _tabix_popdep_lookup,
    ensure_popdep_ref_tabix,
    popdep_chr_ref_path,
)
from genetics.genomics.variant.popdep_region_class import (
    REGION_CLASS_ALPHA,
    REGION_CLASS_COLORS,
    REGION_CLASS_HUE_ORDER,
    REGION_CLASS_LEGEND_LABELS,
    REGION_CLASSES,
    classify_positions,
    load_class_interval_sets,
)
from infra.utils.graph import (
    plot_distribution_with_stats,
    plot_joint_regression,
    plot_scatter_two_categories,
    plot_stacked_hist_counts,
)
from infra.utils.io import load_df_generic, save_df_to_tsv
from infra.wheat.ref_v1 import REF_V1_SUBGENOME_SEGMENT_CHR_IDS, get_chr_id_length

PLINK_CHR_FILE_IDS: tuple[str, ...] = tuple(f"chr{i:03d}" for i in range(45))
TABIX_BATCH_SIZE = 50_000
STREAM_LOOKUP_MIN_SITES = 100_000
TEST_POPDEP_SCATTER_CAP = 2_000_000
TEST_POPDEP_TARGET_N = 100_000_000
TEST_POPDEP_HIST_BINS = 100
TEST_POPDEP_DEPTH_XLIM = (0.0, 40.0)
TEST_POPDEP_REL_XLIM = (0.0, 5.0)
N500_PANEL_N = 500
# June 2026 31run PopDepCrossChr chromosomeLength nTaxa while scanning 500 BAMs
# (missing taxa zero-padded). chr032 used PopDepFull with nTaxa=500.
NTAXA_ZERO_PADDED_BY_SUBGENOME = {
    "A": 7675,
    "B": 7594,
    "D": 7738,
    "Others": 8285,
}
_CHROM_TO_SUBGENOME = {
    int(cid): sg
    for sg, ids in REF_V1_SUBGENOME_SEGMENT_CHR_IDS.items()
    for cid in ids
}
_MEAN_SD_CV_COLS = (
    ("Depth_Mean", "Depth_SD", "Depth_CV"),
    ("RelativeDepth_Mean", "RelativeDepth_SD", "RelativeDepth_CV"),
)


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


def n500_crosschr_ntaxa_used(chrom: int, *, n_new: int = N500_PANEL_N) -> int:
    """Return the nTaxa written into June CrossChr grids for this chromosome.

    chr032 was computed with PopDepFull at the true panel size and is not scaled.
    """
    chrom_i = int(chrom)
    if chrom_i == 32:
        return int(n_new)
    sg = _CHROM_TO_SUBGENOME.get(chrom_i)
    if sg is None:
        raise ValueError(f"Unknown chromosome id: {chrom_i}")
    return NTAXA_ZERO_PADDED_BY_SUBGENOME[sg]


def rescale_zero_padded_mean_sd(
    mean,
    sd,
    *,
    n_old: int,
    n_new: int = N500_PANEL_N,
):
    """Invert mean/SD after n_old-n_new taxa were included as zeros.

    TIGER reports a population SD (divide by n, not n-1). With padded zeros,
    ``sum`` and ``sum of squares`` are unchanged, so::

        mean_new = mean_old * n_old / n_new
        sd_new^2 = (n_old / n_new) * (sd_old^2 + mean_old^2) - mean_new^2
    """
    if n_old <= 0 or n_new <= 0:
        raise ValueError("n_old and n_new must be positive")
    mean_arr = np.asarray(mean, dtype=float)
    sd_arr = np.asarray(sd, dtype=float)
    scalar = np.ndim(mean) == 0 and np.ndim(sd) == 0
    if n_old == n_new:
        out_mean, out_sd = mean_arr.copy(), sd_arr.copy()
    else:
        scale = n_old / n_new
        out_mean = mean_arr * scale
        second = np.square(sd_arr) + np.square(mean_arr)
        var_new = np.maximum(scale * second - np.square(out_mean), 0.0)
        out_sd = np.sqrt(var_new)
    if scalar:
        return float(out_mean), float(out_sd)
    return out_mean, out_sd


def _cv_from_mean_sd(mean, sd):
    mean_arr = np.asarray(mean, dtype=float)
    sd_arr = np.asarray(sd, dtype=float)
    out = np.full(mean_arr.shape, np.nan, dtype=float)
    ok = np.isfinite(mean_arr) & (mean_arr > 0) & np.isfinite(sd_arr)
    out[ok] = sd_arr[ok] / mean_arr[ok]
    return out


def rescale_test_popdep_sites(
    input_tsv: str | Path,
    output_tsv: str | Path,
    *,
    n_new: int = N500_PANEL_N,
) -> pd.DataFrame:
    """Rescale one chromosome test_popdep site table from padded nTaxa to n_new."""
    df = load_df_generic(str(input_tsv))
    if df.empty:
        save_df_to_tsv(df, str(output_tsv))
        return df
    if "CHROM" not in df.columns:
        raise ValueError(f"{input_tsv} has no CHROM column")
    chroms = pd.unique(df["CHROM"].dropna())
    if len(chroms) != 1:
        raise ValueError(f"{input_tsv} expected one CHROM, got {list(chroms)}")
    n_old = n500_crosschr_ntaxa_used(int(chroms[0]), n_new=n_new)
    print(
        f"rescale_test_popdep_sites {input_tsv}: n_rows={len(df)} "
        f"CHROM={int(chroms[0])} n_old={n_old} n_new={n_new}"
    )
    if n_old != n_new:
        for mean_col, sd_col, cv_col in _MEAN_SD_CV_COLS:
            if mean_col not in df.columns or sd_col not in df.columns:
                continue
            mean_new, sd_new = rescale_zero_padded_mean_sd(
                df[mean_col].to_numpy(),
                df[sd_col].to_numpy(),
                n_old=n_old,
                n_new=n_new,
            )
            df[mean_col] = mean_new
            df[sd_col] = sd_new
            if cv_col in df.columns:
                df[cv_col] = _cv_from_mean_sd(mean_new, sd_new)
    save_df_to_tsv(df, str(output_tsv))
    return df


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
        if len(positions) >= STREAM_LOOKUP_MIN_SITES:
            popdep_map = _stream_popdep_lookup(
                popdep_chr_ref_path(popdep_dir, chrom),
                positions,
            )
        else:
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


def sample_chr_sites(chr_file_id: str, n: int, *, seed: int) -> pd.DataFrame:
    """Sample ``n`` distinct 1-based positions on one PLINK segment."""
    chrom = plink_chr_int(chr_file_id)
    length = int(get_chr_id_length(chrom) or 0)
    if n <= 0 or length <= 0:
        return pd.DataFrame(columns=["Chr_id", "CHROM", "Position"])
    take = min(int(n), length)
    rng = random.Random(int(seed) + chrom * 1_000_003)
    if take >= length:
        positions = list(range(1, length + 1))
    else:
        positions = rng.sample(range(1, length + 1), take)
    positions.sort()
    return pd.DataFrame(
        {
            "Chr_id": [chr_file_id] * take,
            "CHROM": [chrom] * take,
            "Position": positions,
        }
    )


def allocate_test_popdep_sites(
    target_bp: int,
    *,
    chr_file_ids: list[str] | None = None,
    chr_exclude: set[int] | None = None,
) -> pd.DataFrame:
    """Proportional site counts per PLINK segment."""
    chr_ids = chr_file_ids or list_plink_chr_file_ids(chr_exclude=chr_exclude)
    lengths: dict[str, int] = {}
    for cid in chr_ids:
        length = int(get_chr_id_length(plink_chr_int(cid)) or 0)
        if length > 0:
            lengths[cid] = length
    alloc = _allocate_proportional(lengths, int(target_bp))
    rows = [
        {"Chr_id": cid, "CHROM": plink_chr_int(cid), "n_sites": int(alloc.get(cid, 0))}
        for cid in chr_ids
    ]
    return pd.DataFrame(rows)


def prepare_test_popdep_chr_panel(
    chr_file_id: str,
    n_sites: int,
    *,
    popdep_dir: str | Path,
    annotation_dir: str | Path,
    output_tsv: str | Path,
    seed: int = 42,
) -> pd.DataFrame:
    """Sample sites on one chromosome, attach popdep, assign region_class."""
    sites = sample_chr_sites(chr_file_id, int(n_sites), seed=seed)
    if sites.empty:
        empty = pd.DataFrame(
            columns=[
                "Chr_id",
                "CHROM",
                "Position",
                "Depth_Mean",
                "Depth_SD",
                "Depth_CV",
                "RelativeDepth_Mean",
                "RelativeDepth_SD",
                "RelativeDepth_CV",
                "region_class",
            ]
        )
        save_df_to_tsv(empty, str(output_tsv))
        return empty

    popdep_sub = extract_popdep_rows_for_chr(chr_file_id, sites, str(popdep_dir))
    class_sets = load_class_interval_sets(annotation_dir)
    popdep_sub["region_class"] = classify_positions(
        plink_chr_int(chr_file_id),
        popdep_sub["Position"].to_numpy(),
        class_sets,
    )
    save_df_to_tsv(popdep_sub, str(output_tsv))
    return popdep_sub


def concat_test_popdep_chr_tables(paths) -> pd.DataFrame:
    frames = []
    for path in paths:
        table = pd.read_csv(path, sep="\t")
        if not table.empty:
            frames.append(table)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def hist_bin_edges(xlim: tuple[float, float], n_bins: int) -> np.ndarray:
    lo, hi = float(xlim[0]), float(xlim[1])
    if hi <= lo:
        raise ValueError(f"xlim must be increasing, got {xlim}")
    return np.linspace(lo, hi, int(n_bins) + 1)


def hist_counts_by_group(
    values,
    groups,
    bin_edges,
    *,
    group_order: tuple[str, ...] = REGION_CLASSES,
) -> dict[str, np.ndarray]:
    """Bin ``values`` per group; clip into ``bin_edges`` so every finite value is counted."""
    edges = np.asarray(bin_edges, dtype=float)
    n_bins = int(edges.size - 1)
    out = {g: np.zeros(n_bins, dtype=np.int64) for g in group_order}
    x = pd.to_numeric(pd.Series(values, dtype="object"), errors="coerce").to_numpy(
        dtype=float, copy=False
    )
    g = pd.Series(groups, dtype="object").astype(str).to_numpy()
    ok = np.isfinite(x)
    if not ok.any():
        return out
    x = np.clip(x[ok], edges[0], edges[-1])
    g = g[ok]
    for label in group_order:
        mask = g == label
        if mask.any():
            out[label] = np.histogram(x[mask], bins=edges)[0].astype(np.int64)
    return out


def _median_from_hist(bin_edges, counts) -> float | None:
    total = int(np.sum(counts))
    if total <= 0:
        return None
    cum = np.cumsum(counts)
    idx = int(np.searchsorted(cum, (total + 1) // 2))
    idx = min(max(idx, 0), len(counts) - 1)
    return float(0.5 * (bin_edges[idx] + bin_edges[idx + 1]))


def _iter_popdep_plot_tables(panel):
    if isinstance(panel, pd.DataFrame):
        yield panel
        return
    if isinstance(panel, (str, Path)):
        yield load_df_generic(str(panel))
        return
    for item in panel:
        if isinstance(item, pd.DataFrame):
            yield item
        else:
            yield load_df_generic(str(item))


def _scatter_quota(n_rows: int, scatter_cap: int) -> int:
    if n_rows <= 0 or scatter_cap <= 0:
        return 0
    take = int(round(n_rows * scatter_cap / float(TEST_POPDEP_TARGET_N)))
    if take < 1:
        take = min(n_rows, int(scatter_cap))
    return min(n_rows, int(scatter_cap), take)


def _stratified_subsample(
    df: pd.DataFrame,
    *,
    group_col: str,
    n: int,
    seed: int = 42,
) -> pd.DataFrame:
    if len(df) <= n:
        return df
    counts = df[group_col].value_counts()
    alloc = _allocate_proportional(counts.to_dict(), n)
    parts = []
    for label, take in alloc.items():
        subset = df[df[group_col] == label]
        if take <= 0 or subset.empty:
            continue
        if take >= len(subset):
            parts.append(subset)
        else:
            parts.append(subset.sample(n=int(take), random_state=seed))
    if not parts:
        return df.sample(n=n, random_state=seed)
    return pd.concat(parts, ignore_index=True)


def _fill_missing_cv(
    scatter: pd.DataFrame,
    mean_col: str,
    sd_col: str,
    cv_col: str,
) -> None:
    """Fill CV from SD/mean when the column is absent or NaN (in-place)."""
    if mean_col not in scatter.columns or sd_col not in scatter.columns:
        return
    if cv_col not in scatter.columns:
        scatter[cv_col] = np.nan
    cv = pd.to_numeric(scatter[cv_col], errors="coerce")
    miss = cv.isna()
    if not bool(miss.any()):
        return
    mean = pd.to_numeric(scatter.loc[miss, mean_col], errors="coerce")
    sd = pd.to_numeric(scatter.loc[miss, sd_col], errors="coerce")
    scatter.loc[miss, cv_col] = np.where(mean > 0, sd / mean, np.nan)


def _plot_popdep_mean_sd_joints(
    scatter: pd.DataFrame,
    *,
    output_prefix: str,
    mean_col: str,
    sd_col: str,
    cv_col: str,
    file_prefix: str,
    x_label: str,
    sd_label: str,
    cv_label: str,
) -> None:
    """JointGrid SD/CV vs mean for one depth scale (absolute or relative)."""
    if mean_col not in scatter.columns:
        return
    for y_col, y_label, tag in (
        (sd_col, sd_label, "sd"),
        (cv_col, cv_label, "cv"),
    ):
        if y_col not in scatter.columns:
            continue
        plot_df = scatter.dropna(subset=[mean_col, y_col])
        if plot_df.empty:
            continue
        y_vals = pd.to_numeric(plot_df[y_col], errors="coerce")
        if y_vals.nunique(dropna=True) < 2:
            continue
        plot_joint_regression(
            plot_df,
            x_col=mean_col,
            y_col=y_col,
            group_col="region_class",
            x_label=x_label,
            y_label=y_label,
            filename=f"{output_prefix}.{file_prefix}.{tag}_vs_mean_by_class.png",
            title=f"{y_label} vs {x_label.lower()} by genomic class (test_popdep)",
            figure_size=(10, 8),
            palette=REGION_CLASS_COLORS,
            hue_order=list(REGION_CLASS_HUE_ORDER),
            draw_regression=False,
            alpha=0.45,
            s=12,
            rasterized=True,
            legend_title="Genomic class",
            bins_x=250,
            bins_y=150,
            min_points=1,
            marginal_multiple="fill",
            stack_last_on_top=True,
            hue_alpha=REGION_CLASS_ALPHA,
            also_stack_counts=True,
            hue_labels=REGION_CLASS_LEGEND_LABELS,
        )
        pos = plot_df[(plot_df[mean_col] > 0) & (plot_df[y_col] > 0)]
        if pos.empty:
            continue
        plot_joint_regression(
            pos,
            x_col=mean_col,
            y_col=y_col,
            group_col="region_class",
            x_label=x_label,
            y_label=y_label,
            filename=f"{output_prefix}.{file_prefix}.{tag}_vs_mean_by_class.log.png",
            title=f"{y_label} vs {x_label.lower()} by genomic class (log, test_popdep)",
            figure_size=(10, 8),
            palette=REGION_CLASS_COLORS,
            hue_order=list(REGION_CLASS_HUE_ORDER),
            log_x=True,
            log_y=True,
            draw_regression=False,
            alpha=0.45,
            s=12,
            rasterized=True,
            legend_title="Genomic class",
            bins_x=250,
            bins_y=150,
            min_points=1,
            marginal_multiple="fill",
            stack_last_on_top=True,
            hue_alpha=REGION_CLASS_ALPHA,
            also_stack_counts=True,
            hue_labels=REGION_CLASS_LEGEND_LABELS,
        )


def ana_test_popdep_class_plots(
    panel: pd.DataFrame | str | Path | list,
    output_prefix: str,
    *,
    scatter_cap: int = TEST_POPDEP_SCATTER_CAP,
    hist_bins: int = TEST_POPDEP_HIST_BINS,
    depth_xlim: tuple[float, float] = TEST_POPDEP_DEPTH_XLIM,
    rel_xlim: tuple[float, float] = TEST_POPDEP_REL_XLIM,
    joints_only: bool = False,
) -> None:
    """Popdep-only QC coloured by gap/repeat/gc/cnv/normal.

    Histograms use pre-binned counts (clip to ``depth_xlim`` / ``rel_xlim``).
    Scatter JointGrids use a stratified subsample. Paths are streamed one
    chromosome table at a time so the full 100M-row panel is never concatenated.
    ``joints_only`` redraws SD/CV JointGrids from an existing scatter table and
    skips the 100M class-mean histograms.
    """
    if joints_only:
        scatter_parts = list(_iter_popdep_plot_tables(panel))
        scatter = (
            pd.concat(scatter_parts, ignore_index=True) if scatter_parts else pd.DataFrame()
        )
        if scatter.empty:
            print("[Error] joints_only: no scatter rows.")
            return
        scatter = scatter.copy()
        if "region_class" in scatter.columns:
            scatter["region_class"] = scatter["region_class"].astype(str)
        _fill_missing_cv(scatter, "Depth_Mean", "Depth_SD", "Depth_CV")
        _fill_missing_cv(scatter, "RelativeDepth_Mean", "RelativeDepth_SD", "RelativeDepth_CV")
        _plot_popdep_mean_sd_joints(
            scatter,
            output_prefix=output_prefix,
            mean_col="Depth_Mean",
            sd_col="Depth_SD",
            cv_col="Depth_CV",
            file_prefix="depth",
            x_label="Mean depth",
            sd_label="Depth SD",
            cv_label="Depth CV",
        )
        _plot_popdep_mean_sd_joints(
            scatter,
            output_prefix=output_prefix,
            mean_col="RelativeDepth_Mean",
            sd_col="RelativeDepth_SD",
            cv_col="RelativeDepth_CV",
            file_prefix="reldepth",
            x_label="Relative mean depth",
            sd_label="Relative depth SD",
            cv_label="Relative depth CV",
        )
        save_df_to_tsv(
            pd.DataFrame([{"n_scatter_sites": int(len(scatter)), "joints_only": 1}]),
            f"{output_prefix}.plot_summary.tsv",
        )
        return

    depth_edges = hist_bin_edges(depth_xlim, hist_bins)
    rel_edges = hist_bin_edges(rel_xlim, hist_bins)
    depth_hist = {cls: np.zeros(hist_bins, dtype=np.int64) for cls in REGION_CLASSES}
    rel_hist = {cls: np.zeros(hist_bins, dtype=np.int64) for cls in REGION_CLASSES}
    class_n = {cls: 0 for cls in REGION_CLASSES}
    n_sites = 0
    depth_sum = 0.0
    rel_sum = 0.0
    n_rel = 0
    scatter_parts: list[pd.DataFrame] = []

    for table in _iter_popdep_plot_tables(panel):
        req = {"Depth_Mean", "Depth_SD", "region_class"}
        if not req.issubset(table.columns):
            missing = req - set(table.columns)
            raise ValueError(f"test_popdep panel missing columns: {sorted(missing)}")
        df = table.replace([np.inf, -np.inf], np.nan)
        df["region_class"] = df["region_class"].astype(str)
        valid = df.dropna(subset=["Depth_Mean", "Depth_SD"])
        if valid.empty:
            continue
        n_sites += int(len(valid))
        depth_sum += float(pd.to_numeric(valid["Depth_Mean"], errors="coerce").sum())
        for cls, n_cls in valid["region_class"].value_counts().items():
            if cls in class_n:
                class_n[cls] += int(n_cls)
        add_d = hist_counts_by_group(
            valid["Depth_Mean"], valid["region_class"], depth_edges
        )
        for cls in REGION_CLASSES:
            depth_hist[cls] += add_d[cls]
        if "RelativeDepth_Mean" in valid.columns:
            rel = valid.dropna(subset=["RelativeDepth_Mean"])
            if not rel.empty:
                n_rel += int(len(rel))
                rel_sum += float(pd.to_numeric(rel["RelativeDepth_Mean"], errors="coerce").sum())
                add_r = hist_counts_by_group(
                    rel["RelativeDepth_Mean"], rel["region_class"], rel_edges
                )
                for cls in REGION_CLASSES:
                    rel_hist[cls] += add_r[cls]
        take = _scatter_quota(len(valid), int(scatter_cap))
        if take > 0:
            if take >= len(valid):
                scatter_parts.append(valid)
            else:
                scatter_parts.append(valid.sample(n=take, random_state=42))
        chrom = int(valid["CHROM"].iloc[0]) if "CHROM" in valid.columns else -1
        print(f"binned CHROM={chrom} n_valid={len(valid)} scatter_take={take} n_sites_cum={n_sites}")

    if n_sites <= 0:
        print("[Error] No valid popdep rows for test_popdep plots.")
        return

    class_counts = pd.DataFrame(
        {"region_class": list(REGION_CLASSES), "n_sites": [class_n[c] for c in REGION_CLASSES]}
    )
    save_df_to_tsv(class_counts, f"{output_prefix}.class_counts.tsv")

    depth_total = sum(depth_hist.values())
    rel_total = sum(rel_hist.values()) if n_rel else None
    plot_stacked_hist_counts(
        depth_edges,
        depth_hist,
        title="Mean depth by genomic class (test_popdep)",
        filename=f"{output_prefix}.depth.mean_dist_by_class.png",
        x_label="Mean depth",
        hue_order=list(REGION_CLASS_HUE_ORDER),
        palette=REGION_CLASS_COLORS,
        hue_labels=REGION_CLASS_LEGEND_LABELS,
        mean_val=depth_sum / n_sites if n_sites else None,
        median_val=_median_from_hist(depth_edges, depth_total),
    )
    if n_rel > 0:
        plot_stacked_hist_counts(
            rel_edges,
            rel_hist,
            title="Relative mean depth by genomic class (test_popdep)",
            filename=f"{output_prefix}.reldepth.mean_dist_by_class.png",
            x_label="Relative mean depth",
            hue_order=list(REGION_CLASS_HUE_ORDER),
            palette=REGION_CLASS_COLORS,
            hue_labels=REGION_CLASS_LEGEND_LABELS,
            mean_val=rel_sum / n_rel,
            median_val=_median_from_hist(rel_edges, rel_total),
        )

    scatter = (
        pd.concat(scatter_parts, ignore_index=True) if scatter_parts else pd.DataFrame()
    )
    if not scatter.empty and len(scatter) > int(scatter_cap):
        scatter = _stratified_subsample(
            scatter, group_col="region_class", n=int(scatter_cap), seed=42
        )
    if scatter.empty:
        print("[Warning] No scatter rows after subsample.")
        save_df_to_tsv(
            pd.DataFrame(
                [{"n_sites": n_sites, "n_scatter_sites": 0, **{f"n_{c}": class_n[c] for c in REGION_CLASSES}}]
            ),
            f"{output_prefix}.plot_summary.tsv",
        )
        return
    scatter = scatter.copy()
    _fill_missing_cv(scatter, "Depth_Mean", "Depth_SD", "Depth_CV")
    _fill_missing_cv(scatter, "RelativeDepth_Mean", "RelativeDepth_SD", "RelativeDepth_CV")
    save_df_to_tsv(scatter, f"{output_prefix}.scatter.sites.tsv")
    _plot_popdep_mean_sd_joints(
        scatter,
        output_prefix=output_prefix,
        mean_col="Depth_Mean",
        sd_col="Depth_SD",
        cv_col="Depth_CV",
        file_prefix="depth",
        x_label="Mean depth",
        sd_label="Depth SD",
        cv_label="Depth CV",
    )
    _plot_popdep_mean_sd_joints(
        scatter,
        output_prefix=output_prefix,
        mean_col="RelativeDepth_Mean",
        sd_col="RelativeDepth_SD",
        cv_col="RelativeDepth_CV",
        file_prefix="reldepth",
        x_label="Relative mean depth",
        sd_label="Relative depth SD",
        cv_label="Relative depth CV",
    )

    save_df_to_tsv(
        pd.DataFrame(
            [
                {
                    "n_sites": int(n_sites),
                    "n_scatter_sites": int(len(scatter)),
                    "hist_bins": int(hist_bins),
                    "depth_xlim_hi": float(depth_xlim[1]),
                    **{f"n_{cls}": int(class_n[cls]) for cls in REGION_CLASSES},
                }
            ]
        ),
        f"{output_prefix}.plot_summary.tsv",
    )
