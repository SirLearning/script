"""Collect per-segment PLINK2 .smiss files into segment index and long tables."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.sample.sample_utils import load_df_from_plink2
from genetics.germplasm.sample.mosdepth_aneuploidy import (
    LOGICAL_AUTOSOMES,
    attach_sample_group,
    load_chr_depth_from_dir,
    normalize_sample_id,
)
from infra.utils.graph import plot_correlation_with_regression
from infra.wheat.ref_v1 import get_ref_v1_chr_name


def parse_segment_id_from_smiss_name(name: str) -> int:
    """Parse ``chr012.smiss`` -> 12."""
    stem = Path(name).name
    if not stem.startswith("chr") or not stem.endswith(".smiss"):
        raise ValueError(f"Unexpected smiss filename: {name}")
    return int(stem[3:].removesuffix(".smiss"))


def collect_segment_smiss_tables(
    smiss_dir: str | Path,
    *,
    thin_job: str,
    thin_mod: str,
    output_dir: str | Path = ".",
) -> dict[str, Path]:
    """
    Merge staged ``chrNNN.smiss`` files into:

    - ``sample_chr_miss.segments.tsv`` — one row per PLINK segment
    - ``sample_chr_miss.long.tsv`` — one row per sample × segment
    """
    smiss_dir = Path(smiss_dir)
    output_dir = Path(output_dir)
    rows: list[dict] = []
    long_parts: list[pd.DataFrame] = []

    for smiss_path in sorted(smiss_dir.glob("chr*.smiss")):
        chr_num = parse_segment_id_from_smiss_name(smiss_path.name)
        logical_chr = get_ref_v1_chr_name(chr_num)
        df = load_df_from_plink2(smiss_path)
        if df is None or df.empty:
            continue
        df = df.copy()
        df["plink_chr"] = chr_num
        df["logical_chr"] = logical_chr
        long_parts.append(df)
        rows.append(
            {
                "plink_chr": chr_num,
                "logical_chr": logical_chr,
                "smiss_file": smiss_path.name,
                "n_samples": len(df),
                "mean_f_miss": float(df["F_MISS"].mean()),
                "source_job": thin_job,
                "source_mod": thin_mod,
            }
        )

    segments = pd.DataFrame(rows).sort_values("plink_chr").reset_index(drop=True)
    long_table = (
        pd.concat(long_parts, ignore_index=True)
        if long_parts
        else pd.DataFrame(columns=["Sample", "F_MISS", "plink_chr", "logical_chr"])
    )

    seg_path = output_dir / "sample_chr_miss.segments.tsv"
    long_path = output_dir / "sample_chr_miss.long.tsv"
    segments.to_csv(seg_path, sep="\t", index=False)
    long_table.to_csv(long_path, sep="\t", index=False)
    return {"segments": seg_path, "long_table": long_path}


def aggregate_logical_chr_f_miss(long_table: pd.DataFrame) -> pd.DataFrame:
    """Weight-merge segment ``F_MISS`` to logical chromosomes (chr1A–chr7D)."""
    df = long_table[long_table["logical_chr"].isin(LOGICAL_AUTOSOMES)].copy()
    if df.empty:
        return pd.DataFrame(columns=["sample", "logical_chr", "F_MISS", "MISSING_CT", "OBS_CT"])
    df["sample"] = df["Sample"].map(normalize_sample_id)
    agg = df.groupby(["sample", "logical_chr"], as_index=False).agg(
        MISSING_CT=("MISSING_CT", "sum"),
        OBS_CT=("OBS_CT", "sum"),
    )
    agg["F_MISS"] = agg["MISSING_CT"] / agg["OBS_CT"].replace(0, np.nan)
    return agg


def merge_chr_depth_with_f_miss(
    chr_all: pd.DataFrame,
    long_table: pd.DataFrame,
) -> pd.DataFrame:
    """Inner-join mosdepth chr depth rows with logical-chromosome ``F_MISS``."""
    if chr_all.empty or long_table.empty:
        return pd.DataFrame()
    miss = aggregate_logical_chr_f_miss(long_table)
    depth = chr_all.copy()
    depth["sample_norm"] = depth["sample"].map(normalize_sample_id)
    return depth.merge(
        miss,
        left_on=["sample_norm", "logical_chr"],
        right_on=["sample", "logical_chr"],
        how="inner",
        suffixes=("", "_miss"),
    )


def _evaluable_depth_miss_rows(sub: pd.DataFrame) -> pd.DataFrame:
    rel = pd.to_numeric(sub["rel_depth"], errors="coerce")
    fmiss = pd.to_numeric(sub["F_MISS"], errors="coerce")
    mask = sub["expected"] & np.isfinite(rel) & np.isfinite(fmiss)
    return sub[mask].copy()


def plot_per_chromosome_rel_depth_vs_fmiss(
    merged: pd.DataFrame,
    plot_dir: str | Path,
    group_file: str | Path | None = None,
    *,
    ylim: tuple[float, float] | None = (0.0, 3.0),
    xlim: tuple[float, float] | None = (0.0, 1.0),
) -> list[Path]:
    """Per logical chromosome: ``F_MISS`` (x) vs relative depth (y), colored by ``Group``."""
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)
    if merged.empty:
        return []

    plot_df = attach_sample_group(merged, group_file)
    written: list[Path] = []
    for logical_chr in LOGICAL_AUTOSOMES:
        sub = _evaluable_depth_miss_rows(plot_df[plot_df["logical_chr"] == logical_chr])
        if len(sub) < 2:
            continue
        n_samples = sub["sample_norm"].nunique()
        path = plot_dir / f"{logical_chr}.rel_depth.vs_fmiss.png"
        plot_correlation_with_regression(
            sub,
            x_col="F_MISS",
            y_col="rel_depth",
            title=f"{logical_chr} relative depth vs genotype missing rate (n={n_samples})",
            filename=str(path),
            x_label="Genotype missing rate (F_MISS, test_thin)",
            y_label="Relative depth (chr depth / genome-wide mean)",
            xlim=xlim,
            ylim=ylim,
            group_col="Group",
        )
        written.append(path)
    return written


def plot_rel_depth_vs_fmiss_from_tables(
    chr_depth_dir: str | Path,
    smiss_long_path: str | Path,
    plot_dir: str | Path,
    group_file: str | Path | None = None,
    *,
    merged_out: str | Path | None = None,
    ylim: tuple[float, float] | None = (0.0, 3.0),
    xlim: tuple[float, float] | None = (0.0, 1.0),
) -> list[Path]:
    """Load depth + smiss tables, optionally write merged TSV, and plot per-chromosome scatters."""
    chr_all = load_chr_depth_from_dir(chr_depth_dir)
    long_table = pd.read_csv(smiss_long_path, sep="\t")
    if chr_all.empty:
        raise FileNotFoundError(f"No chr depth tables under {chr_depth_dir}")
    if long_table.empty:
        raise FileNotFoundError(f"Empty smiss long table: {smiss_long_path}")

    merged = merge_chr_depth_with_f_miss(chr_all, long_table)
    if merged.empty:
        raise ValueError("No overlapping sample × logical_chr rows between depth and F_MISS")

    if merged_out is not None:
        out = Path(merged_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(out, sep="\t", index=False)

    return plot_per_chromosome_rel_depth_vs_fmiss(
        merged,
        plot_dir,
        group_file=group_file,
        ylim=ylim,
        xlim=xlim,
    )
