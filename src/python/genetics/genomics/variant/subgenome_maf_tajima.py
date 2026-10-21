"""Stacked A/B/D MAF and Tajima's D histograms from PLINK2 afreq and vcftools tables."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.variant.maf_miss_group import _maf_from_afreq_df
from infra.utils.graph import (
    SUBGENOME_COLORS,
    SUBGENOME_ORDER,
    plot_faceted_histograms,
)
from infra.utils.io import load_df_from_tsv, save_df_to_tsv


def _strip_hash_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [
        c.replace("#", "") if isinstance(c, str) else c for c in out.columns
    ]
    return out


def _subgenome_from_path(path: str | Path) -> str:
    return Path(path).name.split(".")[0]


def load_maf_from_afreq(path: str | Path, subgenome: str | None = None) -> pd.DataFrame:
    """Load PLINK2 .afreq and return MAF plus subgenome label (no chr remapping)."""
    df = _strip_hash_columns(load_df_from_tsv(path))
    df = _maf_from_afreq_df(df)
    sg = subgenome if subgenome is not None else _subgenome_from_path(path)
    df["subgenome"] = sg
    return df[["MAF", "subgenome"]].copy()


def load_tajima_d_table(
    path: str | Path, subgenome: str | None = None
) -> pd.DataFrame:
    """Load a vcftools ``*.Tajima.D`` table and drop non-finite D values."""
    df = _strip_hash_columns(load_df_from_tsv(path))
    required = {"CHROM", "BIN_START", "N_SNPS", "TajimaD"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{path} missing Tajima.D columns: {sorted(missing)}")
    df["BIN_START"] = pd.to_numeric(df["BIN_START"], errors="coerce")
    df["N_SNPS"] = pd.to_numeric(df["N_SNPS"], errors="coerce")
    df["TajimaD"] = pd.to_numeric(df["TajimaD"], errors="coerce")
    df = df.dropna(subset=["TajimaD"])
    sg = subgenome if subgenome is not None else _subgenome_from_path(path)
    df["subgenome"] = sg
    return df[["CHROM", "BIN_START", "N_SNPS", "TajimaD", "subgenome"]].copy()


def _paths_by_subgenome(paths) -> dict[str, str]:
    mapping = {}
    for path in paths:
        sg = _subgenome_from_path(path)
        if sg in SUBGENOME_ORDER:
            mapping[sg] = str(path)
    missing = [sg for sg in SUBGENOME_ORDER if sg not in mapping]
    if missing:
        raise ValueError(f"missing subgenome inputs: {missing}")
    return mapping


def summarize_numeric_by_subgenome(
    df: pd.DataFrame, value_col: str
) -> pd.DataFrame:
    rows = []
    for sg in SUBGENOME_ORDER:
        vals = pd.to_numeric(
            df.loc[df["subgenome"] == sg, value_col], errors="coerce"
        ).replace([np.inf, -np.inf], np.nan).dropna()
        rows.append(
            {
                "subgenome": sg,
                "n": int(len(vals)),
                "mean": float(vals.mean()) if len(vals) else np.nan,
                "median": float(vals.median()) if len(vals) else np.nan,
                "sd": float(vals.std(ddof=1)) if len(vals) > 1 else np.nan,
                "min": float(vals.min()) if len(vals) else np.nan,
                "max": float(vals.max()) if len(vals) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def plot_subgenome_stacked_distribution(
    df: pd.DataFrame,
    value_col: str,
    output_path: str | Path,
    title: str,
    xlabel: str,
    *,
    log_scale: bool = False,
    bins: int = 100,
    xlim=None,
    sharey: bool = True,
) -> None:
    plot_faceted_histograms(
        data=df,
        x_col=value_col,
        facet_col="subgenome",
        facet_order=list(SUBGENOME_ORDER),
        title=title,
        filename=str(output_path),
        xlabel=xlabel,
        ylabel="Count",
        log_scale=log_scale,
        bins=bins,
        xlim=xlim,
        facet_colors=SUBGENOME_COLORS,
        facet_labels={sg: f"Subgenome {sg}" for sg in SUBGENOME_ORDER},
        sharey=sharey,
    )


def run_subgenome_maf_plots(afreq_paths, output_prefix: str | Path) -> dict[str, Path]:
    """Stacked MAF histograms (linear y and log y) plus a per-subgenome summary TSV."""
    prefix = Path(output_prefix)
    mapping = _paths_by_subgenome(afreq_paths)
    df = pd.concat(
        [load_maf_from_afreq(mapping[sg], sg) for sg in SUBGENOME_ORDER],
        ignore_index=True,
    )
    summary = summarize_numeric_by_subgenome(df, "MAF")
    summary_path = Path(f"{prefix}.subgenome.th.tsv")
    dist_path = Path(f"{prefix}.subgenome.dist.png")
    log_path = Path(f"{prefix}.subgenome.dist.log.png")
    save_df_to_tsv(summary, str(summary_path))
    plot_subgenome_stacked_distribution(
        df,
        "MAF",
        dist_path,
        title="Minor allele frequency by subgenome",
        xlabel="MAF",
        log_scale=False,
        bins=100,
        xlim=(0.0, 0.5),
    )
    plot_subgenome_stacked_distribution(
        df,
        "MAF",
        log_path,
        title="Minor allele frequency by subgenome (log counts)",
        xlabel="MAF",
        log_scale=True,
        bins=100,
        xlim=(0.0, 0.5),
    )
    return {
        "summary": summary_path,
        "dist": dist_path,
        "dist_log": log_path,
    }


def run_subgenome_tajima_plots(
    tajima_paths, output_prefix: str | Path
) -> dict[str, Path]:
    """Stacked Tajima's D histograms plus a per-subgenome summary TSV."""
    prefix = Path(output_prefix)
    mapping = _paths_by_subgenome(tajima_paths)
    df = pd.concat(
        [load_tajima_d_table(mapping[sg], sg) for sg in SUBGENOME_ORDER],
        ignore_index=True,
    )
    summary = summarize_numeric_by_subgenome(df, "TajimaD")
    summary_path = Path(f"{prefix}.subgenome.th.tsv")
    dist_path = Path(f"{prefix}.subgenome.dist.png")
    combined_path = Path(f"{prefix}.subgenome.info.tsv")
    save_df_to_tsv(summary, str(summary_path))
    save_df_to_tsv(df, str(combined_path))
    plot_subgenome_stacked_distribution(
        df,
        "TajimaD",
        dist_path,
        title="Tajima's D by subgenome",
        xlabel="Tajima's D",
        log_scale=False,
        bins=80,
        xlim=None,
    )
    return {
        "summary": summary_path,
        "dist": dist_path,
        "info": combined_path,
    }
