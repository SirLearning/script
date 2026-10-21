"""Join WGS-core and WAP allele frequencies on a shared variant set."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from genetics.genomics.variant.maf_miss_group import _maf_from_afreq_df
from genetics.genomics.variant.variant_utils import load_df_from_plink_variant
from infra.utils.graph import plot_scatter_with_thresholds
from infra.utils.io import save_df_to_tsv, save_thresholds


def _afreq_to_table(path: str | Path, suffix: str) -> pd.DataFrame:
    df = load_df_from_plink_variant(str(path))
    if df is None or df.empty:
        raise ValueError(f"Empty afreq: {path}")
    work = df.copy()
    if "ID" not in work.columns:
        raise ValueError(f"afreq missing ID column: {path}")
    if "ALT_FREQS" in work.columns:
        work["ALT_FREQS"] = pd.to_numeric(work["ALT_FREQS"], errors="coerce")
    if "AF" in work.columns:
        work["AF"] = pd.to_numeric(work["AF"], errors="coerce")
    work = _maf_from_afreq_df(work)
    af_col = "ALT_FREQS" if "ALT_FREQS" in work.columns else "AF"
    out = pd.DataFrame(
        {
            "ID": work["ID"].astype(str),
            f"AF_{suffix}": pd.to_numeric(work[af_col], errors="coerce"),
            f"MAF_{suffix}": pd.to_numeric(work["MAF"], errors="coerce"),
        }
    )
    return out.drop_duplicates(subset=["ID"], keep="first")


def merge_cohort_af(
    wgs_afreq: str | Path,
    wap_afreq: str | Path,
    *,
    rare_max: float = 0.01,
) -> pd.DataFrame:
    """Merge AF/MAF from WGS-core and WAP afreq tables."""

    wgs = _afreq_to_table(wgs_afreq, "wgs")
    wap = _afreq_to_table(wap_afreq, "wap")
    merged = wgs.merge(wap, on="ID", how="inner")
    merged["rare_wgs"] = merged["MAF_wgs"] < rare_max
    merged["rare_wap"] = merged["MAF_wap"] < rare_max
    merged["rare_both"] = merged["rare_wgs"] & merged["rare_wap"]
    merged["rare_discordant"] = merged["rare_wgs"] != merged["rare_wap"]
    return merged


def summarize_dual_af(merged: pd.DataFrame, *, rare_max: float) -> dict[str, int | float]:
    n = len(merged)
    return {
        "n_shared_variants": n,
        "n_rare_wgs": int(merged["rare_wgs"].sum()),
        "n_rare_wap": int(merged["rare_wap"].sum()),
        "n_rare_both": int(merged["rare_both"].sum()),
        "n_rare_discordant": int(merged["rare_discordant"].sum()),
        "rare_max": rare_max,
    }


def plot_dual_af(merged: pd.DataFrame, output_prefix: str | Path, *, rare_max: float) -> None:
    rare_line = {"value": rare_max, "label": f"rare {rare_max}", "color": "#c44e52", "linestyle": "--"}
    plot_scatter_with_thresholds(
        data=merged,
        x_col="MAF_wgs",
        y_col="MAF_wap",
        title="MAF WGS-core vs WAP",
        filename=f"{output_prefix}.maf_wgs_vs_wap.png",
        xlabel="MAF (WGS-core)",
        ylabel="MAF (WAP)",
        thresholds_v=[rare_line],
        thresholds_h=[rare_line],
        identity_line=True,
    )


def run_dual_af(
    wgs_afreq: str | Path,
    wap_afreq: str | Path,
    output_prefix: str | Path,
    *,
    rare_max: float = 0.01,
) -> pd.DataFrame:
    merged = merge_cohort_af(wgs_afreq, wap_afreq, rare_max=rare_max)
    save_df_to_tsv(merged, f"{output_prefix}.dual_af.tsv")
    stats = summarize_dual_af(merged, rare_max=rare_max)
    save_thresholds(stats, f"{output_prefix}.dual_af.th.tsv")
    plot_dual_af(merged, output_prefix, rare_max=rare_max)
    print(f"dual_af: shared={stats['n_shared_variants']} rare_both={stats['n_rare_both']}")
    return merged
