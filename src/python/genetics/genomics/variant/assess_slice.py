"""
PLINK2 slice assess: join .afreq + .vmiss and emit QC tables and standard infra plots.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from infra.utils.graph import (
    plot_bar_chart,
    plot_distribution_with_stats,
    plot_scatter_with_thresholds,
)
from infra.utils.io import load_df_generic, save_df_to_tsv


def _maf_from_alt_freqs(series: pd.Series) -> pd.Series:
    """Minor allele frequency from PLINK2 ALT_FREQS (allele frequency of ALT)."""

    def one(v) -> float:
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return float("nan")
        s = str(v).strip()
        if s in ("", ".", "NA"):
            return float("nan")
        try:
            af = float(s)
        except ValueError:
            return float("nan")
        return af if af <= 0.5 else 1.0 - af

    return series.map(one)


def ana_assess_plink_debug_slice(
    afreq_path: str,
    vmiss_path: str,
    mac_hist_path: str,
    counts_path: str,
    output_prefix: str,
) -> None:
    """
    Read PLINK2 assess slice outputs, save a joined variant table, and write plots
    using infra helpers (same conventions as variant stats modules).

    Parameters
    ----------
    afreq_path
        PLINK2 --freq .afreq file.
    vmiss_path
        PLINK2 --missing .vmiss file.
    mac_hist_path
        Two-column TSV ``MAF_bin\\tn_sites`` from the Nextflow awk helper.
    counts_path
        ``id\\tn_variants\\tn_samples`` TSV from the slice process.
    output_prefix
        Basename prefix for ``*.info.tsv``, ``*.th.tsv``, and ``*.png`` outputs.
    """
    df_a = load_df_generic(afreq_path)
    df_v = load_df_generic(vmiss_path)
    if df_a is None or df_v is None:
        print("[Error] Missing afreq or vmiss input.")
        return

    df_a.columns = [str(c).replace("#", "") for c in df_a.columns]
    df_v.columns = [str(c).replace("#", "") for c in df_v.columns]

    if "ID" not in df_a.columns or "ID" not in df_v.columns:
        print("[Error] Expected ID column in afreq and vmiss.")
        return

    freq_col = "ALT_FREQS" if "ALT_FREQS" in df_a.columns else None
    if freq_col is None:
        print("[Error] ALT_FREQS column not found in afreq.")
        return

    miss_col = "F_MISS" if "F_MISS" in df_v.columns else None
    if miss_col is None:
        print("[Error] F_MISS column not found in vmiss.")
        return

    df_a = df_a.copy()
    df_a["MAF"] = _maf_from_alt_freqs(df_a[freq_col])
    sub_v = df_v[["ID", miss_col]].copy()

    merged = df_a.merge(sub_v, on="ID", how="inner", suffixes=("", "_vmiss"))
    merged = merged.rename(columns={miss_col: "F_MISS"})

    info_path = f"{output_prefix}.maf_miss.info.tsv"
    save_df_to_tsv(merged, info_path)

    df_plot = merged.dropna(subset=["MAF", "F_MISS"]).copy()
    if df_plot.empty:
        print("[Warning] No rows after dropping NA MAF/F_MISS; writing placeholder plot.")
        import seaborn as sns

        sns.set_style("white")
        plt.figure(figsize=(8, 4))
        plt.text(0.5, 0.5, "No plottable variants (MAF/F_MISS)", ha="center", va="center")
        plt.axis("off")
        plt.savefig(f"{output_prefix}.nodata.png", dpi=300, bbox_inches="tight")
        plt.close()
        return

    mean_maf = float(df_plot["MAF"].mean())
    median_maf = float(df_plot["MAF"].median())

    th_path = f"{output_prefix}.maf_miss.th.tsv"
    th = {
        "n_variants": len(df_plot),
        "mean_MAF": mean_maf,
        "median_MAF": median_maf,
        "mean_F_MISS": float(df_plot["F_MISS"].mean()),
        "median_F_MISS": float(df_plot["F_MISS"].median()),
    }
    save_df_to_tsv(pd.DataFrame([th]), th_path)

    plot_scatter_with_thresholds(
        data=df_plot.sample(min(8000, len(df_plot)), random_state=1)
        if len(df_plot) > 8000
        else df_plot,
        x_col="MAF",
        y_col="F_MISS",
        title="Variant MAF vs missing rate (PLINK2 slice)",
        filename=f"{output_prefix}.maf_vs_fmiss.scatter.png",
        xlabel="MAF",
        ylabel="F_MISS",
        alpha=0.35,
        s=12,
        thresholds_h=[{"value": 0.1, "color": "red", "linestyle": "--", "label": "F_MISS 0.1"}],
    )

    plot_distribution_with_stats(
        data=df_plot,
        col="MAF",
        title="MAF distribution (PLINK2 slice)",
        filename=f"{output_prefix}.maf.dist.png",
        x_label="MAF",
        y_label="Variant count",
        bins=80,
        mean_val=mean_maf,
        median_val=median_maf,
        xlim=(0, 0.5),
    )

    mean_fm = float(df_plot["F_MISS"].mean())
    med_fm = float(df_plot["F_MISS"].median())
    plot_distribution_with_stats(
        data=df_plot,
        col="F_MISS",
        title="Variant missing-rate distribution (PLINK2 slice)",
        filename=f"{output_prefix}.fmiss.dist.png",
        x_label="F_MISS",
        y_label="Variant count",
        bins=80,
        mean_val=mean_fm,
        median_val=med_fm,
        thresholds=[{"value": 0.1, "label": "0.1", "color": "red", "linestyle": "--"}],
    )

    hist_df = load_df_generic(mac_hist_path)
    if hist_df is not None and len(hist_df) >= 1:
        hist_df.columns = [c.strip() for c in hist_df.columns]
        if "MAF_bin" in hist_df.columns and "n_sites" in hist_df.columns:
            names = hist_df["MAF_bin"].astype(str).tolist()
            values = [float(x) for x in hist_df["n_sites"].tolist()]
            ymax = max(values) * 1.12 if values else 1.0
            plot_bar_chart(
                names,
                values,
                title="Site counts by MAF bin (from PLINK2 .afreq)",
                ylabel="n sites",
                filename=f"{output_prefix}.maf_bins.bar.png",
                ylim=(0.0, ymax),
                color="steelblue",
                figure_size=(10, 5),
            )

    counts_df = load_df_generic(counts_path)
    if counts_df is not None and not counts_df.empty:
        save_df_to_tsv(counts_df, f"{output_prefix}.counts.echo.tsv")

