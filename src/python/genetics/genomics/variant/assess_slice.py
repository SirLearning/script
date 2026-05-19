"""
PLINK2 slice assess: join .acount + .vmiss and emit QC tables and standard infra plots.
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


def _norm_col(name: str) -> str:
    return str(name).replace("#", "").strip()


def _compute_mac_maf_from_acount(df: pd.DataFrame) -> pd.DataFrame:
    """Derive REF allele count, MAC, and MAF from PLINK2 `--freq counts` (.acount)."""

    out = df.copy()
    if "ALT_CTS" not in out.columns or "OBS_CT" not in out.columns:
        return out.assign(REF_CTS=np.nan, MAC=np.nan, MAF=np.nan)

    alt = pd.to_numeric(out["ALT_CTS"], errors="coerce")
    obs = pd.to_numeric(out["OBS_CT"], errors="coerce")
    ref = obs - alt
    mac = np.minimum(alt, ref)
    with np.errstate(divide="ignore", invalid="ignore"):
        maf = np.where(obs > 0, mac.astype(float) / obs.astype(float), np.nan)
    out["REF_CTS"] = ref
    out["MAC"] = mac
    out["MAF"] = maf
    return out


def _mac_category_bar(mac: pd.Series) -> tuple[list[str], list[float]]:
    """Ordered MAC bucket labels and counts for a bar chart."""

    s = mac.dropna()
    try:
        s = s.astype(np.int64)
    except (ValueError, TypeError):
        s = pd.to_numeric(s, errors="coerce").dropna().astype(np.int64)

    def cnt(lo: int, hi: int) -> float:
        return float(((s >= lo) & (s <= hi)).sum())

    names = ["0", "1", "2", "3-9", "10-99", ">=100"]
    values = [
        cnt(0, 0),
        cnt(1, 1),
        cnt(2, 2),
        cnt(3, 9),
        cnt(10, 99),
        cnt(100, 10**9),
    ]
    return names, values


def ana_assess_plink_debug_slice(
    acount_path: str,
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
    acount_path
        PLINK2 ``--freq counts`` .acount file (ALT_CTS / OBS_CT).
    vmiss_path
        PLINK2 --missing .vmiss file.
    mac_hist_path
        Two-column TSV ``MAF_bin\\tn_sites`` from the Nextflow awk helper.
    counts_path
        ``id\\tn_variants\\tn_samples`` TSV from the slice process.
    output_prefix
        Basename prefix for ``*.info.tsv``, ``*.th.tsv``, ``*.png`` outputs.
    """
    df_a = load_df_generic(acount_path)
    df_v = load_df_generic(vmiss_path)
    if df_a is None or df_v is None:
        print("[Error] Missing acount or vmiss input.")
        return

    df_a.columns = [_norm_col(c) for c in df_a.columns]
    df_v.columns = [_norm_col(c) for c in df_v.columns]

    if "ID" not in df_a.columns or "ID" not in df_v.columns:
        print("[Error] Expected ID column in acount and vmiss.")
        return

    df_a = _compute_mac_maf_from_acount(df_a)
    if df_a["MAF"].isna().all():
        print("[Error] Could not derive MAF from acount (check ALT_CTS / OBS_CT).")
        return

    miss_col = "F_MISS" if "F_MISS" in df_v.columns else None
    if miss_col is None:
        print("[Error] F_MISS column not found in vmiss.")
        return

    sub_v = df_v[["ID", miss_col]].copy()
    merged = df_a.merge(sub_v, on="ID", how="inner", suffixes=("", "_vmiss"))
    merged = merged.rename(columns={miss_col: "F_MISS"})

    info_path = f"{output_prefix}.maf_miss.info.tsv"
    save_df_to_tsv(merged, info_path)

    cat_names, cat_vals = _mac_category_bar(merged["MAC"])
    cat_tbl = pd.DataFrame({"mac_category": cat_names, "n_sites": cat_vals})
    save_df_to_tsv(cat_tbl, f"{output_prefix}.mac_category_counts.tsv")

    n_variants = int(len(merged))
    n_mac1 = int(cat_vals[1]) if len(cat_vals) > 1 else 0
    sing = {
        "n_variants": n_variants,
        "n_mac0": int(cat_vals[0]) if cat_vals else 0,
        "n_mac1": int(cat_vals[1]) if len(cat_vals) > 1 else 0,
        "n_mac2": int(cat_vals[2]) if len(cat_vals) > 2 else 0,
        "n_mac_3_9": int(cat_vals[3]) if len(cat_vals) > 3 else 0,
        "n_mac_10_99": int(cat_vals[4]) if len(cat_vals) > 4 else 0,
        "n_mac_ge_100": int(cat_vals[5]) if len(cat_vals) > 5 else 0,
        "frac_mac1": (n_mac1 / n_variants) if n_variants else 0.0,
    }
    save_df_to_tsv(pd.DataFrame([sing]), f"{output_prefix}.singleton_mac.summary.tsv")

    df_plot = merged.dropna(subset=["MAF", "F_MISS", "MAC"]).copy()
    if df_plot.empty:
        print("[Warning] No rows after dropping NA MAF/F_MISS/MAC; writing placeholder plot.")
        import seaborn as sns

        sns.set_style("white")
        plt.figure(figsize=(8, 4))
        plt.text(0.5, 0.5, "No plottable variants (MAF/F_MISS/MAC)", ha="center", va="center")
        plt.axis("off")
        plt.savefig(f"{output_prefix}.nodata.png", dpi=300, bbox_inches="tight")
        plt.close()
        return

    mean_maf = float(df_plot["MAF"].mean())
    median_maf = float(df_plot["MAF"].median())
    mean_mac = float(df_plot["MAC"].mean())
    median_mac = float(df_plot["MAC"].median())

    th_path = f"{output_prefix}.maf_miss.th.tsv"
    th = {
        "n_variants": len(df_plot),
        "mean_MAF": mean_maf,
        "median_MAF": median_maf,
        "mean_MAC": mean_mac,
        "median_MAC": median_mac,
        "mean_F_MISS": float(df_plot["F_MISS"].mean()),
        "median_F_MISS": float(df_plot["F_MISS"].median()),
        "n_mac1": int((df_plot["MAC"] == 1).sum()),
        "frac_mac1": float((df_plot["MAC"] == 1).sum()) / len(df_plot) if len(df_plot) else 0.0,
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

    ymax = max(cat_vals) * 1.12 if cat_vals else 1.0
    plot_bar_chart(
        cat_names,
        cat_vals,
        title="Variant counts by minor allele count (MAC) bucket",
        ylabel="n variants",
        filename=f"{output_prefix}.mac_category.bar.png",
        ylim=(0.0, ymax),
        color="steelblue",
        figure_size=(10, 5),
    )

    df_macpos = df_plot[df_plot["MAC"] >= 1].copy()
    if not df_macpos.empty:
        mac_max = float(df_macpos["MAC"].max())
        cap = min(500.0, max(20.0, float(df_macpos["MAC"].quantile(0.995))))
        plot_distribution_with_stats(
            data=df_macpos,
            col="MAC",
            title="MAC distribution (variants with MAC >= 1)",
            filename=f"{output_prefix}.mac.dist.png",
            x_label="MAC",
            y_label="Variant count",
            bins=min(80, max(10, int(cap))),
            mean_val=float(df_macpos["MAC"].mean()),
            median_val=float(df_macpos["MAC"].median()),
            xlim=(0.5, cap),
        )

    hist_df = load_df_generic(mac_hist_path)
    if hist_df is not None and len(hist_df) >= 1:
        hist_df.columns = [c.strip() for c in hist_df.columns]
        if "MAF_bin" in hist_df.columns and "n_sites" in hist_df.columns:
            names = hist_df["MAF_bin"].astype(str).tolist()
            values = [float(x) for x in hist_df["n_sites"].tolist()]
            hmax = max(values) * 1.12 if values else 1.0
            plot_bar_chart(
                names,
                values,
                title="Site counts by MAF bin (from PLINK2 .acount)",
                ylabel="n sites",
                filename=f"{output_prefix}.maf_bins.bar.png",
                ylim=(0.0, hmax),
                color="steelblue",
                figure_size=(10, 5),
            )

    counts_df = load_df_generic(counts_path)
    if counts_df is not None and not counts_df.empty:
        save_df_to_tsv(counts_df, f"{output_prefix}.counts.echo.tsv")
