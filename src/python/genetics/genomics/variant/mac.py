from .variant_utils import load_df_from_plink_gcount, load_df_from_plink_variant

import logging
import sys

import numpy as np
import pandas as pd

from infra.utils.errors import DataLoadError, configure_logging, fail
from infra.utils.graph import (
    LEGEND_FONT_SIZE,
    TITLE_FONT_SIZE,
    X_LABEL_FONT_SIZE,
    Y_LABEL_FONT_SIZE,
    TICK_FONT_SIZE,
    plot_heatmap_custom,
    plot_regression_comparison,
    plot_scatter_with_thresholds,
)
from infra.utils.io import save_df_to_tsv, save_thresholds

logger = logging.getLogger(__name__)


def _compute_mac_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Derive allele counts, MAC, and observed heterozygosity among minor-allele genotypes.

    ``Het_Fraction`` = ``HET_REF_ALT_CTS / MAC`` when MAC > 0: among minor-allele copies,
    the fraction carried in heterozygous genotypes (MAC=1 implies 1.0 when segregating).
    """
    required = [
        "HET_REF_ALT_CTS",
        "TWO_ALT_GENO_CTS",
        "HOM_REF_CT",
        "HAP_REF_CT",
        "HAP_ALT_CTS",
    ]
    missing = [col for col in required if col not in df.columns]
    if missing:
        fail(f"Missing required columns for MAC computation: {missing}")

    out = df.copy()
    out["Alt_Count"] = out["HET_REF_ALT_CTS"] + (out["TWO_ALT_GENO_CTS"] * 2) + out["HAP_ALT_CTS"]
    out["Ref_Count"] = out["HET_REF_ALT_CTS"] + (out["HOM_REF_CT"] * 2) + out["HAP_REF_CT"]
    out["Total_Alleles"] = out["Alt_Count"] + out["Ref_Count"]
    out["MAC"] = out[["Alt_Count", "Ref_Count"]].min(axis=1)
    out["MAF"] = np.where(out["Total_Alleles"] > 0, out["MAC"] / out["Total_Alleles"], np.nan)
    out["Total_Samples"] = (
        out["HOM_REF_CT"]
        + out["HET_REF_ALT_CTS"]
        + out["TWO_ALT_GENO_CTS"]
        + out["HAP_REF_CT"]
        + out["HAP_ALT_CTS"]
    )
    out["Het_Fraction"] = np.nan
    mask_mac = out["MAC"] > 0
    out.loc[mask_mac, "Het_Fraction"] = out.loc[mask_mac, "HET_REF_ALT_CTS"] / out.loc[mask_mac, "MAC"]
    return out


def _mac_site_count_table(df: pd.DataFrame) -> pd.DataFrame:
    """All MAC values present in the input, with variant site counts."""
    counts = df["MAC"].value_counts(sort=True).sort_index()
    return pd.DataFrame({"MAC": counts.index.astype(int), "n_sites": counts.values.astype(int)})


def _summarize_mac_zero_sites(df: pd.DataFrame) -> pd.DataFrame:
    """Audit MAC=0 sites (unexpected): classify by allele-count pattern."""
    mac0 = df[df["MAC"] == 0].copy()
    if mac0.empty:
        return pd.DataFrame(
            {
                "category": ["none"],
                "n_sites": [0],
                "description": ["No MAC=0 sites"],
            }
        )

    def _cat(row):
        alt, ref = int(row["Alt_Count"]), int(row["Ref_Count"])
        if alt == 0 and ref == 0:
            return "both_allele_counts_zero"
        if alt == 0 and ref > 0:
            return "alt_absent_hom_ref_only"
        if ref == 0 and alt > 0:
            if int(row["TWO_ALT_GENO_CTS"]) > 0 or int(row["HAP_ALT_CTS"]) > 0:
                return "ref_absent_hom_or_hap_alt_only"
            if int(row["HET_REF_ALT_CTS"]) > 0:
                return "ref_absent_het_only_impossible_mac0"
            return "ref_absent_other"
        return "other_mac0"

    mac0["category"] = mac0.apply(_cat, axis=1)
    summary = (
        mac0.groupby("category", as_index=False)
        .size()
        .rename(columns={"size": "n_sites"})
        .sort_values("n_sites", ascending=False)
    )
    desc = {
        "both_allele_counts_zero": "No observed alleles (all missing or zero genotypes)",
        "alt_absent_hom_ref_only": "Alternate allele absent; hom-ref / hap-ref only",
        "ref_absent_hom_or_hap_alt_only": "Reference allele absent; hom-alt / hap-alt only",
        "ref_absent_het_only_impossible_mac0": "Het genotypes present but MAC=0 (check counts)",
        "ref_absent_other": "Reference absent, other genotype mix",
        "other_mac0": "Unclassified MAC=0 pattern",
    }
    summary["description"] = summary["category"].map(desc)
    summary = pd.concat(
        [
            pd.DataFrame(
                {
                    "category": ["total_mac0"],
                    "n_sites": [len(mac0)],
                    "description": ["All sites with MAC=0"],
                }
            ),
            summary,
        ],
        ignore_index=True,
    )
    return summary


def _plot_mac_site_dist_0_100(
    mac_counts: pd.DataFrame,
    output_prefix: str,
    mac_max: int = 100,
    tick_step: int = 5,
    log_y: bool = False,
) -> None:
    """Bar chart of site counts for MAC 0..mac_max with sparse x tick labels."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    lookup = dict(zip(mac_counts["MAC"].astype(int), mac_counts["n_sites"].astype(int)))
    xs = list(range(0, mac_max + 1))
    values = [float(lookup.get(i, 0)) for i in xs]

    if sum(values) <= 0:
        print(f"[Warning] No variant sites with MAC 0-{mac_max}; writing zero-filled bar chart.")

    suffix = ".dist.0_100.log.png" if log_y else ".dist.0_100.png"
    out_path = f"{output_prefix}{suffix}"

    if log_y and sum(values) <= 0:
        _save_mac_plot_placeholder(
            out_path,
            f"No variant sites with MAC 0-{mac_max} (log scale unavailable)",
        )
        return

    plot_values = values
    if log_y:
        plot_values = [v if v > 0 else np.nan for v in values]
        if not any(np.isfinite(v) and v > 0 for v in plot_values):
            _save_mac_plot_placeholder(
                out_path,
                f"No variant sites with MAC 0-{mac_max} (log scale unavailable)",
            )
            return

    title_scale = " (log y)" if log_y else ""
    y_label = "Number of variant sites (log scale)" if log_y else "Number of variant sites"

    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.bar(xs, plot_values, color="steelblue", alpha=0.85, width=0.9)
    ax.set_xlim(-0.5, mac_max + 0.5)
    ax.set_xticks(list(range(0, mac_max + 1, tick_step)))
    ax.set_xlabel("Minor allele count (MAC)", fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(
        f"Variant site counts by minor allele count (MAC 0-{mac_max}){title_scale}",
        fontsize=TITLE_FONT_SIZE,
    )
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if log_y:
        ax.set_yscale("log")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved bar chart: {out_path}")


def _save_mac_plot_placeholder(filename: str, message: str) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_style("white")
    plt.figure(figsize=(8, 4))
    plt.axis("off")
    plt.text(0.5, 0.5, message, ha="center", va="center", fontsize=12)
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Placeholder saved to: {filename}")


def _plot_mac_hetfrac_heatmap(
    df: pd.DataFrame,
    output_prefix: str,
    mac_max: int = 100,
    het_bin_width: float = 0.05,
    x_tick_step: int = 5,
) -> None:
    """
    Heatmap: x = MAC (0..mac_max), y = Het_Fraction bins among minor-allele genotypes;
    cell value = fraction of sites at that MAC with Het_Fraction in the bin.
    """
    out_file = f"{output_prefix}.heatmap_mac_hetfrac.0_100.png"
    df_sub = df[(df["MAC"] >= 1) & (df["MAC"] <= mac_max)].dropna(subset=["Het_Fraction"])
    if len(df_sub) < 10:
        print("[Warning] Too few rows for MAC vs Het_Fraction heatmap; writing placeholder plot.")
        _save_mac_plot_placeholder(out_file, f"No variants with MAC 1-{mac_max} for heatmap")
        return

    y_bins = np.arange(0.0, 1.0 + het_bin_width, het_bin_width)
    n_y_bins = len(y_bins) - 1
    n_x = mac_max + 1
    matrix = np.zeros((n_y_bins, n_x))

    for mac_val in range(1, mac_max + 1):
        subset = df_sub[df_sub["MAC"] == mac_val]
        if subset.empty:
            continue
        counts, _ = np.histogram(subset["Het_Fraction"], bins=y_bins)
        total = counts.sum()
        if total > 0:
            matrix[:, mac_val] = counts / total

    matrix = np.flipud(matrix)
    x_labels = [str(i) for i in range(0, mac_max + 1)]
    y_labels = [f"{(1.0 - i * het_bin_width):.2f}" for i in range(n_y_bins)]

    plot_heatmap_custom(
        data_matrix=matrix,
        x_labels=x_labels,
        y_labels=y_labels,
        title=f"Observed heterozygosity (minor-allele het fraction) per MAC (0-{mac_max})",
        filename=out_file,
        xlabel="Minor allele count (MAC)",
        ylabel="Het fraction among minor-allele copies",
        cbar_label="Fraction of sites",
        x_tick_step=x_tick_step,
    )


def ana_mac_stats(input_file: str, output_prefix: str = "variant_mac") -> None:
    """
    MAC analytics from PLINK2 ``--geno-counts`` (.gcount).

    Writes:
      - ``{output_prefix}.info.tsv`` — all MAC values with ``n_sites``
      - ``{output_prefix}.mac0.info.tsv`` — MAC=0 site audit
      - ``{output_prefix}.dist.0_100.png`` — site-count bar chart for MAC 0-100 (linear y)
      - ``{output_prefix}.dist.0_100.log.png`` — same distribution with log-scaled y axis
      - ``{output_prefix}.heatmap_mac_hetfrac.0_100.png`` — MAC x het-fraction heatmap
      - ``{output_prefix}.th.tsv`` — summary thresholds
    """
    logger.info("Processing MAC stats: %s", input_file)

    df_raw = load_df_from_plink_gcount(input_file)
    df = _compute_mac_table(df_raw)

    mac_counts = _mac_site_count_table(df)
    save_df_to_tsv(mac_counts, f"{output_prefix}.info.tsv")

    mac0_summary = _summarize_mac_zero_sites(df)
    save_df_to_tsv(mac0_summary, f"{output_prefix}.mac0.info.tsv")
    n_mac0 = int((df["MAC"] == 0).sum())
    if n_mac0:
        print(f"[Warning] Found {n_mac0} MAC=0 sites (see {output_prefix}.mac0.info.tsv)")
    print(mac0_summary.to_string(index=False))

    n_variants = int(len(df))
    n_mac1 = int((df["MAC"] == 1).sum())
    stats_dict = {
        "Total_Variants": n_variants,
        "Distinct_MAC_Values": int(mac_counts.shape[0]),
        "MAC_0_Sites": n_mac0,
        "Frac_MAC_0": (n_mac0 / n_variants) if n_variants else 0.0,
        "MAC_1_Sites": n_mac1,
        "Frac_MAC_1": (n_mac1 / n_variants) if n_variants else 0.0,
        "Max_MAC": int(df["MAC"].max()) if n_variants else 0,
    }
    save_thresholds(stats_dict, f"{output_prefix}.th.tsv")
    print(stats_dict)

    _plot_mac_site_dist_0_100(mac_counts, output_prefix, mac_max=100, tick_step=5, log_y=False)
    _plot_mac_site_dist_0_100(mac_counts, output_prefix, mac_max=100, tick_step=5, log_y=True)
    _plot_mac_hetfrac_heatmap(df, output_prefix, mac_max=100, x_tick_step=5)


def _mac_bin_missing_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate mean/median F_MISS by MAC bucket for compact info.tsv output."""

    def _bucket(mac: float) -> str:
        if pd.isna(mac):
            return "NA"
        m = int(mac)
        if m == 0:
            return "0"
        if m == 1:
            return "1"
        if m == 2:
            return "2"
        if m <= 9:
            return "3-9"
        if m <= 99:
            return "10-99"
        return ">=100"

    work = df[["MAC", "F_MISS"]].dropna().copy()
    if work.empty:
        return pd.DataFrame(columns=["mac_bucket", "n_sites", "mean_F_MISS", "median_F_MISS"])

    work["mac_bucket"] = work["MAC"].map(_bucket)
    order = ["0", "1", "2", "3-9", "10-99", ">=100"]
    summary = (
        work.groupby("mac_bucket", as_index=False)
        .agg(n_sites=("MAC", "size"), mean_F_MISS=("F_MISS", "mean"), median_F_MISS=("F_MISS", "median"))
        .assign(
            mac_bucket=lambda d: pd.Categorical(d["mac_bucket"], categories=order, ordered=True),
        )
        .sort_values("mac_bucket")
        .reset_index(drop=True)
    )
    summary["mac_bucket"] = summary["mac_bucket"].astype(str)
    return summary


def _subsample_for_plot(df: pd.DataFrame, max_points: int, random_seed: int) -> pd.DataFrame:
    if len(df) <= max_points:
        return df
    return df.sample(n=max_points, random_state=random_seed)


_MAC_REG_RANGES = (100, 500, 1000)
_MAC_BIN50_WIDTH = 50
_MAC_BIN50_SAMPLE_PER_BIN = 100
def _mac_bin50_label(mac: int) -> str:
    """MAC bins of width 50: 0-50, 51-100, 101-150, ..."""
    m = int(mac)
    if m <= _MAC_BIN50_WIDTH:
        return f"0-{_MAC_BIN50_WIDTH}"
    lo = _MAC_BIN50_WIDTH + 1 + ((m - (_MAC_BIN50_WIDTH + 1)) // _MAC_BIN50_WIDTH) * _MAC_BIN50_WIDTH
    hi = lo + _MAC_BIN50_WIDTH - 1
    return f"{lo}-{hi}"


def _mac_bin50_label_from_one(mac: int) -> str | None:
    """50-wide MAC bins excluding MAC=0: 1-50, 51-100, 101-150, ..."""
    m = int(mac)
    if m <= 0:
        return None
    lo = 1 + ((m - 1) // _MAC_BIN50_WIDTH) * _MAC_BIN50_WIDTH
    hi = lo + _MAC_BIN50_WIDTH - 1
    return f"{lo}-{hi}"


def _mac_bin50_center(bin_label: str) -> float:
    lo_s, hi_s = bin_label.split("-", 1)
    return (int(lo_s) + int(hi_s)) / 2.0


def _mac_binned_mean_curve(
    df: pd.DataFrame,
    value_col: str,
    mean_col: str,
    mac_min: float = 1.0,
    mac_max: float | None = None,
) -> pd.DataFrame:
    """Mean ``value_col`` per 50-MAC bin (MAC>=1) for overlay curves."""
    work = df.dropna(subset=["MAC", value_col]).copy()
    work = work[work["MAC"] >= mac_min]
    if mac_max is not None:
        work = work[work["MAC"] <= mac_max]
    if work.empty:
        return pd.DataFrame(columns=["MAC", mean_col])

    work["_mac_bin"] = work["MAC"].astype(int).map(_mac_bin50_label_from_one)
    work = work.dropna(subset=["_mac_bin"])
    if work.empty:
        return pd.DataFrame(columns=["MAC", mean_col])

    curve = work.groupby("_mac_bin", as_index=False)[value_col].mean().rename(columns={value_col: mean_col})
    curve["MAC"] = curve["_mac_bin"].map(_mac_bin50_center)
    return curve.sort_values("MAC").reset_index(drop=True)[["MAC", mean_col]]


def _mac_miss_mean_by_mac_bin50(df: pd.DataFrame, mac_min: float = 1.0, mac_max: float | None = None) -> pd.DataFrame:
    """Mean F_MISS per 50-MAC bin (MAC>=1) for overlay curve on regression plots."""
    return _mac_binned_mean_curve(df, "F_MISS", "mean_F_MISS", mac_min=mac_min, mac_max=mac_max)


def _mac_miss_mean_r_by_mac_bin50(df: pd.DataFrame, mac_min: float = 1.0, mac_max: float | None = None) -> pd.DataFrame:
    """Mean R_miss_boundary per 50-MAC bin (MAC>=1) for R_mac plot overlay."""
    return _mac_binned_mean_curve(df, "R_miss_boundary", "mean_R", mac_min=mac_min, mac_max=mac_max)


def _compute_mac_miss_derived_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add per-variant metrics for MAC vs missing / AN plots.

    - ``N``: genotyped cohort size (OBS_CT)
    - ``AN``: allele number among called samples, 2 * N * (1 - F_MISS)
    - ``R_miss_boundary``: F_MISS / (1 - MAC/N), fraction of theoretical max missing
    """
    out = df.copy()
    n = pd.to_numeric(out["N"], errors="coerce").astype(float)
    mac = pd.to_numeric(out["MAC"], errors="coerce").astype(float)
    f_miss = pd.to_numeric(out["F_MISS"], errors="coerce").astype(float)

    boundary_headroom = 1.0 - (mac / n)
    with np.errstate(divide="ignore", invalid="ignore"):
        out["R_miss_boundary"] = np.where(boundary_headroom > 0, f_miss / boundary_headroom, np.nan)
        out["AN"] = 2.0 * n * (1.0 - f_miss)
    return out


def _load_mac_miss_extended_df(gcount_path: str, vmiss_path: str) -> pd.DataFrame:
    """Join gcount + vmiss with MAC, F_MISS, N, AN, and R_miss_boundary."""
    df_raw = load_df_from_plink_gcount(gcount_path)
    df_mac = _compute_mac_table(df_raw)
    if "ID" not in df_mac.columns:
        fail("Expected ID column in gcount input for MAC vs missing regression")

    df_miss = load_df_from_plink_variant(vmiss_path)
    if df_miss is None or "F_MISS" not in df_miss.columns:
        fail("Expected F_MISS column in vmiss input for MAC vs missing regression")
    if "ID" not in df_miss.columns:
        fail("Expected ID column in vmiss input for MAC vs missing regression")

    miss_cols = ["ID", "F_MISS"]
    if "OBS_CT" in df_miss.columns:
        miss_cols.append("OBS_CT")
    merged = df_mac.merge(df_miss[miss_cols], on="ID", how="inner")
    if merged.empty:
        fail("No overlapping variant IDs between gcount and vmiss")

    if "OBS_CT" in merged.columns:
        merged["N"] = pd.to_numeric(merged["OBS_CT"], errors="coerce")
    elif "Total_Samples" in merged.columns:
        merged["N"] = pd.to_numeric(merged["Total_Samples"], errors="coerce")
    else:
        fail("Expected OBS_CT in vmiss or Total_Samples in gcount for MAC vs missing metrics")

    keep = ["MAC", "F_MISS", "N"]
    plot_df = merged[keep].replace([np.inf, -np.inf], np.nan).dropna(subset=["MAC", "F_MISS", "N"])
    if plot_df.empty:
        fail("No plottable variants after dropping NA MAC/F_MISS/N")
    return _compute_mac_miss_derived_metrics(plot_df)


def _stratified_sample_mac_bin50(
    plot_df: pd.DataFrame,
    mac_max: int | None,
    per_bin: int = _MAC_BIN50_SAMPLE_PER_BIN,
    random_seed: int = 1,
) -> pd.DataFrame:
    """Up to ``per_bin`` random variants per 50-wide MAC bin (all if bin is smaller)."""
    work = plot_df.copy()
    if mac_max is not None:
        work = work[(work["MAC"] >= 0) & (work["MAC"] <= mac_max)]
    if work.empty:
        return work

    work = work.assign(_mac_bin50=work["MAC"].astype(int).map(_mac_bin50_label))
    parts: list[pd.DataFrame] = []
    for _, grp in work.groupby("_mac_bin50", sort=False):
        if len(grp) > per_bin:
            parts.append(grp.sample(n=per_bin, random_state=random_seed))
        else:
            parts.append(grp)
    return pd.concat(parts, ignore_index=True).drop(columns="_mac_bin50")


def _plot_mac_miss_regression(
    plot_df: pd.DataFrame,
    output_path: str,
    title: str,
    x_lim: tuple[float, float] | None = None,
    y_lim: tuple[float, float] = (0.0, 1.0),
    mac_max: int | None = None,
    stats_df: pd.DataFrame | None = None,
) -> None:
    """MAC vs F_MISS scatter with OLS/Huber regression, binned mean curve, and boundary."""
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.linear_model import HuberRegressor, LinearRegression

    overlay_df = stats_df if stats_df is not None else plot_df
    work = plot_df.dropna(subset=["MAC", "F_MISS"])
    local_valid = work[["MAC", "F_MISS"]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(local_valid) < 10:
        print(f"Not enough valid data for {output_path}")
        return

    X = local_valid["MAC"].values.reshape(-1, 1)
    y = local_valid["F_MISS"].values

    ols = LinearRegression()
    ols.fit(X, y)
    ols_score = ols.score(X, y)
    ols_eq = f"y = {ols.coef_[0]:.4f}x + {ols.intercept_:.4f}"

    huber = HuberRegressor(epsilon=1.35)
    huber.fit(X, y)
    huber_score = huber.score(X, y)
    huber_eq = f"y = {huber.coef_[0]:.4f}x + {huber.intercept_:.4f}"

    sns.set_style("white")
    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        data=local_valid,
        x="MAC",
        y="F_MISS",
        alpha=0.3,
        s=15,
        color="#1f77b4",
        edgecolor="none",
        label="Variants",
    )

    x_min, x_max = local_valid["MAC"].min(), local_valid["MAC"].max()
    x_span = x_max - x_min
    if x_span == 0:
        x_span = 1
    if x_lim is not None:
        x_line_min, x_line_max = x_lim
    else:
        x_line_min = x_min - 0.05 * x_span
        x_line_max = x_max + 0.05 * x_span
    x_range = np.linspace(x_line_min, x_line_max, 100).reshape(-1, 1)
    plt.plot(
        x_range,
        ols.predict(x_range),
        color="blue",
        linewidth=2,
        linestyle="--",
        label=f"OLS: {ols_eq}, $R^2$={ols_score:.3f}",
    )
    plt.plot(
        x_range,
        huber.predict(x_range),
        color="green",
        linewidth=2,
        label=f"Huber: {huber_eq}, $R^2$={huber_score:.3f}",
    )

    curve_mac_max = x_lim[1] if x_lim is not None else (float(mac_max) if mac_max is not None else None)
    mean_curve = _mac_miss_mean_by_mac_bin50(overlay_df, mac_min=1.0, mac_max=curve_mac_max)
    if len(mean_curve) >= 2:
        plt.plot(
            mean_curve["MAC"],
            mean_curve["mean_F_MISS"],
            color="darkorange",
            linewidth=2,
            linestyle="-",
            label=f"Mean F_MISS ({_MAC_BIN50_WIDTH}-MAC bins, MAC≥1)",
        )

    overlay_work = overlay_df.dropna(subset=["MAC", "F_MISS", "N"]) if "N" in overlay_df.columns else overlay_df
    if "N" in overlay_work.columns and overlay_work["N"].notna().any():
        n_rep = float(overlay_work["N"].median())
        x_bound = np.linspace(0, max(x_line_max, 1.0), 200)
        y_bound = np.clip(1.0 - x_bound / n_rep, 0.0, 1.0)
        plt.plot(
            x_bound,
            y_bound,
            color="red",
            linestyle="--",
            linewidth=2,
            label=f"Boundary: y = 1 - MAC/N (N={n_rep:.0f})",
        )

    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.xlabel("Minor allele count (MAC)", fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel("Variant missing rate (F_MISS)", fontsize=Y_LABEL_FONT_SIZE)
    plt.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if x_lim is not None:
        plt.xlim(x_lim)
    else:
        plt.xlim(x_line_min, x_line_max)
    if y_lim is not None:
        plt.ylim(y_lim)
    plt.legend(loc="upper left", bbox_to_anchor=(0.0, -0.15), fontsize=LEGEND_FONT_SIZE, frameon=False)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved plot to {output_path}")
    plt.close()


def _plot_mac_an_regression(
    plot_df: pd.DataFrame,
    output_path: str,
    title: str,
    x_lim: tuple[float, float] | None = None,
    y_lim: tuple[float, float] | None = None,
) -> None:
    """MAC vs AN scatter with OLS/Huber regression (AN = 2N(1-F_MISS))."""
    reg_df = plot_df.dropna(subset=["MAC", "AN"])
    if len(reg_df) < 10:
        print(f"Not enough valid data for {output_path}")
        return
    plot_regression_comparison(
        reg_df,
        x_col="MAC",
        y_col="AN",
        x_label="Minor allele count (MAC)",
        y_label="Allele number AN = 2N(1 - F_MISS)",
        filename=output_path,
        title=title,
        x_lim=x_lim,
        y_lim=y_lim,
    )


def _plot_mac_missing_reg_by_mac_range(
    plot_df: pd.DataFrame,
    output_prefix: str,
    mac_max: int,
    regression_max_points: int,
    random_seed: int,
) -> dict[str, float | int]:
    """Scatter + OLS/Huber regression for variants with MAC in [0, mac_max]."""
    out_path = f"{output_prefix}.reg.mac0_{mac_max}.png"
    subset = plot_df[(plot_df["MAC"] >= 0) & (plot_df["MAC"] <= mac_max)].copy()
    n_sites = int(len(subset))

    if n_sites < 10:
        _save_mac_plot_placeholder(
            out_path,
            f"Too few variants ({n_sites}) with MAC 0-{mac_max} for regression",
        )
        return {
            f"N_MAC_0_{mac_max}": n_sites,
            f"Pearson_r_MAC_0_{mac_max}": np.nan,
            f"Regression_Points_MAC_0_{mac_max}": 0,
        }

    reg_df = _subsample_for_plot(subset, regression_max_points, random_seed)
    _plot_mac_miss_regression(
        reg_df,
        out_path,
        title=f"Variant MAC vs missing rate (MAC 0-{mac_max})",
        x_lim=(0.0, float(mac_max)),
        y_lim=(0.0, 1.0),
        mac_max=mac_max,
        stats_df=subset,
    )
    pearson_r = float(subset["MAC"].corr(subset["F_MISS"])) if n_sites > 1 else np.nan
    return {
        f"N_MAC_0_{mac_max}": n_sites,
        f"Pearson_r_MAC_0_{mac_max}": pearson_r,
        f"Regression_Points_MAC_0_{mac_max}": int(len(reg_df)),
    }


def _plot_mac_an_reg_by_mac_range(
    plot_df: pd.DataFrame,
    output_prefix: str,
    mac_max: int,
    regression_max_points: int,
    random_seed: int,
) -> dict[str, float | int]:
    """MAC vs AN regression for variants with MAC in [0, mac_max]."""
    out_path = f"{output_prefix}.reg.mac_an.mac0_{mac_max}.png"
    subset = plot_df[(plot_df["MAC"] >= 0) & (plot_df["MAC"] <= mac_max)].copy()
    n_sites = int(len(subset))

    if n_sites < 10:
        _save_mac_plot_placeholder(
            out_path,
            f"Too few variants ({n_sites}) with MAC 0-{mac_max} for MAC vs AN regression",
        )
        return {
            f"N_MAC_AN_0_{mac_max}": n_sites,
            f"Pearson_r_MAC_AN_0_{mac_max}": np.nan,
            f"Regression_Points_MAC_AN_0_{mac_max}": 0,
        }

    reg_df = _subsample_for_plot(subset, regression_max_points, random_seed)
    _plot_mac_an_regression(
        reg_df,
        out_path,
        title=f"Variant MAC vs allele number AN (MAC 0-{mac_max})",
        x_lim=(0.0, float(mac_max)),
    )
    pearson_r = float(subset["MAC"].corr(subset["AN"])) if n_sites > 1 else np.nan
    return {
        f"N_MAC_AN_0_{mac_max}": n_sites,
        f"Pearson_r_MAC_AN_0_{mac_max}": pearson_r,
        f"Regression_Points_MAC_AN_0_{mac_max}": int(len(reg_df)),
    }


def ana_mac_missing_reg(
    gcount_path: str,
    vmiss_path: str,
    output_prefix: str,
    regression_max_points: int = 50000,
    random_seed: int = 1,
) -> None:
    """
    Per-site MAC vs variant missing rate (F_MISS) and MAC vs AN scatter regression.

    Joins ``--geno-counts`` (.gcount) with ``--missing`` variant report (.vmiss), writes
    MAC-bucket missing summaries, regression thresholds, and OLS/Huber plots with the
    boundary line y = 1 - MAC/N (full MAC range plus MAC 0-100, 0-500, and 0-1000 subsets).
    Also emits MAC vs AN (``reg.mac_an*``) and R vs MAC (``R_mac*``) plots on all variants.
    AN = 2 * N * (1 - F_MISS).
    """
    logger.info("Processing MAC vs missing regression: gcount=%s vmiss=%s", gcount_path, vmiss_path)

    plot_df = _load_mac_miss_extended_df(gcount_path, vmiss_path)

    save_df_to_tsv(_mac_bin_missing_summary(plot_df), f"{output_prefix}.info.tsv")

    reg_df = _subsample_for_plot(plot_df, regression_max_points, random_seed)

    n_variants = int(len(plot_df))
    n_reg = int(len(reg_df))
    pearson_r = float(plot_df["MAC"].corr(plot_df["F_MISS"])) if n_variants > 1 else np.nan
    pearson_r_an = float(plot_df["MAC"].corr(plot_df["AN"])) if n_variants > 1 else np.nan

    stats_dict = {
        "Total_Variants": n_variants,
        "Regression_Points": n_reg,
        "Pearson_r_MAC_F_MISS": pearson_r,
        "Pearson_r_MAC_AN": pearson_r_an,
        "Mean_MAC": float(plot_df["MAC"].mean()),
        "Median_MAC": float(plot_df["MAC"].median()),
        "Mean_F_MISS": float(plot_df["F_MISS"].mean()),
        "Median_F_MISS": float(plot_df["F_MISS"].median()),
        "Mean_AN": float(plot_df["AN"].mean()),
        "Median_AN": float(plot_df["AN"].median()),
        "Regression_Subsampled": int(n_reg < n_variants),
    }
    if len(reg_df) < 10:
        _save_mac_plot_placeholder(
            f"{output_prefix}.reg.png",
            f"Too few variants ({len(reg_df)}) for MAC vs F_MISS regression",
        )
        _save_mac_plot_placeholder(
            f"{output_prefix}.reg.mac_an.png",
            f"Too few variants ({len(reg_df)}) for MAC vs AN regression",
        )
    else:
        _plot_mac_miss_regression(
            reg_df,
            f"{output_prefix}.reg.png",
            title="Variant MAC vs missing rate",
            y_lim=(0.0, 1.0),
            stats_df=plot_df,
        )
        _plot_mac_an_regression(
            reg_df,
            f"{output_prefix}.reg.mac_an.png",
            title="Variant MAC vs allele number AN",
        )

    for mac_max in _MAC_REG_RANGES:
        stats_dict.update(
            _plot_mac_missing_reg_by_mac_range(
                plot_df,
                output_prefix,
                mac_max,
                regression_max_points,
                random_seed,
            )
        )
        stats_dict.update(
            _plot_mac_an_reg_by_mac_range(
                plot_df,
                output_prefix,
                mac_max,
                regression_max_points,
                random_seed,
            )
        )

    stats_dict.update(
        _plot_R_mac_by_mac_max(plot_df, output_prefix, mac_max=None, regression_max_points=regression_max_points, random_seed=random_seed)
    )
    for mac_max in _MAC_REG_RANGES:
        stats_dict.update(
            _plot_R_mac_by_mac_max(
                plot_df,
                output_prefix,
                mac_max=mac_max,
                regression_max_points=regression_max_points,
                random_seed=random_seed,
            )
        )

    save_thresholds(stats_dict, f"{output_prefix}.th.tsv")
    print(stats_dict)


def _subset_by_mac_max(plot_df: pd.DataFrame, mac_max: int | None) -> pd.DataFrame:
    if mac_max is None:
        return plot_df.copy()
    return plot_df[(plot_df["MAC"] >= 0) & (plot_df["MAC"] <= mac_max)].copy()


def _plot_R_mac_scatter(
    plot_df: pd.DataFrame,
    output_path: str,
    title: str,
    x_lim: tuple[float, float] | None = None,
    y_lim: tuple[float, float] = (0.0, 1.2),
    stats_df: pd.DataFrame | None = None,
) -> None:
    """R = F_MISS/(1-MAC/N) vs MAC scatter with R=1 threshold and binned mean-R curve."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    overlay_df = stats_df if stats_df is not None else plot_df
    scatter = plot_df.dropna(subset=["MAC", "R_miss_boundary"])
    if len(scatter) < 10:
        print(f"Not enough valid data for {output_path}")
        return

    sns.set_style("white")
    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        data=scatter,
        x="MAC",
        y="R_miss_boundary",
        alpha=0.35,
        s=12,
        color="royalblue",
        edgecolor="w",
        label="Variants",
    )
    plt.axhline(y=1.0, color="red", linestyle="--", linewidth=2, label="R=1 (at boundary)")

    curve_mac_max = x_lim[1] if x_lim is not None else None
    mean_curve = _mac_miss_mean_r_by_mac_bin50(overlay_df, mac_min=1.0, mac_max=curve_mac_max)
    if len(mean_curve) >= 2:
        plt.plot(
            mean_curve["MAC"],
            mean_curve["mean_R"],
            color="darkorange",
            linewidth=2,
            linestyle="-",
            label=f"Mean R ({_MAC_BIN50_WIDTH}-MAC bins, MAC≥1)",
        )

    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.xlabel("Minor allele count (MAC)", fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel("R = F_MISS / (1 - MAC/N)", fontsize=Y_LABEL_FONT_SIZE)
    plt.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if x_lim is not None:
        plt.xlim(x_lim)
    if y_lim is not None:
        plt.ylim(y_lim)
    plt.legend(loc="upper left", bbox_to_anchor=(0.0, -0.15), fontsize=LEGEND_FONT_SIZE, frameon=False)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved scatter plot to {output_path}")
    plt.close()


def _plot_R_mac_by_mac_max(
    plot_df: pd.DataFrame,
    output_prefix: str,
    mac_max: int | None,
    regression_max_points: int,
    random_seed: int,
) -> dict[str, float | int]:
    """R = F_MISS / (1 - MAC/N) vs MAC scatter using all variants in range."""
    if mac_max is None:
        out_path = f"{output_prefix}.R_mac.png"
        key = "ALL"
        title_suffix = "all MAC"
    else:
        out_path = f"{output_prefix}.R_mac.mac0_{mac_max}.png"
        key = f"MAC0_{mac_max}"
        title_suffix = f"MAC 0-{mac_max}"

    in_range = _subset_by_mac_max(plot_df, mac_max)
    n_in_range = int(len(in_range))
    plot_pts = _subsample_for_plot(in_range, regression_max_points, random_seed)
    n_plotted = int(len(plot_pts))
    stats: dict[str, float | int] = {
        f"N_in_range_R_mac_{key}": n_in_range,
        f"N_plotted_R_mac_{key}": n_plotted,
        f"Plot_subsampled_R_mac_{key}": int(n_plotted < n_in_range),
    }

    if n_plotted < 10:
        _save_mac_plot_placeholder(
            out_path,
            f"Too few variants ({n_plotted}) for R vs MAC ({title_suffix})",
        )
        return stats

    _plot_R_mac_scatter(
        plot_pts,
        out_path,
        title=f"R = F_MISS/(1-MAC/N) vs MAC (all variants, {title_suffix})",
        x_lim=None if mac_max is None else (0.0, float(mac_max)),
        y_lim=(0.0, 1.2),
        stats_df=in_range,
    )

    near_bound = in_range["R_miss_boundary"].dropna()
    if len(near_bound):
        stats[f"Frac_R_ge_0p9_R_mac_{key}"] = float((near_bound >= 0.9).mean())
        stats[f"Median_R_mac_{key}"] = float(near_bound.median())
    return stats


def _plot_mac_missing_reg_bin50_sample(
    plot_df: pd.DataFrame,
    output_prefix: str,
    random_seed: int,
) -> dict[str, float | int]:
    """MAC vs F_MISS regression on stratified 50-wide MAC bin samples (full MAC range)."""
    out_path = f"{output_prefix}.reg.bin50s.png"
    n_in_range = int(len(plot_df))
    sample_df = _stratified_sample_mac_bin50(plot_df, mac_max=None, random_seed=random_seed)
    n_sampled = int(len(sample_df))
    n_bins = int(plot_df["MAC"].astype(int).map(_mac_bin50_label).nunique()) if n_in_range else 0

    if n_sampled < 10:
        _save_mac_plot_placeholder(
            out_path,
            f"Too few stratified sample points ({n_sampled}) for MAC vs F_MISS (all MAC)",
        )
        pearson_r = np.nan
    else:
        _plot_mac_miss_regression(
            sample_df,
            out_path,
            title=f"MAC vs missing (50-MAC bin sample, n<={_MAC_BIN50_SAMPLE_PER_BIN}/bin, all MAC)",
            y_lim=(0.0, 1.0),
            stats_df=plot_df,
        )
        pearson_r = float(sample_df["MAC"].corr(sample_df["F_MISS"])) if n_sampled > 1 else np.nan

    return {
        "N_in_range_BIN50S": n_in_range,
        "N_mac_bin50_bins_BIN50S": n_bins,
        "N_sampled_BIN50S": n_sampled,
        "Pearson_r_sampled_BIN50S": pearson_r,
    }


def _plot_mac_an_reg_bin50_sample(
    plot_df: pd.DataFrame,
    output_prefix: str,
    random_seed: int,
) -> dict[str, float | int]:
    """MAC vs AN regression on stratified 50-MAC-bin samples (full MAC range)."""
    out_path = f"{output_prefix}.reg.mac_an.bin50s.png"
    n_in_range = int(len(plot_df))
    sample_df = _stratified_sample_mac_bin50(plot_df, mac_max=None, random_seed=random_seed)
    n_sampled = int(len(sample_df))

    if n_sampled < 10:
        _save_mac_plot_placeholder(
            out_path,
            f"Too few stratified sample points ({n_sampled}) for MAC vs AN (all MAC)",
        )
        pearson_r = np.nan
    else:
        _plot_mac_an_regression(
            sample_df,
            out_path,
            title=f"MAC vs AN (50-MAC bin sample, n<={_MAC_BIN50_SAMPLE_PER_BIN}/bin, all MAC)",
        )
        pearson_r = float(sample_df["MAC"].corr(sample_df["AN"])) if n_sampled > 1 else np.nan

    return {
        "N_in_range_MAC_AN_BIN50S": n_in_range,
        "N_sampled_MAC_AN_BIN50S": n_sampled,
        "Pearson_r_MAC_AN_sampled_BIN50S": pearson_r,
    }


def _plot_R_mac_bin50_sample(
    plot_df: pd.DataFrame,
    output_prefix: str,
    random_seed: int,
) -> dict[str, float | int]:
    """R vs MAC plot on stratified 50-MAC-bin samples (full MAC range)."""
    out_path = f"{output_prefix}.R_mac.bin50s.png"
    sample_df = _stratified_sample_mac_bin50(plot_df, mac_max=None, random_seed=random_seed)
    n_sampled = int(len(sample_df))
    stats: dict[str, float | int] = {"N_plotted_BIN50S": n_sampled}

    if n_sampled < 10:
        _save_mac_plot_placeholder(
            out_path,
            f"Too few stratified sample points ({n_sampled}) for R vs MAC (all MAC)",
        )
        return stats

    _plot_R_mac_scatter(
        sample_df,
        out_path,
        title="R = F_MISS/(1-MAC/N) vs MAC (50-MAC bin sample, n<=100/bin, all MAC)",
        y_lim=(0.0, 1.2),
        stats_df=plot_df,
    )

    near_bound = sample_df["R_miss_boundary"].dropna()
    if len(near_bound):
        stats["Frac_R_ge_0p9_BIN50S"] = float((near_bound >= 0.9).mean())
        stats["Median_R_BIN50S"] = float(near_bound.median())
    stats["N_bin50s_sampled_BIN50S"] = n_sampled
    return stats


def ana_mac_missing_reg_bin50_sample(
    gcount_path: str,
    vmiss_path: str,
    output_prefix: str,
    random_seed: int = 1,
) -> None:
    """
    MAC vs F_MISS / MAC vs AN / R vs MAC plots from stratified 50-MAC-bin samples.

    Draws up to 100 variants per 50-wide MAC bin across all MAC values (full MAC range
    only). MAC vs F_MISS plots include the boundary line y = 1 - MAC/N.
    """
    logger.info("MAC vs missing bin50-sample plots: gcount=%s vmiss=%s", gcount_path, vmiss_path)
    plot_df = _load_mac_miss_extended_df(gcount_path, vmiss_path)

    stats_dict: dict[str, float | int] = {
        "Total_Variants": int(len(plot_df)),
        "N_median": float(plot_df["N"].median()),
        "MAC_bin_width": _MAC_BIN50_WIDTH,
        "Sample_per_bin_cap": _MAC_BIN50_SAMPLE_PER_BIN,
        "Random_seed": random_seed,
    }

    stats_dict.update(_plot_mac_missing_reg_bin50_sample(plot_df, output_prefix, random_seed))
    stats_dict.update(_plot_mac_an_reg_bin50_sample(plot_df, output_prefix, random_seed))
    stats_dict.update(_plot_R_mac_bin50_sample(plot_df, output_prefix, random_seed))

    save_thresholds(stats_dict, f"{output_prefix}.bin50sample.th.tsv")
    print(stats_dict)


def _mac_bin_maf_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate mean/median MAF by MAC bucket for compact info.tsv output."""

    def _bucket(mac: float) -> str:
        if pd.isna(mac):
            return "NA"
        m = int(mac)
        if m == 0:
            return "0"
        if m == 1:
            return "1"
        if m == 2:
            return "2"
        if m <= 9:
            return "3-9"
        if m <= 99:
            return "10-99"
        return ">=100"

    work = df[["MAC", "MAF"]].dropna().copy()
    if work.empty:
        return pd.DataFrame(columns=["mac_bucket", "n_sites", "mean_MAF", "median_MAF"])

    work["mac_bucket"] = work["MAC"].map(_bucket)
    order = ["0", "1", "2", "3-9", "10-99", ">=100"]
    summary = (
        work.groupby("mac_bucket", as_index=False)
        .agg(n_sites=("MAC", "size"), mean_MAF=("MAF", "mean"), median_MAF=("MAF", "median"))
        .assign(
            mac_bucket=lambda d: pd.Categorical(d["mac_bucket"], categories=order, ordered=True),
        )
        .sort_values("mac_bucket")
        .reset_index(drop=True)
    )
    summary["mac_bucket"] = summary["mac_bucket"].astype(str)
    return summary


def _plot_mac_maf_reg_by_mac_range(
    plot_df: pd.DataFrame,
    output_prefix: str,
    mac_max: int,
    regression_max_points: int,
    random_seed: int,
) -> dict[str, float | int]:
    """Scatter + OLS/Huber regression for MAC vs MAF within MAC [0, mac_max]."""
    out_path = f"{output_prefix}.reg.mac0_{mac_max}.png"
    subset = plot_df[(plot_df["MAC"] >= 0) & (plot_df["MAC"] <= mac_max)].copy()
    n_sites = int(len(subset))

    if n_sites < 10:
        _save_mac_plot_placeholder(
            out_path,
            f"Too few variants ({n_sites}) with MAC 0-{mac_max} for MAC vs MAF regression",
        )
        return {
            f"N_MAC_0_{mac_max}": n_sites,
            f"Pearson_r_MAC_MAF_0_{mac_max}": np.nan,
            f"Regression_Points_MAC_0_{mac_max}": 0,
        }

    reg_df = _subsample_for_plot(subset, regression_max_points, random_seed)
    plot_regression_comparison(
        reg_df,
        x_col="MAC",
        y_col="MAF",
        x_label="Minor allele count (MAC)",
        y_label="Minor allele frequency (MAF)",
        filename=out_path,
        title=f"Variant MAC vs MAF (MAC 0-{mac_max})",
        x_lim=(0.0, float(mac_max)),
        y_lim=(0.0, 0.5),
    )
    pearson_r = float(subset["MAC"].corr(subset["MAF"])) if n_sites > 1 else np.nan
    return {
        f"N_MAC_0_{mac_max}": n_sites,
        f"Pearson_r_MAC_MAF_0_{mac_max}": pearson_r,
        f"Regression_Points_MAC_0_{mac_max}": int(len(reg_df)),
    }


def ana_mac_maf_reg(
    gcount_path: str,
    output_prefix: str,
    regression_max_points: int = 50000,
    random_seed: int = 1,
) -> None:
    """
    Per-site MAC vs MAF scatter regression from PLINK2 ``--geno-counts`` (.gcount).

    Writes MAC-bucket MAF summaries, regression thresholds, and OLS/Huber plots
    (full MAC range plus MAC 0-100, 0-500, and 0-1000 subsets).
    """
    logger.info("Processing MAC vs MAF regression: gcount=%s", gcount_path)

    df_raw = load_df_from_plink_gcount(gcount_path)
    df_mac = _compute_mac_table(df_raw)

    plot_df = df_mac[["MAC", "MAF"]].replace([np.inf, -np.inf], np.nan).dropna()
    if plot_df.empty:
        fail("No plottable variants after dropping NA MAC/MAF")

    save_df_to_tsv(_mac_bin_maf_summary(plot_df), f"{output_prefix}.info.tsv")

    reg_df = _subsample_for_plot(plot_df, regression_max_points, random_seed)
    n_variants = int(len(plot_df))
    n_reg = int(len(reg_df))
    pearson_r = float(plot_df["MAC"].corr(plot_df["MAF"])) if n_variants > 1 else np.nan

    stats_dict = {
        "Total_Variants": n_variants,
        "Regression_Points": n_reg,
        "Pearson_r_MAC_MAF": pearson_r,
        "Mean_MAC": float(plot_df["MAC"].mean()),
        "Median_MAC": float(plot_df["MAC"].median()),
        "Mean_MAF": float(plot_df["MAF"].mean()),
        "Median_MAF": float(plot_df["MAF"].median()),
        "Regression_Subsampled": int(n_reg < n_variants),
    }

    if len(reg_df) < 10:
        _save_mac_plot_placeholder(
            f"{output_prefix}.reg.png",
            f"Too few variants ({len(reg_df)}) for MAC vs MAF regression",
        )
    else:
        plot_regression_comparison(
            reg_df,
            x_col="MAC",
            y_col="MAF",
            x_label="Minor allele count (MAC)",
            y_label="Minor allele frequency (MAF)",
            filename=f"{output_prefix}.reg.png",
            title="Variant MAC vs MAF",
            y_lim=(0.0, 0.5),
        )

    for mac_max in _MAC_REG_RANGES:
        stats_dict.update(
            _plot_mac_maf_reg_by_mac_range(
                plot_df,
                output_prefix,
                mac_max,
                regression_max_points,
                random_seed,
            )
        )

    save_thresholds(stats_dict, f"{output_prefix}.th.tsv")
    print(stats_dict)


def redraw_mac_dist_log_from_info(
    info_path: str,
    output_prefix: str,
    mac_max: int = 100,
    tick_step: int = 5,
) -> None:
    """Regenerate ``*.dist.0_100.log.png`` from an existing ``*.variant.mac.info.tsv``."""
    from infra.utils.io import load_df_from_tsv

    logger.info("Redraw MAC log distribution from: %s", info_path)
    mac_counts = load_df_from_tsv(info_path)
    if mac_counts.empty:
        fail(f"No MAC site counts in {info_path}")
    if not {"MAC", "n_sites"}.issubset(mac_counts.columns):
        fail(f"Expected columns MAC, n_sites in {info_path}")
    _plot_mac_site_dist_0_100(
        mac_counts,
        output_prefix,
        mac_max=mac_max,
        tick_step=tick_step,
        log_y=True,
    )


if __name__ == "__main__":
    configure_logging()
    if len(sys.argv) > 1:
        try:
            ana_mac_stats(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "mac_out")
        except (DataLoadError, FileNotFoundError, ValueError) as exc:
            logger.error("%s", exc)
            sys.exit(1)
    else:
        print("Usage: python mac.py <input.gcount> <output_prefix>")
        sys.exit(1)
