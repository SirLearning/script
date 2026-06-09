from .variant_utils import load_df_from_plink_gcount

import logging
import sys

import numpy as np
import pandas as pd

from infra.utils.errors import DataLoadError, configure_logging, fail
from infra.utils.graph import (
    TITLE_FONT_SIZE,
    X_LABEL_FONT_SIZE,
    Y_LABEL_FONT_SIZE,
    TICK_FONT_SIZE,
    plot_heatmap_custom,
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
