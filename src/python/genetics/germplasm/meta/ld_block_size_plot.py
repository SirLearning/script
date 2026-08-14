"""Plot LD block size distributions from Wu et al. 2025 supplementary tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from infra.utils.graph import plot_distribution_with_stats, plot_scatter_with_thresholds
from infra.utils.io import load_df_generic, save_df_to_tsv

BLOCK_SIZE_COL = "Block size (Mb)"
SUMMARY_COLUMNS = (
    "block_number",
    "chromosome",
    "start_bp",
    "end_bp",
    "n_genes",
    "block_size_mb",
)


def load_wu_ld_block_table(csv_path: str | Path) -> pd.DataFrame:
    """
    Load Wu 2025 Table S26 LD block listing.

    The export includes a title row before the column header; it is skipped automatically.
    """
    raw = load_df_generic(str(csv_path), header=1)
    rename = {
        "Block number": "block_number",
        "Chromosome": "chromosome",
        "Start (V2.1)": "start_bp",
        "End (V2.1)": "end_bp",
        "Number of Gene": "n_genes",
        BLOCK_SIZE_COL: "block_size_mb",
    }
    missing = [src for src in rename if src not in raw.columns]
    if missing:
        raise ValueError(f"Missing expected columns in LD block CSV: {missing}")

    out = raw.rename(columns=rename)[list(SUMMARY_COLUMNS)].copy()
    out["block_size_mb"] = pd.to_numeric(out["block_size_mb"], errors="coerce")
    out["start_bp"] = pd.to_numeric(out["start_bp"], errors="coerce")
    out["end_bp"] = pd.to_numeric(out["end_bp"], errors="coerce")
    out["n_genes"] = pd.to_numeric(out["n_genes"], errors="coerce")
    return out.dropna(subset=["block_size_mb", "n_genes"]).reset_index(drop=True)


def summarize_ld_block_sizes(block_df: pd.DataFrame) -> pd.DataFrame:
    """Return one-row summary statistics for block size (Mb)."""
    sizes = block_df["block_size_mb"]
    return pd.DataFrame(
        [
            {
                "n_blocks": int(sizes.shape[0]),
                "mean_mb": float(sizes.mean()),
                "median_mb": float(sizes.median()),
                "std_mb": float(sizes.std(ddof=1)),
                "min_mb": float(sizes.min()),
                "max_mb": float(sizes.max()),
                "q25_mb": float(sizes.quantile(0.25)),
                "q75_mb": float(sizes.quantile(0.75)),
                "q95_mb": float(sizes.quantile(0.95)),
                "q99_mb": float(sizes.quantile(0.99)),
            }
        ]
    )


def plot_ld_block_size_distribution(
    block_df: pd.DataFrame,
    output_png: str | Path,
    *,
    bins: int = 80,
) -> Path:
    """Histogram of LD block length (Mb) with mean and median markers."""
    sizes = block_df["block_size_mb"]
    plot_distribution_with_stats(
        block_df,
        col="block_size_mb",
        title="LD block size distribution (Wu et al. 2025)",
        filename=str(output_png),
        mean_val=float(sizes.mean()),
        median_val=float(sizes.median()),
        std_val=float(sizes.std(ddof=1)),
        x_label="Block size (Mb)",
        y_label="Number of blocks",
        bins=bins,
        color="#4c72b0",
    )
    return Path(output_png)


def plot_ld_block_size_distribution_log_y(
    block_df: pd.DataFrame,
    output_png: str | Path,
    *,
    bins: int = 80,
) -> Path:
    """Histogram of LD block length (Mb) with log-scaled block counts on y-axis."""
    sizes = block_df["block_size_mb"]
    plot_distribution_with_stats(
        block_df,
        col="block_size_mb",
        title="LD block size distribution (Wu et al. 2025, log count)",
        filename=str(output_png),
        mean_val=float(sizes.mean()),
        median_val=float(sizes.median()),
        std_val=float(sizes.std(ddof=1)),
        x_label="Block size (Mb)",
        y_label="Number of blocks",
        bins=bins,
        color="#4c72b0",
        log_scale=True,
    )
    return Path(output_png)


def plot_ld_block_size_vs_n_genes(
    block_df: pd.DataFrame,
    output_png: str | Path,
) -> Path:
    """Scatter of block size (Mb) vs number of genes with mean/median vertical guides."""
    sizes = block_df["block_size_mb"]
    mean_mb = float(sizes.mean())
    median_mb = float(sizes.median())
    std_mb = float(sizes.std(ddof=1))
    mean_label = f"Mean: {mean_mb:.4f}\nSD: {std_mb:.4f}"
    plot_scatter_with_thresholds(
        block_df,
        x_col="block_size_mb",
        y_col="n_genes",
        title="LD block size vs number of genes (Wu et al. 2025)",
        filename=str(output_png),
        xlabel="Block size (Mb)",
        ylabel="Number of genes",
        color="#4c72b0",
        alpha=0.35,
        s=18,
        figure_size=(10, 6),
        thresholds_v=[
            {
                "value": mean_mb,
                "color": "red",
                "linestyle": "--",
                "label": mean_label,
            },
            {
                "value": median_mb,
                "color": "orange",
                "linestyle": "-",
                "label": f"Median: {median_mb:.4f}",
            },
        ],
    )
    return Path(output_png)


def run_ld_block_size_plot(csv_path: str | Path, output_dir: str | Path) -> dict[str, Path]:
    """Parse Wu 2025 LD block table, write summary TSV, and emit distribution/scatter PNGs."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    block_df = load_wu_ld_block_table(csv_path)
    summary_df = summarize_ld_block_sizes(block_df)

    blocks_tsv = out_dir / "wu2025_ld_block_sizes.tsv"
    summary_tsv = out_dir / "wu2025_ld_block_size_summary.tsv"
    save_df_to_tsv(block_df, str(blocks_tsv))
    save_df_to_tsv(summary_df, str(summary_tsv))

    plot_path = plot_ld_block_size_distribution(block_df, out_dir / "wu2025_ld_block_size_distribution.png")
    log_y_plot_path = plot_ld_block_size_distribution_log_y(
        block_df, out_dir / "wu2025_ld_block_size_distribution_logy.png"
    )
    scatter_plot_path = plot_ld_block_size_vs_n_genes(
        block_df, out_dir / "wu2025_ld_block_size_vs_n_genes.png"
    )
    return {
        "blocks_tsv": blocks_tsv,
        "summary_tsv": summary_tsv,
        "distribution_plot": plot_path,
        "distribution_logy_plot": log_y_plot_path,
        "size_vs_genes_plot": scatter_plot_path,
    }
