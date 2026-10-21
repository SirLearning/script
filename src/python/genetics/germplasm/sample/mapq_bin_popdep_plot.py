"""QC scatter plots for 1 Mb population depth grids."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from genetics.germplasm.sample.mapq_bin_popdep import DepthMode, normalize_depth_mode
from infra.utils.graph import plot_joint_regression

# Marginal histogram bin counts (match popdep het-filter scatters).
POPDEP_CNV_BINS_X = 250
POPDEP_CNV_BINS_Y = 150

CNV_FRAC_BIN_LABELS = (
    "0.0-0.2",
    "0.2-0.4",
    "0.4-0.6",
    "0.6-0.8",
    "0.8-1.0",
)

CNV_FRAC_BIN_COLORS = {
    "0.0-0.2": "#ffffff",
    "0.2-0.4": "#fee0d2",
    "0.4-0.6": "#fc9272",
    "0.6-0.8": "#cb181d",
    "0.8-1.0": "#67000d",
}

CNV_FRAC_BIN_ALPHA = {
    "0.0-0.2": 0.15,
    "0.2-0.4": 0.40,
    "0.4-0.6": 0.55,
    "0.6-0.8": 0.70,
    "0.8-1.0": 0.85,
}


def _scatter_specs(depth_mode: DepthMode) -> tuple[dict, ...]:
    mode = normalize_depth_mode(depth_mode)
    if mode == "all":
        depth_label = "Mean depth"
        depth_sd = "Depth SD"
        depth_cv = "Depth CV"
    else:
        depth_label = "Mean MQ30+ depth"
        depth_sd = "MQ30+ depth SD"
        depth_cv = "MQ30+ depth CV"
    return (
        {
            "tag": "depth",
            "mean": "depth_mean",
            "sd": "depth_sd",
            "cv": "depth_cv",
            "log_y_sd": "depth_sd",
            "log_y_cv": "depth_cv",
            "label_mean": depth_label,
            "label_sd": depth_sd,
            "label_cv": depth_cv,
        },
        {
            "tag": "reldepth",
            "mean": "rel_mean",
            "sd": "rel_sd",
            "cv": "rel_cv",
            "log_y_sd": "rel_var",
            "log_y_cv": "rel_cv",
            "label_mean": "Relative mean depth",
            "label_sd": "Relative depth SD",
            "label_cv": "Relative depth CV",
            "label_log_y_sd": "Relative depth variance",
        },
    )


# Backward-compatible alias for tests / imports
MAPQ_BIN_POPDEP_SCATTER_SPECS = _scatter_specs("all")


def assign_cnv_frac_bin(
    cnv_frac: pd.Series,
    *,
    labels: tuple[str, ...] = CNV_FRAC_BIN_LABELS,
) -> pd.Series:
    """Bucket cohort CNV fraction into labelled intervals for hue / marginals."""
    frac = pd.to_numeric(cnv_frac, errors="coerce").clip(lower=0.0, upper=1.0)
    edges = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0 + 1e-9]
    return pd.cut(
        frac,
        bins=edges,
        labels=list(labels),
        include_lowest=True,
        right=False,
    )


def _add_cv_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    depth_mean = pd.to_numeric(out["depth_mean"], errors="coerce")
    depth_sd = pd.to_numeric(out["depth_sd"], errors="coerce")
    rel_mean = pd.to_numeric(out["rel_mean"], errors="coerce")
    rel_sd = pd.to_numeric(out["rel_sd"], errors="coerce")
    out["depth_cv"] = np.where(depth_mean > 0, depth_sd / depth_mean, np.nan)
    out["rel_cv"] = np.where(rel_mean > 0, rel_sd / rel_mean, np.nan)
    out["depth_var"] = depth_sd**2
    out["rel_var"] = rel_sd**2
    return out


def _prepare_cnv_joint_frame(
    grid: pd.DataFrame,
    *,
    x_col: str,
    y_col: str,
    color_col: str,
    log_axes: bool,
) -> pd.DataFrame:
    plot_df = grid[[x_col, y_col, color_col]].replace([np.inf, -np.inf], np.nan).dropna()
    if log_axes:
        plot_df = plot_df[(plot_df[x_col] > 0) & (plot_df[y_col] > 0)]
    if plot_df.empty:
        return plot_df
    plot_df = plot_df.copy()
    plot_df["cnv_frac_bin"] = assign_cnv_frac_bin(plot_df[color_col])
    plot_df = plot_df.dropna(subset=["cnv_frac_bin"])
    plot_df["cnv_frac_bin"] = plot_df["cnv_frac_bin"].astype(str)
    return plot_df


def plot_mapq_bin_popdep_joint_scatter(
    grid: pd.DataFrame,
    *,
    x_col: str,
    y_col: str,
    color_col: str,
    title: str,
    filename: str | Path,
    x_label: str,
    y_label: str,
    log_axes: bool = False,
) -> bool:
    """
    JointGrid scatter with 100% fill x/y marginals stacked by CNV fraction bin.

    Point colour and marginal stacks use the same five CNV cohort-fraction
    intervals (0.0-0.2 … 0.8-1.0) with increasing red intensity.
    """
    plot_df = _prepare_cnv_joint_frame(
        grid,
        x_col=x_col,
        y_col=y_col,
        color_col=color_col,
        log_axes=log_axes,
    )
    if plot_df.empty:
        print(f"[Warning] plot_mapq_bin_popdep_joint_scatter: no data for {filename}")
        return False

    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    plot_joint_regression(
        plot_df,
        x_col=x_col,
        y_col=y_col,
        group_col="cnv_frac_bin",
        x_label=x_label,
        y_label=y_label,
        filename=str(filename),
        title=title,
        figure_size=(10, 8),
        palette=CNV_FRAC_BIN_COLORS,
        hue_order=list(CNV_FRAC_BIN_LABELS),
        log_x=log_axes,
        log_y=log_axes,
        draw_regression=False,
        alpha=0.55,
        s=18,
        rasterized=True,
        legend_title="CNV sample fraction",
        bins_x=POPDEP_CNV_BINS_X,
        bins_y=POPDEP_CNV_BINS_Y,
        min_points=1,
        marginal_multiple="fill",
        stack_last_on_top=True,
        hue_alpha=CNV_FRAC_BIN_ALPHA,
    )
    return Path(filename).exists()


def plot_mapq_bin_popdep_qc_scatters(
    grid: pd.DataFrame,
    output_prefix: str | Path,
    *,
    color_col: str = "cnv_frac",
    depth_mode: DepthMode = "all",
) -> list[Path]:
    """
    Write linear and log-log mean-vs-dispersion JointGrid scatters.

    Eight PNGs: SD vs mean and CV vs mean (linear + log) for depth and reldepth.
    Marginal x/y histograms show stacked CNV fraction intervals (fill, 100%).
    """
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    if color_col not in grid.columns:
        raise ValueError(f"Grid missing color column {color_col!r} for CNV overlay.")

    work = _add_cv_columns(grid)
    scatter_specs = _scatter_specs(depth_mode)

    for spec in scatter_specs:
        mean_col = spec["mean"]
        tag = spec["tag"]

        for disp, y_key, log_y_key, label_key in (
            ("sd", "sd", "log_y_sd", "label_sd"),
            ("cv", "cv", "log_y_cv", "label_cv"),
        ):
            y_col = spec[y_key]
            log_y_col = spec[log_y_key]
            label_y = spec[label_key]
            log_label_y = spec.get(f"label_{log_y_key}", label_y)

            linear_path = (
                output_prefix.parent
                / f"{output_prefix.name}.{tag}.scatter_{disp}_vs_mean.png"
            )
            if plot_mapq_bin_popdep_joint_scatter(
                work,
                x_col=mean_col,
                y_col=y_col,
                color_col=color_col,
                title=f"{label_y} vs {spec['label_mean']} ({tag})",
                filename=linear_path,
                x_label=spec["label_mean"],
                y_label=label_y,
                log_axes=False,
            ):
                written.append(linear_path)

            log_path = (
                output_prefix.parent
                / f"{output_prefix.name}.{tag}.scatter_log_{disp}_vs_log_mean.png"
            )
            if plot_mapq_bin_popdep_joint_scatter(
                work,
                x_col=mean_col,
                y_col=log_y_col,
                color_col=color_col,
                title=f"Log {log_label_y} vs log {spec['label_mean']} ({tag})",
                filename=log_path,
                x_label=f"Log {spec['label_mean']}",
                y_label=f"Log {log_label_y}",
                log_axes=True,
            ):
                written.append(log_path)

    return written
