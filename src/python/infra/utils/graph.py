import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
import seaborn as sns
import numpy as np
import pandas as pd
import sys
from sklearn.linear_model import HuberRegressor, LinearRegression
from scipy.ndimage import gaussian_filter

TITLE_FONT_SIZE = 18
X_LABEL_FONT_SIZE = 14
Y_LABEL_FONT_SIZE = 14
TICK_FONT_SIZE = 12
LEGEND_FONT_SIZE = 12
# Heatmaps use figsize (14, 10); line-plot sizes look small at full-frame view.
HEATMAP_TITLE_FONT_SIZE = 26
HEATMAP_AXIS_FONT_SIZE = 20
HEATMAP_TICK_FONT_SIZE = 16
HEATMAP_LEGEND_FONT_SIZE = 16

# Two-category contrast: blue vs red (not red vs green — readable with deuteranopia).
TWO_CATEGORY_BLUE = "#4c72b0"
TWO_CATEGORY_RED = "#c44e52"
TWO_CATEGORY_COLORS = (TWO_CATEGORY_BLUE, TWO_CATEGORY_RED)

# A/B/D panels: blue / red / purple (no green — same deuteranopia rule as two-category).
SUBGENOME_ORDER = ("A", "B", "D")
SUBGENOME_PURPLE = "#8172b3"
SUBGENOME_COLORS = {
    "A": TWO_CATEGORY_BLUE,
    "B": TWO_CATEGORY_RED,
    "D": SUBGENOME_PURPLE,
}

# Missing-rate heatmaps: red = low missing (better), blue = high missing (worse).
MISSING_RATE_CMAP = "coolwarm_r"

# Sequential intensity: colorless (white) to red. Use with vmin/vmax, not TwoSlopeNorm.
SEQUENTIAL_WHITE_RED = mcolors.LinearSegmentedColormap.from_list(
    "sequential_white_red",
    [
        (0.00, "#ffffff"),
        (0.20, "#fee0d2"),
        (0.45, "#fc9272"),
        (0.70, "#cb181d"),
        (1.00, "#67000d"),
    ],
)

# Signed GAM residual scatter: vivid blue / gray at 0 / vivid red (used with TwoSlopeNorm).
RESIDUAL_SIGNED_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "residual_signed_vivid",
    [
        (0.00, "#08306b"),
        (0.30, "#2171b5"),
        (0.46, "#6baed6"),
        (0.50, "#bdbdbd"),
        (0.54, "#fc9272"),
        (0.70, "#cb181d"),
        (1.00, "#67000d"),
    ],
)

# Canonical sample_group.txt categories — fixed hue order and colors across all stats plots.
SAMPLE_GROUP_ORDER = (
    "A",
    "AB",
    "ABD",
    "D",
    "HZNU",
    "Nature",
    "S",
    "WAP",
    "Watkins",
    "Unknown",
)
SAMPLE_GROUP_UNKNOWN_COLOR = (0.5, 0.5, 0.5)

DEFAULT_AXIS_PADDING_FRACTION = 0.05


def padded_axis_limits(values, fraction=DEFAULT_AXIS_PADDING_FRACTION):
    """
    Return (lo, hi) spanning finite data with proportional padding on each side.

    When all values are identical, use a small absolute span so the axis is still visible.
    """
    s = pd.to_numeric(pd.Series(values), errors="coerce").replace(
        [np.inf, -np.inf], np.nan
    ).dropna()
    if s.empty:
        return None
    lo, hi = float(s.min()), float(s.max())
    span = hi - lo
    if span == 0:
        span = max(abs(lo), abs(hi), 1.0) * 0.01
    pad = span * fraction
    return lo - pad, hi + pad


def padded_log_axis_limits(values, fraction=DEFAULT_AXIS_PADDING_FRACTION):
    """Return (lo, hi) in log10 units with proportional padding; ignores non-positive values."""
    s = pd.to_numeric(pd.Series(values), errors="coerce").replace(
        [np.inf, -np.inf], np.nan
    ).dropna()
    s = s[s > 0]
    if s.empty:
        return None
    return padded_axis_limits(np.log10(s.to_numpy(dtype=np.float64)), fraction=fraction)


def adaptive_axis_max(
    values,
    *,
    min_points: int = 10,
    gap_multiplier: float = 5.0,
) -> float | None:
    """
    Upper display bound by peeling isolated high values from the maximum downward.

    From the sorted maximum, forms the top cluster (consecutive points separated by
    at most ``gap_multiplier`` times the typical step). If that cluster is a single
    point and the gap to the rest of the data is much larger than typical spacing,
    exclude it and repeat. Dense multi-point tails are kept even when separated from
    the bulk.

    The minimum is never adjusted. Returns ``None`` when trimming does not apply.
    """
    s = pd.to_numeric(pd.Series(values), errors="coerce").replace(
        [np.inf, -np.inf], np.nan
    ).dropna()
    if len(s) < min_points:
        return None

    vals = np.sort(s.values.astype(float))
    raw_max = float(vals[-1])
    trimmed_max = raw_max

    while len(vals) >= 2:
        diffs = np.diff(vals)
        if len(diffs) == 0:
            break

        sorted_diffs = np.sort(diffs)
        half = max(1, len(sorted_diffs) // 2)
        baseline = float(np.median(sorted_diffs[:half]))
        if baseline <= 0:
            baseline = max(float(np.std(vals)), 1e-12)

        top_idx = len(vals) - 1
        cluster_start = top_idx
        while cluster_start > 0:
            step = vals[cluster_start] - vals[cluster_start - 1]
            if step <= gap_multiplier * baseline:
                cluster_start -= 1
                continue
            break

        if cluster_start == 0:
            break

        gap_below = vals[cluster_start] - vals[cluster_start - 1]
        cluster_size = top_idx - cluster_start + 1
        if gap_below <= gap_multiplier * baseline or cluster_size >= 2:
            break

        trimmed_max = float(vals[cluster_start - 1])
        vals = vals[:cluster_start]

    tol = 1e-12 * max(abs(raw_max), 1.0)
    if trimmed_max >= raw_max - tol:
        return None
    return trimmed_max


def axis_limits_with_adaptive_upper(
    values,
    *,
    min_points: int = 10,
    gap_multiplier: float = 5.0,
    padding_fraction: float = DEFAULT_AXIS_PADDING_FRACTION,
    lower_override: float | None = None,
) -> tuple[float, float] | None:
    """
    Matplotlib axis limits with adaptive upper trim and padded lower bound.

    The lower bound defaults to the data minimum minus padding; pass
    ``lower_override`` (e.g. ``0.0`` for non-negative CV axes) to fix the floor.
    Returns ``None`` when no upper trim is needed.
    """
    s = pd.to_numeric(pd.Series(values), errors="coerce").replace(
        [np.inf, -np.inf], np.nan
    ).dropna()
    if len(s) < min_points:
        return None

    lo_raw = float(s.min())
    hi_raw = float(s.max())
    hi_trim = adaptive_axis_max(
        s,
        min_points=min_points,
        gap_multiplier=gap_multiplier,
    )
    if hi_trim is None:
        return None

    span_full = hi_raw - lo_raw
    if span_full <= 0:
        span_full = max(abs(lo_raw), abs(hi_raw), 1.0) * 0.01
    lo = lo_raw - span_full * padding_fraction
    if lower_override is not None:
        lo = lower_override
    hi = hi_trim + (hi_trim - lo_raw) * padding_fraction
    return lo, hi


def sample_group_palette() -> dict[str, tuple]:
    """Fixed Group → color map (seaborn ``deep``, same family as default categorical plots)."""
    colors = sns.color_palette("deep", len(SAMPLE_GROUP_ORDER))
    return dict(zip(SAMPLE_GROUP_ORDER, colors))


def sample_group_hue_config(values) -> tuple[list[str], dict[str, tuple]]:
    """
    Return ``(hue_order, palette)`` for a ``Group`` column.

    Groups in ``SAMPLE_GROUP_ORDER`` keep stable colors; any extra labels get husl fallbacks.
    """
    base = sample_group_palette()
    present = {str(v) for v in pd.Series(values).dropna().unique()}
    hue_order = [g for g in SAMPLE_GROUP_ORDER if g in present]
    extras = sorted(g for g in present if g not in SAMPLE_GROUP_ORDER)
    if extras:
        for g, c in zip(extras, sns.color_palette("husl", len(extras))):
            base[g] = c
        hue_order.extend(extras)
    palette_map = {g: base.get(g, SAMPLE_GROUP_UNKNOWN_COLOR) for g in hue_order}
    return hue_order, palette_map


def _hue_kwargs_for_sample_group(df: pd.DataFrame, group_col: str) -> dict:
    """Seaborn ``hue_order`` / ``palette`` when ``group_col`` is the sample ``Group`` column."""
    if group_col != "Group" or group_col not in df.columns:
        return {}
    hue_order, palette_map = sample_group_hue_config(df[group_col])
    return {"hue_order": hue_order, "palette": palette_map}


def draw_x_interval_backgrounds(
    ax,
    intervals: list[dict],
    *,
    xlim: tuple[float, float] | None = None,
    zorder: int = 0,
) -> list:
    """
    Shade vertical intervals on ``ax`` (e.g. copy-number classes on rel_depth).

    Each interval dict: ``label``, ``lo``, ``hi``, optional ``color``, ``alpha``.
    Returns legend handles (matplotlib Patch objects).
    """
    import matplotlib.patches as mpatches

    handles = []
    x_hi_cap = xlim[1] if xlim is not None else None
    x_lo_cap = xlim[0] if xlim is not None else None
    for spec in intervals:
        lo = float(spec["lo"])
        hi = float(spec["hi"])
        if x_lo_cap is not None:
            lo = max(lo, x_lo_cap)
        if x_hi_cap is not None and np.isfinite(hi):
            hi = min(hi, x_hi_cap)
        if not np.isfinite(hi):
            hi = x_hi_cap if x_hi_cap is not None else lo + 1.0
        if hi <= lo:
            continue
        color = spec.get("color", "#eeeeee")
        alpha = spec.get("alpha", 0.35)
        ax.axvspan(lo, hi, color=color, alpha=alpha, linewidth=0, zorder=zorder)
        handles.append(
            mpatches.Patch(
                facecolor=color,
                alpha=alpha,
                edgecolor="0.5",
                linewidth=0.5,
                label=str(spec.get("label", f"{lo}-{hi}")),
            )
        )
    return handles


def _finite_for_vline(val) -> bool:
    """True if val can be drawn as a vertical line on a plot axis."""
    if val is None:
        return False
    try:
        return bool(np.isfinite(float(val)))
    except (TypeError, ValueError):
        return False


def combine_plots(
    images,
    output_file="combined_plot.png",
    orientation='h',
):
    """
    Combines multiple images into a single figure.
    orientation: 'h' or 'v'
    """
    # 1. Pre-load images and calculate aspect ratio
    loaded_imgs = []
    aspect_ratio = 1.2 # Default backup ratio (width / height)
    ratio_determined = False
    
    for img_path in images:
        try:
            img = mpimg.imread(img_path)
            loaded_imgs.append(img)
            
            # Use the first valid image to determine uniform aspect ratio
            if not ratio_determined:
                h, w = img.shape[:2]
                aspect_ratio = w / h
                ratio_determined = True
        except Exception as e:
            print(f"Error loading {img_path}: {e}")
            loaded_imgs.append(None)

    if not any(loaded_imgs):
        print("No valid images to combine.")
        return

    n_images = len(images)
    
    # 2. Dynamic Figure Configuration
    if orientation == 'h':
        rows = 1
        cols = n_images
        # Fix height, scale width to maintain aspect ratio
        plot_height = 5
        plot_width = plot_height * aspect_ratio * n_images
    else: # vertical
        rows = n_images
        cols = 1
        # Fix width, scale height to maintain aspect ratio
        plot_width = 6
        # w / h = aspect_ratio  =>  h = w / aspect_ratio
        img_height = plot_width / aspect_ratio
        plot_height = img_height * n_images

    # Create figure
    fig, axes = plt.subplots(rows, cols, figsize=(plot_width, plot_height))

    # Ensure axes is iterable even for a single image
    if n_images == 1:
        axes = [axes]
    elif isinstance(axes, np.ndarray):
        axes = axes.flatten()

    for ax, img in zip(axes, loaded_imgs):
        if img is not None:
            ax.imshow(img)
            ax.axis('off')  # Hide axis
        else:
            ax.text(0.5, 0.5, 'Image Not Found', ha='center', va='center')
            ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close(fig)
    print(f"Combined image saved to {output_file}")


def plot_multi_line_series(
    data,
    x_col,
    y_specs,
    title,
    filename,
    x_label=None,
    y_label=None,
    figure_size=(10, 6),
    xlim=None,
    ylim=None,
    rotate_xlabels=None,
):
    """
    Plot multiple numeric y-series against one x column (line plot).

    y_specs: list of dicts, each with required key ``y_col`` and optional
    ``label``, ``color``, ``linestyle``, ``linewidth``, ``alpha``, ``zorder``
    (defaults match common decay-style plots).
    """
    sns.set_style("white")
    plt.figure(figsize=figure_size)

    if x_label is None:
        x_label = x_col
    if y_label is None:
        y_label = ""

    for spec in y_specs:
        yc = spec["y_col"]
        if yc not in data.columns:
            continue
        lbl = spec.get("label", yc)
        clr = spec.get("color", "steelblue")
        ls = spec.get("linestyle", "-")
        lw = float(spec.get("linewidth", 1.5))
        marker = spec.get("marker")
        plot_kwargs = dict(color=clr, linestyle=ls, linewidth=lw, label=lbl)
        if marker:
            plot_kwargs["marker"] = marker
        if spec.get("zorder") is not None:
            plot_kwargs["zorder"] = spec["zorder"]
        if "alpha" in spec and spec.get("alpha") is not None:
            plot_kwargs["alpha"] = float(spec["alpha"])
        x = data[x_col]
        y = data[yc]
        mask = pd.notna(x) & pd.notna(y)
        plt.plot(x[mask], y[mask], **plot_kwargs)

    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    plt.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    plt.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, -0.15),
        fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    if rotate_xlabels is not None:
        plt.xticks(rotation=rotate_xlabels, ha="right")

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved plot to {filename}")
    plt.close()


def plot_dual_axis_line_series(
    data,
    x_col,
    left_spec,
    right_spec,
    title,
    filename,
    x_label=None,
    figure_size=(10, 6),
    rotate_xlabels=None,
):
    """
    Line plot with independent left/right y-axes (twinx).

    left_spec / right_spec: dicts with required ``y_col`` and optional
    ``label``, ``color``, ``linestyle``, ``linewidth``, ``marker``, ``y_label``.
    """
    sns.set_style("white")
    fig, ax_left = plt.subplots(figsize=figure_size)
    ax_right = ax_left.twinx()

    def _plot_axis(ax, spec, axis_name):
        y_col = spec["y_col"]
        color = spec.get("color", "steelblue" if axis_name == "left" else "darkorange")
        label = spec.get("label", y_col)
        ls = spec.get("linestyle", "-")
        lw = float(spec.get("linewidth", 1.8))
        marker = spec.get("marker")
        plot_kwargs = dict(color=color, linestyle=ls, linewidth=lw, label=label)
        if marker:
            plot_kwargs["marker"] = marker
        ax.plot(data[x_col], data[y_col], **plot_kwargs)
        y_label = spec.get("y_label", y_col)
        ax.set_ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE, color=color)
        ax.tick_params(axis="y", labelcolor=color, labelsize=TICK_FONT_SIZE)
        return ax

    _plot_axis(ax_left, left_spec, "left")
    _plot_axis(ax_right, right_spec, "right")

    if x_label is None:
        x_label = x_col
    ax_left.set_xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    ax_left.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax_left.tick_params(axis="x", labelsize=TICK_FONT_SIZE)
    if rotate_xlabels is not None:
        ax_left.tick_params(axis="x", rotation=rotate_xlabels)
        for tick in ax_left.get_xticklabels():
            tick.set_ha("right")

    lines_left, labels_left = ax_left.get_legend_handles_labels()
    lines_right, labels_right = ax_right.get_legend_handles_labels()
    ax_left.legend(
        lines_left + lines_right,
        labels_left + labels_right,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.15),
        fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved plot to {filename}")
    plt.close()


def plot_distribution_with_stats(
    data, 
    col, 
    title, 
    filename, 
    mean_val=None, 
    median_val=None,
    std_val=None, 
    thresholds=None,
    x_label=None, 
    y_label="Count", 
    log_scale=False,
    color='steelblue',
    bins=100,
    xlim=None,
    figure_size=(10,6)
):
    """
    Plots a histogram distribution with mean, median, and optional threshold lines.
    If xlim is provided, bins will be calculated strictly within that range.
    
    Args:
        thresholds (list of dict): Optional list of thresholds to draw vertical lines for.
                                   Format: [{'value': float, 'label': str, 'color': str, 'linestyle': str}]
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    
    sns.set_style("white") # No background pattern

    if x_label is None:
        x_label = col

    plt.figure(figsize=figure_size)
    
    # Determine bin edges if xlim is provided to ensure fixed bin count within the visible range
    plot_bins = bins
    if xlim is not None and isinstance(bins, int):
        plot_bins = np.linspace(xlim[0], xlim[1], bins + 1)

    # Histogram
    sns.histplot(data=data, x=col, bins=plot_bins, color=color, edgecolor='white', linewidth=0.1)
    
    # Statistics Lines
    if _finite_for_vline(mean_val):
        label_text = f'Mean: {float(mean_val):.4f}'
        if _finite_for_vline(std_val):
            label_text += f'\nSD: {float(std_val):.4f}'
        plt.axvline(
            x=float(mean_val),
            color=TWO_CATEGORY_RED,
            linestyle='--',
            linewidth=1.5,
            label=label_text,
        )

    if _finite_for_vline(median_val):
        plt.axvline(
            x=float(median_val),
            color=TWO_CATEGORY_BLUE,
            linestyle='-',
            linewidth=1.5,
            label=f'Median: {float(median_val):.4f}',
        )

    # Additional Threshold Lines
    if thresholds:
        for thr in thresholds:
            val = thr.get('value')
            if not _finite_for_vline(val):
                continue
            v = float(val)
            lbl = thr.get('label', f'{v:.4f}')
            clr = thr.get('color', 'black')
            ls = thr.get('linestyle', '--')
            lw = thr.get('linewidth', 1.5)
            plt.axvline(x=v, color=clr, linestyle=ls, linewidth=lw, label=lbl)
    
    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    
    plt.tick_params(axis='both', which='major', labelsize=TICK_FONT_SIZE)

    if log_scale:
        plt.yscale('log')
        plt.ylabel(f"{y_label} (Log Scale)", fontsize=Y_LABEL_FONT_SIZE)
    else:
        plt.ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
        
    if xlim is not None:
        plt.xlim(xlim)

    # Legend below left
    plt.legend(loc='upper left', bbox_to_anchor=(0.0, -0.15), fontsize=LEGEND_FONT_SIZE, frameon=False)
    
    # plt.tight_layout() # Removed to prevent shrinking of axes
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved plot to {filename}")
    plt.close()


def plot_faceted_histograms(
    data,
    x_col,
    facet_col,
    facet_order,
    title,
    filename,
    xlabel=None,
    ylabel="Count",
    log_scale=False,
    bins=100,
    xlim=None,
    ylim=None,
    facet_colors=None,
    facet_labels=None,
    show_mean=True,
    show_median=True,
    sharey=True,
    figure_size=None,
    thresholds=None,
):
    """One histogram panel per facet, stacked vertically, optional log y-axis."""
    sns.set_style("white")
    n_facet = max(1, len(facet_order))
    if figure_size is None:
        figure_size = (10, 3.2 * n_facet)
    if facet_labels is None:
        facet_labels = {k: k for k in facet_order}
    if xlabel is None:
        xlabel = x_col

    plot_bins = bins
    if xlim is not None and isinstance(bins, int):
        plot_bins = np.linspace(xlim[0], xlim[1], bins + 1)

    fig, axes = plt.subplots(
        n_facet, 1, figsize=figure_size, sharex=True, sharey=sharey
    )
    if n_facet == 1:
        axes = [axes]

    y_text = f"{ylabel} (Log Scale)" if log_scale else ylabel
    plot_df = data.copy() if data is not None else pd.DataFrame()

    for ax, facet in zip(axes, facet_order):
        sub = (
            plot_df[plot_df[facet_col] == facet]
            if (not plot_df.empty and facet_col in plot_df.columns)
            else pd.DataFrame()
        )
        color = "steelblue"
        if facet_colors:
            color = facet_colors.get(facet, color)
        n = 0
        if not sub.empty and x_col in sub.columns:
            x = pd.to_numeric(sub[x_col], errors="coerce").replace(
                [np.inf, -np.inf], np.nan
            ).dropna()
            n = int(len(x))
            if n:
                sns.histplot(
                    x=x,
                    bins=plot_bins,
                    color=color,
                    edgecolor="white",
                    linewidth=0.1,
                    ax=ax,
                )
                mean_val = float(x.mean())
                median_val = float(x.median())
                if show_mean and _finite_for_vline(mean_val):
                    ax.axvline(
                        x=mean_val,
                        color=TWO_CATEGORY_RED,
                        linestyle="--",
                        linewidth=1.5,
                        label=f"Mean: {mean_val:.4f}",
                    )
                if show_median and _finite_for_vline(median_val):
                    ax.axvline(
                        x=median_val,
                        color=TWO_CATEGORY_BLUE,
                        linestyle="-",
                        linewidth=1.5,
                        label=f"Median: {median_val:.4f}",
                    )
        if thresholds:
            for thr in thresholds:
                val = thr.get("value")
                if not _finite_for_vline(val):
                    continue
                ax.axvline(
                    x=float(val),
                    color=thr.get("color", "black"),
                    linestyle=thr.get("linestyle", "--"),
                    linewidth=thr.get("linewidth", 1.5),
                    label=thr.get("label", f"{float(val):.4f}"),
                )
        ax.set_title(facet_labels.get(facet, facet), fontsize=TITLE_FONT_SIZE)
        ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
        ax.set_ylabel(y_text, fontsize=Y_LABEL_FONT_SIZE)
        if log_scale:
            ax.set_yscale("log")
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        handles, _labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(
                loc="upper right",
                fontsize=LEGEND_FONT_SIZE,
                frameon=False,
            )

    axes[-1].set_xlabel(xlabel, fontsize=X_LABEL_FONT_SIZE)
    fig.suptitle(title, fontsize=TITLE_FONT_SIZE, y=1.02)
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved {filename}")
    plt.close(fig)


def plot_regression_comparison(
    df, 
    x_col, 
    y_col, 
    x_label, 
    y_label,
    filename,
    x_lim=None,
    y_lim=None,
    title=None,
    figure_size=(10,6)
):
    """
    Plots a scatter plot with OLS and Huber regression lines.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn.linear_model import HuberRegressor, LinearRegression
    
    sns.set_style("white")

    print(f"Running regressions for {y_label} vs {x_label}...")
    
    # Filter valid data
    local_valid = df[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    
    if len(local_valid) < 10:
        print(f"Not enough valid data for {filename}")
        return

    X = local_valid[x_col].values.reshape(-1, 1)
    y = local_valid[y_col].values
    
    # 1. OLS Regression
    ols = LinearRegression()
    ols.fit(X, y)
    ols_score = ols.score(X, y)
    ols_eq = f"y = {ols.coef_[0]:.4f}x + {ols.intercept_:.4f}"
    
    # 2. Huber Regression
    huber = HuberRegressor(epsilon=1.35)
    huber.fit(X, y)
    huber_score = huber.score(X, y) 
    huber_eq = f"y = {huber.coef_[0]:.4f}x + {huber.intercept_:.4f}"
    
    # Plotting
    plt.figure(figsize=figure_size)
    
    # Scatter plot
    sns.scatterplot(data=local_valid, x=x_col, y=y_col, 
                    alpha=0.3, s=15, color='#1f77b4', edgecolor='none', label='Samples')
    
    # Generate line points
    x_min, x_max = local_valid[x_col].min(), local_valid[x_col].max()
    x_span = x_max - x_min
    if x_span == 0: x_span = 1
    x_range = np.linspace(x_min - 0.05*x_span, x_max + 0.05*x_span, 100).reshape(-1, 1)
    
    y_ols = ols.predict(x_range)
    y_huber = huber.predict(x_range)
    
    # Plot lines
    plt.plot(x_range, y_ols, color=TWO_CATEGORY_BLUE, linewidth=2, linestyle='--', 
            label=f'OLS: {ols_eq}, $R^2$={ols_score:.3f}')
    plt.plot(x_range, y_huber, color=TWO_CATEGORY_RED, linewidth=2, 
            label=f'Huber: {huber_eq}, $R^2$={huber_score:.3f}')
    
    if title:
        plt.title(title, fontsize=TITLE_FONT_SIZE)
    else:
        plt.title(f'{y_label} vs {x_label}', fontsize=TITLE_FONT_SIZE)

    plt.xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    plt.tick_params(axis='both', which='major', labelsize=TICK_FONT_SIZE)
    if x_lim is not None:
        x_kwargs = {}
        if x_lim[0] is not None:
            x_kwargs["left"] = x_lim[0]
        if len(x_lim) > 1 and x_lim[1] is not None:
            x_kwargs["right"] = x_lim[1]
        if x_kwargs:
            plt.xlim(**x_kwargs)
    else:
        plt.xlim(x_min - 0.05*x_span, x_max + 0.05*x_span)
    if y_lim is not None:
        y_kwargs = {}
        if y_lim[0] is not None:
            y_kwargs["bottom"] = y_lim[0]
        if len(y_lim) > 1 and y_lim[1] is not None:
            y_kwargs["top"] = y_lim[1]
        if y_kwargs:
            plt.ylim(**y_kwargs)
    
    plt.legend(loc='upper left', bbox_to_anchor=(0.0, -0.15), fontsize=LEGEND_FONT_SIZE, frameon=False)
    
    # plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved plot to {filename}")
    plt.close()


def plot_stacked_distribution(
    df, 
    col, 
    group_col,
    title, 
    filename, 
    mean_val=None, 
    median_val=None,
    x_label=None, 
    y_label="Count", 
    log_scale=False,
    bins=100,
    figure_size=(10,6),
    xlim=None,
    auto_xlim=True,
    background_intervals: list[dict] | None = None,
):
    """
    Plots a stacked histogram distribution split by a group column.
    Includes Mean and Median vertical lines.
    Optional ``background_intervals`` shades the x-axis (drawn beneath bars).
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    sns.set_style("white")

    if x_label is None:
        x_label = col

    if df is None or getattr(df, "empty", True) or col not in df.columns:
        print(f"[Warning] Skipping stacked distribution (empty or missing column): {filename}")
        return
    plot_df = df[[c for c in (col, group_col) if c in df.columns]].replace(
        [np.inf, -np.inf], np.nan
    ).dropna(subset=[col])
    if plot_df.empty:
        print(f"[Warning] Skipping stacked distribution (no valid rows): {filename}")
        return
    if group_col in plot_df.columns and plot_df[group_col].dropna().empty:
        print(f"[Warning] Skipping stacked distribution (no group levels): {filename}")
        return

    fig, ax = plt.subplots(figsize=figure_size)
    bg_handles: list = []
    if background_intervals:
        bg_handles = draw_x_interval_backgrounds(ax, background_intervals, xlim=xlim, zorder=0)

    plot_bins = bins
    binrange = None
    if xlim is not None:
        binrange = xlim
        if isinstance(bins, int):
            plot_bins = np.linspace(xlim[0], xlim[1], bins + 1)

    # Histogram with Stacked Groups
    hue_kw = _hue_kwargs_for_sample_group(plot_df, group_col)
    ax = sns.histplot(
        data=plot_df,
        x=col,
        hue=group_col if group_col in plot_df.columns else None,
        multiple='stack',
        bins=plot_bins,
        binrange=binrange,
        linewidth=0.1,
        ax=ax,
        **hue_kw,
    )
    
    # Statistics Lines
    handles, labels = list(bg_handles), [h.get_label() for h in bg_handles]
    
    # Extract existing legend handles (Groups) if present
    if ax.get_legend():
        try:
            handles.extend(ax.get_legend().legend_handles)
            texts = ax.get_legend().get_texts()
            labels.extend(t.get_text() for t in texts)
            ax.get_legend().remove()
        except Exception:
            pass
    
    # Add Stats Lines
    if mean_val is not None:
        l1 = ax.axvline(x=mean_val, color=TWO_CATEGORY_RED, linestyle='--', linewidth=1, label=f'Mean: {mean_val:.4f}')
        handles.append(l1)
        labels.append(l1.get_label())
        
    if median_val is not None:
        l2 = ax.axvline(x=median_val, color=TWO_CATEGORY_BLUE, linestyle=':', linewidth=1.5, label=f'Median: {median_val:.4f}')
        handles.append(l2)
        labels.append(l2.get_label())
    
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.set_xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    ax.tick_params(axis='both', which='major', labelsize=TICK_FONT_SIZE)
    
    if log_scale:
        ax.set_yscale('log')
        ax.set_ylabel(f"{y_label} (Log Scale)", fontsize=Y_LABEL_FONT_SIZE)
    else:
        ax.set_ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)

    if xlim is not None:
        ax.set_xlim(xlim)
    elif auto_xlim and col in df.columns:
        lim = padded_axis_limits(df[col])
        if lim is not None:
            ax.set_xlim(lim)
        
    if handles:
        legend_title = 'CN class & Groups' if background_intervals else f'{group_col} & Stats'
        ax.legend(
            handles=handles,
            labels=labels,
            title=legend_title,
            loc='upper left',
            bbox_to_anchor=(0.0, -0.18),
            fontsize=LEGEND_FONT_SIZE,
            title_fontsize=LEGEND_FONT_SIZE,
            frameon=False,
            ncol=2 if background_intervals else 1,
        )
    
    fig.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved plot to {filename}")
    plt.close(fig)


def plot_stacked_hist_counts(
    bin_edges,
    counts_by_group: dict,
    title,
    filename,
    *,
    x_label,
    y_label="Count",
    hue_order=None,
    palette=None,
    hue_labels=None,
    mean_val=None,
    median_val=None,
    log_scale=False,
    figure_size=(10, 6),
):
    """Stacked histogram from pre-binned counts (no per-row seaborn scan)."""
    edges = np.asarray(bin_edges, dtype=float)
    if edges.size < 2:
        print(f"[Warning] Skipping stacked hist counts (need 2+ edges): {filename}")
        return
    widths = np.diff(edges)
    if np.any(widths <= 0):
        print(f"[Warning] Skipping stacked hist counts (non-increasing edges): {filename}")
        return

    order = list(hue_order) if hue_order is not None else list(counts_by_group.keys())
    n_bins = edges.size - 1
    series = []
    for key in order:
        raw = counts_by_group.get(key)
        if raw is None:
            continue
        arr = np.asarray(raw, dtype=float)
        if arr.size != n_bins:
            print(f"[Warning] Skipping group {key} (counts length {arr.size} != {n_bins})")
            continue
        series.append((str(key), arr))
    if not series:
        print(f"[Warning] Skipping stacked hist counts (no groups): {filename}")
        return

    sns.set_style("white")
    fig, ax = plt.subplots(figsize=figure_size)
    bottom = np.zeros(n_bins, dtype=float)
    handles = []
    labels = []
    for key, arr in series:
        color = None
        if isinstance(palette, dict):
            color = palette.get(key)
        bar = ax.bar(
            edges[:-1],
            arr,
            width=widths,
            align="edge",
            bottom=bottom,
            color=color,
            linewidth=0.1,
            edgecolor="white",
            label=key,
        )
        handles.append(bar)
        labels.append(
            hue_labels.get(key, key) if hue_labels is not None else key
        )
        bottom = bottom + arr

    if mean_val is not None:
        l1 = ax.axvline(
            x=mean_val,
            color=TWO_CATEGORY_RED,
            linestyle="--",
            linewidth=1,
            label=f"Mean: {mean_val:.4f}",
        )
        handles.append(l1)
        labels.append(l1.get_label())
    if median_val is not None:
        l2 = ax.axvline(
            x=median_val,
            color=TWO_CATEGORY_BLUE,
            linestyle=":",
            linewidth=1.5,
            label=f"Median: {median_val:.4f}",
        )
        handles.append(l2)
        labels.append(l2.get_label())

    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.set_xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if log_scale:
        ax.set_yscale("log")
        ax.set_ylabel(f"{y_label} (Log Scale)", fontsize=Y_LABEL_FONT_SIZE)
    else:
        ax.set_ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_xlim(float(edges[0]), float(edges[-1]))
    ax.legend(
        handles=handles,
        labels=labels,
        title="Groups & Stats" if (mean_val is not None or median_val is not None) else None,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.18),
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved plot to {filename}")
    plt.close(fig)


def hist_counts_by_label(
    values,
    labels,
    bin_edges,
    label_order,
) -> dict[str, np.ndarray]:
    """Bin ``values`` per label; clip into ``bin_edges`` so every finite value is counted."""
    edges = np.asarray(bin_edges, dtype=float)
    n_bins = int(edges.size - 1)
    order = [str(g) for g in label_order]
    out = {g: np.zeros(n_bins, dtype=np.int64) for g in order}
    x = pd.to_numeric(pd.Series(values, dtype="object"), errors="coerce").to_numpy(
        dtype=float, copy=False
    )
    g = pd.Series(labels, dtype="object").astype(str).to_numpy()
    ok = np.isfinite(x)
    if not ok.any():
        return out
    x = np.clip(x[ok], edges[0], edges[-1])
    g = g[ok]
    for label in order:
        mask = g == label
        if mask.any():
            out[label] = np.histogram(x[mask], bins=edges)[0].astype(np.int64)
    return out


def joint_hist_bin_edges(lim, n_bins: int, *, log_scale: bool = False) -> np.ndarray:
    """Bin edges spanning ``lim``; log-spaced when ``log_scale``."""
    lo, hi = float(lim[0]), float(lim[1])
    n = max(int(n_bins), 1)
    if log_scale:
        lo = max(lo, np.finfo(float).tiny)
        if hi <= lo:
            hi = lo * 10.0
        return np.logspace(np.log10(lo), np.log10(hi), n + 1)
    if hi <= lo:
        hi = lo + 1.0
    return np.linspace(lo, hi, n + 1)


def draw_stacked_hist_counts_on_ax(
    ax,
    bin_edges,
    counts_by_group: dict,
    *,
    orientation: str = "x",
    multiple: str = "stack",
    hue_order=None,
    palette=None,
) -> None:
    """Draw pre-binned stacked or filled bars on an existing axes (no seaborn scan)."""
    edges = np.asarray(bin_edges, dtype=float)
    if edges.size < 2:
        return
    widths = np.diff(edges)
    if np.any(widths <= 0):
        return
    order = list(hue_order) if hue_order is not None else list(counts_by_group.keys())
    n_bins = edges.size - 1
    series = []
    for key in order:
        raw = counts_by_group.get(key)
        if raw is None:
            raw = counts_by_group.get(str(key))
        if raw is None:
            continue
        arr = np.asarray(raw, dtype=float)
        if arr.size != n_bins:
            continue
        series.append((str(key), arr))
    if not series:
        return
    if multiple == "fill":
        tot = np.sum([arr for _, arr in series], axis=0)
        tot = np.where(tot == 0, 1.0, tot)
        series = [(key, arr / tot) for key, arr in series]
    left = np.zeros(n_bins, dtype=float)
    for key, arr in series:
        color = palette.get(key) if isinstance(palette, dict) else None
        if orientation == "y":
            ax.barh(
                edges[:-1],
                arr,
                height=widths,
                align="edge",
                left=left,
                color=color,
                linewidth=0.1,
                edgecolor="white",
            )
        else:
            ax.bar(
                edges[:-1],
                arr,
                width=widths,
                align="edge",
                bottom=left,
                color=color,
                linewidth=0.1,
                edgecolor="white",
            )
        left = left + arr


def plot_stacked_proportion_bar(
    prop_df: pd.DataFrame,
    title: str,
    filename: str,
    *,
    x_label: str | None = None,
    y_label: str = "Fraction of variants",
    stack_col_title: str = "Group",
    palette_map: dict[str, str] | None = None,
    annotate: bool = False,
    percent_axis: bool = False,
    figure_size=(10, 6),
) -> None:
    """
    100% stacked bar chart: rows in ``prop_df`` are x categories; columns are stack segments.

    Values should be non-negative and sum to ~1 along each row (within float tolerance).
    When ``stack_col_title`` is ``Group``, uses :func:`sample_group_hue_config` colors.
    ``palette_map`` overrides those colors. ``annotate`` writes each segment as a percent.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_style("white")
    if prop_df.empty:
        print(f"[Warning] plot_stacked_proportion_bar: empty data for {filename}")
        return

    work = prop_df.fillna(0.0).copy()
    row_sums = work.sum(axis=1)
    if not np.allclose(row_sums, 1.0, atol=1e-6, equal_nan=True):
        print(f"[Warning] Row sums not ~1 in stacked proportion bar; renormalizing for {filename}")
        work = work.div(row_sums.replace(0, np.nan), axis=0).fillna(0.0)

    x_labels = [str(x) for x in work.index.tolist()]
    stack_cols = list(work.columns)
    if palette_map is not None:
        colors = [palette_map.get(str(c), "#cccccc") for c in stack_cols]
    elif stack_col_title == "Group":
        _, group_palette = sample_group_hue_config(stack_cols)
        colors = [group_palette.get(str(c), SAMPLE_GROUP_UNKNOWN_COLOR) for c in stack_cols]
    else:
        colors = sns.color_palette("deep", len(stack_cols))

    fig, ax = plt.subplots(figsize=figure_size)
    bottom = np.zeros(len(work))
    x_pos = np.arange(len(work))

    for col, color in zip(stack_cols, colors):
        vals = work[col].to_numpy(dtype=float)
        ax.bar(x_pos, vals, bottom=bottom, label=str(col), color=color, width=0.65, edgecolor="white", linewidth=0.4)
        if annotate:
            for xpos, val, base in zip(x_pos, vals, bottom):
                if val <= 0:
                    continue
                text = f"{val * 100:.1f}%"
                if val >= 0.08:
                    ax.text(xpos, base + val / 2.0, text, ha="center", va="center", color="white", fontsize=TICK_FONT_SIZE)
                else:
                    ax.text(xpos + 0.38, base + val / 2.0, text, ha="left", va="center", color="#222222", fontsize=TICK_FONT_SIZE)
        bottom = bottom + vals

    ax.set_ylim(0.0, 1.05 if annotate else 1.0)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels)
    if percent_axis:
        ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    if x_label:
        ax.set_xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    ax.legend(
        title=stack_col_title,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.18),
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved plot to {filename}")
    plt.close(fig)


def plot_stacked_proportion_barh(
    prop_df: pd.DataFrame,
    title: str,
    filename: str,
    *,
    y_label: str | None = None,
    x_label: str = "Fraction of variants",
    stack_col_title: str = "Group",
    palette_map: dict[str, str] | None = None,
    figure_size=(8, 3.5),
    row_order: list | None = None,
) -> None:
    """100% horizontal stacked bars: rows in ``prop_df`` are y categories; columns are stack segments.

    ``row_order`` sets bottom-to-top bar order (matplotlib ``barh`` places the first row at the bottom).
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_style("white")
    if prop_df.empty:
        print(f"[Warning] plot_stacked_proportion_barh: empty data for {filename}")
        return

    work = prop_df.fillna(0.0).copy()
    row_sums = work.sum(axis=1)
    if not np.allclose(row_sums, 1.0, atol=1e-6, equal_nan=True):
        print(f"[Warning] Row sums not ~1 in horizontal stacked proportion bar; renormalizing for {filename}")
        work = work.div(row_sums.replace(0, np.nan), axis=0).fillna(0.0)

    if row_order is not None:
        work = work.reindex([str(r) for r in row_order])

    y_labels = [str(y) for y in work.index.tolist()]
    stack_cols = list(work.columns)
    if palette_map is not None:
        colors = [palette_map.get(str(c), "#999999") for c in stack_cols]
    elif stack_col_title == "Group":
        _, palette = sample_group_hue_config(stack_cols)
        colors = [palette.get(str(c), SAMPLE_GROUP_UNKNOWN_COLOR) for c in stack_cols]
    else:
        colors = sns.color_palette("deep", len(stack_cols))

    fig, ax = plt.subplots(figsize=figure_size)
    y_pos = np.arange(len(work))
    left = np.zeros(len(work))

    for col, color in zip(stack_cols, colors):
        vals = work[col].to_numpy(dtype=float)
        ax.barh(y_pos, vals, left=left, label=str(col), color=color, height=0.55, edgecolor="white", linewidth=0.4)
        left = left + vals

    ax.set_xlim(0.0, 1.0)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels)
    ax.set_xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    if y_label:
        ax.set_ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _pos: f"{x:.0%}"))
    ax.legend(
        title=stack_col_title,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.22),
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved plot to {filename}")
    plt.close(fig)


def plot_faceted_multi_line_series(
    data,
    x_col,
    facet_col,
    facet_order,
    y_specs,
    title,
    filename,
    xlabel=None,
    ylabel=None,
    xlim=None,
    ylim=None,
    facet_labels=None,
    figure_size=None,
    sharey=True,
):
    """One subplot per facet; each panel draws the same y-series specs vs ``x_col``."""
    sns.set_style("white")
    n_facet = max(1, len(facet_order))
    if figure_size is None:
        figure_size = (10, 3.2 * n_facet)
    if facet_labels is None:
        facet_labels = {k: k for k in facet_order}
    fig, axes = plt.subplots(n_facet, 1, figsize=figure_size, sharex=True, sharey=sharey)
    if n_facet == 1:
        axes = [axes]
    plot_df = data.copy()
    for ax, facet in zip(axes, facet_order):
        sub = plot_df[plot_df[facet_col] == facet]
        for spec in y_specs:
            yc = spec["y_col"]
            if sub.empty or yc not in sub.columns:
                continue
            ax.plot(
                sub[x_col],
                sub[yc],
                color=spec.get("color", "steelblue"),
                linestyle=spec.get("linestyle", "-"),
                linewidth=float(spec.get("linewidth", 1.5)),
                marker=spec.get("marker"),
                label=spec.get("label", yc),
            )
        ax.set_title(facet_labels.get(facet, facet), fontsize=TITLE_FONT_SIZE)
        ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        ax.grid(True, which="major", axis="y", alpha=0.3)
    axes[-1].set_xlabel(xlabel if xlabel is not None else x_col, fontsize=X_LABEL_FONT_SIZE)
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="upper left",
            bbox_to_anchor=(0.0, 0.02),
            fontsize=LEGEND_FONT_SIZE,
            frameon=False,
        )
    fig.suptitle(title, fontsize=TITLE_FONT_SIZE, y=1.02)
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved {filename}")
    plt.close(fig)


def plot_faceted_stacked_bars(
    data,
    x_col,
    facet_col,
    facet_order,
    stack_cols,
    title,
    filename,
    xlabel=None,
    ylabel="Proportion",
    stack_colors=None,
    stack_labels=None,
    facet_labels=None,
    bar_width=None,
    xlim=None,
    ylim=(0.0, 1.0),
    figure_size=None,
):
    """Numeric-x 100% stacked bars, one panel per facet."""
    sns.set_style("white")
    n_facet = max(1, len(facet_order))
    if figure_size is None:
        figure_size = (10, 3.2 * n_facet)
    if facet_labels is None:
        facet_labels = {k: k for k in facet_order}
    if stack_labels is None:
        stack_labels = {c: c for c in stack_cols}
    fig, axes = plt.subplots(n_facet, 1, figsize=figure_size, sharex=True, sharey=True)
    if n_facet == 1:
        axes = [axes]
    x_all = data[x_col].to_numpy(dtype=np.float64) if x_col in data.columns and not data.empty else np.array([])
    if bar_width is None and len(x_all) >= 2:
        diffs = np.diff(np.unique(x_all))
        bar_width = float(np.median(diffs)) if len(diffs) else 1.0
    if bar_width is None:
        bar_width = 1.0
    for ax, facet in zip(axes, facet_order):
        sub = data[data[facet_col] == facet].sort_values(x_col) if not data.empty else data
        x = sub[x_col].to_numpy(dtype=np.float64) if not sub.empty else np.array([])
        bottom = np.zeros(len(x), dtype=np.float64)
        for col in stack_cols:
            vals = (
                sub[col].to_numpy(dtype=np.float64)
                if (not sub.empty and col in sub.columns)
                else np.array([])
            )
            color = "#cccccc"
            if stack_colors:
                color = stack_colors.get(col, color)
            ax.bar(
                x,
                vals,
                bottom=bottom,
                width=float(bar_width),
                align="center",
                color=color,
                edgecolor="white",
                linewidth=0.2,
                label=stack_labels.get(col, col),
            )
            if len(vals):
                bottom = bottom + np.nan_to_num(vals, nan=0.0)
        ax.set_title(facet_labels.get(facet, facet), fontsize=TITLE_FONT_SIZE)
        ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
        ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
    axes[-1].set_xlabel(xlabel if xlabel is not None else x_col, fontsize=X_LABEL_FONT_SIZE)
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="upper left",
            bbox_to_anchor=(0.0, 0.02),
            fontsize=LEGEND_FONT_SIZE,
            frameon=False,
        )
    fig.suptitle(title, fontsize=TITLE_FONT_SIZE, y=1.02)
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved {filename}")
    plt.close(fig)


def plot_scatter_facet_grid(
    data,
    x_col,
    y_col,
    row_col,
    col_col,
    row_order,
    col_order,
    title,
    filename,
    xlabel=None,
    ylabel=None,
    xlim=None,
    ylim=None,
    identity_line=False,
    identity_label=r"$y=x$",
    size_col=None,
    color=None,
    connect=True,
    row_labels=None,
    col_labels=None,
    figure_size=None,
    sharex=True,
    sharey=True,
    row_xlim=None,
    row_ylim=None,
):
    """Scatter grid with optional identity line; one panel per row/column facet."""
    sns.set_style("white")
    n_row = max(1, len(row_order))
    n_col = max(1, len(col_order))
    if figure_size is None:
        figure_size = (4.2 * n_col, 3.4 * n_row)
    if row_labels is None:
        row_labels = {k: k for k in row_order}
    if col_labels is None:
        col_labels = {k: k for k in col_order}
    if color is None:
        color = TWO_CATEGORY_BLUE
    fig, axes = plt.subplots(
        n_row,
        n_col,
        figsize=figure_size,
        sharex=sharex,
        sharey=sharey,
        squeeze=False,
    )
    plot_df = data.copy() if data is not None else pd.DataFrame()
    for i, row_key in enumerate(row_order):
        for j, col_key in enumerate(col_order):
            ax = axes[i][j]
            if plot_df.empty or row_col not in plot_df.columns or col_col not in plot_df.columns:
                sub = plot_df.iloc[0:0]
            else:
                sub = plot_df[
                    (plot_df[row_col] == row_key) & (plot_df[col_col] == col_key)
                ]
            row_txt = row_labels.get(row_key, row_key)
            col_txt = col_labels.get(col_key, col_key)
            panel_title = str(row_txt) if not col_txt else f"{row_txt} | {col_txt}"
            panel_xlim = xlim
            panel_ylim = ylim
            if row_xlim:
                panel_xlim = row_xlim.get(row_key, xlim)
            if row_ylim:
                panel_ylim = row_ylim.get(row_key, ylim)
            plot_scatter_with_thresholds(
                data=sub if not sub.empty else pd.DataFrame({x_col: [], y_col: []}),
                x_col=x_col,
                y_col=y_col,
                title=panel_title,
                filename=filename,
                xlabel=xlabel if i == n_row - 1 else "",
                ylabel=ylabel if j == 0 else "",
                color=color,
                alpha=0.9,
                s=42,
                xlim=panel_xlim,
                ylim=panel_ylim,
                connect=connect and (not sub.empty),
                size_col=size_col if (size_col and not sub.empty) else None,
                identity_line=identity_line,
                identity_label=identity_label,
                ax=ax,
                save=False,
                show_legend=False,
                legend_outside=False,
            )
    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="upper left",
            bbox_to_anchor=(0.0, 0.02),
            fontsize=LEGEND_FONT_SIZE,
            frameon=False,
        )
    fig.suptitle(title, fontsize=TITLE_FONT_SIZE, y=1.02)
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved {filename}")
    plt.close(fig)


def plot_stacked_bar_chart(
    count_df: pd.DataFrame,
    title: str,
    filename: str,
    *,
    x_label: str | None = None,
    y_label: str = "Count",
    stack_col_title: str = "Category",
    stack_colors: dict[str, str] | None = None,
    figure_size=(10, 6),
    rotate_xlabels: int | None = 45,
) -> None:
    """
    Absolute-count stacked bar chart: rows in ``count_df`` are x categories; columns are stack segments.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_style("white")
    if count_df.empty:
        print(f"[Warning] plot_stacked_bar_chart: empty data for {filename}")
        return

    work = count_df.fillna(0).copy()
    stack_cols = list(work.columns)
    x_labels = [str(x) for x in work.index.tolist()]
    if stack_colors is None:
        colors = sns.color_palette("deep", len(stack_cols))
    else:
        colors = [stack_colors.get(str(col), "#cccccc") for col in stack_cols]

    fig, ax = plt.subplots(figsize=figure_size)
    bottom = np.zeros(len(work))
    x_pos = np.arange(len(work))

    for col, color in zip(stack_cols, colors):
        vals = work[col].to_numpy(dtype=float)
        ax.bar(
            x_pos,
            vals,
            bottom=bottom,
            label=str(col),
            color=color,
            width=0.75,
            edgecolor="white",
            linewidth=0.4,
        )
        bottom = bottom + vals

    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels)
    if x_label:
        ax.set_xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if rotate_xlabels is not None:
        ax.tick_params(axis="x", rotation=rotate_xlabels)
        for label in ax.get_xticklabels():
            label.set_ha("right")
    ax.legend(
        title=stack_col_title,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.18),
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved stacked bar chart to {filename}")
    plt.close(fig)


def opaque_legend_handles(handles):
    """Return legend proxies at alpha 1.0 without changing plotted artists."""
    from copy import copy

    from matplotlib.collections import PathCollection
    from matplotlib.lines import Line2D

    out = []
    for handle in handles:
        if isinstance(handle, PathCollection):
            facecolors = np.asarray(handle.get_facecolor())
            if facecolors.size == 0:
                out.append(handle)
                continue
            row = facecolors[0] if facecolors.ndim > 1 else facecolors
            rgb = tuple(float(c) for c in row[:3])
            sizes = handle.get_sizes()
            markersize = 8.0
            if sizes is not None and len(sizes):
                markersize = float(np.sqrt(sizes[0]))
            out.append(
                Line2D(
                    [0],
                    [0],
                    linestyle="none",
                    marker="o",
                    markersize=markersize,
                    markerfacecolor=rgb,
                    markeredgecolor="none",
                    alpha=1.0,
                )
            )
            continue
        proxy = copy(handle)
        if hasattr(proxy, "set_alpha"):
            proxy.set_alpha(1.0)
        out.append(proxy)
    return out


def plot_joint_regression(
    df, 
    x_col, 
    y_col, 
    group_col,
    x_label, 
    y_label, 
    filename,
    title=None,
    x_lim=None,
    y_lim=None,
    figure_size=(10,6),
    *,
    palette=None,
    hue_order=None,
    log_x=False,
    log_y=False,
    draw_regression=True,
    alpha=0.6,
    s=20,
    rasterized=False,
    legend_title=None,
    bins_x=100,
    bins_y=50,
    min_points=10,
    marginal_multiple="stack",
    stack_last_on_top=False,
    hue_alpha=None,
    also_stack_counts=False,
    hue_labels=None,
):
    """
    Plots a JointGrid (Scatter + stacked marginal histograms), optionally with
    OLS and Huber regression lines. Supports grouping/hue in scatter and histograms.

    Seaborn ``histplot`` draws hue levels in reverse, so the first ``hue_order``
    level lands on top of a stack. Set ``stack_last_on_top=True`` to put the last
    ``hue_order`` level on top (same z-order as the scatter). ``marginal_multiple``
    is passed to ``histplot`` (``stack`` or ``fill``). When ``hue_order`` and a
    dict ``palette`` are set, scatter is drawn one hue at a time (first bottom,
    last top) so rasterized PNG output keeps painter order. ``hue_alpha`` maps
    hue label to point alpha; otherwise ``alpha`` is used for every hue.
    ``also_stack_counts`` adds an outer pair of filled (0–1) marginals for
    class share. The inner strips (next to the scatter) show stacked counts.
    Outer boxes use 0–1 ticks with no Count/Frac axis title. Marginal
    histograms are pre-binned with numpy (same approach as
    ``plot_stacked_hist_counts``); they do not scan rows with seaborn.
    Legend markers use the palette colour at alpha 1.0 even when scatter
    points are translucent.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn.linear_model import HuberRegressor, LinearRegression
    
    sns.set_style("white")

    print(f"Running joint plot for {y_label} vs {x_label}...")
    
    # Filter valid data
    local_valid = df[[x_col, y_col, group_col]].replace([np.inf, -np.inf], np.nan).dropna()
    
    if len(local_valid) < int(min_points):
        print(f"Not enough valid data for {filename}")
        return

    # Plotting: JointGrid
    g = sns.JointGrid(data=local_valid, x=x_col, y=y_col, height=figure_size[1], ratio=5)
    hue_kw = _hue_kwargs_for_sample_group(local_valid, group_col)
    if palette is not None:
        hue_kw["palette"] = palette
    if hue_order is not None:
        hue_kw["hue_order"] = hue_order

    # 1. Main Scatter Plot
    # sns.scatterplot(hue=...) can flatten z-order when rasterized; draw each
    # hue_order level as its own collection so the last level stays on top.
    palette_map = palette if isinstance(palette, dict) else hue_kw.get("palette")
    scatter_order = list(hue_order) if hue_order is not None else None
    if scatter_order is not None and isinstance(palette_map, dict):
        for i, level in enumerate(scatter_order):
            sub = local_valid[local_valid[group_col] == level]
            if sub.empty:
                continue
            level_alpha = alpha
            if hue_alpha is not None:
                level_alpha = hue_alpha.get(level, alpha)
            g.ax_joint.scatter(
                sub[x_col],
                sub[y_col],
                s=s,
                alpha=level_alpha,
                c=palette_map[level],
                edgecolors="none",
                rasterized=rasterized,
                zorder=3 + i,
                label=(
                    hue_labels.get(level, level) if hue_labels is not None else level
                ),
            )
    else:
        sns.scatterplot(
            data=local_valid,
            x=x_col,
            y=y_col,
            hue=group_col,
            alpha=alpha,
            s=s,
            edgecolor="none",
            ax=g.ax_joint,
            legend="full",
            rasterized=rasterized,
            **hue_kw,
        )

    if draw_regression:
        X = local_valid[x_col].values.reshape(-1, 1)
        y = local_valid[y_col].values
        ols = LinearRegression()
        ols.fit(X, y)
        ols_score = ols.score(X, y)
        ols_eq = f"y = {ols.coef_[0]:.4f}x + {ols.intercept_:.4f}"
        huber = HuberRegressor(epsilon=1.35)
        huber.fit(X, y)
        huber_score = huber.score(X, y)
        huber_eq = f"y = {huber.coef_[0]:.4f}x + {huber.intercept_:.4f}"
        x_min, x_max = local_valid[x_col].min(), local_valid[x_col].max()
        x_span = x_max - x_min
        if x_span == 0:
            x_span = 0.01
        x_range = np.linspace(x_min - 0.05 * x_span, x_max + 0.05 * x_span, 100).reshape(-1, 1)
        g.ax_joint.plot(
            x_range,
            ols.predict(x_range),
            color=TWO_CATEGORY_BLUE,
            linewidth=2,
            linestyle="--",
            label=f"OLS: {ols_eq}\n$R^2$={ols_score:.3f}",
        )
        g.ax_joint.plot(
            x_range,
            huber.predict(x_range),
            color=TWO_CATEGORY_RED,
            linewidth=2,
            label=f"Huber: {huber_eq}\n$R^2$={huber_score:.3f}",
        )
    
    # --- Legend Management ---
    
    # Remove default legend inside (if any)
    if g.ax_joint.get_legend():
        g.ax_joint.get_legend().remove()

    # Retrieve all handles (Groups + Regression Lines); legend uses opaque palette colours.
    final_handles, final_labels = g.ax_joint.get_legend_handles_labels()
    final_handles = opaque_legend_handles(final_handles)
    if legend_title is None:
        legend_title = f"{group_col} & Regression" if draw_regression else None
    g.ax_joint.legend(
        handles=final_handles,
        labels=final_labels,
        title=legend_title,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.15),
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )
    
    # 3. Axis limits first, then pre-binned marginals (numpy histogram + bar).
    g.ax_joint.set_xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    g.ax_joint.set_ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    g.ax_joint.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if log_x:
        g.ax_joint.set_xscale("log")
    if log_y:
        g.ax_joint.set_yscale("log")

    if x_lim is not None:
        g.ax_joint.set_xlim(x_lim)
    elif not log_x:
        auto_x = padded_axis_limits(local_valid[x_col])
        if auto_x is not None:
            g.ax_joint.set_xlim(auto_x)
    if y_lim is not None:
        g.ax_joint.set_ylim(y_lim)
    elif not log_y:
        auto_y = padded_axis_limits(local_valid[y_col])
        if auto_y is not None:
            g.ax_joint.set_ylim(auto_y)

    g.ax_marg_x.set_xlim(g.ax_joint.get_xlim())
    g.ax_marg_y.set_ylim(g.ax_joint.get_ylim())

    draw_order = list(hue_order) if hue_order is not None else sorted(
        local_valid[group_col].dropna().astype(str).unique()
    )
    if not stack_last_on_top:
        draw_order = list(reversed(draw_order))
    palette_for_hist = palette if isinstance(palette, dict) else hue_kw.get("palette")
    x_edges = joint_hist_bin_edges(g.ax_joint.get_xlim(), bins_x, log_scale=log_x)
    y_edges = joint_hist_bin_edges(g.ax_joint.get_ylim(), bins_y, log_scale=log_y)
    x_counts = hist_counts_by_label(
        local_valid[x_col], local_valid[group_col], x_edges, draw_order
    )
    y_counts = hist_counts_by_label(
        local_valid[y_col], local_valid[group_col], y_edges, draw_order
    )
    inner_multiple = "stack" if also_stack_counts else marginal_multiple
    draw_stacked_hist_counts_on_ax(
        g.ax_marg_x,
        x_edges,
        x_counts,
        orientation="x",
        multiple=inner_multiple,
        hue_order=draw_order,
        palette=palette_for_hist,
    )
    draw_stacked_hist_counts_on_ax(
        g.ax_marg_y,
        y_edges,
        y_counts,
        orientation="y",
        multiple=inner_multiple,
        hue_order=draw_order,
        palette=palette_for_hist,
    )

    if also_stack_counts:
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        div_x = make_axes_locatable(g.ax_marg_x)
        ax_frac_x = div_x.append_axes("top", size="100%", pad=0.12)
        draw_stacked_hist_counts_on_ax(
            ax_frac_x,
            x_edges,
            x_counts,
            orientation="x",
            multiple="fill",
            hue_order=draw_order,
            palette=palette_for_hist,
        )
        ax_frac_x.set_xlim(g.ax_joint.get_xlim())
        ax_frac_x.set_ylim(0.0, 1.0)
        ax_frac_x.set_yticks([0, 1])
        ax_frac_x.set_yticklabels(["0", "1"])
        if log_x:
            ax_frac_x.set_xscale("log")
        ax_frac_x.tick_params(axis="x", labelbottom=False)
        ax_frac_x.set_xlabel("")
        ax_frac_x.set_ylabel("")
        ax_frac_x.tick_params(axis="y", labelsize=TICK_FONT_SIZE)
        div_y = make_axes_locatable(g.ax_marg_y)
        ax_frac_y = div_y.append_axes("right", size="100%", pad=0.12)
        draw_stacked_hist_counts_on_ax(
            ax_frac_y,
            y_edges,
            y_counts,
            orientation="y",
            multiple="fill",
            hue_order=draw_order,
            palette=palette_for_hist,
        )
        ax_frac_y.set_ylim(g.ax_joint.get_ylim())
        ax_frac_y.set_xlim(0.0, 1.0)
        ax_frac_y.set_xticks([0, 1])
        ax_frac_y.set_xticklabels(["0", "1"])
        if log_y:
            ax_frac_y.set_yscale("log")
        ax_frac_y.tick_params(axis="y", labelleft=False)
        ax_frac_y.set_ylabel("")
        ax_frac_y.set_xlabel("")
        ax_frac_y.tick_params(axis="x", labelsize=TICK_FONT_SIZE)

    if title:
        plt.suptitle(title, y=0.98, fontsize=TITLE_FONT_SIZE)
    else:
        plt.suptitle(f'{y_label} vs {x_label}', y=0.98, fontsize=TITLE_FONT_SIZE)
    
    g.fig.subplots_adjust(top=0.9)
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved plot to {filename}")
    plt.close()


def plot_mean_variance_fit(
    x_data,
    y_data,
    x_line,
    y_line_poisson,
    y_line_nb,
    phi_est,
    filename,
    title='Mean-Variance Relationship',
    x_label='Mean Depth (log)',
    y_label='Depth Variance (log)',
    sample_size=50000,
    figure_size=(10,6)
):
    """
    Plots Mean vs Variance with fitted Poisson and Negative Binomial lines.
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns
    sns.set_style("white")

    plt.figure(figsize=figure_size)
    
    # Downsample for scatter
    if len(x_data) > sample_size:
        # Combine to sample rows
        temp_df = pd.DataFrame({'Mean': x_data, 'Var': y_data})
        temp_df = temp_df.sample(n=sample_size, random_state=42)
        plt.scatter(temp_df['Mean'], temp_df['Var'], 
                    alpha=0.2, s=10, color='gray', label='Observed (Sampled)')
    else:
        plt.scatter(x_data, y_data, 
                    alpha=0.2, s=10, color='gray', label='Observed')
    
    # Lines
    plt.plot(x_line, y_line_poisson, 'b--', linewidth=2, label='Poisson')
    plt.plot(x_line, y_line_nb, 'r-', linewidth=2, label=f'NB Fit (phi={phi_est:.4f})')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.tick_params(axis='both', which='major', labelsize=TICK_FONT_SIZE)
    
    plt.legend(loc='upper left', bbox_to_anchor=(0.0, -0.15), fontsize=LEGEND_FONT_SIZE, frameon=False)
    # plt.grid(True, alpha=0.3)
    
    # plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved fit plot: {filename}")
    plt.close()

def plot_qq_residuals(
    residuals,
    filename,
    title="QQ Plot of Standardized Residuals",
    figure_size=(8,6)
):
    """
    Plots a QQ plot of residuals using statsmodels.
    """
    import matplotlib.pyplot as plt
    import statsmodels.api as sm
    import seaborn as sns
    sns.set_style("white")

    plt.figure(figsize=figure_size)
    sm.qqplot(residuals, line='45', fit=True)
    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.xlabel("Theoretical Quantiles", fontsize=X_LABEL_FONT_SIZE) # sm.qqplot usually adds labels, but we can override or let it be
    plt.ylabel("Sample Quantiles", fontsize=Y_LABEL_FONT_SIZE)
    plt.tick_params(axis='both', which='major', labelsize=TICK_FONT_SIZE)
    
    # plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved QQ plot: {filename}")
    plt.close()


def plot_bar_chart(
    names,
    values,
    title,
    ylabel,
    filename,
    ylim=(0.0, 1.05),
    color='steelblue',
    figure_size=(10, 5),
    dpi=300,
    rotate_xlabels=None,
):
    """
    Plots a simple bar chart with value labels on top.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os
    
    sns.set_style("white")

    # Create directory if needed (defensive)
    folder = os.path.dirname(filename)
    if folder:
        os.makedirs(folder, exist_ok=True)

    fig, ax = plt.subplots(figsize=figure_size)
    ax.bar(names, values, color=color, alpha=0.8)
    ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.tick_params(axis='both', which='major', labelsize=TICK_FONT_SIZE)
    if ylim:
        ax.set_ylim(ylim)

    ymax = max(values) if values else 0.0
    label_pad = 0.01 if ymax <= 1.0 else ymax * 0.02

    # Add text labels
    for i, v in enumerate(values):
        if ymax <= 1.0:
            lbl = f"{v:.3f}"
        else:
            lbl = f"{v:.0f}" if v >= 10 else f"{v:.2f}"
        ax.text(i, v + label_pad, lbl, ha='center', va='bottom', fontsize=TICK_FONT_SIZE)

    if rotate_xlabels is not None:
        ax.tick_params(axis='x', rotation=rotate_xlabels)
        for label in ax.get_xticklabels():
            label.set_ha('right')

    # fig.tight_layout()
    fig.savefig(filename, dpi=dpi, bbox_inches='tight')
    print(f"Saved bar chart: {filename}")
    plt.close(fig)


def plot_categorical_bar_chart(
    categories: list[str],
    values: list[float],
    title: str,
    ylabel: str,
    filename: str,
    *,
    colors: list | None = None,
    ylim: tuple[float, float] | None = None,
    value_formatter=None,
    figure_size=(10, 5),
    rotate_xlabels: int | None = 45,
):
    """Bar chart with per-category colors and value labels (counts or percentages)."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_style("white")
    fig, ax = plt.subplots(figsize=figure_size)
    bar_colors = colors if colors is not None else ["steelblue"] * len(categories)
    ax.bar(categories, values, color=bar_colors, alpha=0.85, edgecolor="white", linewidth=0.5)
    ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if ylim is not None:
        ax.set_ylim(ylim)
    ymax = max(values) if values else 0.0
    label_pad = ymax * 0.02 if ymax > 1.0 else 0.02
    fmt = value_formatter or (lambda v: f"{v:.0f}" if v >= 10 else f"{v:.1f}")
    for i, v in enumerate(values):
        ax.text(i, v + label_pad, fmt(v), ha="center", va="bottom", fontsize=TICK_FONT_SIZE)
    if rotate_xlabels is not None:
        ax.tick_params(axis="x", rotation=rotate_xlabels)
        for label in ax.get_xticklabels():
            label.set_ha("right")
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved bar chart: {filename}")
    plt.close(fig)


def plot_gwas_qq(
    expected,
    observed,
    title,
    filename,
    xlabel="Expected -log10(P)",
    ylabel="Observed -log10(P)",
    figure_size=(6, 6)
):
    """
    Plots a GWAS QQ plot (Expected vs Observed -log10 P-values).
    """
    import matplotlib.pyplot as plt
    import os
    import seaborn as sns
    sns.set_style("white")
    
     # Create directory if needed
    folder = os.path.dirname(filename)
    if folder:
        os.makedirs(folder, exist_ok=True)

    fig, ax = plt.subplots(figsize=figure_size)
    ax.scatter(expected, observed, s=6, alpha=0.6)
    
    # Diagonal line
    # Determine limit based on max value in both arrays
    lim = max(expected.max(), observed.max())
    ax.plot([0, lim], [0, lim], color='red', lw=1)
    
    ax.set_xlabel(xlabel, fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.tick_params(axis='both', which='major', labelsize=TICK_FONT_SIZE)
    
    # fig.tight_layout()
    fig.savefig(filename, dpi=160, bbox_inches='tight')
    print(f"Saved GWAS QQ plot: {filename}")
    plt.close(fig)


def plot_manhattan(
    chromosomes,
    positions,
    values,
    title,
    filename,
    xlabel="Chromosome",
    ylabel="-log10(P)",
    thresholds=None,
    figure_size=(14, 6),
    band_colors=("#4c72b0", "#95a5a6"),
):
    """
    Plots values against cumulative genomic position, banded by chromosome.

    Chromosome offsets come from the observed maximum position per chromosome,
    so this works for any chromosome coding (segment ids or linear ids) without
    an external length table. ``thresholds`` is a mapping of label to y value.
    """
    import matplotlib.pyplot as plt
    import os
    import numpy as np
    import pandas as pd
    import seaborn as sns
    sns.set_style("white")

    folder = os.path.dirname(filename)
    if folder:
        os.makedirs(folder, exist_ok=True)

    frame = pd.DataFrame({
        "chrom": np.asarray(chromosomes),
        "pos": pd.to_numeric(pd.Series(positions), errors="coerce").to_numpy(),
        "value": pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(),
    }).dropna()
    if frame.empty:
        print(f"[Warning] No finite points to plot: {filename}")
        return

    order = sorted(frame["chrom"].unique(), key=lambda c: (len(str(c)), str(c)))
    offset = 0.0
    ticks, tick_labels = [], []
    fig, ax = plt.subplots(figsize=figure_size)
    for index, chrom in enumerate(order):
        block = frame[frame["chrom"] == chrom]
        span = float(block["pos"].max())
        ax.scatter(
            block["pos"] + offset,
            block["value"],
            s=3,
            alpha=0.7,
            color=band_colors[index % len(band_colors)],
            linewidths=0,
        )
        ticks.append(offset + span / 2.0)
        tick_labels.append(str(chrom))
        offset += span

    for label, level in (thresholds or {}).items():
        ax.axhline(y=level, color="red", linestyle="--", lw=1, label=label)
    if thresholds:
        ax.legend(fontsize=LEGEND_FONT_SIZE, loc="upper right")

    ax.set_xlim(0, offset)
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, fontsize=TICK_FONT_SIZE)
    ax.set_xlabel(xlabel, fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.tick_params(axis="y", which="major", labelsize=TICK_FONT_SIZE)

    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved Manhattan plot: {filename}")
    plt.close(fig)


def plot_correlation_with_regression(
    data, x_col, y_col, 
    title, filename, 
    x_label=None, y_label=None,
    color='gray', line_color='red',
    figure_size=(10,6),
    xlim=None, ylim=None,
    group_col=None,
):
    """
    Plots a scatter plot with linear regression line and statistics (slope, intercept, R^2).

    When ``group_col`` is set (e.g. ``Group``), points are colored by that column using
    the canonical sample_group palette; the regression line is fit on all points.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    
    sns.set_style("white")

    if x_label is None: x_label = x_col
    if y_label is None: y_label = y_col

    use_hue = group_col and group_col in data.columns
    cols = [x_col, y_col] + ([group_col] if use_hue else [])
    clean_data = data[cols].replace([np.inf, -np.inf], np.nan).dropna()
    x = clean_data[x_col]
    y = clean_data[y_col]

    if len(x) < 2:
        print(f"[Warning] Not enough data points to plot regression for {x_col} vs {y_col}")
        return

    slope, intercept = np.polyfit(x, y, 1)
    r_squared = np.corrcoef(x, y)[0, 1] ** 2

    plt.figure(figsize=figure_size)
    if use_hue:
        hue_kw = _hue_kwargs_for_sample_group(clean_data, group_col)
        sns.scatterplot(
            data=clean_data, x=x_col, y=y_col, hue=group_col,
            alpha=0.45, s=12, edgecolor='w', **hue_kw,
        )
        sns.regplot(
            data=clean_data, x=x_col, y=y_col, scatter=False,
            line_kws={'color': line_color, 'label': 'Linear Regression'},
        )
    else:
        sns.regplot(
            x=x_col, y=y_col, data=clean_data,
            scatter_kws={'alpha':0.4, 's':10, 'color': color, 'label': 'Data Points'},
            line_kws={'color': line_color, 'label': 'Linear Regression'}
        )

    stats_text = f'$y = {slope:.4f}x + {intercept:.4f}$\n$R^2 = {r_squared:.4f}$'
    plt.text(
        0.05, 0.9, stats_text,
        transform=plt.gca().transAxes, 
        fontsize=LEGEND_FONT_SIZE, color='darkred',
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray', boxstyle='round')
    )

    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    plt.tick_params(axis='both', which='major', labelsize=TICK_FONT_SIZE)
    if ylim is not None:
        plt.ylim(*ylim)
    if xlim is not None:
        plt.xlim(*xlim)
    
    plt.legend(
        loc='upper left', bbox_to_anchor=(0.0, -0.15),
        fontsize=LEGEND_FONT_SIZE, frameon=False,
        title=group_col if use_hue else None,
    )
    
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {filename}")
    plt.close()


def plot_dual_regression(
    df_valid, df_plot, x_col, y_col, x_label, y_label, filename,
    xlim=None, ylim=None
):
    """
    Plots a scatter plot with OLS and Huber regression lines.
    """
    print(f"Running regressions for {y_label} vs {x_label}...")
    
    # Filter for valid data
    local_valid = df_valid[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    
    if len(local_valid) < 10:
        print(f"Not enough valid data for {filename}")
        return

    X = local_valid[x_col].values.reshape(-1, 1)
    y = local_valid[y_col].values
    
    # 1. OLS Regression
    ols = LinearRegression()
    ols.fit(X, y)
    ols_slope = ols.coef_[0]
    ols_intercept = ols.intercept_
    ols_score = ols.score(X, y)
    ols_eq = f"y = {ols_slope:.4f}x + {ols_intercept:.4f}"
    
    # 2. Huber Regression
    huber = HuberRegressor(epsilon=1.35)
    huber.fit(X, y)
    huber_slope = huber.coef_[0]
    huber_intercept = huber.intercept_
    huber_score = huber.score(X, y) 
    huber_eq = f"y = {huber_slope:.4f}x + {huber_intercept:.4f}"
    
    # Plotting
    plt.figure(figsize=(10, 8))
    
    # Calculate limits if not provided
    if xlim is None:
        x_max = local_valid[x_col].max()
        if x_max > 0.5: x_max = 0.55 # specific adjustment from original script
        xlim = (0, x_max)
    
    if ylim is None:
        y_max = local_valid[y_col].max() * 1.1
        ylim = (0, y_max)

    # Scatter plot (downsampled)
    local_plot = df_plot[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    sns.scatterplot(data=local_plot, x=x_col, y=y_col, 
                    alpha=0.1, s=10, color='gray', edgecolor=None, label='Data (Sampled)')
    
    # Generate line points
    x_range = np.linspace(xlim[0], xlim[1], 100).reshape(-1, 1)
    y_ols = ols.predict(x_range)
    y_huber = huber.predict(x_range)
    
    # OLS Line
    plt.plot(x_range, y_ols, color=TWO_CATEGORY_BLUE, linewidth=2, linestyle='--', 
             label=f'OLS: {ols_eq}, $R^2$={ols_score:.3f}')
    # Huber Line
    plt.plot(x_range, y_huber, color=TWO_CATEGORY_RED, linewidth=2, 
             label=f'Huber: {huber_eq}, $R^2$={huber_score:.3f}')
    
    if "Heterozygosity" in y_label and "Frequency" in x_label:
         y_hwe = 2 * x_range * (1 - x_range)
         plt.plot(x_range, y_hwe, color="#7f7f7f", linewidth=2, linestyle=":",
             label='Expected (HWE): $2p(1-p)$')

    plt.title(f'{y_label} vs {x_label}')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    plt.close()


def plot_heatmap_custom(
    data_matrix, x_labels, y_labels,
    title, filename,
    xlabel=None, ylabel=None,
    cmap='viridis', cbar_label='',
    x_tick_step=None,
    vmin=None,
    vmax=None,
    gamma=None,
    overlays=None,
    vlines=None,
):
    """
    Plots a heatmap from a matrix.

    ``gamma`` < 1 (PowerNorm) stretches low values so small fractions stay
    visible. ``overlays`` / ``vlines`` are drawn in heatmap coordinates
    (cell *i* spans ``[i, i+1]``). Each overlay/vline dict may include
    ``x``/``y`` (overlays), ``x`` (vlines), plus matplotlib line kwargs.
    """
    sns.set_style("white")
    plt.figure(figsize=(14, 10))
    heatmap_kw = dict(cmap=cmap, cbar_kws={'label': cbar_label})
    if gamma is not None:
        lo = 0.0 if vmin is None else float(vmin)
        hi = 1.0 if vmax is None else float(vmax)
        heatmap_kw['norm'] = mcolors.PowerNorm(
            gamma=float(gamma), vmin=lo, vmax=hi, clip=True
        )
    else:
        if vmin is not None:
            heatmap_kw['vmin'] = vmin
        if vmax is not None:
            heatmap_kw['vmax'] = vmax
    ax = sns.heatmap(data_matrix, **heatmap_kw)

    if xlabel:
        ax.set_xlabel(xlabel, fontsize=HEATMAP_AXIS_FONT_SIZE)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=HEATMAP_AXIS_FONT_SIZE)
    ax.set_title(title, fontsize=HEATMAP_TITLE_FONT_SIZE, pad=12)

    n_x = len(x_labels)
    step_x = x_tick_step if x_tick_step is not None else max(1, n_x // 10)
    xticks = np.arange(0, n_x, step_x)
    ax.set_xticks(xticks + 0.5)
    ax.set_xticklabels([x_labels[i] for i in xticks], rotation=0, fontsize=HEATMAP_TICK_FONT_SIZE)

    n_y = len(y_labels)
    step_y = max(1, n_y // 10)
    yticks = np.arange(0, n_y, step_y)
    ax.set_yticks(yticks + 0.5)
    ax.set_yticklabels([y_labels[i] for i in yticks], rotation=0, fontsize=HEATMAP_TICK_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=HEATMAP_TICK_FONT_SIZE)

    cbar = ax.collections[0].colorbar
    if cbar is not None:
        cbar.set_label(cbar_label, fontsize=HEATMAP_AXIS_FONT_SIZE)
        cbar.ax.tick_params(labelsize=HEATMAP_TICK_FONT_SIZE)

    legend_handles = False
    if overlays:
        for spec in overlays:
            plot_kw = {k: v for k, v in spec.items() if k not in ("x", "y")}
            plot_kw.setdefault("zorder", 10)
            ax.plot(spec["x"], spec["y"], **plot_kw)
            if spec.get("label"):
                legend_handles = True
    if vlines:
        for spec in vlines:
            line_kw = {k: v for k, v in spec.items() if k != "x"}
            line_kw.setdefault("zorder", 10)
            ax.axvline(spec["x"], **line_kw)
            if spec.get("label"):
                legend_handles = True
    if legend_handles:
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(0.0, -0.15),
            fontsize=HEATMAP_LEGEND_FONT_SIZE,
            frameon=False,
        )

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Heatmap saved to: {filename}")
    plt.close()


def plot_binned_mean_heatmap(
    df,
    x_col,
    y_col,
    value_col,
    x_label,
    y_label,
    cbar_label,
    title,
    filename,
    n_x_bins=50,
    n_y_bins=50,
    x_edges=None,
    y_edges=None,
    cmap=MISSING_RATE_CMAP,
    vmin=None,
    vmax=None,
    figure_size=(10, 8),
    x_tick_step=None,
    y_tick_step=None,
    x_label_fmt="{:.2f}",
    y_label_fmt="{:.0f}",
    y_log10=False,
):
    """
    2D heatmap of mean ``value_col`` over ``(x_col, y_col)`` bins.

    Rows are ``y_col`` (high values at the top). ``cmap`` defaults to ``MISSING_RATE_CMAP``
    (reversed coolwarm): red = low values, blue = high values (lower missing rate is better).

    When ``y_log10`` is True, samples with ``y_col`` <= 0 are dropped and binning
    uses log10(y) with padded edges in log space.
    """
    sns.set_style("white")

    work = df[[x_col, y_col, value_col]].replace([np.inf, -np.inf], np.nan).dropna()
    if y_log10:
        work = work[work[y_col] > 0].copy()
        work["_y_plot"] = np.log10(work[y_col].astype(float))
        y_bin_source = "_y_plot"
    else:
        work = work.copy()
        y_bin_source = y_col

    if len(work) < 10:
        print(f"Not enough valid data for binned heatmap: {filename}")
        return

    x_lo, x_hi = float(work[x_col].min()), float(work[x_col].max())
    y_lo, y_hi = float(work[y_bin_source].min()), float(work[y_bin_source].max())
    x_span = x_hi - x_lo or max(abs(x_lo), 1.0) * 0.01
    y_span = y_hi - y_lo or max(abs(y_lo), 1.0) * 0.01
    x_pad = x_span * DEFAULT_AXIS_PADDING_FRACTION
    y_pad = y_span * DEFAULT_AXIS_PADDING_FRACTION

    if x_edges is None:
        x_edges = np.linspace(x_lo - x_pad, x_hi + x_pad, n_x_bins + 1)
    if y_edges is None:
        y_edges = np.linspace(y_lo - y_pad, y_hi + y_pad, n_y_bins + 1)

    binned = work.copy()
    binned["_x_bin"] = pd.cut(binned[x_col], bins=x_edges, include_lowest=True)
    binned["_y_bin"] = pd.cut(binned[y_bin_source], bins=y_edges, include_lowest=True)
    matrix_df = binned.pivot_table(
        index="_y_bin", columns="_x_bin", values=value_col, aggfunc="mean", observed=True
    )
    if matrix_df.empty:
        print(f"No binned cells for heatmap: {filename}")
        return

    matrix_df = matrix_df.sort_index(ascending=False)
    matrix = matrix_df.to_numpy(dtype=float)
    x_labels = [x_label_fmt.format(interval.mid) for interval in matrix_df.columns]
    y_labels = [y_label_fmt.format(interval.mid) for interval in matrix_df.index]

    if vmin is None:
        vmin = float(np.nanmin(matrix))
    if vmax is None:
        vmax = float(np.nanmax(matrix))
    if vmin == vmax:
        vmax = vmin + 1e-6

    plt.figure(figsize=figure_size)
    ax = sns.heatmap(
        matrix,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        mask=np.isnan(matrix),
        cbar_kws={"label": cbar_label},
        linewidths=0,
    )
    ax.set_xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)

    n_x = len(x_labels)
    step_x = x_tick_step if x_tick_step is not None else max(1, n_x // 10)
    xticks = np.arange(0, n_x, step_x)
    ax.set_xticks(xticks + 0.5)
    ax.set_xticklabels([x_labels[i] for i in xticks], rotation=0, fontsize=TICK_FONT_SIZE)

    n_y = len(y_labels)
    step_y = y_tick_step if y_tick_step is not None else max(1, n_y // 10)
    yticks = np.arange(0, n_y, step_y)
    ax.set_yticks(yticks + 0.5)
    ax.set_yticklabels([y_labels[i] for i in yticks], rotation=0, fontsize=TICK_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved binned mean heatmap to {filename}")
    plt.close()


def _fit_loess_curve(x, y, frac=0.2, n_grid=200):
    """Return sorted (x_line, y_line) from statsmodels LOWESS."""
    from statsmodels.nonparametric.smoothers_lowess import lowess

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 10:
        return None, None
    order = np.argsort(x)
    x, y = x[order], y[order]
    smoothed = lowess(y, x, frac=frac, return_sorted=True)
    x_line = np.linspace(float(x.min()), float(x.max()), n_grid)
    y_line = np.interp(x_line, smoothed[:, 0], smoothed[:, 1])
    return x_line, y_line


def plot_loess_scatter(
    df,
    x_col,
    y_col,
    x_label,
    y_label,
    filename,
    title=None,
    y_lim=(0, 1),
    x_lim=None,
    frac=0.2,
    scatter_alpha=0.15,
    scatter_size=10,
    line_color="#c0392b",
    line_width=2.5,
    figure_size=(10, 6),
):
    """
    Scatter with a single LOESS smooth (statsmodels LOWESS).
    """
    sns.set_style("white")

    work = df[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(work) < 10:
        print(f"Not enough valid data for LOESS plot: {filename}")
        return

    x_line, y_line = _fit_loess_curve(work[x_col], work[y_col], frac=frac)
    if x_line is None:
        print(f"LOESS fit failed for {filename}")
        return

    plt.figure(figsize=figure_size)
    plt.scatter(
        work[x_col],
        work[y_col],
        alpha=scatter_alpha,
        s=scatter_size,
        color="#1f77b4",
        edgecolors="none",
        label="Samples",
    )
    plt.plot(x_line, y_line, color=line_color, linewidth=line_width, label="LOESS")

    if title:
        plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    plt.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if x_lim is not None:
        plt.xlim(x_lim)
    else:
        lim = padded_axis_limits(work[x_col])
        if lim is not None:
            plt.xlim(lim)
    if y_lim is not None:
        plt.ylim(y_lim)

    plt.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, -0.12),
        fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved LOESS plot to {filename}")
    plt.close()


def plot_loess_scatter_strata(
    df,
    x_col,
    y_col,
    strata,
    x_label,
    y_label,
    filename,
    title=None,
    y_lim=(0, 1),
    x_lim=None,
    frac=0.2,
    scatter_alpha=0.08,
    scatter_size=8,
    line_width=2.5,
    figure_size=(10, 6),
):
    """
    Overlay LOESS curves for multiple row subsets of ``df``.

    ``strata`` is a list of dicts with keys ``label``, ``mask`` (boolean Series/array),
    and optional ``color``.
    """
    sns.set_style("white")

    base = df.replace([np.inf, -np.inf], np.nan)
    if base[[x_col, y_col]].dropna().shape[0] < 10:
        print(f"Not enough valid data for stratified LOESS plot: {filename}")
        return

    palette = sns.color_palette("deep", max(len(strata), 3))
    plt.figure(figsize=figure_size)

    plotted = 0
    for i, spec in enumerate(strata):
        label = spec["label"]
        mask = spec["mask"]
        color = spec.get("color", palette[i % len(palette)])
        sub = base.loc[mask, [x_col, y_col]].dropna()
        if len(sub) < 10:
            print(f"[Warning] Skipping stratum '{label}': only {len(sub)} samples")
            continue
        plt.scatter(
            sub[x_col],
            sub[y_col],
            alpha=scatter_alpha,
            s=scatter_size,
            color=color,
            edgecolors="none",
        )
        x_line, y_line = _fit_loess_curve(sub[x_col], sub[y_col], frac=frac)
        if x_line is None:
            continue
        plt.plot(x_line, y_line, color=color, linewidth=line_width, label=label)
        plotted += 1

    if plotted == 0:
        print(f"No strata had enough samples for LOESS: {filename}")
        plt.close()
        return

    if title:
        plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    plt.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if x_lim is not None:
        plt.xlim(x_lim)
    else:
        lim = padded_axis_limits(base[x_col].dropna())
        if lim is not None:
            plt.xlim(lim)
    if y_lim is not None:
        plt.ylim(y_lim)

    plt.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, -0.12),
        fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved stratified LOESS plot to {filename}")
    plt.close()


def plot_gam_tensor_surface(
    x_grid,
    y_grid,
    z_grid,
    x_label,
    y_label,
    cbar_label,
    title,
    filename,
    cmap=MISSING_RATE_CMAP,
    vmin=0.0,
    vmax=1.0,
    figure_size=(10, 8),
):
    """Contour/heatmap of a bivariate GAM tensor-product surface."""
    sns.set_style("white")

    plt.figure(figsize=figure_size)
    ax = plt.gca()
    cf = ax.contourf(x_grid, y_grid, z_grid, levels=40, cmap=cmap, vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(cf, ax=ax)
    cbar.set_label(cbar_label, fontsize=Y_LABEL_FONT_SIZE)
    cbar.ax.tick_params(labelsize=TICK_FONT_SIZE)
    ax.set_xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved GAM tensor surface to {filename}")
    plt.close()


def plot_gam_partial_curve(
    x_line,
    y_line,
    x_label,
    y_label,
    filename,
    title=None,
    y_lim=(0, 1),
    line_color="#c0392b",
    line_width=2.5,
    figure_size=(10, 6),
    reference_note=None,
):
    """1D GAM partial effect curve (optionally annotated with conditioning note)."""
    sns.set_style("white")

    plt.figure(figsize=figure_size)
    plt.plot(x_line, y_line, color=line_color, linewidth=line_width, label="GAM partial")
    if title:
        plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    plt.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if y_lim is not None:
        plt.ylim(y_lim)
    lim = padded_axis_limits(x_line)
    if lim is not None:
        plt.xlim(lim)
    if reference_note:
        plt.text(
            0.02,
            0.02,
            reference_note,
            transform=plt.gca().transAxes,
            fontsize=TICK_FONT_SIZE,
            va="bottom",
            ha="left",
        )
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved GAM partial curve to {filename}")
    plt.close()


def _finalize_panel_figure(fig, suptitle=None, **adjust_kwargs):
    """Apply subplot margins and place suptitle just above subplot row titles."""
    top = adjust_kwargs.get("top", 0.87)
    fig.subplots_adjust(**adjust_kwargs)
    if suptitle:
        fig.suptitle(suptitle, fontsize=TITLE_FONT_SIZE, y=min(top + 0.085, 0.98))


def plot_contourf_panels(
    panels,
    x_label,
    y_label,
    cbar_label,
    filename,
    suptitle=None,
    cmap=MISSING_RATE_CMAP,
    vmin=0.0,
    vmax=1.0,
    figure_size=(18, 5.5),
):
    """
    Side-by-side contour panels sharing one color scale.

    Each entry in ``panels`` is a dict with keys ``x_grid``, ``y_grid``, ``z_grid``, ``title``.
    """
    sns.set_style("white")
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=figure_size)
    if n == 1:
        axes = [axes]

    mappable = None
    for ax, spec in zip(axes, panels):
        cf = ax.contourf(
            spec["x_grid"],
            spec["y_grid"],
            spec["z_grid"],
            levels=40,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        mappable = cf
        ax.set_title(spec.get("title", ""), fontsize=LEGEND_FONT_SIZE, pad=6)
        ax.set_xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
        if ax is axes[0]:
            ax.set_ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
        ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)

    panel_adjust = dict(left=0.08, right=0.79, top=0.87, bottom=0.13, wspace=0.16)
    _finalize_panel_figure(fig, suptitle=suptitle, **panel_adjust)

    if mappable is not None:
        cbar_ax = fig.add_axes([0.825, 0.13, 0.012, 0.68])
        cbar = fig.colorbar(mappable, cax=cbar_ax)
        cbar.set_label(cbar_label, fontsize=Y_LABEL_FONT_SIZE, labelpad=10)
        cbar.ax.tick_params(labelsize=TICK_FONT_SIZE)

    plt.savefig(filename, dpi=300, pad_inches=0.06)
    print(f"Saved contour panel figure to {filename}")
    plt.close()


def plot_line_panels(
    panels,
    x_label,
    y_label,
    filename,
    suptitle=None,
    y_lim=(0, 1),
    figure_size=(18, 5),
    line_width=2.5,
):
    """Side-by-side line panels (e.g. GAM partial effects per subgenome)."""
    sns.set_style("white")
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=figure_size)
    if n == 1:
        axes = [axes]

    palette = sns.color_palette("deep", n)
    for ax, spec, color in zip(axes, panels, palette):
        ax.plot(spec["x"], spec["y"], color=color, linewidth=line_width)
        ax.set_title(spec.get("title", ""), fontsize=LEGEND_FONT_SIZE, pad=6)
        ax.set_xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
        if ax is axes[0]:
            ax.set_ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
        ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
        if y_lim is not None:
            ax.set_ylim(y_lim)
        lim = padded_axis_limits(spec["x"])
        if lim is not None:
            ax.set_xlim(lim)
        note = spec.get("note")
        if note:
            ax.text(
                0.02,
                0.02,
                note,
                transform=ax.transAxes,
                fontsize=TICK_FONT_SIZE,
                va="bottom",
                ha="left",
            )

    panel_adjust = dict(left=0.08, right=0.98, top=0.87, bottom=0.13, wspace=0.16)
    _finalize_panel_figure(fig, suptitle=suptitle, **panel_adjust)
    plt.savefig(filename, dpi=300, pad_inches=0.04)
    print(f"Saved line panel figure to {filename}")
    plt.close()


def plot_clustermap(
    data_matrix,
    row_colors=None,
    col_colors=None,
    title="Clustermap",
    filename="clustermap.png",
    figsize=(14, 14),
    cmap='viridis',
    xticklabels=False,
    yticklabels=False
):
    """
    Plots a clustermap (heatmap with hierarchical clustering).
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    # Increase recursion limit for large matrices
    sys.setrecursionlimit(20000)
    
    try:
        g = sns.clustermap(
            data_matrix,
            metric="euclidean",
            method="average",
            cmap=cmap,
            figsize=figsize,
            row_colors=row_colors,
            col_colors=col_colors,
            xticklabels=xticklabels,
            yticklabels=yticklabels,
            dendrogram_ratio=(0.15, 0.15),
            cbar_pos=(0.02, 0.8, 0.03, 0.15) # left, bottom, width, height
        )
        
        # Add Title
        g.fig.suptitle(title, fontsize=TITLE_FONT_SIZE, y=0.98)
        
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Saved clustermap to {filename}")
        plt.close()
    except Exception as e:
        print(f"Error plotting clustermap: {e}")


def plot_stacked_histogram(
    data, x_col, hue_col, 
    hue_order, palette_map,
    title, filename,
    xlabel, ylabel,
    mean_val=None, median_val=None,
    xlim=None, log_scale=False,
    legend_title="Categories",
):
    """
    Plots a stacked histogram using seaborn.
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Bins logic
    bins = 100
    binrange = xlim if xlim else None

    sns.histplot(
        data=data,
        x=x_col,
        hue=hue_col,
        hue_order=hue_order,
        bins=bins,
        binrange=binrange,
        multiple='stack',
        palette=palette_map,
        edgecolor='white',
        alpha=0.8,
        kde=False,
        legend=False,
        ax=ax,
    )

    # seaborn histplot legend handles are unreliable via get_legend_handles_labels()
    # in some matplotlib/seaborn builds — build hue patches explicitly.
    hue_handles = [
        Patch(
            facecolor=palette_map[label],
            edgecolor='white',
            linewidth=0.5,
            label=label,
        )
        for label in hue_order
    ]
    hue_labels = list(hue_order)

    stat_handles: list[Line2D] = []
    stat_labels: list[str] = []
    if mean_val is not None:
        ax.axvline(mean_val, color=TWO_CATEGORY_RED, linestyle='--', linewidth=1.5)
        stat_handles.append(Line2D([0], [0], color=TWO_CATEGORY_RED, linestyle='--', linewidth=1.5))
        stat_labels.append(f'Mean: {mean_val:.4f}')
    if median_val is not None:
        ax.axvline(median_val, color=TWO_CATEGORY_BLUE, linestyle='-', linewidth=1.5)
        stat_handles.append(Line2D([0], [0], color=TWO_CATEGORY_BLUE, linestyle='-', linewidth=1.5))
        stat_labels.append(f'Median: {median_val:.4f}')

    ax.legend(
        handles=hue_handles + stat_handles,
        labels=hue_labels + stat_labels,
        title=legend_title,
        loc='upper left',
        bbox_to_anchor=(0.0, -0.15),
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )

    if log_scale:
        ax.set_yscale('log')
        ax.set_ylabel(f"{ylabel} (Log Scale)", fontsize=Y_LABEL_FONT_SIZE)
    else:
        ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)

    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.set_xlabel(xlabel, fontsize=X_LABEL_FONT_SIZE)
    ax.tick_params(axis='both', which='major', labelsize=TICK_FONT_SIZE)
    if xlim:
        ax.set_xlim(xlim)
    ax.grid(True, which="major", axis='y', alpha=0.3)
    fig.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved {filename}")
    plt.close(fig)


def log10_offset(values, offset=1.0):
    """Return ``log10(values + offset)`` so zeros stay defined and finite."""
    return np.log10(np.asarray(values, dtype=np.float64) + float(offset))


def plot_grouped_density_lines(
    data,
    x_col,
    hue_col,
    hue_order,
    palette_map,
    title,
    filename,
    xlabel=None,
    ylabel="Density",
    xlim=None,
    n_bins=400,
    density=True,
    log_scale=False,
    log_offset=1.0,
    hue_labels=None,
    figure_size=(10, 6),
    linewidth=2.0,
):
    """Frequency-polygon density: many bins, line through bin tops, one curve per hue.

    ``n_bins`` large approximates a continuous empirical density. Each group is
    drawn as its own line (not stacked). ``density=True`` normalises each group
    so curves remain comparable when counts differ.

    When ``log_scale`` is true, plot ``log10(density + log_offset)`` on a linear
    axis (default offset 1) so empty bins stay at 0 instead of breaking the
    line. The transform is labelled on the axes and in-plot.
    """
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=figure_size)
    if xlabel is None:
        xlabel = x_col
    if hue_labels is None:
        hue_labels = {k: k for k in hue_order}
    offset = float(log_offset)
    if offset <= 0.0:
        raise ValueError("log_offset must be positive so log10(density + offset) is defined")

    plot_df = data[[x_col, hue_col]].replace([np.inf, -np.inf], np.nan).dropna()
    x_lo, x_hi = (0.0, 1.0) if xlim is None else (float(xlim[0]), float(xlim[1]))
    if xlim is None and not plot_df.empty:
        x_lo = float(np.nanmin(plot_df[x_col].to_numpy(dtype=np.float64)))
        x_hi = float(np.nanmax(plot_df[x_col].to_numpy(dtype=np.float64)))
        if not np.isfinite(x_lo) or not np.isfinite(x_hi) or x_hi <= x_lo:
            x_lo, x_hi = 0.0, 1.0
    n_bins = max(8, int(n_bins))
    edges = np.linspace(x_lo, x_hi, n_bins + 1)
    mids = 0.5 * (edges[:-1] + edges[1:])

    for hue in hue_order:
        vals = plot_df.loc[plot_df[hue_col] == hue, x_col].to_numpy(dtype=np.float64)
        if len(vals) == 0:
            continue
        heights, _ = np.histogram(vals, bins=edges, density=bool(density))
        y = heights.astype(np.float64)
        if log_scale:
            y = log10_offset(y, offset=offset)
        ax.plot(
            mids,
            y,
            color=palette_map.get(hue, "steelblue"),
            linestyle="-",
            linewidth=float(linewidth),
            label=hue_labels.get(hue, hue),
        )

    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.set_xlabel(xlabel, fontsize=X_LABEL_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if log_scale:
        offset_label = f"{offset:g}"
        log_ylabel = rf"$\log_{{10}}(\mathrm{{density}} + {offset_label})$"
        ax.set_ylabel(log_ylabel, fontsize=Y_LABEL_FONT_SIZE)
        ax.text(
            0.04,
            0.92,
            rf"$y = \log_{{10}}(\mathrm{{density}} + {offset_label})$",
            transform=ax.transAxes,
            fontsize=LEGEND_FONT_SIZE,
            va="top",
            ha="left",
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="gray", boxstyle="round"),
        )
        ax.set_ylim(bottom=0.0)
    else:
        ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_xlim(x_lo, x_hi)
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, -0.15),
        fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved {filename}")
    plt.close(fig)


def plot_kde_density_lines(
    data,
    x_col,
    hue_col,
    hue_order,
    palette_map,
    title,
    filename,
    xlabel=None,
    ylabel="Density",
    xlim=None,
    bw_method=0.15,
    hue_labels=None,
    figure_size=(10, 6),
    linewidth=2.0,
    fill=True,
    fill_alpha=0.16,
    x_points=600,
):
    """Gaussian-KDE density curves, one per hue, with optional fill.

    Use this instead of :func:`plot_grouped_density_lines` when a group has
    only a few hundred values: the frequency polygon becomes jagged at that
    size, while a KDE stays smooth. ``bw_method`` follows
    :class:`scipy.stats.gaussian_kde` (float, ``"scott"`` or ``"silverman"``).
    """
    from scipy.stats import gaussian_kde

    sns.set_style("white")
    fig, ax = plt.subplots(figsize=figure_size)
    if xlabel is None:
        xlabel = x_col
    if hue_labels is None:
        hue_labels = {k: k for k in hue_order}

    plot_df = data[[x_col, hue_col]].replace([np.inf, -np.inf], np.nan).dropna()
    x_lo, x_hi = (0.0, 1.0) if xlim is None else (float(xlim[0]), float(xlim[1]))
    if xlim is None and not plot_df.empty:
        x_lo = float(np.nanmin(plot_df[x_col].to_numpy(dtype=np.float64)))
        x_hi = float(np.nanmax(plot_df[x_col].to_numpy(dtype=np.float64)))
        if not np.isfinite(x_lo) or not np.isfinite(x_hi) or x_hi <= x_lo:
            x_lo, x_hi = 0.0, 1.0
    xs = np.linspace(x_lo, x_hi, max(50, int(x_points)))

    for hue in hue_order:
        vals = plot_df.loc[plot_df[hue_col] == hue, x_col].to_numpy(dtype=np.float64)
        if len(vals) < 2 or float(np.nanstd(vals)) <= 0.0:
            continue
        density = gaussian_kde(vals, bw_method=bw_method)(xs)
        color = palette_map.get(hue, "steelblue")
        if fill:
            ax.fill_between(xs, 0.0, density, color=color, alpha=fill_alpha, linewidth=0)
        ax.plot(
            xs,
            density,
            color=color,
            linestyle="-",
            linewidth=float(linewidth),
            label=hue_labels.get(hue, hue),
        )

    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.set_xlabel(xlabel, fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    ax.set_xlim(x_lo, x_hi)
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, -0.15),
        fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved {filename}")
    plt.close(fig)


def plot_grouped_histogram(
    data,
    x_col,
    hue_col,
    hue_order,
    palette_map,
    title,
    filename,
    xlabel=None,
    ylabel="Density",
    xlim=None,
    bins=80,
    density=True,
    hue_labels=None,
    figure_size=(10, 6),
    thresholds_v=None,
):
    """Grouped histogram bars, one bar series per hue.

    This is the bar-chart counterpart of :func:`plot_kde_density_lines`.
    Each hue gets ``1 / len(hue_order)`` of every bin, so the bars stay
    side-by-side and density-normalised when group sizes differ.
    ``thresholds_v`` draws vertical reference lines, each a dict with
    ``value`` plus optional ``color``, ``linestyle`` and ``label``.
    """
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=figure_size)
    if xlabel is None:
        xlabel = x_col
    if hue_labels is None:
        hue_labels = {k: k for k in hue_order}
    plot_df = data[[x_col, hue_col]].replace([np.inf, -np.inf], np.nan).dropna()
    x_lo, x_hi = (0.0, 1.0) if xlim is None else (float(xlim[0]), float(xlim[1]))
    if xlim is None and not plot_df.empty:
        x_lo = float(np.nanmin(plot_df[x_col].to_numpy(dtype=np.float64)))
        x_hi = float(np.nanmax(plot_df[x_col].to_numpy(dtype=np.float64)))
        if not np.isfinite(x_lo) or not np.isfinite(x_hi) or x_hi <= x_lo:
            x_lo, x_hi = 0.0, 1.0
    bins = max(4, int(bins))
    edges = np.linspace(x_lo, x_hi, bins + 1)
    mids = 0.5 * (edges[:-1] + edges[1:])
    n_hues = max(1, len(hue_order))
    width = (edges[1] - edges[0]) / n_hues * 0.9

    for position, hue in enumerate(hue_order):
        vals = plot_df.loc[plot_df[hue_col] == hue, x_col].to_numpy(dtype=np.float64)
        if len(vals) == 0:
            continue
        heights, _ = np.histogram(vals, bins=edges, density=bool(density))
        offset = (position - (n_hues - 1) / 2.0) * width
        ax.bar(
            mids + offset,
            heights.astype(np.float64),
            width=width,
            color=palette_map.get(hue, "steelblue"),
            label=hue_labels.get(hue, hue),
            linewidth=0,
        )

    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.set_xlabel(xlabel, fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    ax.set_xlim(x_lo, x_hi)
    if thresholds_v:
        for thr in thresholds_v:
            value = thr.get("value")
            if value is None:
                continue
            ax.axvline(
                x=value,
                color=thr.get("color", "black"),
                linestyle=thr.get("linestyle", "--"),
                linewidth=thr.get("linewidth", 1.0),
                label=thr.get("label", f"{value:.4f}"),
                zorder=4,
            )
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, -0.15),
        fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved {filename}")
    plt.close(fig)


def plot_connected_scatter(
    points,
    edges,
    title,
    filename,
    x_col="PC1",
    y_col="PC2",
    id_col="id",
    source_col="source",
    target_col="target",
    category_col="category",
    edge_colors=None,
    label_categories=None,
    label_font_size=6,
    background_color="#cccccc",
    background_size=6,
    point_size=26,
    point_color="#444444",
    xlabel=None,
    ylabel=None,
    figure_size=(10, 8),
    focus=False,
    focus_padding=0.08,
    focus_percentiles=(0.01, 0.99),
    xlim=None,
    ylim=None,
):
    """Scatter with group edges drawn between its points.

    ``points`` is one row per item (``id``, ``x_col``, ``y_col``); every point is
    drawn as a light background. ``edges`` is one row per connection
    (``source``, ``target``, ``category``); each edge is a line coloured by its
    category, and the endpoints of edges belonging to ``label_categories`` are
    annotated with their id. Edges pointing outside ``points`` are dropped.
    ``focus`` crops the axes to the ``focus_percentiles`` band of the whole
    point cloud plus ``focus_padding``: that keeps a few outlier samples from
    flattening the figure without cropping away the population the edges sit in
    (cropping to the connected points alone hides most of a small group).
    ``xlim`` / ``ylim`` override that crop with explicit limits.
    """
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D

    sns.set_style("white")
    fig, ax = plt.subplots(figsize=figure_size)
    if xlabel is None:
        xlabel = x_col
    if ylabel is None:
        ylabel = y_col
    if edge_colors is None:
        edge_colors = {}
    label_categories = list(label_categories or [])

    coords = points[[id_col, x_col, y_col]].dropna()
    position = coords.set_index(id_col)[[x_col, y_col]]

    edge_df = edges.dropna(subset=[source_col, target_col]).copy()
    edge_df = edge_df.loc[
        edge_df[source_col].isin(position.index)
        & edge_df[target_col].isin(position.index)
    ]
    categories = list(edge_colors)
    for category in edge_df[category_col].unique():
        if category not in categories:
            categories.append(category)

    ax.scatter(
        coords[x_col],
        coords[y_col],
        s=background_size,
        color=background_color,
        linewidths=0,
        zorder=1,
    )

    for category in categories:
        sub = edge_df.loc[edge_df[category_col] == category]
        if sub.empty:
            continue
        segments = [
            [
                position.loc[source].to_numpy(dtype=float),
                position.loc[target].to_numpy(dtype=float),
            ]
            for source, target in zip(sub[source_col], sub[target_col])
        ]
        ax.add_collection(
            LineCollection(
                segments,
                colors=edge_colors.get(category, "grey"),
                linewidths=1.2,
                alpha=0.7,
                zorder=2,
            )
        )

    involved = sorted(
        set(edge_df[source_col]) | set(edge_df[target_col])
    )
    if involved:
        ax.scatter(
            position.loc[involved, x_col],
            position.loc[involved, y_col],
            s=point_size,
            color=point_color,
            linewidths=0,
            zorder=3,
        )

    labelled = sorted(
        set(edge_df.loc[edge_df[category_col].isin(label_categories), source_col])
        | set(edge_df.loc[edge_df[category_col].isin(label_categories), target_col])
    )
    for item in labelled:
        ax.annotate(
            str(item),
            (position.loc[item, x_col], position.loc[item, y_col]),
            fontsize=label_font_size,
            xytext=(3, 3),
            textcoords="offset points",
            color="#333333",
            zorder=4,
        )

    handles = [
        Line2D(
            [0],
            [0],
            color=edge_colors.get(category, "grey"),
            linewidth=1.5,
            label=str(category),
        )
        for category in categories
        if not edge_df.loc[edge_df[category_col] == category].empty
    ]
    if handles:
        ax.legend(
            handles=handles,
            loc="upper right",
            fontsize=LEGEND_FONT_SIZE,
            frameon=False,
        )

    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.set_xlabel(xlabel, fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if focus and not coords.empty:
        lower_q, upper_q = focus_percentiles
        for axis, column in ((ax.set_xlim, x_col), (ax.set_ylim, y_col)):
            values = coords[column].to_numpy(dtype=float)
            low = float(np.nanquantile(values, lower_q))
            high = float(np.nanquantile(values, upper_q))
            span = high - low
            if span <= 0:
                span = abs(high) or 1.0
            pad = span * focus_padding
            axis(low - pad, high + pad)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved {filename}")
    plt.close(fig)


def plot_scatter_with_thresholds(
    data, x_col, y_col,
    title, filename,
    thresholds_h=None, 
    thresholds_v=None,
    xlabel=None, ylabel=None,
    color='royalblue',
    alpha=0.5,
    s=25,
    figure_size=(10,8),
    xlim=None,
    ylim=None,
    group_col=None,
    hue_palette=None,
    highlight_ids=None,
    id_col=None,
    highlight_color='crimson',
    highlight_marker='*',
    highlight_size=160,
    x_spans=None,
    connect=False,
    connect_flag_col=None,
    connect_alt_color='crimson',
    connect_alt_linestyle='--',
    connect_flag_true_label=None,
    connect_flag_false_label=None,
    identity_line=False,
    identity_label=r"$y=x$",
    identity_color="#7f7f7f",
    size_col=None,
    line_overlays=None,
    ax=None,
    save=True,
    show_legend=True,
    legend_outside=True,
    log_x=False,
    log_y=False,
):
    """
    Plots a scatter plot with optional horizontal and vertical threshold lines.
    
    Args:
        thresholds_h (list of dict): e.g. [{'value': 0.95, 'color': 'red', 'linestyle': '--', 'label': 'Cutoff'}]
        thresholds_v (list of dict): e.g. [{'value': 0.5, 'color': 'orange', 'linestyle': '--', 'label': 'Limit'}]
        group_col: Optional column name for hue (e.g. ``Group`` from ``anno_group``).
        hue_palette: Optional seaborn palette (dict or list) used when
            ``group_col`` is set but is not the canonical ``Group`` column.
        highlight_ids: Sample IDs to mark with a distinct marker and text label.
        id_col: Column holding sample IDs for ``highlight_ids`` (defaults to Sample, IID, #IID).
        x_spans: Optional interval shading on x, each dict with start/end and optional color/alpha/label.
        connect: If True, draw a line through points sorted by ``x_col``.
        connect_flag_col: Optional boolean column. True points use ``color`` and a
            solid connector; False points use ``connect_alt_color`` and dashed
            segments to neighbors (including True neighbors).
        connect_alt_color: Color for False-flag points and mixed/False segments.
        connect_alt_linestyle: Linestyle for mixed/False segments (default dashed).
        connect_flag_true_label: Legend label for True-flag points.
        connect_flag_false_label: Legend label for False-flag points.
        identity_line: If True, draw y=x (two-category contrast uses blue vs this gray).
        size_col: Optional numeric column mapped to point area.
        line_overlays: Optional list of dicts with ``x`` / ``y`` arrays and
            optional ``color``, ``linestyle``, ``linewidth``, ``marker``,
            ``markersize``, ``label``, and ``zorder`` (same shape as
            ``plot_scatter_continuous_color``). Use it for a fitted trend or a
            binned summary drawn over the points.
        ax: Optional existing axes. When set, draw on it instead of creating a figure.
        save: If False, skip ``savefig`` (use with ``ax`` for multi-panel figures).
        show_legend: If False, skip the legend even when labels exist.
        legend_outside: If True (default), place the legend below the axes.
        log_x: If True, use a log x-axis.
        log_y: If True, use a log y-axis.
    """
    sns.set_style("white")

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figure_size)
    else:
        fig = ax.get_figure()

    if x_spans:
        labeled_spans = set()
        for span in x_spans:
            start = span.get("start")
            end = span.get("end")
            if start is None or end is None:
                continue
            lbl = span.get("label")
            if lbl in labeled_spans:
                lbl = None
            elif lbl:
                labeled_spans.add(lbl)
            ax.axvspan(
                float(start),
                float(end),
                color=span.get("color", "#2ca02c"),
                alpha=span.get("alpha", 0.18),
                linewidth=0,
                zorder=0,
                label=lbl,
            )

    use_flag = bool(connect_flag_col) and connect_flag_col in data.columns
    if connect:
        line_cols = [x_col, y_col]
        if use_flag:
            line_cols.append(connect_flag_col)
        line_df = data[line_cols].replace([np.inf, -np.inf], np.nan).dropna().sort_values(x_col)
        if not line_df.empty and use_flag:
            xs = line_df[x_col].to_numpy(dtype=np.float64)
            ys = line_df[y_col].to_numpy(dtype=np.float64)
            flags = line_df[connect_flag_col].to_numpy(dtype=bool)
            for i in range(len(line_df) - 1):
                alt = not (bool(flags[i]) and bool(flags[i + 1]))
                ax.plot(
                    [xs[i], xs[i + 1]],
                    [ys[i], ys[i + 1]],
                    color=connect_alt_color if alt else color,
                    linestyle=connect_alt_linestyle if alt else "-",
                    linewidth=1.2,
                    zorder=1,
                    label="_nolegend_",
                )
        elif not line_df.empty:
            ax.plot(
                line_df[x_col],
                line_df[y_col],
                color=color,
                linewidth=1.2,
                zorder=1,
            )

    use_hue = group_col and group_col in data.columns
    if use_flag:
        scatter_df = data[[x_col, y_col, connect_flag_col]].replace(
            [np.inf, -np.inf], np.nan
        ).dropna()
        flag = scatter_df[connect_flag_col].astype(bool)
        pri = scatter_df.loc[flag]
        alt = scatter_df.loc[~flag]
        if not pri.empty:
            ax.scatter(
                pri[x_col],
                pri[y_col],
                c=color,
                alpha=alpha,
                s=s,
                edgecolors="w",
                zorder=2,
                label=connect_flag_true_label,
            )
        if not alt.empty:
            ax.scatter(
                alt[x_col],
                alt[y_col],
                c=connect_alt_color,
                alpha=alpha,
                s=max(float(s) * 1.6, 44.0),
                edgecolors="darkred",
                linewidths=0.5,
                zorder=3,
                label=connect_flag_false_label,
            )
    elif use_hue:
        plot_data = data[[x_col, y_col, group_col]].replace([np.inf, -np.inf], np.nan).dropna()
        hue_kw = _hue_kwargs_for_sample_group(plot_data, group_col)
        if hue_palette is not None and "palette" not in hue_kw:
            hue_kw["palette"] = hue_palette
        sns.scatterplot(
            data=plot_data, x=x_col, y=y_col,
            hue=group_col, alpha=alpha, s=s, edgecolor='w', zorder=2, ax=ax, **hue_kw,
        )
    else:
        if data is None or getattr(data, "empty", True) or x_col not in data.columns:
            plot_data = pd.DataFrame(columns=[x_col, y_col])
        else:
            cols = [x_col, y_col]
            if size_col and size_col in data.columns:
                cols.append(size_col)
            plot_data = data[cols].replace([np.inf, -np.inf], np.nan).dropna(subset=[x_col, y_col])
        obs_label = r"Observed $P_{+}$" if identity_line else None
        if not plot_data.empty and size_col and size_col in plot_data.columns:
            nvals = pd.to_numeric(plot_data[size_col], errors="coerce").to_numpy(dtype=np.float64)
            finite = np.isfinite(nvals)
            nmax = float(np.nanmax(nvals)) if finite.any() else 1.0
            nmax = max(nmax, 1.0)
            sizes = 40.0 + 140.0 * np.sqrt(np.clip(np.where(finite, nvals, 0.0), 0.0, None) / nmax)
            ax.scatter(
                plot_data[x_col],
                plot_data[y_col],
                s=sizes,
                alpha=alpha,
                c=color,
                edgecolors="w",
                zorder=2,
                label=obs_label,
            )
        elif not plot_data.empty:
            sns.scatterplot(
                data=plot_data, x=x_col, y=y_col,
                alpha=alpha, s=s, edgecolor='w', color=color, zorder=2, ax=ax,
            )

    if highlight_ids:
        hi_ids = {str(v) for v in highlight_ids if str(v).strip()}
        if hi_ids:
            id_candidates = [id_col] if id_col else []
            id_candidates.extend(c for c in ('Sample', 'IID', '#IID', 'FID') if c not in id_candidates)
            resolved_id_col = next((c for c in id_candidates if c in data.columns), None)
            if resolved_id_col is None:
                print(f"[Warning] highlight_ids requested but no ID column found in {list(data.columns)}")
            else:
                hi = data[data[resolved_id_col].astype(str).isin(hi_ids)].copy()
                hi = hi[[x_col, y_col, resolved_id_col]].replace([np.inf, -np.inf], np.nan).dropna()
                if hi.empty:
                    print(f"[Warning] highlight_ids not found in column {resolved_id_col}: {sorted(hi_ids)}")
                else:
                    ax.scatter(
                        hi[x_col],
                        hi[y_col],
                        c=highlight_color,
                        s=highlight_size,
                        marker=highlight_marker,
                        edgecolors='black',
                        linewidths=0.6,
                        zorder=5,
                        label='Highlighted',
                    )
                    for _, row in hi.iterrows():
                        ax.annotate(
                            str(row[resolved_id_col]),
                            (row[x_col], row[y_col]),
                            xytext=(6, 6),
                            textcoords='offset points',
                            fontsize=TICK_FONT_SIZE,
                            color=highlight_color,
                            fontweight='bold',
                            zorder=6,
                        )
    if line_overlays:
        for overlay in line_overlays:
            xs = overlay.get("x")
            ys = overlay.get("y")
            if xs is None or ys is None:
                continue
            ax.plot(
                xs,
                ys,
                color=overlay.get("color", "black"),
                linestyle=overlay.get("linestyle", "-"),
                linewidth=overlay.get("linewidth", 2.0),
                marker=overlay.get("marker"),
                markersize=overlay.get("markersize", 5),
                zorder=overlay.get("zorder", 4),
                label=overlay.get("label"),
            )

    # Horizontal Thresholds
    if thresholds_h:
        for thr in thresholds_h:
            val = thr.get('value')
            if val is None: continue
            lbl = thr.get('label', f'{val:.4f}')
            clr = thr.get('color', 'red')
            ls = thr.get('linestyle', '--')
            ax.axhline(y=val, color=clr, linestyle=ls, label=lbl)

    # Vertical Thresholds
    if thresholds_v:
        for thr in thresholds_v:
            val = thr.get('value')
            if val is None: continue
            lbl = thr.get('label', f'{val:.4f}')
            clr = thr.get('color', 'orange')
            ls = thr.get('linestyle', '--')
            ax.axvline(x=val, color=clr, linestyle=ls, label=lbl)

    if xlabel is None: xlabel = x_col
    if ylabel is None: ylabel = y_col
    
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.set_xlabel(xlabel, fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    ax.tick_params(axis='both', which='major', labelsize=TICK_FONT_SIZE)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    if identity_line:
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        lo = max(min(x0, x1), min(y0, y1))
        hi = min(max(x0, x1), max(y0, y1))
        ax.plot(
            [lo, hi],
            [lo, hi],
            color=identity_color,
            linestyle="--",
            linewidth=1.5,
            zorder=1,
            label=identity_label,
        )

    span_labeled = bool(x_spans) and any(span.get("label") for span in x_spans)
    flag_labeled = use_flag and bool(connect_flag_true_label or connect_flag_false_label)
    line_labeled = bool(line_overlays) and any(
        overlay.get("label") for overlay in line_overlays
    )
    if show_legend and (use_hue or thresholds_h or thresholds_v or highlight_ids or span_labeled or flag_labeled or line_labeled or identity_line):
        legend_kw = dict(
            fontsize=LEGEND_FONT_SIZE,
            frameon=False,
            title=group_col if use_hue else None,
        )
        if legend_outside:
            ax.legend(loc='upper left', bbox_to_anchor=(0.0, -0.15), **legend_kw)
        else:
            ax.legend(loc='upper right', **legend_kw)
    
    if save and filename:
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Saved scatter plot to {filename}")
    if own_fig and save:
        plt.close(fig)
    return ax


def plot_slope_chart(
    df, col1, col2,
    xlabel1, xlabel2,
    title, filename,
    line_color='gray',
    highlight_color='red',
    figure_size=(6, 8),
    group_col=None,
    ylabel='Value',
    ylim=None,
    line_alpha=None,
):
    """
    Plots a slope chart (parallel coordinates for 2 variables).

    When ``group_col`` is set (typically ``Group``), each sample line uses
    :func:`sample_group_palette` colors; otherwise all lines share ``line_color``.
    """
    from matplotlib.lines import Line2D

    sns.set_style("white")
    plt.figure(figsize=figure_size)

    plot_df = df[[col1, col2]].replace([np.inf, -np.inf], np.nan).dropna()
    if group_col and group_col in df.columns:
        plot_df = df[[col1, col2, group_col]].replace([np.inf, -np.inf], np.nan).dropna()
        hue_order, palette_map = sample_group_hue_config(plot_df[group_col])
        alpha = line_alpha if line_alpha is not None else 0.35
        for grp in hue_order:
            sub = plot_df[plot_df[group_col] == grp]
            color = palette_map.get(grp, SAMPLE_GROUP_UNKNOWN_COLOR)
            for _, row in sub.iterrows():
                plt.plot(
                    [1, 2], [row[col1], row[col2]],
                    marker='o', alpha=alpha, color=color, linewidth=1,
                )
        group_handles = [
            Line2D([0], [0], color=palette_map[g], lw=2, label=g)
            for g in hue_order
        ]
    else:
        alpha = line_alpha if line_alpha is not None else 0.1
        for _, row in plot_df.iterrows():
            plt.plot(
                [1, 2], [row[col1], row[col2]],
                marker='o', alpha=alpha, color=line_color,
            )
        group_handles = []

    avg1 = plot_df[col1].mean()
    avg2 = plot_df[col2].mean()
    plt.plot(
        [1, 2], [avg1, avg2], marker='o', color=highlight_color,
        linewidth=3, label='Population Mean',
    )

    plt.xticks([1, 2], [xlabel1, xlabel2], fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    if ylim is not None:
        plt.ylim(ylim)

    handles = group_handles + [
        Line2D([0], [0], color=highlight_color, lw=3, label='Population Mean'),
    ]
    plt.legend(
        handles=handles,
        loc='upper left', bbox_to_anchor=(0.0, -0.12),
        fontsize=LEGEND_FONT_SIZE, frameon=False,
        title=group_col if group_handles else None,
    )

    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved slope chart to {filename}")
    plt.close()


def plot_parallel_coordinates(
    data,
    value_cols,
    title,
    filename,
    series_col,
    palette_map=None,
    x_labels=None,
    ylabel=None,
    invert_yaxis=False,
    yticks=None,
    annotate=False,
    annotate_template="{series}{axis}",
    marker="o",
    markersize=8,
    linewidth=2.0,
    figure_size=(7, 8),
):
    """Parallel-coordinate / bump chart: one polyline per series across axes.

    ``value_cols`` are y-values at equally spaced vertical axes (x = 1, 2, …).
    ``series_col`` colours each polyline. When ``invert_yaxis`` is true, smaller
    y sits at the top (rank charts with 1 = highest).
    """
    from matplotlib.lines import Line2D

    sns.set_style("white")
    plt.figure(figsize=figure_size)
    cols = list(value_cols)
    if not cols:
        raise ValueError("value_cols must be non-empty")
    xs = list(range(1, len(cols) + 1))
    if x_labels is None:
        x_labels = cols
    if ylabel is None:
        ylabel = ""
    needed = [series_col] + cols
    plot_df = data[needed].copy()
    for col in cols:
        plot_df[col] = pd.to_numeric(plot_df[col], errors="coerce")
    series_order = [s for s in plot_df[series_col].drop_duplicates().tolist()]
    if palette_map is None:
        palette = sns.color_palette("tab10", max(len(series_order), 1))
        palette_map = {s: palette[i % len(palette)] for i, s in enumerate(series_order)}

    handles = []
    for series in series_order:
        row = plot_df[plot_df[series_col] == series].iloc[0]
        ys = [row[col] for col in cols]
        mask = [pd.notna(y) for y in ys]
        if not any(mask):
            continue
        x_draw = [x for x, keep in zip(xs, mask) if keep]
        y_draw = [y for y, keep in zip(ys, mask) if keep]
        color = palette_map.get(series, palette_map.get(int(series) if str(series).isdigit() else series, "#444444"))
        label = f"chr{series}" if str(series).isdigit() else str(series)
        plt.plot(
            x_draw,
            y_draw,
            color=color,
            linestyle="-",
            linewidth=linewidth,
            marker=marker,
            markersize=markersize,
            label=label,
            zorder=3,
        )
        if annotate:
            for x, y, axis in zip(x_draw, y_draw, [lab for lab, keep in zip(x_labels, mask) if keep]):
                text = annotate_template.format(series=series, axis=axis)
                ha = "right" if x == xs[0] else ("left" if x == xs[-1] else "center")
                dx = -0.06 if ha == "right" else (0.06 if ha == "left" else 0.0)
                plt.text(
                    x + dx,
                    y,
                    text,
                    color=color,
                    fontsize=TICK_FONT_SIZE,
                    ha=ha,
                    va="center",
                    zorder=4,
                )
        handles.append(Line2D([0], [0], color=color, marker=marker, lw=linewidth, label=label))

    plt.xticks(xs, x_labels, fontsize=X_LABEL_FONT_SIZE)
    plt.xlim(xs[0] - 0.35, xs[-1] + 0.35)
    plt.ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    for x in xs:
        plt.axvline(x, color="#dddddd", linewidth=1.0, zorder=0)
    if yticks is not None:
        plt.yticks(yticks)
    if invert_yaxis:
        plt.gca().invert_yaxis()
    plt.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.12),
        fontsize=LEGEND_FONT_SIZE,
        frameon=False,
        ncol=min(4, max(1, len(handles))),
    )
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved plot to {filename}")
    plt.close()


def plot_scatter_two_categories(
    data,
    x_col,
    y_col,
    category_col,
    *,
    category_labels,
    category_colors,
    title,
    filename,
    xlabel=None,
    ylabel=None,
    alpha=0.5,
    s=25,
    figure_size=(10, 8),
    dpi=300,
    xlim=None,
    ylim=None,
    legend_below=False,
    rasterized=False,
    log_x=False,
    log_y=False,
):
    """Scatter plot with two categorical point colors (e.g. variant vs non-variant sites)."""
    import os

    sns.set_style("white")
    plt.figure(figsize=figure_size)

    categories = sorted(category_labels.keys(), key=lambda k: (str(k), repr(k)))
    for cat in categories:
        subset = data[data[category_col] == cat]
        if subset.empty:
            continue
        plt.scatter(
            subset[x_col],
            subset[y_col],
            s=s,
            alpha=alpha,
            label=category_labels[cat],
            color=category_colors.get(cat, "steelblue"),
            rasterized=rasterized,
        )

    plt.xlabel(xlabel if xlabel else x_col, fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel(ylabel if ylabel else y_col, fontsize=Y_LABEL_FONT_SIZE)
    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if legend_below:
        plt.legend(
            loc="upper left",
            bbox_to_anchor=(0.0, -0.15),
            fontsize=LEGEND_FONT_SIZE,
            frameon=False,
        )
    else:
        plt.legend(loc="best", fontsize=LEGEND_FONT_SIZE)
    if log_x:
        plt.xscale("log")
    if log_y:
        plt.yscale("log")
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)

    parent_dir = os.path.dirname(filename)
    if parent_dir and not os.path.exists(parent_dir):
        os.makedirs(parent_dir, exist_ok=True)

    plt.savefig(filename, dpi=dpi, bbox_inches="tight")
    plt.close()
    print(f"Saved two-category scatter plot to {filename}")


def _bbox_overlap_wh(a, b, pad=2.0):
    """Return (overlap_width, overlap_height) in pixels; negatives mean a gap."""
    return (
        min(a.x1, b.x1) - max(a.x0, b.x0) + pad,
        min(a.y1, b.y1) - max(a.y0, b.y0) + pad,
    )


def place_nonoverlapping_labels(
    ax,
    x,
    y,
    labels,
    *,
    fontsize=9,
    italic=False,
    color="#222222",
    max_iter=500,
    seed=0,
):
    """
    Place text labels at data points and nudge them in display space so boxes do not overlap.

    Returns the matplotlib Annotation artists. Call after limits and scales are set.
    """
    xs = np.asarray(x, dtype=np.float64)
    ys = np.asarray(y, dtype=np.float64)
    n = len(labels)
    if n == 0:
        return []
    if xs.shape[0] != n or ys.shape[0] != n:
        raise ValueError("x, y, and labels must have the same length")

    rng = np.random.default_rng(seed)
    texts = []
    for i, (xi, yi, lab) in enumerate(zip(xs, ys, labels)):
        angle = 2.0 * np.pi * (i / max(n, 1)) + 0.35
        texts.append(
            ax.annotate(
                str(lab),
                (xi, yi),
                xytext=(12.0 * np.cos(angle), 12.0 * np.sin(angle)),
                textcoords="offset points",
                fontsize=fontsize,
                fontstyle="italic" if italic else "normal",
                color=color,
                ha="center",
                va="center",
                clip_on=False,
                zorder=6,
                arrowprops={
                    "arrowstyle": "-",
                    "color": "#8a8a8a",
                    "lw": 0.5,
                    "shrinkA": 0,
                    "shrinkB": 3,
                },
            )
        )

    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    sizes = np.array(
        [[t.get_window_extent(renderer).width, t.get_window_extent(renderer).height] for t in texts],
        dtype=np.float64,
    )
    pts = ax.transData.transform(np.column_stack([xs, ys]))
    offsets_pt = np.array([t.get_position() for t in texts], dtype=np.float64)
    px_per_pt = fig.dpi / 72.0
    max_offset_pt = 180.0
    half = sizes * 0.55

    for _ in range(max_iter):
        centers = pts + offsets_pt * px_per_pt
        x0 = centers[:, 0] - half[:, 0]
        x1 = centers[:, 0] + half[:, 0]
        y0 = centers[:, 1] - half[:, 1]
        y1 = centers[:, 1] + half[:, 1]
        forces = np.zeros((n, 2), dtype=np.float64)
        still_overlap = False

        for i in range(n):
            for j in range(i + 1, n):
                ow = min(x1[i], x1[j]) - max(x0[i], x0[j]) + 3.0
                oh = min(y1[i], y1[j]) - max(y0[i], y0[j]) + 3.0
                if ow <= 0 or oh <= 0:
                    continue
                still_overlap = True
                if ow <= oh:
                    sign = 1.0 if centers[i, 0] >= centers[j, 0] else -1.0
                    if centers[i, 0] == centers[j, 0]:
                        sign = 1.0 if rng.random() > 0.5 else -1.0
                    push = sign * (ow * 0.5 + 2.0)
                    forces[i, 0] += push
                    forces[j, 0] -= push
                else:
                    sign = 1.0 if centers[i, 1] >= centers[j, 1] else -1.0
                    if centers[i, 1] == centers[j, 1]:
                        sign = 1.0 if rng.random() > 0.5 else -1.0
                    push = sign * (oh * 0.5 + 2.0)
                    forces[i, 1] += push
                    forces[j, 1] -= push

        for i in range(n):
            for k, (px, py) in enumerate(pts):
                if k == i:
                    continue
                if x0[i] - 5 <= px <= x1[i] + 5 and y0[i] - 5 <= py <= y1[i] + 5:
                    still_overlap = True
                    delta = np.array([centers[i, 0] - px, centers[i, 1] - py], dtype=np.float64)
                    dist = float(np.linalg.norm(delta))
                    if dist < 1e-6:
                        delta = rng.normal(size=2)
                        dist = float(np.linalg.norm(delta))
                    forces[i] += (delta / dist) * 6.0

        if not still_overlap:
            break

        offsets_pt += forces / px_per_pt
        hyp = np.hypot(offsets_pt[:, 0], offsets_pt[:, 1])
        too_far = hyp > max_offset_pt
        if np.any(too_far):
            scale = max_offset_pt / np.maximum(hyp, 1e-9)
            offsets_pt[too_far] *= scale[too_far, None]

    for text, (ox, oy) in zip(texts, offsets_pt):
        text.set_position((float(ox), float(oy)))
    return texts


def plot_scatter_labeled(
    data,
    x_col,
    y_col,
    *,
    hue_col,
    title,
    filename,
    label_col=None,
    hue_order=None,
    hue_labels=None,
    hue_colors=None,
    xlabel=None,
    ylabel=None,
    alpha=0.9,
    s=140,
    figure_size=(14, 10),
    dpi=300,
    log_x=False,
    log_y=False,
    italic_labels=False,
    label_fontsize=11,
    legend_below=True,
    xlim=None,
    ylim=None,
):
    """Scatter plot with categorical colors and optional non-overlapping point labels."""
    import os

    needed = [x_col, y_col, hue_col]
    if label_col:
        needed.append(label_col)
    missing = [c for c in needed if c not in data.columns]
    if missing:
        raise ValueError(f"Missing columns for labeled scatter: {missing}")

    plot_df = data[needed].replace([np.inf, -np.inf], np.nan).dropna().copy()
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=figure_size)

    if hue_order is None:
        categories = list(dict.fromkeys(plot_df[hue_col].tolist()))
    else:
        categories = [c for c in hue_order if (plot_df[hue_col] == c).any()]
        extras = [c for c in dict.fromkeys(plot_df[hue_col].tolist()) if c not in categories]
        categories.extend(extras)

    if hue_colors is None:
        if len(categories) <= 2:
            palette = list(TWO_CATEGORY_COLORS)
            color_map = {cat: palette[i] for i, cat in enumerate(categories)}
        else:
            palette = sns.color_palette("deep", len(categories))
            color_map = {cat: palette[i] for i, cat in enumerate(categories)}
    else:
        color_map = {cat: hue_colors.get(cat, TWO_CATEGORY_BLUE) for cat in categories}

    if hue_labels is None:
        hue_labels = {cat: str(cat) for cat in categories}

    for cat in categories:
        subset = plot_df[plot_df[hue_col] == cat]
        ax.scatter(
            subset[x_col],
            subset[y_col],
            s=s,
            alpha=alpha,
            color=color_map[cat],
            edgecolors="white",
            linewidths=0.6,
            zorder=3,
            label=hue_labels.get(cat, str(cat)),
        )

    ax.set_xlabel(xlabel if xlabel else x_col, fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(ylabel if ylabel else y_col, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if log_x:
        ax.set_xscale("log")
        if xlim is None:
            log_lim = padded_log_axis_limits(plot_df[x_col], fraction=0.06)
            if log_lim is not None:
                xlim = (10 ** log_lim[0], 10 ** log_lim[1])
    elif xlim is None:
        xlim = padded_axis_limits(plot_df[x_col], fraction=0.06)
    if log_y:
        ax.set_yscale("log")
        if ylim is None:
            log_lim = padded_log_axis_limits(plot_df[y_col], fraction=0.06)
            if log_lim is not None:
                ylim = (10 ** log_lim[0], 10 ** log_lim[1])
    elif ylim is None:
        ylim = padded_axis_limits(plot_df[y_col], fraction=0.06)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if label_col:
        place_nonoverlapping_labels(
            ax,
            plot_df[x_col],
            plot_df[y_col],
            plot_df[label_col],
            fontsize=label_fontsize,
            italic=italic_labels,
        )

    if legend_below:
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(0.0, -0.12),
            fontsize=LEGEND_FONT_SIZE,
            frameon=False,
        )
    else:
        ax.legend(loc="best", fontsize=LEGEND_FONT_SIZE, frameon=False)

    parent_dir = os.path.dirname(filename)
    if parent_dir and not os.path.exists(parent_dir):
        os.makedirs(parent_dir, exist_ok=True)
    fig.savefig(filename, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved labeled scatter plot to {filename}")


def plot_scatter_with_outliers(
    data, x_col, y_col, outlier_col,
    title, filename,
    xlabel=None, ylabel=None,
    color_normal=TWO_CATEGORY_BLUE, color_outlier=TWO_CATEGORY_RED,
    alpha=0.5, s=25,
    figure_size=(10,8)
):
    """
    Plots a scatter plot separating normal and outlier points.
    
    Args:
        data: pandas DataFrame
        x_col: string, column name for X axis
        y_col: string, column name for Y axis
        outlier_col: string, boolean column indicating outliers (True=outlier)
        title: string, plot title
        filename: string, output file path
        xlabel: string, label for X axis (defaults to x_col)
        ylabel: string, label for Y axis (defaults to y_col)
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_style("white")
    
    plt.figure(figsize=figure_size)
    
    normal_pts = data[~data[outlier_col]]
    outlier_pts = data[data[outlier_col]]

    plt.scatter(normal_pts[x_col], normal_pts[y_col], s=s, alpha=alpha, label='Normal', color=color_normal)
    plt.scatter(outlier_pts[x_col], outlier_pts[y_col], s=s, alpha=alpha, color=color_outlier, label='Outlier')
    
    plt.xlabel(xlabel if xlabel else x_col, fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel(ylabel if ylabel else y_col, fontsize=Y_LABEL_FONT_SIZE)
    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.legend(loc='best', fontsize=LEGEND_FONT_SIZE)
    
    import os
    parent_dir = os.path.dirname(filename)
    if parent_dir and not os.path.exists(parent_dir):
        os.makedirs(parent_dir, exist_ok=True)
        
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved scatter plot with outliers to {filename}")


def local_xy_density(x, y, *, n_bins=None):
    """
    Per-point local density in standardized x/y space.

    Each point is scored by a smoothed 2D histogram of neighboring points
    (higher = more points in that region). Axes are scaled by a robust MAD
    estimate so a long tail on one axis does not dominate the bins.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.full(x.shape, np.nan, dtype=np.float64)
    finite = np.isfinite(x) & np.isfinite(y)
    xf = x[finite]
    yf = y[finite]
    n = int(xf.size)
    if n == 0:
        return z
    if n == 1:
        z[finite] = 1.0
        return z

    def _scale(arr):
        med = float(np.median(arr))
        mad = float(np.median(np.abs(arr - med)))
        if mad > 0.0:
            return 1.4826 * mad
        std = float(np.std(arr))
        return std if std > 0.0 else 1.0

    xs = (xf - float(np.median(xf))) / _scale(xf)
    ys = (yf - float(np.median(yf))) / _scale(yf)
    bins = int(n_bins) if n_bins is not None else max(12, min(60, int(np.sqrt(n))))
    counts, xedges, yedges = np.histogram2d(xs, ys, bins=bins)
    density = gaussian_filter(counts.astype(np.float64), sigma=1.0)
    ix = np.clip(np.digitize(xs, xedges) - 1, 0, density.shape[0] - 1)
    iy = np.clip(np.digitize(ys, yedges) - 1, 0, density.shape[1] - 1)
    z[finite] = density[ix, iy]
    return z


def plot_scatter_continuous_color(
    data,
    x_col,
    y_col,
    color_col,
    title,
    filename,
    xlabel=None,
    ylabel=None,
    cbar_label=None,
    cmap="coolwarm",
    center=0.0,
    color_pad_fraction=DEFAULT_AXIS_PADDING_FRACTION,
    center_fade_alpha=False,
    center_min_alpha=0.08,
    center_fade_power=0.55,
    alpha=0.55,
    s=22,
    figure_size=(10, 8),
    dpi=300,
    vmin=None,
    vmax=None,
    xlim=None,
    ylim=None,
    aspect=None,
    line_overlays=None,
    thresholds_h=None,
    thresholds_v=None,
    sort_by_color=False,
    log_color=False,
    legend_outside=True,
    bbox_inches="tight",
):
    """
    Scatter plot with point color mapped to a continuous column (diverging cmap by default).

    Color limits default to data min/max plus ``color_pad_fraction`` padding on each side.
    Pass ``vmin`` / ``vmax`` to fix the color scale (e.g. LD r2 from 0 to 1).
    When ``center`` falls inside that range, a two-slope norm keeps the neutral color at
    ``center`` (e.g. zero residual).

    With ``center_fade_alpha=True``, marker alpha increases with distance from ``center`` so
    values near zero appear gray and transparent while extremes stay vivid.

    ``line_overlays`` is an optional list of dicts with ``x`` / ``y`` arrays and optional
    ``color``, ``linestyle``, ``linewidth``, ``label``, and ``zorder``.
    ``thresholds_h`` / ``thresholds_v`` are optional lists of dicts with ``value`` and
    optional ``color``, ``linestyle``, and ``label`` (same shape as
    ``plot_scatter_with_thresholds``).
    When ``sort_by_color`` is True, points are drawn in ascending color order so the
    highest values sit on top.
    ``log_color`` uses a log color scale (``LogNorm``) so a heavy high-density tail
    does not crush the rest of the points to the dark end of the cmap.
    ``bbox_inches=None`` keeps the saved canvas equal to ``figure_size`` (use for a
    fixed aspect ratio).
    """
    import matplotlib.pyplot as plt
    import os

    sns.set_style("white")

    plot_df = data[[x_col, y_col, color_col]].dropna()
    if plot_df.empty:
        print(f"[Warning] plot_scatter_continuous_color: no data for {filename}")
        return
    if sort_by_color:
        plot_df = plot_df.sort_values(color_col)

    cvals = plot_df[color_col].astype(float).to_numpy()
    cmin, cmax = float(cvals.min()), float(cvals.max())
    span = cmax - cmin
    pad = span * color_pad_fraction if span > 0 else max(abs(cmin), abs(cmax), 1e-9) * color_pad_fraction
    if vmin is None:
        vmin = cmin - pad
        if center is None and cmin >= 0:
            vmin = max(0.0, vmin)
    else:
        vmin = float(vmin)
    if vmax is None:
        vmax = cmax + pad
    else:
        vmax = float(vmax)

    if log_color:
        positive = cvals[cvals > 0]
        if vmin is None or vmin <= 0:
            vmin = float(positive.min()) if positive.size else 1e-6
        else:
            vmin = float(vmin)
        vmin = max(vmin, 1e-12)
        if vmax is None:
            vmax = max(cmax, vmin * 1.0001)
        else:
            vmax = float(vmax)
        if vmax <= vmin:
            vmax = vmin * 10.0
        norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
    elif center is not None and vmin < center < vmax:
        norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)
    else:
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    if isinstance(cmap, str):
        colormap = plt.get_cmap(cmap)
    else:
        colormap = cmap

    folder = os.path.dirname(filename)
    if folder:
        os.makedirs(folder, exist_ok=True)

    fig, ax = plt.subplots(
        figsize=figure_size,
        layout="constrained" if bbox_inches is None else None,
    )

    if center_fade_alpha and center is not None:
        neg_span = center - vmin
        pos_span = vmax - center
        frac = np.where(
            cvals >= center,
            np.divide(cvals - center, pos_span, out=np.zeros_like(cvals), where=pos_span > 0),
            np.divide(center - cvals, neg_span, out=np.zeros_like(cvals), where=neg_span > 0),
        )
        frac = np.clip(frac, 0.0, 1.0)
        point_alpha = center_min_alpha + (1.0 - center_min_alpha) * np.power(frac, center_fade_power)
        rgba = colormap(norm(cvals))
        rgba[:, 3] = point_alpha
        ax.scatter(
            plot_df[x_col],
            plot_df[y_col],
            c=rgba,
            s=s,
            linewidths=0,
        )
        sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax)
    else:
        sc = ax.scatter(
            plot_df[x_col],
            plot_df[y_col],
            c=cvals,
            cmap=colormap,
            norm=norm,
            s=s,
            alpha=alpha,
            linewidths=0,
        )
        cbar = fig.colorbar(sc, ax=ax)

    cbar.set_label(
        cbar_label if cbar_label else color_col,
        fontsize=Y_LABEL_FONT_SIZE,
    )
    cbar.ax.tick_params(labelsize=TICK_FONT_SIZE)

    ax.set_xlabel(xlabel if xlabel else x_col, fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(ylabel if ylabel else y_col, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if aspect is not None:
        ax.set_aspect(aspect)

    if thresholds_h:
        for thr in thresholds_h:
            val = thr.get("value")
            if val is None:
                continue
            ax.axhline(
                y=val,
                color=thr.get("color", "red"),
                linestyle=thr.get("linestyle", "--"),
                label=thr.get("label", f"{val:.4f}"),
                zorder=4,
            )
    if thresholds_v:
        for thr in thresholds_v:
            val = thr.get("value")
            if val is None:
                continue
            ax.axvline(
                x=val,
                color=thr.get("color", "orange"),
                linestyle=thr.get("linestyle", "--"),
                label=thr.get("label", f"{val:.4f}"),
                zorder=4,
            )

    if line_overlays:
        for spec in line_overlays:
            plot_kw = {
                "color": spec.get("color", "black"),
                "linestyle": spec.get("linestyle", "-"),
                "linewidth": spec.get("linewidth", 1.2),
                "zorder": spec.get("zorder", 5),
            }
            if spec.get("label"):
                plot_kw["label"] = spec["label"]
            ax.plot(spec["x"], spec["y"], **plot_kw)

    _handles, labels = ax.get_legend_handles_labels()
    if labels:
        legend_kw = dict(fontsize=LEGEND_FONT_SIZE, frameon=False)
        if legend_outside:
            if bbox_inches is None:
                fig.legend(_handles, labels, loc="outside lower left", **legend_kw)
            else:
                ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.15), **legend_kw)
        else:
            ax.legend(loc="lower right", **legend_kw)

    save_kw = {"dpi": dpi}
    if bbox_inches is not None:
        save_kw["bbox_inches"] = bbox_inches
    fig.savefig(filename, **save_kw)
    plt.close(fig)
    print(f"Saved continuous-color scatter to {filename}")


def plot_scatter_density(
    data,
    x_col,
    y_col,
    title,
    filename,
    xlabel=None,
    ylabel=None,
    cbar_label="Local point density",
    cmap="inferno",
    alpha=0.9,
    s=56,
    figure_size=(12, 6),
    thresholds_h=None,
    thresholds_v=None,
    line_overlays=None,
    xlim=None,
    ylim=None,
    aspect=None,
    dpi=300,
):
    """
    Scatter plot colored as a density heatmap: brighter points have more neighbors.

    Canvas default is width:height = 2:1 and is kept on save (no tight crop).
    Color uses a log scale so the dense core does not hide mid-density structure.
    Dense regions are drawn last so they sit on top of sparse outliers.
    Threshold legends sit below the axes (same placement as other stats plots).
    """
    plot_df = data[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna().copy()
    plot_df["_local_density"] = local_xy_density(plot_df[x_col], plot_df[y_col])
    plot_scatter_continuous_color(
        plot_df,
        x_col,
        y_col,
        "_local_density",
        title,
        filename,
        xlabel=xlabel,
        ylabel=ylabel,
        cbar_label=cbar_label,
        cmap=cmap,
        center=None,
        color_pad_fraction=0.0,
        vmin=None,
        log_color=True,
        alpha=alpha,
        s=s,
        figure_size=figure_size,
        dpi=dpi,
        xlim=xlim,
        ylim=ylim,
        aspect=aspect,
        thresholds_h=thresholds_h,
        thresholds_v=thresholds_v,
        line_overlays=line_overlays,
        sort_by_color=True,
        legend_outside=True,
        bbox_inches=None,
    )


def plot_scatter_with_outliers_sized(
    data,
    x_col,
    y_col,
    outlier_col,
    size_col,
    title,
    filename,
    xlabel=None,
    ylabel=None,
    color_normal=TWO_CATEGORY_BLUE,
    color_outlier=TWO_CATEGORY_RED,
    alpha=0.45,
    size_range=(10, 140),
    figure_size=(10, 8),
    dpi=300,
):
    """
    Scatter with normal/outlier colors (like ``plot_scatter_with_outliers``) and marker
    size scaled by ``size_col`` (typically |residual|).
    """
    import matplotlib.pyplot as plt
    import os

    sns.set_style("white")

    plot_df = data[[x_col, y_col, outlier_col, size_col]].dropna()
    if plot_df.empty:
        print(f"[Warning] plot_scatter_with_outliers_sized: no data for {filename}")
        return

    sizes_raw = plot_df[size_col].astype(float)
    lo, hi = float(sizes_raw.min()), float(sizes_raw.max())
    if hi <= lo:
        point_sizes = np.full(len(plot_df), np.mean(size_range))
    else:
        scale = (sizes_raw - lo) / (hi - lo)
        point_sizes = size_range[0] + scale * (size_range[1] - size_range[0])

    folder = os.path.dirname(filename)
    if folder:
        os.makedirs(folder, exist_ok=True)

    fig, ax = plt.subplots(figsize=figure_size)
    for is_out, color, label in (
        (False, color_normal, "Normal"),
        (True, color_outlier, "Outlier"),
    ):
        mask = plot_df[outlier_col].astype(bool) == is_out
        if not mask.any():
            continue
        ax.scatter(
            plot_df.loc[mask, x_col],
            plot_df.loc[mask, y_col],
            s=point_sizes[mask],
            alpha=alpha,
            color=color,
            label=label,
            linewidths=0,
        )

    ax.set_xlabel(xlabel if xlabel else x_col, fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(ylabel if ylabel else y_col, fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.legend(loc="best", fontsize=LEGEND_FONT_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)

    fig.savefig(filename, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved outlier-sized scatter to {filename}")


def _consecutive_ref_name_blocks(data, ref_name_col='ref_name', x_col='bin_index'):
    """Merge consecutive rows with the same ref_name into axis annotation blocks."""
    blocks = []
    if data is None or data.empty:
        return blocks
    cur_ref = data.iloc[0][ref_name_col]
    start = data.iloc[0][x_col]
    prev_x = start
    for _, row in data.iloc[1:].iterrows():
        x = row[x_col]
        ref = row[ref_name_col]
        if ref != cur_ref:
            blocks.append({'ref_name': cur_ref, 'start_bin': start, 'end_bin': prev_x})
            cur_ref = ref
            start = x
        prev_x = x
    blocks.append({'ref_name': cur_ref, 'start_bin': start, 'end_bin': prev_x})
    return blocks


def _decorate_genome_ref_block_axis(ax, blocks, n_bins, show_labels=True, block_alpha=0.06):
    """Shade ref-chromosome blocks, draw boundaries, and label ref names (no numeric x ticks)."""
    for i, block in enumerate(blocks):
        left = block['start_bin'] - 0.5
        right = block['end_bin'] + 0.5
        if i % 2 == 0:
            ax.axvspan(left, right, color='gray', alpha=block_alpha, linewidth=0, zorder=0)
        if i > 0:
            ax.axvline(left, color='0.75', linewidth=0.6, zorder=0)
    if show_labels:
        tick_pos = [(b['start_bin'] + b['end_bin']) / 2.0 for b in blocks]
        tick_lbl = [b['ref_name'] for b in blocks]
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(tick_lbl, rotation=45, ha='right', fontsize=TICK_FONT_SIZE)
    else:
        ax.set_xticklabels([])
    ax.set_xlim(-0.5, max(n_bins - 0.5, 0.5))


def _shade_centromere_bin_spans(ax, spans, alpha=0.18, color='#2ca02c', zorder=1):
    """Highlight centromere intervals on a binned genome axis."""
    for span in spans:
        left = span['start_bin'] - 0.5
        right = span['end_bin'] + 0.5
        ax.axvspan(left, right, color=color, alpha=alpha, linewidth=0, zorder=zorder)


def plot_genome_binned_density_panels(
    data,
    x_col,
    ref_name_col,
    panel_specs,
    filename,
    suptitle=None,
    x_label='Chromosome',
    y_label='Variants per Mb',
    bin_size_mb=5,
    figure_size=(18, 8),
    block_alpha=0.06,
    linewidth=0.8,
    show_centromere=True,
    centromere_spans=None,
    centromere_alpha=0.18,
):
    """
    Stack one line panel per series along a binned genome axis with ref-chromosome blocks.

    ``data`` must contain ``x_col`` (bin index), ``ref_name_col`` (e.g. chr1A), and one
    numeric column per entry in ``panel_specs`` (keys: ``y_col``, ``title``, optional ``color``).
    X-axis shows ref names at block centres; numeric bin / bp positions are not drawn.

    When ``show_centromere`` is True, vmap4 centromere intervals are shaded in green
    (auto-loaded from ``infra.wheat.ref_v1`` unless ``centromere_spans`` is supplied).
    """
    import os
    from matplotlib.patches import Patch

    sns.set_style('white')
    blocks = _consecutive_ref_name_blocks(data, ref_name_col=ref_name_col, x_col=x_col)
    n_bins = int(data[x_col].max()) + 1 if not data.empty else 0

    if show_centromere and centromere_spans is None:
        from infra.wheat.ref_v1 import get_ref_v1_centromere_bin_spans

        centromere_spans = get_ref_v1_centromere_bin_spans(int(bin_size_mb * 1_000_000))

    y_max = 1.0
    for spec in panel_specs:
        y_max = max(y_max, float(data[spec['y_col']].max()))
    y_max *= 1.05

    n_panels = len(panel_specs)
    fig, axes = plt.subplots(n_panels, 1, figsize=figure_size, sharex=True)
    if n_panels == 1:
        axes = [axes]

    centromere_color = '#2ca02c'
    for i, (ax, spec) in enumerate(zip(axes, panel_specs)):
        y_col = spec['y_col']
        color = spec.get('color', 'steelblue')
        _decorate_genome_ref_block_axis(
            ax,
            blocks,
            n_bins,
            show_labels=(i == n_panels - 1),
            block_alpha=block_alpha,
        )
        if show_centromere and centromere_spans:
            _shade_centromere_bin_spans(
                ax,
                centromere_spans,
                alpha=centromere_alpha,
                color=centromere_color,
            )
        ax.plot(data[x_col], data[y_col], color=color, linewidth=linewidth, zorder=2)
        ax.set_ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
        ax.set_title(spec.get('title', y_col), fontsize=TITLE_FONT_SIZE)
        ax.set_ylim(0, y_max)
        ax.tick_params(axis='y', which='major', labelsize=TICK_FONT_SIZE)
        if show_centromere and centromere_spans and i == n_panels - 1:
            centromere_patch = Patch(
                facecolor=centromere_color,
                alpha=centromere_alpha,
                edgecolor='none',
                label='Centromere',
            )
            ax.legend(
                handles=[centromere_patch],
                loc='upper right',
                fontsize=LEGEND_FONT_SIZE,
                framealpha=0.85,
            )

    axes[-1].set_xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)

    if suptitle:
        fig.suptitle(suptitle, fontsize=TITLE_FONT_SIZE, y=1.01)
    fig.tight_layout()

    parent_dir = os.path.dirname(filename)
    if parent_dir:
        os.makedirs(parent_dir, exist_ok=True)
    fig.savefig(filename, dpi=300, bbox_inches='tight')
    print(f'Saved genome density panels to {filename}')
    plt.close(fig)


def plot_circos_chromosome_tracks(
    chrom_lengths,
    tracks,
    filename,
    chrom_order=None,
    title=None,
    gap_frac=0.008,
    wrap_gap_frac=0.0,
    figure_size=(12, 12),
    dpi=300,
    tick_bp=200_000_000,
    legend_fontsize=None,
    legend_loc="center",
    chrom_label_pad=0.05,
):
    """Polar karyotype with concentric histogram, heatmap, and interval tracks.

    ``chrom_lengths`` maps chromosome name to length in bp. Each track is a
    dict with ``kind`` ``histogram`` (columns chrom / pos / value), ``heatmap``
    (chrom / bin_start / bin_end / value), or ``karyotype`` (optional
    ``intervals`` with chrom / start_pos / end_pos, and optional per-row
    ``color`` / ``label``), plus ``color``, ``r_inner``, ``r_outer``, and
    ``label``. ``legend_loc`` is ``center`` (hole), ``outside`` (below the
    ring), ``wrap_gap`` (labels in the 7D–1A gap at each track radius), or
    ``none``. ``wrap_gap_frac`` is an extra fraction of the circle reserved
    between the last and first chromosome. ``chrom_label_pad`` is radial
    space beyond the outermost track for chromosome names.
    ``legend_fontsize`` overrides the module legend size.
    """
    import os
    from matplotlib.patches import Patch

    sns.set_style("white")
    if not chrom_lengths:
        raise ValueError("chrom_lengths is empty")
    order = list(chrom_order) if chrom_order else list(chrom_lengths.keys())
    order = [c for c in order if c in chrom_lengths and float(chrom_lengths[c]) > 0]
    if not order:
        raise ValueError("No chromosomes with positive length to plot")

    total_len = float(sum(float(chrom_lengths[c]) for c in order))
    n_chrom = len(order)
    inter_gap = 2.0 * np.pi * float(gap_frac)
    wrap_gap = 2.0 * np.pi * float(wrap_gap_frac)
    if wrap_gap > 0:
        usable = max(2.0 * np.pi - inter_gap * max(n_chrom - 1, 0) - wrap_gap, np.pi)
        cursor = wrap_gap / 2.0
    else:
        usable = max(2.0 * np.pi - inter_gap * n_chrom, np.pi)
        cursor = 0.0
    theta_start = {}
    theta_span = {}
    for i, chrom in enumerate(order):
        span = usable * (float(chrom_lengths[chrom]) / total_len)
        theta_start[chrom] = cursor
        theta_span[chrom] = span
        cursor += span
        if wrap_gap > 0:
            if i < n_chrom - 1:
                cursor += inter_gap
        else:
            cursor += inter_gap

    def pos_to_theta(chrom, pos):
        length = float(chrom_lengths[chrom])
        frac = np.clip(float(pos) / length, 0.0, 1.0) if length > 0 else 0.0
        return theta_start[chrom] + theta_span[chrom] * frac

    max_outer = max((float(track["r_outer"]) for track in tracks), default=1.0)
    label_r = max_outer + float(chrom_label_pad)
    fig = plt.figure(figsize=figure_size)
    ax = fig.add_subplot(111, polar=True)
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2.0)
    ax.set_ylim(0.0, label_r + 0.06)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines["polar"].set_visible(False)
    ax.grid(False)

    legend_handles = []
    heatmap_mappables = []
    gap_labels = []
    min_theta_width = 0.003
    for track in tracks:
        kind = track.get("kind")
        r0 = float(track["r_inner"])
        r1 = float(track["r_outer"])
        color = track.get("color", "steelblue")
        label = track.get("label")
        if kind == "histogram":
            data = track.get("data")
            if data is None or data.empty:
                continue
            values = pd.to_numeric(data["value"], errors="coerce").fillna(0.0)
            vmax = float(values.max()) if len(values) else 1.0
            if vmax <= 0:
                vmax = 1.0
            height_scale = (r1 - r0) / vmax
            has_bin = "bin_start" in data.columns and "bin_end" in data.columns
            for row in data.to_dict("records"):
                chrom = str(row["chrom"])
                if chrom not in theta_start:
                    continue
                value = float(row["value"])
                if value <= 0:
                    continue
                width = max(theta_span[chrom] * 0.004, min_theta_width)
                if has_bin:
                    t0 = pos_to_theta(chrom, row["bin_start"])
                    t1 = pos_to_theta(chrom, row["bin_end"])
                    width = max(abs(t1 - t0), min_theta_width)
                    theta = (t0 + t1) / 2.0
                else:
                    theta = pos_to_theta(chrom, row["pos"])
                ax.bar(
                    theta,
                    value * height_scale,
                    width=width,
                    bottom=r0,
                    color=color,
                    edgecolor="none",
                    linewidth=0,
                    align="center",
                )
            if label:
                legend_handles.append(Patch(facecolor=color, edgecolor="none", label=label))
                gap_labels.append(
                    {"theta": 0.0, "r": (r0 + r1) / 2.0, "text": label, "color": color}
                )
        elif kind == "heatmap":
            data = track.get("data")
            if data is None or data.empty:
                continue
            values = pd.to_numeric(data["value"], errors="coerce").fillna(0.0)
            vmin = track.get("vmin")
            vmax = track.get("vmax")
            lo = float(values.min()) if vmin is None else float(vmin)
            hi = float(values.max()) if vmax is None else float(vmax)
            if hi <= lo:
                hi = lo + 1.0
            cmap = plt.get_cmap(track.get("cmap", "Blues"))
            norm = mcolors.Normalize(vmin=lo, vmax=hi, clip=True)
            for row in data.to_dict("records"):
                chrom = str(row["chrom"])
                if chrom not in theta_start:
                    continue
                t0 = pos_to_theta(chrom, row["bin_start"])
                t1 = pos_to_theta(chrom, row["bin_end"])
                width = max(abs(t1 - t0), min_theta_width)
                ax.bar(
                    (t0 + t1) / 2.0,
                    r1 - r0,
                    width=width,
                    bottom=r0,
                    color=cmap(norm(float(row["value"]))),
                    edgecolor="none",
                    linewidth=0,
                    align="center",
                )
            heatmap_mappables.append(
                (
                    mcm.ScalarMappable(norm=norm, cmap=cmap),
                    track.get("cbar_label") or label or "Genes",
                )
            )
        elif kind == "karyotype":
            base_color = track.get("base_color", "#1a1a1a")
            interval_color = track.get("interval_color", color)
            show_legend = track.get("show_legend", True)
            for chrom in order:
                mid = theta_start[chrom] + theta_span[chrom] / 2.0
                ax.bar(
                    mid,
                    r1 - r0,
                    width=max(theta_span[chrom], min_theta_width),
                    bottom=r0,
                    color=base_color,
                    edgecolor="none",
                    linewidth=0,
                    align="center",
                    zorder=1,
                )
            intervals = track.get("intervals")
            seen_interval_labels = set()
            track_gap_items = []
            if intervals is not None and not intervals.empty:
                for row in intervals.to_dict("records"):
                    chrom = str(row["chrom"])
                    if chrom not in theta_start:
                        continue
                    t0 = pos_to_theta(chrom, row["start_pos"])
                    t1 = pos_to_theta(chrom, row["end_pos"])
                    width = abs(t1 - t0)
                    if width <= 0:
                        continue
                    raw_color = row.get("color")
                    row_color = raw_color if isinstance(raw_color, str) and raw_color else interval_color
                    ax.bar(
                        (t0 + t1) / 2.0,
                        r1 - r0,
                        width=width,
                        bottom=r0,
                        color=row_color,
                        edgecolor="none",
                        linewidth=0,
                        align="center",
                        zorder=2,
                    )
                    row_label = row.get("label")
                    if not isinstance(row_label, str) or not row_label:
                        row_label = None
                    if show_legend and row_label and row_label not in seen_interval_labels:
                        legend_handles.append(
                            Patch(facecolor=row_color, edgecolor="none", label=row_label)
                        )
                        track_gap_items.append((row_label, row_color))
                        seen_interval_labels.add(row_label)
            if show_legend and label and not seen_interval_labels:
                legend_handles.append(
                    Patch(facecolor=interval_color, edgecolor="none", label=label)
                )
                gap_labels.append(
                    {
                        "theta": 0.0,
                        "r": (r0 + r1) / 2.0,
                        "text": label,
                        "color": interval_color,
                    }
                )
            elif track_gap_items:
                for i, (name, face) in enumerate(track_gap_items):
                    gap_labels.append(
                        {
                            "theta": 0.0,
                            "r": r1 + 0.028 + i * 0.040,
                            "text": name,
                            "color": face,
                        }
                    )
            if track.get("show_ticks", True):
                tick_r = r1 + 0.012
                for chrom in order:
                    length = float(chrom_lengths[chrom])
                    tick = 0.0
                    while tick <= length + 1.0:
                        theta = pos_to_theta(chrom, tick)
                        ax.plot(
                            [theta, theta],
                            [r1, r1 + 0.010],
                            color="#333333",
                            linewidth=0.4,
                            zorder=3,
                        )
                        if tick > 0:
                            ax.text(
                                theta,
                                tick_r + 0.012,
                                f"{int(tick / 1_000_000)}",
                                ha="center",
                                va="center",
                                fontsize=max(TICK_FONT_SIZE - 5, 6),
                                color="#333333",
                                rotation=0,
                            )
                        tick += float(tick_bp)

    for chrom in order:
        mid = theta_start[chrom] + theta_span[chrom] / 2.0
        ax.text(
            mid,
            label_r,
            chrom,
            ha="center",
            va="center",
            fontsize=X_LABEL_FONT_SIZE,
            fontweight="bold",
        )

    if title:
        ax.set_title(title, fontsize=TITLE_FONT_SIZE, pad=18)
    legend_kw = dict(
        handles=legend_handles,
        fontsize=LEGEND_FONT_SIZE if legend_fontsize is None else legend_fontsize,
        frameon=False,
        handlelength=1.2,
        handletextpad=0.4,
        labelspacing=0.35,
        borderpad=0.2,
    )
    if legend_handles and legend_loc == "center":
        ax.legend(
            loc="center",
            bbox_to_anchor=(0.5, 0.5),
            bbox_transform=ax.transAxes,
            ncol=1,
            **legend_kw,
        )
    elif legend_handles and legend_loc == "outside":
        fig.legend(
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            ncol=1,
            **legend_kw,
        )
    elif legend_loc == "wrap_gap" and gap_labels:
        gap_fs = LEGEND_FONT_SIZE if legend_fontsize is None else legend_fontsize
        for item in gap_labels:
            face = item.get("color", "#333333")
            ax.plot(
                item["theta"] - 0.055,
                item["r"],
                marker="s",
                markersize=5.5,
                color=face,
                markeredgecolor="#333333",
                markeredgewidth=0.3,
                zorder=5,
            )
            ax.text(
                item["theta"] + 0.01,
                item["r"],
                item["text"],
                ha="left",
                va="center",
                fontsize=gap_fs,
                color="#222222",
                zorder=6,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 0.12},
            )
    for mappable, cbar_label in heatmap_mappables:
        cbar = fig.colorbar(mappable, ax=ax, pad=0.08, fraction=0.035, shrink=0.45)
        cbar.set_label(cbar_label, fontsize=Y_LABEL_FONT_SIZE)
        cbar.ax.tick_params(labelsize=TICK_FONT_SIZE)

    parent_dir = os.path.dirname(filename)
    if parent_dir:
        os.makedirs(parent_dir, exist_ok=True)
    fig.savefig(filename, dpi=dpi, bbox_inches="tight")
    print(f"Saved Circos chromosome tracks to {filename}")
    plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) >= 4:
        combine_plots([sys.argv[1], sys.argv[2]], sys.argv[3])
    else:
        print("Usage: python graph.py <image1> <image2> <output_file>")
    
