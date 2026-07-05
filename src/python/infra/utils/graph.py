import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as mcolors
import seaborn as sns
import numpy as np
import pandas as pd
import sys
from sklearn.linear_model import HuberRegressor, LinearRegression

TITLE_FONT_SIZE = 18
X_LABEL_FONT_SIZE = 14
Y_LABEL_FONT_SIZE = 14
TICK_FONT_SIZE = 12
LEGEND_FONT_SIZE = 12

# Missing-rate heatmaps: red = low missing (better), blue = high missing (worse).
MISSING_RATE_CMAP = "coolwarm_r"

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
    ``label``, ``color``, ``linestyle``, ``linewidth`` (defaults match common decay-style plots).
    """
    sns.set_style("white")
    plt.figure(figsize=figure_size)

    if x_label is None:
        x_label = x_col
    if y_label is None:
        y_label = ""

    for spec in y_specs:
        yc = spec["y_col"]
        lbl = spec.get("label", yc)
        clr = spec.get("color", "steelblue")
        ls = spec.get("linestyle", "-")
        lw = float(spec.get("linewidth", 1.5))
        plt.plot(data[x_col], data[yc], color=clr, linestyle=ls, linewidth=lw, label=lbl)

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
        plt.axvline(x=float(mean_val), color='red', linestyle='--', linewidth=1.5, label=label_text)

    if _finite_for_vline(median_val):
        plt.axvline(
            x=float(median_val),
            color='orange',
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
    plt.plot(x_range, y_ols, color='blue', linewidth=2, linestyle='--', 
            label=f'OLS: {ols_eq}, $R^2$={ols_score:.3f}')
    plt.plot(x_range, y_huber, color='red', linewidth=2, 
            label=f'Huber: {huber_eq}, $R^2$={huber_score:.3f}')
    
    if title:
        plt.title(title, fontsize=TITLE_FONT_SIZE)
    else:
        plt.title(f'{y_label} vs {x_label}', fontsize=TITLE_FONT_SIZE)

    plt.xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    plt.tick_params(axis='both', which='major', labelsize=TICK_FONT_SIZE)
    if x_lim is not None:
        plt.xlim(x_lim)
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
    hue_kw = _hue_kwargs_for_sample_group(df, group_col)
    ax = sns.histplot(
        data=df,
        x=col,
        hue=group_col,
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
        l1 = ax.axvline(x=mean_val, color='blue', linestyle='--', linewidth=1, label=f'Mean: {mean_val:.4f}')
        handles.append(l1)
        labels.append(l1.get_label())
        
    if median_val is not None:
        l2 = ax.axvline(x=median_val, color='orange', linestyle=':', linewidth=1.5, label=f'Median: {median_val:.4f}')
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
    figure_size=(10,6)
):
    """
    Plots a JointGrid (Scatter + Marginal Histograms) with OLS and Huber regression lines.
    Supports grouping/hue in scatter and histograms.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn.linear_model import HuberRegressor, LinearRegression
    
    sns.set_style("white")

    print(f"Running regressions for {y_label} vs {x_label}...")
    
    # Filter valid data
    local_valid = df[[x_col, y_col, group_col]].replace([np.inf, -np.inf], np.nan).dropna()
    
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
    
    # Plotting: JointGrid
    g = sns.JointGrid(data=local_valid, x=x_col, y=y_col, height=figure_size[1], ratio=5)
    hue_kw = _hue_kwargs_for_sample_group(local_valid, group_col)

    # 1. Main Scatter Plot
    sns.scatterplot(data=local_valid, x=x_col, y=y_col, hue=group_col,
                    alpha=0.6, s=20, edgecolor='none', ax=g.ax_joint, legend='full', **hue_kw)
    
    # 2. Add Regression Lines
    x_min, x_max = local_valid[x_col].min(), local_valid[x_col].max()
    x_span = x_max - x_min
    if x_span == 0: x_span = 0.01
    x_range = np.linspace(x_min - 0.05*x_span, x_max + 0.05*x_span, 100).reshape(-1, 1)
    
    y_ols = ols.predict(x_range)
    y_huber = huber.predict(x_range)
    
    line1, = g.ax_joint.plot(x_range, y_ols, color='blue', linewidth=2, linestyle='--', 
                    label=f'OLS: {ols_eq}\n$R^2$={ols_score:.3f}')
    line2, = g.ax_joint.plot(x_range, y_huber, color='red', linewidth=2, 
                    label=f'Huber: {huber_eq}\n$R^2$={huber_score:.3f}')
    
    # --- Legend Management ---
    
    # Remove default legend inside (if any)
    if g.ax_joint.get_legend():
        g.ax_joint.get_legend().remove()

    # Retrieve all handles (Groups + Regression Lines)
    final_handles, final_labels = g.ax_joint.get_legend_handles_labels()
        
    # Create unified legend below left
    g.ax_joint.legend(handles=final_handles, labels=final_labels, title=f"{group_col} & Regression",
                      loc='upper left', bbox_to_anchor=(0.0, -0.15), 
                      fontsize=LEGEND_FONT_SIZE, title_fontsize=LEGEND_FONT_SIZE, frameon=False)
    
    # 3. Marginals (Stacked Histograms)
    sns.histplot(data=local_valid, x=x_col, hue=group_col, multiple='stack', 
                 bins=100, linewidth=0.1, ax=g.ax_marg_x, legend=False, **hue_kw)
    sns.histplot(data=local_valid, y=y_col, hue=group_col, multiple='stack', 
                 bins=50, linewidth=0.1, ax=g.ax_marg_y, legend=False, **hue_kw)

    g.ax_joint.set_xlabel(x_label, fontsize=X_LABEL_FONT_SIZE)
    g.ax_joint.set_ylabel(y_label, fontsize=Y_LABEL_FONT_SIZE)
    g.ax_joint.tick_params(axis='both', which='major', labelsize=TICK_FONT_SIZE)

    if x_lim is not None:
        g.ax_joint.set_xlim(x_lim)
    else:
        x_lim = padded_axis_limits(local_valid[x_col])
        if x_lim is not None:
            g.ax_joint.set_xlim(x_lim)
    if y_lim is not None:
        g.ax_joint.set_ylim(y_lim)
    else:
        y_lim = padded_axis_limits(local_valid[y_col])
        if y_lim is not None:
            g.ax_joint.set_ylim(y_lim)

    g.ax_marg_x.set_xlim(g.ax_joint.get_xlim())
    g.ax_marg_y.set_ylim(g.ax_joint.get_ylim())
    
    if title:
        plt.suptitle(title, y=0.98, fontsize=TITLE_FONT_SIZE)
    else:
        plt.suptitle(f'{y_label} vs {x_label}', y=0.98, fontsize=TITLE_FONT_SIZE)
    
    g.fig.subplots_adjust(top=0.9)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
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
    plt.plot(x_range, y_ols, color='blue', linewidth=2, linestyle='--', 
             label=f'OLS: {ols_eq}, $R^2$={ols_score:.3f}')
    # Huber Line
    plt.plot(x_range, y_huber, color='red', linewidth=2, 
             label=f'Huber: {huber_eq}, $R^2$={huber_score:.3f}')
    
    if "Heterozygosity" in y_label and "Frequency" in x_label:
         y_hwe = 2 * x_range * (1 - x_range)
         plt.plot(x_range, y_hwe, color='green', linewidth=2, linestyle='-.', 
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
):
    """
    Plots a heatmap from a matrix.
    """
    plt.figure(figsize=(14, 10))
    ax = sns.heatmap(data_matrix, cmap=cmap, 
                     cbar_kws={'label': cbar_label})
    
    # Set labels
    if xlabel: ax.set_xlabel(xlabel, fontsize=X_LABEL_FONT_SIZE)
    if ylabel: ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    
    # Adjust ticks
    n_x = len(x_labels)
    step_x = x_tick_step if x_tick_step is not None else max(1, n_x // 10)
    xticks = np.arange(0, n_x, step_x)
    ax.set_xticks(xticks + 0.5)
    ax.set_xticklabels([x_labels[i] for i in xticks], rotation=0)
    
    n_y = len(y_labels)
    step_y = max(1, n_y // 10)
    yticks = np.arange(0, n_y, step_y)
    ax.set_yticks(yticks + 0.5)
    ax.set_yticklabels([y_labels[i] for i in yticks], rotation=0)

    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.savefig(filename, dpi=300)
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
    xlim=None, log_scale=False
):
    """
    Plots a stacked histogram using seaborn.
    """
    plt.figure(figsize=(10, 6))
    
    # Bins logic
    bins = 100
    binrange = xlim if xlim else None

    sns.histplot(data=data, x=x_col, hue=hue_col, hue_order=hue_order,
                 bins=bins, binrange=binrange, multiple='stack',
                 palette=palette_map, edgecolor='white', alpha=0.8, kde=False)

    if mean_val is not None:
        plt.axvline(mean_val, color='red', linestyle='--', label=f'Mean: {mean_val:.4f}')
    if median_val is not None:
        plt.axvline(median_val, color='green', linestyle='--', label=f'Median: {median_val:.4f}')

    plt.legend(loc='upper right', title='Categories')

    if log_scale:
        plt.yscale('log')
        plt.ylabel(f"{ylabel} (Log Scale)", fontsize=Y_LABEL_FONT_SIZE)
    else:
        plt.ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)

    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.xlabel(xlabel, fontsize=X_LABEL_FONT_SIZE)
    if xlim:
        plt.xlim(xlim)
    plt.grid(True, which="major", axis='y', alpha=0.3)
    # plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved {filename}")
    plt.close()



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
):
    """
    Plots a scatter plot with optional horizontal and vertical threshold lines.
    
    Args:
        thresholds_h (list of dict): e.g. [{'value': 0.95, 'color': 'red', 'linestyle': '--', 'label': 'Cutoff'}]
        thresholds_v (list of dict): e.g. [{'value': 0.5, 'color': 'orange', 'linestyle': '--', 'label': 'Limit'}]
        group_col: Optional column name for hue (e.g. ``Group`` from ``anno_group``).
    """
    sns.set_style("white")
    
    plt.figure(figsize=figure_size)

    use_hue = group_col and group_col in data.columns
    if use_hue:
        plot_data = data[[x_col, y_col, group_col]].replace([np.inf, -np.inf], np.nan).dropna()
        hue_kw = _hue_kwargs_for_sample_group(plot_data, group_col)
        sns.scatterplot(
            data=plot_data, x=x_col, y=y_col,
            hue=group_col, alpha=alpha, s=s, edgecolor='w', **hue_kw,
        )
    else:
        plot_data = data
        sns.scatterplot(
            data=plot_data, x=x_col, y=y_col,
            alpha=alpha, s=s, edgecolor='w', color=color,
        )
    
    # Horizontal Thresholds
    if thresholds_h:
        for thr in thresholds_h:
            val = thr.get('value')
            if val is None: continue
            lbl = thr.get('label', f'{val:.4f}')
            clr = thr.get('color', 'red')
            ls = thr.get('linestyle', '--')
            plt.axhline(y=val, color=clr, linestyle=ls, label=lbl)

    # Vertical Thresholds
    if thresholds_v:
        for thr in thresholds_v:
            val = thr.get('value')
            if val is None: continue
            lbl = thr.get('label', f'{val:.4f}')
            clr = thr.get('color', 'orange')
            ls = thr.get('linestyle', '--')
            plt.axvline(x=val, color=clr, linestyle=ls, label=lbl)

    if xlabel is None: xlabel = x_col
    if ylabel is None: ylabel = y_col
    
    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.xlabel(xlabel, fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    plt.tick_params(axis='both', which='major', labelsize=TICK_FONT_SIZE)
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)

    if use_hue or thresholds_h or thresholds_v:
        plt.legend(
            loc='upper left', bbox_to_anchor=(0.0, -0.15),
            fontsize=LEGEND_FONT_SIZE, frameon=False,
            title=group_col if use_hue else None,
        )
    
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved scatter plot to {filename}")
    plt.close()


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


def plot_scatter_with_outliers(
    data, x_col, y_col, outlier_col,
    title, filename,
    xlabel=None, ylabel=None,
    color_normal='royalblue', color_outlier='red',
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
):
    """
    Scatter plot with point color mapped to a continuous column (diverging cmap by default).

    Color limits default to data min/max plus ``color_pad_fraction`` padding on each side.
    When ``center`` falls inside that range, a two-slope norm keeps the neutral color at
    ``center`` (e.g. zero residual).

    With ``center_fade_alpha=True``, marker alpha increases with distance from ``center`` so
    values near zero appear gray and transparent while extremes stay vivid.
    """
    import matplotlib.pyplot as plt
    import os

    sns.set_style("white")

    plot_df = data[[x_col, y_col, color_col]].dropna()
    if plot_df.empty:
        print(f"[Warning] plot_scatter_continuous_color: no data for {filename}")
        return

    cvals = plot_df[color_col].astype(float).to_numpy()
    cmin, cmax = float(cvals.min()), float(cvals.max())
    span = cmax - cmin
    pad = span * color_pad_fraction if span > 0 else max(abs(cmin), abs(cmax), 1e-9) * color_pad_fraction
    vmin = cmin - pad
    vmax = cmax + pad

    if center is not None and vmin < center < vmax:
        norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)
    else:
        if center is None and cmin >= 0:
            vmin = max(0.0, cmin - pad)
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    if isinstance(cmap, str):
        colormap = plt.get_cmap(cmap)
    else:
        colormap = cmap

    folder = os.path.dirname(filename)
    if folder:
        os.makedirs(folder, exist_ok=True)

    fig, ax = plt.subplots(figsize=figure_size)

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

    fig.savefig(filename, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved continuous-color scatter to {filename}")


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
    color_normal="royalblue",
    color_outlier="red",
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


if __name__ == "__main__":
    if len(sys.argv) >= 4:
        combine_plots([sys.argv[1], sys.argv[2]], sys.argv[3])
    else:
        print("Usage: python graph.py <image1> <image2> <output_file>")
    
