import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
import numpy as np
import pandas as pd
import sys


def combine_plots(
    images,
    output_file = "/data/home/tusr1/git/script/out/combined_plot.png",
    orientation = 'horizontal'
):
    """
    Combines multiple images into a single figure.
    orientation: 'horizontal' or 'vertical'
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
    if orientation == 'horizontal':
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
    xlim=None
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
    
    if x_label is None:
        x_label = col

    plt.figure(figsize=(8, 6))
    
    # Determine bin edges if xlim is provided to ensure fixed bin count within the visible range
    plot_bins = bins
    if xlim is not None and isinstance(bins, int):
        plot_bins = np.linspace(xlim[0], xlim[1], bins + 1)

    # Histogram
    sns.histplot(data=data, x=col, bins=plot_bins, color=color, edgecolor='white', linewidth=0.1)
    
    # Statistics Lines
    if mean_val is not None:
        label_text = f'Mean: {mean_val:.4f}'
        if std_val is not None:
            label_text += f'\nSD: {std_val:.4f}'
        plt.axvline(x=mean_val, color='red', linestyle='--', linewidth=1.5, label=label_text)
        
    if median_val is not None:
        plt.axvline(x=median_val, color='orange', linestyle='-', linewidth=1.5, label=f'Median: {median_val:.4f}')
        
    # Additional Threshold Lines
    if thresholds:
        for thr in thresholds:
            val = thr.get('value')
            if val is None: continue
            lbl = thr.get('label', f'{val:.4f}')
            clr = thr.get('color', 'black')
            ls = thr.get('linestyle', '--')
            lw = thr.get('linewidth', 1.5)
            plt.axvline(x=val, color=clr, linestyle=ls, linewidth=lw, label=lbl)
    
    plt.title(title, fontsize=15)
    plt.xlabel(x_label, fontsize=12)
    
    if log_scale:
        plt.yscale('log')
        plt.ylabel(f"{y_label} (Log Scale)", fontsize=12)
    else:
        plt.ylabel(y_label, fontsize=12)
        
    if xlim is not None:
        plt.xlim(xlim)

    plt.legend(loc='upper right')
    plt.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"Saved plot to {filename}")
    plt.close()


def plot_regression_comparison(
    df, 
    x_col, 
    y_col, 
    x_label, 
    y_label, 
    filename,
    title=None
):
    """
    Plots a scatter plot with OLS and Huber regression lines.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn.linear_model import HuberRegressor, LinearRegression
    
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
    plt.figure(figsize=(10, 8))
    
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
        plt.title(title)
    else:
        plt.title(f'{y_label} vs {x_label}')

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(x_min - 0.05*x_span, x_max + 0.05*x_span)
    
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(filename, dpi=300)
    print(f"Saved plot to {filename}")
    plt.close()


def plot_stacked_distribution(
    data, 
    col, 
    group_col,
    title, 
    filename, 
    mean_val=None, 
    median_val=None,
    x_label=None, 
    y_label="Count", 
    log_scale=False,
    bins=100
):
    """
    Plots a stacked histogram distribution split by a group column.
    Includes Mean and Median vertical lines.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    if x_label is None:
        x_label = col

    plt.figure(figsize=(10, 6))
    
    # Histogram with Stacked Groups
    ax = sns.histplot(data=data, x=col, hue=group_col, multiple='stack', 
                      bins=bins, linewidth=0.1)
    
    # Statistics Lines
    handles, labels = [], []
    
    # Extract existing legend handles (Groups) if present
    if ax.get_legend():
        try:
            handles = ax.get_legend().legend_handles
            texts = ax.get_legend().get_texts()
            labels = [t.get_text() for t in texts]
            ax.get_legend().remove()
        except:
            pass # Fallback if legend structure is unexpected
    
    # Add Stats Lines
    if mean_val is not None:
        l1 = plt.axvline(x=mean_val, color='blue', linestyle='--', linewidth=1, label=f'Mean: {mean_val:.4f}')
        handles.append(l1)
        labels.append(l1.get_label())
        
    if median_val is not None:
        l2 = plt.axvline(x=median_val, color='orange', linestyle=':', linewidth=1.5, label=f'Median: {median_val:.4f}')
        handles.append(l2)
        labels.append(l2.get_label())
    
    plt.title(title, fontsize=15)
    plt.xlabel(x_label, fontsize=12)
    
    if log_scale:
        plt.yscale('log')
        plt.ylabel(f"{y_label} (Log Scale)", fontsize=12)
    else:
        plt.ylabel(y_label, fontsize=12)
        
    if handles:
        plt.legend(handles=handles, labels=labels, title=f'{group_col} & Stats', loc='upper right')
    plt.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"Saved plot to {filename}")
    plt.close()


def plot_joint_regression(
    df, 
    x_col, 
    y_col, 
    group_col,
    x_label, 
    y_label, 
    filename,
    title=None
):
    """
    Plots a JointGrid (Scatter + Marginal Histograms) with OLS and Huber regression lines.
    Supports grouping/hue in scatter and histograms.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn.linear_model import HuberRegressor, LinearRegression
    
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
    g = sns.JointGrid(data=local_valid, x=x_col, y=y_col, height=10, ratio=5)

    # 1. Main Scatter Plot
    sns.scatterplot(data=local_valid, x=x_col, y=y_col, hue=group_col,
                    alpha=0.6, s=20, edgecolor='none', ax=g.ax_joint, legend='full')
    
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
    # Retrieve Group Legend
    group_legend = g.ax_joint.get_legend()
    if group_legend:
        group_legend.set_bbox_to_anchor((1.25, 1.0))
        group_legend.set_loc('upper left')
        group_legend.set_title(group_col)
        group_legend.set_frame_on(False)
    
    # Add Regression Legend (this replaces the standard legend slot)
    g.ax_joint.legend(handles=[line1, line2], title="Regression Models", 
                      loc='upper left', bbox_to_anchor=(1.25, 0.4), frameon=False)
    
    # Add Group Legend back manually
    if group_legend:
        g.ax_joint.add_artist(group_legend)
    
    # 3. Marginals (Stacked Histograms)
    sns.histplot(data=local_valid, x=x_col, hue=group_col, multiple='stack', 
                 bins=100, linewidth=0.1, ax=g.ax_marg_x, legend=False)
    sns.histplot(data=local_valid, y=y_col, hue=group_col, multiple='stack', 
                 bins=50, linewidth=0.1, ax=g.ax_marg_y, legend=False)

    g.ax_joint.set_xlabel(x_label)
    g.ax_joint.set_ylabel(y_label)
    
    if title:
        plt.suptitle(title, y=1.02, fontsize=16)
    else:
        plt.suptitle(f'{y_label} vs {x_label}', y=1.02, fontsize=16)
    
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
    sample_size=50000
):
    """
    Plots Mean vs Variance with fitted Poisson and Negative Binomial lines.
    """
    import matplotlib.pyplot as plt
    import pandas as pd

    plt.figure(figsize=(10, 8))
    
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
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"Saved fit plot: {filename}")
    plt.close()

def plot_qq_residuals(
    residuals,
    filename,
    title="QQ Plot of Standardized Residuals"
):
    """
    Plots a QQ plot of residuals using statsmodels.
    """
    import matplotlib.pyplot as plt
    import statsmodels.api as sm

    plt.figure(figsize=(8, 6))
    sm.qqplot(residuals, line='45', fit=True)
    plt.title(title)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"Saved QQ plot: {filename}")
    plt.close()


def plot_bar_chart(
    names,
    values,
    title,
    ylabel,
    filename,
    ylim=(0.0, 1.05),
    color="#1f77b4",
    figsize=(10, 5)
):
    """
    Plots a simple bar chart with value labels on top.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os
    
    # Create directory if needed (defensive)
    folder = os.path.dirname(filename)
    if folder:
        os.makedirs(folder, exist_ok=True)

    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(names, values, color=color, alpha=0.8)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if ylim:
        ax.set_ylim(ylim)
        
    # Add text labels
    for i, v in enumerate(values):
        ax.text(i, v + 0.01, f"{v:.3f}", ha='center', va='bottom', fontsize=9)
        
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    print(f"Saved bar chart: {filename}")
    plt.close(fig)


def plot_gwas_qq(
    expected,
    observed,
    title,
    filename,
    xlabel="Expected -log10(P)",
    ylabel="Observed -log10(P)",
    figsize=(6, 6)
):
    """
    Plots a GWAS QQ plot (Expected vs Observed -log10 P-values).
    """
    import matplotlib.pyplot as plt
    import os
    
     # Create directory if needed
    folder = os.path.dirname(filename)
    if folder:
        os.makedirs(folder, exist_ok=True)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(expected, observed, s=6, alpha=0.6)
    
    # Diagonal line
    # Determine limit based on max value in both arrays
    lim = max(expected.max(), observed.max())
    ax.plot([0, lim], [0, lim], color='red', lw=1)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    fig.tight_layout()
    fig.savefig(filename, dpi=160)
    print(f"Saved GWAS QQ plot: {filename}")
    plt.close(fig)


def plot_correlation_with_regression(
    data, x_col, y_col, 
    title, filename, 
    x_label=None, y_label=None,
    color='gray', line_color='red'
):
    """
    Plots a scatter plot with linear regression line and statistics (slope, intercept, R^2).
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np

    if x_label is None: x_label = x_col
    if y_label is None: y_label = y_col

    clean_data = data[[x_col, y_col]].dropna()
    x = clean_data[x_col]
    y = clean_data[y_col]

    if len(x) < 2:
        print(f"[Warning] Not enough data points to plot regression for {x_col} vs {y_col}")
        return

    slope, intercept = np.polyfit(x, y, 1)
    r_squared = np.corrcoef(x, y)[0, 1] ** 2

    plt.figure(figsize=(10, 6))
    sns.regplot(
        x=x_col, y=y_col, data=data,
        scatter_kws={'alpha':0.4, 's':10, 'color': color},
        line_kws={'color': line_color, 'label': 'Linear Regression'}
    )

    stats_text = f'$y = {slope:.4f}x + {intercept:.4f}$\n$R^2 = {r_squared:.4f}$'
    plt.text(
        0.5, 0.9, stats_text,
        transform=plt.gca().transAxes, 
        fontsize=12, color='darkred',
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray', boxstyle='round')
    )

    plt.title(title, fontsize=15)
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel(y_label, fontsize=12)
    plt.ylim(0, 1)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"Plot saved to {filename}")


if __name__ == "__main__":
    if len(sys.argv) >= 4:
        combine_plots([sys.argv[1], sys.argv[2]], sys.argv[3])
    else:
        print("Usage: python graph.py <image1> <image2> <output_file>")
    
