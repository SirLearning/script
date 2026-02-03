from ast import arg
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def combine_plots(
    image1 = "/data/home/tusr1/git/script/out/site_depth_variant_reg_log_mean_vs_log_var.png",
    image2 = "/data/home/tusr1/git/script/out/site_depth_variant_reg_log_mean_vs_log_cv.png",
    output_file = "/data/home/tusr1/git/script/out/site_depth_variant_reg_log_mean_vs_log_combined.png"
):
    # Settings
    images = [
        image1,
        image2
    ]

    # Create a figure with subplots
    fig, axes = plt.subplots(1, len(images), figsize=(6 * len(images), 5))

    for ax, img_path in zip(axes, images):
        try:
            img = mpimg.imread(img_path)
            ax.imshow(img)
            ax.axis('off')  # Hide axis
        except FileNotFoundError:
            print(f"Error: File not found {img_path}")
            ax.text(0.5, 0.5, 'Image Not Found', ha='center', va='center')
            ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
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


if __name__ == "__main__":
    combine_plots(arg[1], arg[2], arg[3])
    
