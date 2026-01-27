import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
from sklearn.linear_model import HuberRegressor, LinearRegression

# Global settings
PLOT_SAMPLE_SIZE = 100000 
OUTPUT_PREFIX = "site_depth_variant_"

print("Processing Depth Data Analysis...")
# 1. Read Depth Data
try:
    # Adjust path if necessary based on previous context
    input_file = "/data/home/tusr1/01projects/vmap4/05reliable.lib/01test/chr002_popdep.txt" 
    print(f"Reading {input_file}...")
    df_popdep = pd.read_csv(input_file, sep=r'\s+', usecols=['Position', 'Depth_SD', 'Depth_Mean'])
except Exception as e:
    print(f"Failed to read file: {e}")
    sys.exit(1)

# 2. Calculate Derived Metrics
print("Calculating Variance and CV...")
df_popdep['Depth_Var'] = df_popdep['Depth_SD'] ** 2
# Handle division by zero for CV
df_popdep['Depth_CV'] = np.where(df_popdep['Depth_Mean'] > 0, 
                                 df_popdep['Depth_SD'] / df_popdep['Depth_Mean'], 
                                 np.nan)

# Drop NaNs created (where Mean=0)
initial_len = len(df_popdep)
df_valid = df_popdep.dropna(subset=['Depth_CV']).copy()
print(f"Dropped {initial_len - len(df_valid)} sites with 0 mean depth.")

# Calculate Log metrics (handle 0s by replacing with NaN for SD/Var/CV)
# Depth_Mean is already > 0 due to dropna(subset=['Depth_CV'])
df_valid['Log_Depth_Mean'] = np.log(df_valid['Depth_Mean'])
df_valid['Log_Depth_SD'] = np.log(df_valid['Depth_SD'].replace(0, np.nan))
df_valid['Log_Depth_Var'] = np.log(df_valid['Depth_Var'].replace(0, np.nan))
df_valid['Log_Depth_CV'] = np.log(df_valid['Depth_CV'].replace(0, np.nan))

# --- NEW: Read Highlight Positions ---
print("Reading Position Highlight File...")
df_highlight = pd.DataFrame()
try:
    pos_file = "/data/home/tusr1/01projects/vmap4/05reliable.lib/01test/2_1_122798052.pos.txt"
    # No header, columns: Chrom, Position
    df_pos_target = pd.read_csv(pos_file, sep=r'\s+', header=None, names=['Chrom', 'Position'])
    target_positions = set(df_pos_target['Position'])
    
    # Create Highlight Subset
    df_highlight = df_valid[df_valid['Position'].isin(target_positions)].copy()
    print(f"Found {len(df_highlight)} sites to highlight out of {len(df_pos_target)} provided positions.")
    
    # Downsample Highlight for plotting (avoid slowness with huge highlight sets)
    if len(df_highlight) > PLOT_SAMPLE_SIZE:
        df_highlight = df_highlight.sample(n=PLOT_SAMPLE_SIZE, random_state=42)
        print(f"Downsampled Highlight to {PLOT_SAMPLE_SIZE} points for plotting.")

except Exception as e:
    print(f"Warning: Failed to read/process position file: {e}")

# Downsample for plotting if data is too large
if len(df_valid) > PLOT_SAMPLE_SIZE:
    df_plot = df_valid.sample(n=PLOT_SAMPLE_SIZE, random_state=42)
    print(f"Downsampled to {PLOT_SAMPLE_SIZE} points for scatter plots.")
else:
    df_plot = df_valid

# --- Helper Function for Distribution Plots ---
def plot_distribution(col, label, color, filename):
    plt.figure(figsize=(10, 6))
    
    # 1. 先计算关注的范围上限 (99% 分位数) - 计算分位数很快，可以用全量数据
    x_limit = df_valid[col].quantile(0.99) * 1.1
    
    # 2. 关键修改: 使用 df_plot (降采样数据) 进行绘图
    #kde=True 在大数据量下(百万级)是极慢的，改用 5万个采样点瞬间即可完成，且分布形状几乎无损
    sns.histplot(df_plot[col], bins=500, binrange=(0, x_limit), kde=True, color=color, stat="density", label='Background')
    
    # Highlight Overlay
    if not df_highlight.empty:
        sns.histplot(df_highlight[col], bins=500, binrange=(0, x_limit), kde=True, color='red', stat="density", label='Highlight', alpha=0.5)
        # Add HL Mean
        hl_mean = df_highlight[col].mean()
        plt.axvline(hl_mean, color='magenta', linestyle=':', linewidth=2, label=f'HL Mean: {hl_mean:.2f}')

    # Add stats - 统计值(均值/中位数)仍然使用全量数据 df_valid 以保持高精度
    mean_val = df_valid[col].mean()
    median_val = df_valid[col].median()
    min_val = df_valid[col].min()
    max_val = df_valid[col].max()
    
    plt.axvline(mean_val, color='red', linestyle='--', label=f'Mean: {mean_val:.2f}')
    plt.axvline(median_val, color='green', linestyle='--', label=f'Median: {median_val:.2f}')
    # Add Min/Max to legend (using invisible lines so they appear in legend)
    plt.plot([], [], ' ', label=f'Min: {min_val:.2f}')
    plt.plot([], [], ' ', label=f'Max: {max_val:.2f}')
    
    plt.title(f'Distribution of {label}')
    plt.xlabel(label)
    plt.xlim(0, x_limit)  # 设置由于 binrange 已经限制了范围，这里主要是为了对齐视觉
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    plt.close()

# --- Helper Function for Regression Plots ---
def plot_dual_regression(x_col, y_col, x_label, y_label, filename):
    print(f"Running regressions for {y_label} vs {x_label}...")
    
    # Filter for valid data for this specific pair (handling NaNs/Infs from Log transform)
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
    
    # Calculate limits from local_valid (compatible with negative Log values)
    # Use 0.1% and 99.9% quantiles to ignore extreme outliers
    x_low = local_valid[x_col].quantile(0.001)
    x_high = local_valid[x_col].quantile(0.999)
    y_low = local_valid[y_col].quantile(0.001)
    y_high = local_valid[y_col].quantile(0.999)
    
    # Pad limits slightly
    x_span = x_high - x_low
    y_span = y_high - y_low
    x_lims = (x_low - 0.05*x_span, x_high + 0.05*x_span)
    y_lims = (y_low - 0.05*y_span, y_high + 0.05*y_span)

    # Scatter plot (downsampled) - also filter for valid values
    local_plot = df_plot[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    sns.scatterplot(data=local_plot, x=x_col, y=y_col, 
                    alpha=0.1, s=10, color='gray', edgecolor=None, label='Data (Sampled)')
    
    # Highlight Points Overlay
    if not df_highlight.empty:
        local_hl = df_highlight[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
        if not local_hl.empty:
            sns.scatterplot(data=local_hl, x=x_col, y=y_col, 
                            alpha=0.9, s=30, color='red', marker='X', edgecolor='white', linewidth=0.5, label='Highlight')

    # Generate line points
    x_range = np.linspace(x_lims[0], x_lims[1], 100).reshape(-1, 1)
    y_ols = ols.predict(x_range)
    y_huber = huber.predict(x_range)
    
    # Plot lines
    plt.plot(x_range, y_ols, color='blue', linewidth=2, linestyle='--', 
             label=f'OLS: {ols_eq}, $R^2$={ols_score:.3f}')
    plt.plot(x_range, y_huber, color='red', linewidth=2, 
             label=f'Huber: {huber_eq}, $R^2$={huber_score:.3f}')
    
    plt.title(f'{y_label} vs {x_label}')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(x_lims)
    plt.ylim(y_lims)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    plt.close()

# --- Execution ---

# 1. Distributions
plot_distribution('Depth_Mean', 'Mean Depth', 'teal', f"{OUTPUT_PREFIX}dist_mean.png")
plot_distribution('Depth_SD', 'Depth SD', 'orange', f"{OUTPUT_PREFIX}dist_sd.png")
plot_distribution('Depth_Var', 'Depth Variance', 'purple', f"{OUTPUT_PREFIX}dist_var.png")
plot_distribution('Depth_CV', 'Depth CV', 'green', f"{OUTPUT_PREFIX}dist_cv.png")

# 2. Regressions (Mean as independent variable)
# Mean vs SD
plot_dual_regression('Depth_Mean', 'Depth_SD', 'Mean Depth', 'Standard Deviation', 
                     f"{OUTPUT_PREFIX}reg_mean_vs_sd.png")

# Mean vs Variance
plot_dual_regression('Depth_Mean', 'Depth_Var', 'Mean Depth', 'Variance', 
                     f"{OUTPUT_PREFIX}reg_mean_vs_var.png")

# Mean vs CV
plot_dual_regression('Depth_Mean', 'Depth_CV', 'Mean Depth', 'Coefficient of Variation', 
                     f"{OUTPUT_PREFIX}reg_mean_vs_cv.png")

# 3. Log-Log Regressions
plot_dual_regression('Log_Depth_Mean', 'Log_Depth_SD', 'Log Mean Depth', 'Log SD', 
                     f"{OUTPUT_PREFIX}reg_log_mean_vs_log_sd.png")

plot_dual_regression('Log_Depth_Mean', 'Log_Depth_Var', 'Log Mean Depth', 'Log Variance', 
                     f"{OUTPUT_PREFIX}reg_log_mean_vs_log_var.png")

plot_dual_regression('Log_Depth_Mean', 'Log_Depth_CV', 'Log Mean Depth', 'Log CV', 
                     f"{OUTPUT_PREFIX}reg_log_mean_vs_log_cv.png")

print("Analysis Complete.")