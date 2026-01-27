import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
from sklearn.linear_model import HuberRegressor, LinearRegression

# Global settings
PLOT_SAMPLE_SIZE = 100000 
OUTPUT_PREFIX = "missing_analysis_"

print("1. Reading Data...")
# 1. Read vmiss
try:
    df_vmiss = pd.read_csv("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.missing.vmiss", sep='\s+')
    df_vmiss.columns = [c.replace('#', '') for c in df_vmiss.columns]
    # Extract Pos from ID (e.g. '2-952' -> 952)
    df_vmiss['Pos_key'] = df_vmiss['ID'].str.split('-').str[1].astype(int)
except Exception as e:
    print(f"Error reading vmiss: {e}")
    sys.exit(1)

# 2. Read popdep
try:
    df_popdep = pd.read_csv("/data/home/tusr1/01projects/vmap4/05reliable.lib/01test/chr002_popdep.txt", 
                            sep='\s+', 
                            usecols=['Position', 'Depth_SD', 'Depth_Mean'])
except Exception as e:
    print(f"Error reading popdep: {e}")
    sys.exit(1)

# 3. Merge
print("Merging data...")
merged = pd.merge(df_vmiss[['Pos_key', 'F_MISS']], 
                  df_popdep, 
                  left_on='Pos_key', 
                  right_on='Position', 
                  how='inner')

if merged.empty:
    print("Error: No overlapping positions.")
    sys.exit(1)

print(f"Merged {len(merged)} sites.")

# 4. Calculate Derived Metrics
print("Calculating metrics...")
# Variance
merged['Depth_Var'] = merged['Depth_SD'] ** 2
# CV (SD / Mean)
merged['Depth_CV'] = np.where(merged['Depth_Mean'] > 0, 
                              merged['Depth_SD'] / merged['Depth_Mean'], 
                              np.nan)

# Drop NaNs for basic metrics
merged = merged.dropna(subset=['Depth_CV', 'F_MISS'])

# 5. Calculate Log Metrics
# Handle 0s by replacing with NaN before log to avoid -inf columns everywhere
# For F_MISS, 0 is common (no missing). log(0) is -inf.
# Depending on goal, we either drop 0s or use log1p. 
# Prompt says "log of missing values", usually implies treating 0s as excluded or undefined.
merged['Log_Depth_Mean'] = np.log(merged['Depth_Mean'].replace(0, np.nan))
merged['Log_Depth_SD'] = np.log(merged['Depth_SD'].replace(0, np.nan))
merged['Log_Depth_Var'] = np.log(merged['Depth_Var'].replace(0, np.nan))
merged['Log_Depth_CV'] = np.log(merged['Depth_CV'].replace(0, np.nan))
merged['Log_F_MISS'] = np.log(merged['F_MISS'].replace(0, np.nan))

# --- Plotting Functions ---

def plot_missing_distribution(df, col, output_file):
    print(f"Plotting distribution for {col}...")
    plt.figure(figsize=(10, 6))
    
    data = df[col].dropna()
    if len(data) == 0:
        print(f"No data for {col}")
        plt.close()
        return

    # Histogram
    sns.histplot(data, bins=50, kde=False, color='skyblue', edgecolor='black', stat='count')
    
    # Calculate Stats
    min_val = data.min()
    max_val = data.max()
    mean_val = data.mean()
    median_val = data.median()
    
    # Annotate stats on the plot
    stats_text = (f"Min: {min_val:.5f}\n"
                  f"Max: {max_val:.5f}\n"
                  f"Mean: {mean_val:.4f}\n"
                  f"Median: {median_val:.4f}")
    
    plt.text(0.95, 0.95, stats_text, 
             transform=plt.gca().transAxes, 
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
             fontsize=12)

    # Mark Mean/Median on Axis
    plt.axvline(mean_val, color='red', linestyle='--', linewidth=1.5, label=f'Mean: {mean_val:.4f}')
    plt.axvline(median_val, color='orange', linestyle='-.', linewidth=1.5, label=f'Median: {median_val:.4f}')
    plt.legend()

    plt.title(f'Distribution of Missing Rate')
    plt.xlabel('Missing Rate (F_MISS)')
    plt.ylabel('Count')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"   Saved {output_file}")
    plt.close()

def plot_dual_regression(df, x_col, y_col, x_label, y_label, output_file):
    print(f"Plotting Regression: {y_col} vs {x_col}")
    
    # Clean data: drop NaNs and Infs (especially from log)
    df_clean = df[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    
    if len(df_clean) < 10:
        print(f"   Skipping {output_file}: Not enough valid data points ({len(df_clean)}).")
        return

    X = df_clean[x_col].values.reshape(-1, 1)
    y = df_clean[y_col].values
    
    # 1. OLS Regression
    ols = LinearRegression()
    ols.fit(X, y)
    ols_r2 = ols.score(X, y)
    ols_eq = f"y = {ols.coef_[0]:.4f}x + {ols.intercept_:.4f}"

    # 2. Huber Regression (Robust)
    huber = HuberRegressor(epsilon=1.35, max_iter=200)
    huber.fit(X, y)
    huber_r2 = huber.score(X, y)
    huber_eq = f"y = {huber.coef_[0]:.4f}x + {huber.intercept_:.4f}"

    # Prepare Plot
    plt.figure(figsize=(10, 8))
    
    # Downsample for Scatter Plot only (calculations used full data)
    if len(df_clean) > PLOT_SAMPLE_SIZE:
        df_plot = df_clean.sample(n=PLOT_SAMPLE_SIZE, random_state=42)
    else:
        df_plot = df_clean
        
    sns.scatterplot(x=df_plot[x_col], y=df_plot[y_col], 
                    alpha=0.2, s=10, color='gray', edgecolor=None, label='Data points (subsampled)')

    # Generate regression lines
    x_min, x_max = df_clean[x_col].min(), df_clean[x_col].max()
    x_vec = np.linspace(x_min, x_max, 100).reshape(-1, 1)
    
    plt.plot(x_vec, ols.predict(x_vec), 'b--', linewidth=2, label=f'OLS: {ols_eq} ($R^2$={ols_r2:.3f})')
    plt.plot(x_vec, huber.predict(x_vec), 'r-', linewidth=2, label=f'Huber: {huber_eq} ($R^2$={huber_r2:.3f})')
    
    # Set Limits
    # X Axis: Use 0.99 quantile for upper limit to ignore extreme outliers
    x_min_val = df_clean[x_col].min()
    x_max_val = df_clean[x_col].max()
    x_99_val = df_clean[x_col].quantile(0.99)
    
    # If 99% quantile is significantly smaller than max, clip view to slightly above 99%
    # Otherwise use max.
    if x_99_val < x_max_val:
        x_upper_lim = x_99_val * 1.05 # Add 5% buffer
    else:
        x_upper_lim = x_max_val + (x_max_val - x_min_val) * 0.05
    
    if x_upper_lim <= x_min_val: x_upper_lim = x_min_val + 1.0 # fallback

    # Add small buffer to min as well
    x_range_view = x_upper_lim - x_min_val
    plt.xlim(x_min_val - 0.02 * x_range_view, x_upper_lim)
    
    # Y Axis: 0-1 for Missing Rate, Auto for others (Log)
    if y_col == 'F_MISS':
        plt.ylim(-0.02, 1.02)
    else:
        y_min, y_max = df_clean[y_col].min(), df_clean[y_col].max()
        y_range_val = y_max - y_min
        if y_range_val == 0: y_range_val = 1.0
        plt.ylim(y_min - 0.05 * y_range_val, y_max + 0.05 * y_range_val)

    plt.title(f'{y_label} vs {x_label}')
    plt.xlabel(f'{x_label} ({x_col})')
    plt.ylabel(f'{y_label} ({y_col})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"   Saved {output_file}")
    plt.close()

# --- Execution ---

# # 1. Distribution of Missing Rate
# plot_missing_distribution(merged, 'F_MISS', f"{OUTPUT_PREFIX}dist_fmiss.png")

# # 2. Standard Regressions (F_MISS vs Depth Metrics)
# metrics = [
#     ('Depth_Mean', 'Mean Depth'),
#     ('Depth_SD', 'Depth SD'),
#     ('Depth_Var', 'Depth Variance'),
#     ('Depth_CV', 'Depth CV')
# ]

# for col, label in metrics:
#     plot_dual_regression(merged, col, 'F_MISS', label, 'Missing Rate', 
#                          f"{OUTPUT_PREFIX}reg_fmiss_vs_{col.lower()}.png")

# 3. Log-Log Regressions (Log F_MISS vs Log Depth Metrics)
# Note: Log F_MISS implies filtering out F_MISS=0 (sites with complete data).
# This focuses on the relationship magnitude where missingness exists.
log_metrics = [
    ('Log_Depth_Mean', 'Log Mean Depth'),
    ('Log_Depth_SD', 'Log Depth SD'),
    ('Log_Depth_Var', 'Log Depth Variance'),
    ('Log_Depth_CV', 'Log Depth CV')
]

for col, label in log_metrics:
    plot_dual_regression(merged, col, 'F_MISS', label, 'Missing Rate', 
                         f"{OUTPUT_PREFIX}reg_fmiss_vs_{col.lower()}.png")

print("All analyses completed.")