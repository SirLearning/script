import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
from sklearn.linear_model import HuberRegressor, LinearRegression

# Settings
SCOUNT_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.sample_stats.scount"
MISSING_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.missing.smiss"
OUTPUT_PREFIX = "sample_het_analysis"

print("Processing Sample Heterozygosity Analysis...")

# 1. Read Heterozygosity Data (.scount)
print("Reading .scount data...")
try:
    df_het = pd.read_csv(SCOUNT_FILE, sep=r'\s+')
except Exception as e:
    print(f"Failed to read scount file: {e}")
    sys.exit(1)

# Check scount columns
required_cols_het = ['#IID', 'HOM_REF_CT', 'HET_SNP_CT', 'HOM_ALT_SNP_CT']
if not all(col in df_het.columns for col in required_cols_het):
    print(f"Error: Missing columns in scount. Found: {df_het.columns}")
    sys.exit(1)

# Calculate Heterozygosity Rate
df_het['Total_Called'] = df_het['HOM_REF_CT'] + df_het['HET_SNP_CT'] + df_het['HOM_ALT_SNP_CT']
df_het['Het_Rate'] = df_het['HET_SNP_CT'] / df_het['Total_Called']
df_het = df_het[['#IID', 'Het_Rate']].rename(columns={'#IID': 'Sample'})

# Mean/Median Heterozygosity
mean_het = df_het['Het_Rate'].mean()
median_het = df_het['Het_Rate'].median()
print(f"Heterozygosity - Mean: {mean_het:.4f}, Median: {median_het:.4f}")

# 2. Read Missing Data (.smiss)
print("Reading .smiss data...")
try:
    df_mid = pd.read_csv(MISSING_FILE, sep=r'\s+')
except Exception as e:
    print(f"Failed to read smiss file: {e}")
    sys.exit(1)

# Check smiss columns
# 通常包含 #IID (或 IID) 和 F_MISS
if '#IID' in df_mid.columns:
    iid_col = '#IID'
elif 'IID' in df_mid.columns:
    iid_col = 'IID'
else:
    print(f"Error: Missing IID column in smiss. Found: {df_mid.columns}")
    sys.exit(1)

missing_col = 'F_MISS'
if missing_col not in df_mid.columns:
     print(f"Error: Missing '{missing_col}' in smiss. Found: {df_mid.columns}")
     sys.exit(1)
     
df_mid = df_mid[[iid_col, missing_col]].rename(columns={iid_col: 'Sample', missing_col: 'Missing_Rate'})

# 3. Merge Data
print("Merging Heterozygosity and Missing Rate data...")
df_merged = pd.merge(df_het, df_mid, on='Sample', how='inner')
print(f"Merged samples: {len(df_merged)}")


# Set Seaborn Style
sns.set_theme(style="ticks")

# ==========================================
# Plot 1: Heterozygosity Distribution
# ==========================================
def plot_distribution(data, col, mean_val, median_val, title, filename, log_scale=False):
    plt.figure(figsize=(8, 6))
    
    # Histogram
    sns.histplot(data[col], bins=100, color='forestgreen', edgecolor='white', linewidth=0.1)
    
    # Statistics Lines
    plt.axvline(x=mean_val, color='blue', linestyle='--', linewidth=1, label=f'Mean: {mean_val:.4f}')
    plt.axvline(x=median_val, color='orange', linestyle=':', linewidth=1.5, label=f'Median: {median_val:.4f}')
    
    plt.title(title, fontsize=15)
    plt.xlabel("Sample Heterozygosity Rate", fontsize=12)
    
    if log_scale:
        plt.yscale('log')
        plt.ylabel("Count (Log Scale)", fontsize=12)
    else:
        plt.ylabel("Count", fontsize=12)
        
    plt.legend(loc='upper right')
    plt.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    plt.close()

print("Generating Het Distribution Plots...")
plot_distribution(df_merged, 'Het_Rate', mean_het, median_het, 
                  "Distribution of Sample Heterozygosity", f"{OUTPUT_PREFIX}_dist_het.png")

plot_distribution(df_merged, 'Het_Rate', mean_het, median_het, 
                  "Distribution of Sample Heterozygosity (Log Scale)", f"{OUTPUT_PREFIX}_dist_het_log.png", 
                  log_scale=True)


# ==========================================
# Plot 2: Heterozygosity vs Missing Rate Regression
# ==========================================
def plot_dual_regression(df, x_col, y_col, x_label, y_label, filename):
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
    
    # Scatter plot
    sns.scatterplot(data=local_valid, x=x_col, y=y_col, 
                    alpha=0.3, s=15, color='#1f77b4', edgecolor='none', label='Samples')
    
    # Generate line points
    x_min, x_max = local_valid[x_col].min(), local_valid[x_col].max()
    x_span = x_max - x_min
    x_range = np.linspace(x_min - 0.05*x_span, x_max + 0.05*x_span, 100).reshape(-1, 1)
    
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
    plt.xlim(x_min - 0.05*x_span, x_max + 0.05*x_span)
    # plt.ylim(...) # Let matplotlib autoscale y usually works well
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    plt.close()

print("Generating Regression Plot (Het vs Missing)...")
# X: Missing Rate, Y: Heterozygosity (Usually Het depends on Missing via artifacts, or vice versa)
# User asked "Heterozygosity and Missing Rate scatter", let's put Missing on X as it's often a quality metric predicting others
plot_dual_regression(df_merged, 'Missing_Rate', 'Het_Rate', 
                     'Missing Rate', 'Heterozygosity Rate', 
                     f"{OUTPUT_PREFIX}_reg_miss_vs_het.png")

print("Analysis Complete.")
