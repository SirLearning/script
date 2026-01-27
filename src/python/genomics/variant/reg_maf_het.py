import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import numpy as np
from sklearn.linear_model import HuberRegressor, LinearRegression

# Settings
PLOT_SAMPLE_SIZE = 1000000
INPUT_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.hardy.hardy"
OUTPUT_REG_FILE = "reg_maf_vs_het.png"

print("Processing Variant Stats (Regression MAF vs Heterozygosity)...")

# 1. Read Data
try:
    print(f"Reading {INPUT_FILE}...")
    df = pd.read_csv(INPUT_FILE, sep=r'\s+')
except Exception as e:
    print(f"Failed to read file: {e}")
    sys.exit(1)

# Check columns
required_cols = ['HOM_A1_CT', 'HET_A1_CT', 'TWO_AX_CT']
if not all(col in df.columns for col in required_cols):
    print(f"Error: Missing columns. Found: {df.columns}")
    sys.exit(1)

# ==========================================
# Data Processing
# ==========================================
print("\n--- Processing Data ---")
# Calculate Metrics & MAF 
df['Total_Samples'] = df['HOM_A1_CT'] + df['HET_A1_CT'] + df['TWO_AX_CT']

# Filter valid
df_valid = df[df['Total_Samples'] > 0].copy()

# Hobs (Observed Heterozygosity)
df_valid['Hobs'] = df_valid['HET_A1_CT'] / df_valid['Total_Samples']
# p (Freq of A1)
df_valid['p'] = (df_valid['HOM_A1_CT'] * 2 + df_valid['HET_A1_CT']) / (df_valid['Total_Samples'] * 2)
# Hexp
df_valid['Hexp'] = 2 * df_valid['p'] * (1 - df_valid['p'])

# Filter Hexp > 0 (Consistent with original script Part 2 logic)
df_f = df_valid[df_valid['Hexp'] > 0].copy()

# MAF
df_f['MAF'] = df_f['p'].apply(lambda x: min(x, 1-x))

# Stats
print(f"Valid Sites: {len(df_f)}")

# Downsample for plotting
if len(df_f) > PLOT_SAMPLE_SIZE:
    df_plot_f = df_f.sample(n=PLOT_SAMPLE_SIZE, random_state=42)
else:
    df_plot_f = df_f

# ==========================================
# Regression Plotting Function
# ==========================================
def plot_dual_regression(df_valid, df_plot, x_col, y_col, x_label, y_label, filename):
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
    
    # Calculate limits
    x_limit = 0.55 # MAF max is 0.5
    y_limit = local_valid[y_col].max() * 1.1

    # Scatter plot (downsampled)
    local_plot = df_plot[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    sns.scatterplot(data=local_plot, x=x_col, y=y_col, 
                    alpha=0.1, s=10, color='gray', edgecolor=None, label='Data (Sampled)')
    
    # Generate line points
    x_range = np.linspace(0, 0.5, 100).reshape(-1, 1)
    y_ols = ols.predict(x_range)
    y_huber = huber.predict(x_range)
    y_hwe = 2 * x_range * (1 - x_range)
    
    # Plot lines
    plt.plot(x_range, y_ols, color='blue', linewidth=2, linestyle='--', 
             label=f'OLS: {ols_eq}, $R^2$={ols_score:.3f}')
    plt.plot(x_range, y_huber, color='red', linewidth=2, 
             label=f'Huber: {huber_eq}, $R^2$={huber_score:.3f}')
    plt.plot(x_range, y_hwe, color='green', linewidth=2, linestyle='-.', 
             label='Expected (HWE): $2p(1-p)$')
    
    plt.title(f'{y_label} vs {x_label}')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(0, 0.55)
    plt.ylim(0, y_limit)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    plt.close()

# Use df_f and df_plot_f (if defined above)
plot_dual_regression(df_f, df_plot_f, 'MAF', 'Hobs', 'Minor Allele Frequency (MAF)', 'Observed Heterozygosity', OUTPUT_REG_FILE)
