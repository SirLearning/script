import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import numpy as np
from sklearn.linear_model import HuberRegressor, LinearRegression

# Settings
PLOT_SAMPLE_SIZE = 1000000
INPUT_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.hardy.hardy"
OUTPUT_HET_FILE = "variant_het_dist.png"
OUTPUT_F_FILE = "variant_f_dist.png"

print("Processing Variant Stats (Heterozygosity & Inbreeding Coefficient)...")

# 1. Read Data
try:
    print(f"Reading {INPUT_FILE}...")
    df = pd.read_csv(INPUT_FILE, sep=r'\s+')
except Exception as e:
    print(f"Failed to read file: {e}")
    sys.exit(1)

# Check columns
# PLINK2 headers usually: #CHROM ID A1 AX HOM_A1_CT HET_A1_CT TWO_AX_CT O(HET_A1) E(HET_A1) P
required_cols = ['HOM_A1_CT', 'HET_A1_CT', 'TWO_AX_CT']
if not all(col in df.columns for col in required_cols):
    print(f"Error: Missing columns. Found: {df.columns}")
    sys.exit(1)

# ==========================================
# Part 1: Observed Heterozygosity Distribution
# ==========================================
print("\n--- 1. Generating Heterozygosity Distribution ---")
# Calculate Metrics & MAF for Grouping
df['Total_Samples'] = df['HOM_A1_CT'] + df['HET_A1_CT'] + df['TWO_AX_CT']
het_col = 'O(HET_A1)'
if het_col not in df.columns:
    df[het_col] = df['HET_A1_CT'] / df['Total_Samples']

# Calculate MAF (for grouping)
df['p'] = (df['HOM_A1_CT'] * 2 + df['HET_A1_CT']) / (df['Total_Samples'] * 2)
df['MAF'] = df['p'].apply(lambda x: min(x, 1-x))

# MAF Groups
conditions = [
    (df['MAF'] <= 0.001),
    (df['MAF'] > 0.001) & (df['MAF'] <= 0.05),
    (df['MAF'] > 0.05)
]
choices = ['MAF <= 0.001', '0.001 < MAF <= 0.05', 'MAF > 0.05']
df['MAF_Group'] = np.select(conditions, choices, default='Unknown')

# Stats Calculation
df_valid_het = df.dropna(subset=[het_col]).copy()
mean_val = df_valid_het[het_col].mean()
median_val = df_valid_het[het_col].median()
min_val = df_valid_het[het_col].min()
max_val = df_valid_het[het_col].max()

# Downsample for plotting
if len(df_valid_het) > PLOT_SAMPLE_SIZE:
    df_plot_het = df_valid_het.sample(n=PLOT_SAMPLE_SIZE, random_state=42)
else:
    df_plot_het = df_valid_het

# Plot Het
plt.figure(figsize=(10, 6))
# X limit logic (99% quantile)
x_limit = df_valid_het[het_col].quantile(0.99)
if x_limit < 0.01: x_limit = max(0.01, df_valid_het[het_col].max() * 0.1) 
x_limit = min(x_limit * 1.2, 1.0) 

# Stacked Histogram by MAF
hue_order = ['MAF <= 0.001', '0.001 < MAF <= 0.05', 'MAF > 0.05']
palette_map = {'MAF <= 0.001': 'lightgray', 
               '0.001 < MAF <= 0.05': 'teal', 
               'MAF > 0.05': 'purple'}

sns.histplot(data=df_plot_het, x=het_col, hue='MAF_Group', hue_order=hue_order,
             bins=100, binrange=(0, x_limit), multiple='stack',
             palette=palette_map, edgecolor='white', alpha=0.8, kde=False)

plt.axvline(mean_val, color='red', linestyle='--', label=f'Mean: {mean_val:.4f}')
plt.axvline(median_val, color='green', linestyle='--', label=f'Median: {median_val:.4f}')

# Custom Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=palette_map[label], label=label) for label in hue_order]
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles=legend_elements + handles, loc='upper right', title='MAF Categories & Stats')

plt.yscale('log')
plt.title(f'Distribution of Variant Heterozygosity - Stacked by MAF')
plt.xlabel('Observed Heterozygosity')
plt.ylabel('Number of Variants (Log Scale)')
plt.xlim(0, x_limit)
plt.grid(True, which="major", axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(OUTPUT_HET_FILE, dpi=300)
print(f"Saved {OUTPUT_HET_FILE}")

# --- 1b. Heterozygosity Distribution (Linear Scale) ---
print("Generating Heterozygosity Distribution (Linear Scale)...")
OUTPUT_HET_LINEAR_FILE = "variant_het_dist_linear.png"

plt.figure(figsize=(10, 6))

sns.histplot(data=df_plot_het, x=het_col, hue='MAF_Group', hue_order=hue_order,
             bins=100, binrange=(0, x_limit), multiple='stack',
             palette=palette_map, edgecolor='white', alpha=0.8, kde=False)

plt.axvline(mean_val, color='red', linestyle='--', label=f'Mean: {mean_val:.4f}')
plt.axvline(median_val, color='green', linestyle='--', label=f'Median: {median_val:.4f}')

# Custom Legend
# (legend_elements is already defined in previous block)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles=legend_elements + handles, loc='upper right', title='MAF Categories & Stats')

# No LOG scale
plt.title(f'Distribution of Variant Heterozygosity - Stacked by MAF (Linear Scale)')
plt.xlabel('Observed Heterozygosity')
plt.ylabel('Number of Variants')
plt.xlim(0, x_limit)
plt.grid(True, axis='y', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig(OUTPUT_HET_LINEAR_FILE, dpi=300)
print(f"Saved {OUTPUT_HET_LINEAR_FILE}")

# ==========================================
# Part 2: Inbreeding Coefficient (F) Distribution
# ==========================================
print("\n--- 2. Generating Inbreeding Coefficient (F) Distribution ---")

# Calculate metrics from scratch to ensure F logic
df['Total_Samples'] = df['HOM_A1_CT'] + df['HET_A1_CT'] + df['TWO_AX_CT']
df_f = df[df['Total_Samples'] > 0].copy()

# Hobs
df_f['Hobs'] = df_f['HET_A1_CT'] / df_f['Total_Samples']
# p (Freq of A1)
df_f['p'] = (df_f['HOM_A1_CT'] * 2 + df_f['HET_A1_CT']) / (df_f['Total_Samples'] * 2)
# Hexp
df_f['Hexp'] = 2 * df_f['p'] * (1 - df_f['p'])
# Filter Hexp > 0
df_f = df_f[df_f['Hexp'] > 0].copy()
# F
df_f['F'] = 1 - (df_f['Hobs'] / df_f['Hexp'])
# MAF
df_f['MAF'] = df_f['p'].apply(lambda x: min(x, 1-x))

# --- MAF Grouping for Stacked Plot ---
conditions = [
    (df_f['MAF'] <= 0.001),
    (df_f['MAF'] > 0.001) & (df_f['MAF'] <= 0.05),
    (df_f['MAF'] > 0.05)
]
choices = ['MAF <= 0.001', '0.001 < MAF <= 0.05', 'MAF > 0.05']
df_f['MAF_Group'] = np.select(conditions, choices, default='Unknown')

# Stats
min_f = df_f['F'].min()
max_f = df_f['F'].max()
mean_f = df_f['F'].mean()
median_f_all = df_f['F'].median() # Global median (no filter)

# Reference Median (Used in papers)
f_filtered_series = df_f[(df_f['F'] > 0) & (df_f['MAF'] > 0.05)]['F']
if len(f_filtered_series) > 0:
    f_median_filtered = f_filtered_series.median()
    f_label = f'Filtered Median (F>0 & MAF>0.05): {f_median_filtered:.3f}'
else:
    f_median_filtered = None
    f_label = 'Filtered Median: N/A'
print(f"Valid Sites for F: {len(df_f)}")
print(f"F Stats: Min={min_f:.4f}, Max={max_f:.4f}, Mean={mean_f:.4f}, Median={median_f_all:.4f}")

# Plot F (Stacked by MAF Group)
if len(df_f) > PLOT_SAMPLE_SIZE:
    df_plot_f = df_f.sample(n=PLOT_SAMPLE_SIZE, random_state=42)
else:
    df_plot_f = df_f

plt.figure(figsize=(10, 6))

# Stacked Histogram
hue_order = ['MAF <= 0.001', '0.001 < MAF <= 0.05', 'MAF > 0.05']
palette_map = {'MAF <= 0.001': 'lightgray', 
               '0.001 < MAF <= 0.05': 'teal', 
               'MAF > 0.05': 'purple'}

sns.histplot(data=df_plot_f, x='F', hue='MAF_Group', hue_order=hue_order,
             bins=100, multiple='stack', kde=False, 
             palette=palette_map,
             binrange=(-1.5, 1.1), edgecolor='white', alpha=0.8)

# Keep original reference lines
plt.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5, label='F=0')
if 'median_f_all' in locals():
    plt.axvline(median_f_all, color='orange', linestyle='--', linewidth=2, label=f'Global Median: {median_f_all:.3f}')
if 'f_median_filtered' in locals() and f_median_filtered is not None:
    plt.axvline(f_median_filtered, color='red', linestyle='--', linewidth=2, label=f_label)

# Custom Legend to explain colors
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=palette_map[label], label=label) for label in hue_order]

# Combine with existing line handles (from axvline)
handles, labels = plt.gca().get_legend_handles_labels()
# Note: handles/labels from get_legend_handles_labels() contains the lines we just added.
# We prepend our color patches.
plt.legend(handles=legend_elements + handles, loc='upper left', title='MAF Categories & Stats')

plt.title('Distribution of Inbreeding Coefficient (F) - Stacked by MAF')
plt.xlabel('Inbreeding Coefficient (F = 1 - Hobs/Hexp)')
plt.ylabel('Count of Variants')
plt.xlim(-1.5, 1.1)
# plt.legend(...) handled above
plt.grid(True, axis='y', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig(OUTPUT_F_FILE, dpi=300)
print(f"Saved {OUTPUT_F_FILE}")

# --- 2b. Inbreeding Coefficient (F) Distribution (Log Scale) ---
print("Generating Inbreeding Coefficient (F) Distribution (Log Scale)...")
OUTPUT_F_LOG_FILE = "variant_f_dist_log.png"

plt.figure(figsize=(10, 6))

sns.histplot(data=df_plot_f, x='F', hue='MAF_Group', hue_order=hue_order,
             bins=100, multiple='stack', kde=False, 
             palette=palette_map,
             binrange=(-1.5, 1.1), edgecolor='white', alpha=0.8)

# Keep original reference lines
plt.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5, label='F=0')
if 'median_f_all' in locals():
    plt.axvline(median_f_all, color='orange', linestyle='--', linewidth=2, label=f'Global Median: {median_f_all:.3f}')
if 'f_median_filtered' in locals() and f_median_filtered is not None:
    plt.axvline(f_median_filtered, color='red', linestyle='--', linewidth=2, label=f_label)

# Custom Legend
# Reuse legend_elements from previous block
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles=legend_elements + handles, loc='upper left', title='MAF Categories & Stats')

plt.yscale('log')
plt.title('Distribution of Inbreeding Coefficient (F) - Stacked by MAF (Log Scale)')
plt.xlabel('Inbreeding Coefficient (F = 1 - Hobs/Hexp)')
plt.ylabel('Count of Variants (Log Scale)')
plt.xlim(-1.5, 1.1)

plt.grid(True, axis='y', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig(OUTPUT_F_LOG_FILE, dpi=300)
print(f"Saved {OUTPUT_F_LOG_FILE}")

# --- Part 5: Regression MAF vs Heterozygosity ---
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
    
    # Plot lines
    plt.plot(x_range, y_ols, color='blue', linewidth=2, linestyle='--', 
             label=f'OLS: {ols_eq}, $R^2$={ols_score:.3f}')
    plt.plot(x_range, y_huber, color='red', linewidth=2, 
             label=f'Huber: {huber_eq}, $R^2$={huber_score:.3f}')
    
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
if 'df_plot_f' in locals():
    plot_dual_regression(df_f, df_plot_f, 'MAF', 'Hobs', 'Minor Allele Frequency (MAF)', 'Observed Heterozygosity', 'reg_maf_vs_het.png')

# ==========================================
# Part 3: Inbreeding Coefficient for MAF > 0.001
# ==========================================
print("\n--- 3. Generating F Distribution (MAF > 0.001) ---")
OUTPUT_F_001_FILE = "variant_f_dist_maf0.001.png"

# Filter by MAF
df_f001 = df_f[df_f['MAF'] > 0.001].copy()
if len(df_f001) > 0:
    # Stats
    min_f001 = df_f001['F'].min()
    max_f001 = df_f001['F'].max()
    mean_f001 = df_f001['F'].mean()
    median_f001 = df_f001['F'].median()
    
    print(f"MAF>0.001 Sites: {len(df_f001)}")
    print(f"F Stats (MAF>0.001): Min={min_f001:.4f}, Max={max_f001:.4f}, Mean={mean_f001:.4f}, Median={median_f001:.4f}")

    # Downsample
    if len(df_f001) > PLOT_SAMPLE_SIZE:
        df_plot_f001 = df_f001.sample(n=PLOT_SAMPLE_SIZE, random_state=42)
    else:
        df_plot_f001 = df_f001

    plt.figure(figsize=(10, 6))
    
    # Histogram
    sns.histplot(df_plot_f001['F'], bins=100, kde=True, color='teal', edgecolor='white', alpha=0.7, binrange=(-1.5, 1.1))

    # Lines
    plt.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5, label='F=0')
    plt.axvline(median_f001, color='orange', linestyle='--', linewidth=2, label=f'Median: {median_f001:.3f}')

    # Add Min/Max to Legend
    plt.plot([], [], ' ', label=f'Min: {min_f001:.2f}')
    plt.plot([], [], ' ', label=f'Max: {max_f001:.2f}')
    plt.plot([], [], ' ', label=f'Mean: {mean_f001:.2f}')
    plt.plot([], [], ' ', label=f'N: {len(df_f001)}')

    plt.title('Distribution of Inbreeding Coefficient (MAF > 0.001)')
    plt.xlabel('Inbreeding Coefficient')
    plt.ylabel('Count')
    plt.xlim(-1.5, 1.1)
    plt.legend(loc='upper left')
    plt.grid(True, axis='y', linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTPUT_F_001_FILE, dpi=300)
    print(f"Saved {OUTPUT_F_001_FILE}")
else:
    print("No sites with MAF > 0.001 found.")


# ==========================================
# Part 4: Inbreeding Coefficient for MAF > 0.05
# ==========================================
print("\n--- 4. Generating F Distribution (MAF > 0.05) ---")
OUTPUT_F_05_FILE = "variant_f_dist_maf0.05.png"

# Filter by MAF
df_f05 = df_f[df_f['MAF'] > 0.05].copy()
if len(df_f05) > 0:
    # Stats
    min_f05 = df_f05['F'].min()
    max_f05 = df_f05['F'].max()
    mean_f05 = df_f05['F'].mean()
    median_f05 = df_f05['F'].median()
    
    print(f"MAF>0.05 Sites: {len(df_f05)}")
    print(f"F Stats (MAF>0.05): Min={min_f05:.4f}, Max={max_f05:.4f}, Mean={mean_f05:.4f}, Median={median_f05:.4f}")

    # Downsample
    if len(df_f05) > PLOT_SAMPLE_SIZE:
        df_plot_f05 = df_f05.sample(n=PLOT_SAMPLE_SIZE, random_state=42)
    else:
        df_plot_f05 = df_f05

    plt.figure(figsize=(10, 6))
    
    # Histogram
    sns.histplot(df_plot_f05['F'], bins=100, kde=True, color='purple', edgecolor='white', alpha=0.7, binrange=(-1.5, 1.1))

    # Lines
    plt.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5, label='F=0')
    plt.axvline(median_f05, color='orange', linestyle='--', linewidth=2, label=f'Median: {median_f05:.3f}')

    # Add Min/Max to Legend
    plt.plot([], [], ' ', label=f'Min: {min_f05:.2f}')
    plt.plot([], [], ' ', label=f'Max: {max_f05:.2f}')
    plt.plot([], [], ' ', label=f'Mean: {mean_f05:.2f}')
    plt.plot([], [], ' ', label=f'N: {len(df_f05)}')

    plt.title('Distribution of Inbreeding Coefficient (MAF > 0.05)')
    plt.xlabel('Inbreeding Coefficient')
    plt.ylabel('Count')
    plt.xlim(-1.5, 1.1)
    plt.legend(loc='upper left')
    plt.grid(True, axis='y', linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTPUT_F_05_FILE, dpi=300)
    print(f"Saved {OUTPUT_F_05_FILE}")
else:
    print("No sites with MAF > 0.05 found.")