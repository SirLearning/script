import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

# Settings
INPUT_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.king.kin0"
OUTPUT_PREFIX = "derived_het"

print("Processing Derived Heterozygosity (Hi+Hj) Analysis...")

# 1. Read Data
try:
    print(f"Reading {INPUT_FILE}...")
    df = pd.read_csv(INPUT_FILE, sep=r'\s+')
except Exception as e:
    print(f"Failed to read file: {e}")
    sys.exit(1)

# Check columns
required_cols = ['HETHET', 'IBS0', 'KINSHIP']
if not all(col in df.columns for col in required_cols):
    print(f"Error: Missing required columns {required_cols}. Found: {df.columns}")
    sys.exit(1)

# 2. Calculate Hi + Hj = (HETHET - 2*IBS0) / KINSHIP
# Handle KINSHIP = 0 to avoid division by zero (though KING usually non-zero, but let's be safe)
# Filter out rows where Kinship is extremely close to 0 or results in inf/nan?
# Formula from user: (HETHET - 2*IBS0) / KINSHIP
print("Calculating Hi + Hj...")

# Filter valid data for calculation
# Avoid division by zero
epsilon = 1e-9
df_calc = df[abs(df['KINSHIP']) > epsilon].copy()
dropped = len(df) - len(df_calc)
if dropped > 0:
    print(f"Dropped {dropped} rows with KINSHIP ~ 0")

df_calc['Hi_Hj'] = (df_calc['HETHET'] - 2 * df_calc['IBS0']) / df_calc['KINSHIP']

# Inspect values - there might be extreme outliers if Kinship is very small
print(f"Hi_Hj - Min: {df_calc['Hi_Hj'].min()}, Max: {df_calc['Hi_Hj'].max()}")

# Filter extremely unreasonable values? 
# Theoretically H is a rate 0-1, so Hi+Hj should be roughly 0-2 (or slightly more if noise, or negative?)
# If KINSHIP is negative (unrelated), the formula might behave weirdly if meant for relateds?
# User prompt implies calculating it directly. Assuming it's valid for the dataset.
# However, if Kinship is negative (common for unrelated in KING), the denominator is negative. 
# HETHET and IBS0 are positive counts/fractions. 
# NOTE: IBS0 is fraction. HETHET is usually count or fraction? 
# Let's check typical kin0 format.
# PLINK2 .kin0: IBS0 is fraction. KINSHIP is coefficient. 
# HETHET: "Proportion of sites that are HET in both samples" (fraction).
# So all are fractions. Result should be reasonable.

# Stats
mean_val = df_calc['Hi_Hj'].mean()
median_val = df_calc['Hi_Hj'].median()
print(f"Hi_Hj - Mean: {mean_val:.4f}, Median: {median_val:.4f}")


# Set Seaborn Style
sns.set_theme(style="ticks")

# ==========================================
# Plot Distribution Function
# ==========================================
def plot_distribution(data_series, title, filename_suffix, log_scale=False, bins=100, xlim=None, binrange=None):
    plt.figure(figsize=(8, 6))
    
    # Histogram
    sns.histplot(data_series, bins=bins, binrange=binrange, color='forestgreen', edgecolor='white', linewidth=0.1)
    
    # Statistics Lines
    plt.axvline(x=mean_val, color='blue', linestyle='--', linewidth=1, label=f'Mean: {mean_val:.4f}')
    plt.axvline(x=median_val, color='orange', linestyle=':', linewidth=1.5, label=f'Median: {median_val:.4f}')
    
    # Custom Annotation Lines for specific plots (e.g. 0-0.1 range)
    # Check if we are plotting the zoom range where these make sense
    if xlim == (0, 0.1):
        custom_lines = [0.00133, 0.00333, 0.00476, 0.008] # 0.333 is outside 0.1 range
        colors = ['red', 'purple', 'brown', 'cyan', 'magenta']
        for i, val in enumerate(custom_lines):
            if xlim[0] <= val <= xlim[1]:
                 plt.axvline(x=val, color=colors[i % len(colors)], linestyle='-.', linewidth=1, alpha=0.7, label=f'{val}')
    
    # Add Min/Max to Legend
    current_min = data_series.min()
    current_max = data_series.max()
    plt.plot([], [], ' ', label=f'Min: {current_min:.4f}')
    plt.plot([], [], ' ', label=f'Max: {current_max:.4f}')
    
    # Settings
    plt.title(title, fontsize=15)
    plt.xlabel("Derived $H_i + H_j$ = (HetHet - 2*IBS0) / Kinship", fontsize=12)
    
    if log_scale:
        plt.yscale('log')
        plt.ylabel("Count (Log Scale)", fontsize=12)
    else:
        plt.ylabel("Count", fontsize=12)
        
    if xlim:
        plt.xlim(xlim)
        
    plt.legend(loc='upper right')
    plt.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_PREFIX}_{filename_suffix}.png", dpi=300)
    print(f"Saved {OUTPUT_PREFIX}_{filename_suffix}.png")
    plt.close()

# ==========================================
# 3. Plots
# ==========================================

# # 1. Linear Scale
# print("Generating Distribution (Linear)...")
# plot_distribution(df_calc['Hi_Hj'], "Distribution of Derived $H_i + H_j$", "dist_linear")

# # 2. Log Scale
# print("Generating Distribution (Log Count)...")
# plot_distribution(df_calc['Hi_Hj'], "Distribution of Derived $H_i + H_j$ (Log Scale)", "dist_log", log_scale=True)

# # 3. Zoomed Range -2 to 2 (Linear)
# print("Generating Distribution (Zoom -2 to 2)...")
# df_zoom = df_calc[(df_calc['Hi_Hj'] >= -2) & (df_calc['Hi_Hj'] <= 2)]
# plot_distribution(df_zoom['Hi_Hj'], "Distribution of Derived $H_i + H_j$ (-2 to 2)", "dist_zoom_linear", 
#                   xlim=(-2, 2), bins=100, binrange=(-2, 2))

# # 4. Zoomed Range -2 to 2 (Log Scale)
# print("Generating Distribution (Zoom -2 to 2, Log)...")
# plot_distribution(df_zoom['Hi_Hj'], "Distribution of Derived $H_i + H_j$ (-2 to 2, Log Scale)", "dist_zoom_log", 
#                   log_scale=True, xlim=(-2, 2), bins=100, binrange=(-2, 2))

# # 5. Zoomed Range 0 to 0.5 (Linear)
# print("Generating Distribution (Zoom 0 to 0.5)...")
# df_zoom = df_calc[(df_calc['Hi_Hj'] >= 0) & (df_calc['Hi_Hj'] <= 0.5)]
# plot_distribution(df_zoom['Hi_Hj'], "Distribution of Derived $H_i + H_j$ (0 to 0.5)", "dist_zoom_0_0.5_linear", 
#                   xlim=(0, 0.5), bins=500, binrange=(0, 0.5))

# # 6. Zoomed Range 0 to 0.5 (Log Scale)
# print("Generating Distribution (Zoom 0 to 0.5, Log)...")
# plot_distribution(df_zoom['Hi_Hj'], "Distribution of Derived $H_i + H_j$ (0 to 0.5, Log Scale)", "dist_zoom_0_0.5_log", 
#                   log_scale=True, xlim=(0, 0.5), bins=500, binrange=(0, 0.5))

# 7. Zoomed Range 0 to 0.1 (Linear)
print("Generating Distribution (Zoom 0 to 0.1)...")
df_zoom = df_calc[(df_calc['Hi_Hj'] >= 0) & (df_calc['Hi_Hj'] <= 0.1)]
plot_distribution(df_zoom['Hi_Hj'], "Distribution of Derived $H_i + H_j$ (0 to 0.1)", "dist_zoom_0_0.1_linear", 
                  xlim=(0, 0.1), bins=200, binrange=(0, 0.1))

# 8. Zoomed Range 0 to 0.1 (Log Scale)
print("Generating Distribution (Zoom 0 to 0.1, Log)...")
plot_distribution(df_zoom['Hi_Hj'], "Distribution of Derived $H_i + H_j$ (0 to 0.1, Log Scale)", "dist_zoom_0_0.1_log", 
                  log_scale=True, xlim=(0, 0.1), bins=200, binrange=(0, 0.1))

print("Analysis Complete.")
