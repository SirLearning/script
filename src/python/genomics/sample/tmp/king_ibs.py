import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

# Settings
INPUT_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.king.kin0"
OUTPUT_PREFIX = "king_ibs"

print("Processing KING IBS/Kinship Analysis...")

# 1. Read Data
try:
    print(f"Reading {INPUT_FILE}...")
    df = pd.read_csv(INPUT_FILE, sep=r'\s+')
except Exception as e:
    print(f"Failed to read file: {e}")
    sys.exit(1)

# Check columns (Expected: IBS0, KINSHIP)
if 'IBS0' not in df.columns or 'KINSHIP' not in df.columns:
    print(f"Error: Missing required columns (IBS0, KINSHIP). Found: {df.columns}")
    sys.exit(1)

# Calculate Stats
mean_val = df['KINSHIP'].mean()
median_val = df['KINSHIP'].median()
print(f"Mean: {mean_val:.4f}, Median: {median_val:.4f}")

# Set Seaborn Style
sns.set_theme(style="ticks")

# ==========================================
# 1. Scatter Plot: Kinship vs. IBS0 (Base)
# ==========================================
print("Generating Scatter Plot (Base)...")
plt.figure(figsize=(8, 6))

sns.scatterplot(data=df, x='IBS0', y='KINSHIP', 
                color='#1f77b4', alpha=0.4, s=10, edgecolor=None)

# Reference Lines
lines = [0.0442, 0.0884, 0.177, 0.354]
for y in lines:
    plt.axhline(y=y, color='red', linestyle='--', linewidth=1, alpha=0.5)

plt.text(0.005, 0.36, "Duplicate/MZ Twin", color='red', fontsize=9, ha='left')

plt.title("Kinship vs. IBS0 Scatter Plot", fontsize=15)
plt.xlabel("Fraction of Zero IBS Sites (IBS0)", fontsize=12)
plt.ylabel("KING Kinship Coefficient", fontsize=12)
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTPUT_PREFIX}_plot_reform_py.png", dpi=300)
print(f"Saved {OUTPUT_PREFIX}_plot_reform_py.png")
plt.close()

# ==========================================
# 2. Scatter Plot: Cutoff (-1 to 0.5)
# ==========================================
print("Generating Scatter Plot (Cutoff)...")
plt.figure(figsize=(8, 6))

# Split Data
df_main = df[df['KINSHIP'] >= -1].copy()
df_low = df[df['KINSHIP'] < -1].copy()

# Plot Main Points
sns.scatterplot(data=df_main, x='IBS0', y='KINSHIP', 
                color='#1f77b4', alpha=0.4, s=10, edgecolor=None) # muted blue

# Plot Low Points (Collapsed to -1.05)
if not df_low.empty:
    plt.scatter(df_low['IBS0'], [-1.05]*len(df_low), 
                color='#8b0000', alpha=0.1, s=5, marker='|') # dark red

# Reference Lines
for y in lines:
    plt.axhline(y=y, color='red', linestyle='--', linewidth=1, alpha=0.5)

# Boundary Line
plt.axhline(y=-1, color='gray', linestyle='-', linewidth=1)

# Annotations
plt.text(0.002, 0.4, ">0.354: Duplicate", color='red', fontsize=8, ha='left')
plt.text(0.002, 0.26, "0.177-0.354: 1st Degree", color='red', fontsize=8, ha='left')
plt.text(0.002, 0.13, "0.0884-0.177: 2nd Degree", color='red', fontsize=8, ha='left')
plt.text(0.002, 0.06, "0.0442-0.0884: 3rd Degree", color='red', fontsize=8, ha='left')

plt.title("Kinship vs. IBS0 (Zoomed -1 to 0.5)\nIntervals annotated; points < -1 shown at y=-1.05", fontsize=15, pad=15)
# Removed suptitle

plt.xlabel("Fraction of Zero IBS Sites (IBS0)", fontsize=12)
plt.ylabel("KING Kinship Coefficient", fontsize=12)
plt.ylim(-1.1, 0.5)
plt.yticks([-1, -0.5, 0, 0.0442, 0.0884, 0.177, 0.354, 0.5])
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTPUT_PREFIX}_plot_cutoff_py.png", dpi=300)
print(f"Saved {OUTPUT_PREFIX}_plot_cutoff_py.png")
plt.close()

# ==========================================
# Distribution Plots Helper
# ==========================================
def plot_distribution(data_series, title, filename_suffix, log_scale=False, xlim=None, bins=100):
    plt.figure(figsize=(8, 6))
    
    # Histogram
    sns.histplot(data_series, bins=bins, color='forestgreen', edgecolor='white', linewidth=0.1)
    
    # Statistics Lines (Global mean/median)
    plt.axvline(x=mean_val, color='blue', linestyle='--', linewidth=1, label=f'Mean: {mean_val:.4f}')
    plt.axvline(x=median_val, color='orange', linestyle=':', linewidth=1.5, label=f'Median: {median_val:.4f}')
    
    # Settings
    plt.title(title, fontsize=15)
    plt.xlabel("KING Kinship Coefficient", fontsize=12)
    
    if log_scale:
        plt.yscale('log')
        plt.ylabel("Count (Log Scale)", fontsize=12)
    else:
        plt.ylabel("Count", fontsize=12)
        
    if xlim:
        plt.xlim(xlim)
        
    plt.legend(loc='upper left')
    plt.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_PREFIX}_{filename_suffix}.png", dpi=300)
    print(f"Saved {OUTPUT_PREFIX}_{filename_suffix}.png")
    plt.close()

# # ==========================================
# # 3. Distribution: Full Range
# # ==========================================
# print("Generating Distribution (Full Range)...")
# plot_distribution(df['KINSHIP'], "Distribution of KING Kinship (Full Range)", "dist_full_py")

# # ==========================================
# # 4. Distribution: Full Range (Log Scale)
# # ==========================================
# print("Generating Distribution (Full Range, Log)...")
# plot_distribution(df['KINSHIP'], "Distribution of KING Kinship (Full Range, Log Scale)", "dist_full_log_py", log_scale=True)

# # ==========================================
# # 5. Distribution: Zoom (-1 to 0.5)
# # ==========================================
# print("Generating Distribution (Zoom -1 to 0.5)...")
# # Filter data for histplot to ensure bins are calculated correctly for the range
# df_zoom = df[(df['KINSHIP'] >= -1) & (df['KINSHIP'] <= 0.5)]
# plot_distribution(df_zoom['KINSHIP'], "Distribution of KING Kinship (-1 to 0.5)", "dist_zoom_py", xlim=(-1, 0.5), bins=100)

# # ==========================================
# # 6. Distribution: Neg12 (-12 to 0.5)
# # ==========================================
# print("Generating Distribution (Zoom -12 to 0.5)...")
# df_neg12 = df[(df['KINSHIP'] >= -12) & (df['KINSHIP'] <= 0.5)]
# plot_distribution(df_neg12['KINSHIP'], "Distribution of KING Kinship (-12 to 0.5)", "dist_neg12_py", xlim=(-12, 0.5), bins=200)

print("Done.")
