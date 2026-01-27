import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
import os
from sklearn.linear_model import HuberRegressor, LinearRegression

# Settings
SCOUNT_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.sample_stats.scount"
MISSING_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.missing.smiss"
GROUP_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/sample_group.txt"
OUTPUT_PREFIX = "ref_ibs_analysis"

print("Processing Reference IBS Analysis...")

# 1. Read .scount data
print("Reading .scount data...")
try:
    df_scount = pd.read_csv(SCOUNT_FILE, sep=r'\s+')
except Exception as e:
    print(f"Failed to read scount file: {e}")
    sys.exit(1)

# Check columns
required_cols = ['#IID', 'HOM_REF_CT', 'HET_SNP_CT', 'HOM_ALT_SNP_CT']
if not all(col in df_scount.columns for col in required_cols):
    print(f"Error: Missing columns in scount. Found: {df_scount.columns}")
    sys.exit(1)

# Calculate IBS_ref
# Total alleles = 2 * Total_Sites (assuming diploid)
# Ref Alleles = 2 * HOM_REF + 1 * HET
# IBS_ref = Ref Alleles / Total Alleles
# Formula: (2*HOM_REF_CT + 1*HET_SNP_CT) / (2 * (HOM_REF_CT + HOM_ALT_SNP_CT + HET_SNP_CT))

df_scount['Total_Sites'] = df_scount['HOM_REF_CT'] + df_scount['HOM_ALT_SNP_CT'] + df_scount['HET_SNP_CT']
# Filter out 0 sites if any
df_scount = df_scount[df_scount['Total_Sites'] > 0].copy()

df_scount['Total_Alleles'] = 2 * df_scount['Total_Sites']
df_scount['Ref_Alleles'] = 2 * df_scount['HOM_REF_CT'] + df_scount['HET_SNP_CT']
df_scount['IBS_Ref'] = df_scount['Ref_Alleles'] / df_scount['Total_Alleles']

df_ibs = df_scount[['#IID', 'IBS_Ref']].rename(columns={'#IID': 'Sample'})

# Stats
mean_ibs = df_ibs['IBS_Ref'].mean()
median_ibs = df_ibs['IBS_Ref'].median()
print(f"IBS_Ref - Mean: {mean_ibs:.4f}, Median: {median_ibs:.4f}")

# 2. Read .smiss data
print("Reading .smiss data...")
try:
    df_missing = pd.read_csv(MISSING_FILE, sep=r'\s+')
except Exception as e:
    print(f"Failed to read smiss file: {e}")
    sys.exit(1)

# Handle IID header variations
if '#IID' in df_missing.columns:
    iid_col = '#IID'
elif 'IID' in df_missing.columns:
    iid_col = 'IID'
else:
    print(f"Error: Missing IID column in smiss. Found: {df_missing.columns}")
    sys.exit(1)
    
missing_col = 'F_MISS'
if missing_col not in df_missing.columns:
     print(f"Error: Missing '{missing_col}' in smiss. Found: {df_missing.columns}")
     sys.exit(1)

df_missing = df_missing[[iid_col, missing_col]].rename(columns={iid_col: 'Sample', missing_col: 'Missing_Rate'})

# 3. Merge
print("Merging IBS and Missing Rate data...")
df_merged = pd.merge(df_ibs, df_missing, on='Sample', how='inner')
print(f"Merged samples: {len(df_merged)}")

# 4. Integrate Group Info
print("Reading and integrating sample group info...")
if os.path.exists(GROUP_FILE):
    try:
        df_group = pd.read_csv(GROUP_FILE, sep=r'\s+', header=None, names=['Sample', 'Group'])
        # Deduplicate
        df_group = df_group.drop_duplicates(subset=['Sample'])
        
        # Merge
        df_merged = pd.merge(df_merged, df_group, on='Sample', how='left')
        df_merged['Group'] = df_merged['Group'].fillna('Unknown')
    except Exception as e:
        print(f"Warning: Failed to process group file: {e}")
        df_merged['Group'] = 'Unknown'
else:
    print(f"Warning: Group file not found at {GROUP_FILE}")
    df_merged['Group'] = 'Unknown'

print(f"Groups found: {df_merged['Group'].unique()}")

# Set Seaborn Style
sns.set_theme(style="ticks")

# ==========================================
# Plot 1: Distribution
# ==========================================
def plot_distribution(data, col, mean_val, median_val, title, filename, log_scale=False):
    plt.figure(figsize=(10, 6))
    
    # Histogram with Stacked Groups
    sns.histplot(data=data, x=col, hue='Group', multiple='stack', 
                 bins=100, linewidth=0.1)
    
    # Statistics Lines
    plt.axvline(x=mean_val, color='blue', linestyle='--', linewidth=1, label=f'Mean: {mean_val:.4f}')
    plt.axvline(x=median_val, color='orange', linestyle=':', linewidth=1.5, label=f'Median: {median_val:.4f}')
    
    plt.title(title, fontsize=15)
    plt.xlabel("IBS with Reference Genome", fontsize=12)
    
    if log_scale:
        plt.yscale('log')
        plt.ylabel("Count (Log Scale)", fontsize=12)
    else:
        plt.ylabel("Count", fontsize=12)
        
    plt.legend(handles=lines+labels, loc='upper right', title='Groups & Stats') 
    
    plt.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    plt.close()

print("Generating Distribution Plots...")
# Note: seaborn histplot automatically adds a legend for 'hue'. 
# We just need to ensure our axvline stats labels are also included or don't conflict.
# seaborn's legend usually overwrites or is overwritten.
# Improved strategy: Draw histplot, get legend, draw lines, add separate legend or combine.
# seaborn histplot(multiple='stack') might put legend outside.

# Modified plot_distribution to handle legend better
def plot_distribution_refined(data, col, mean_val, median_val, title, filename, log_scale=False):
    plt.figure(figsize=(10, 6))
    
    # Histogram with Stacked Groups
    ax = sns.histplot(data=data, x=col, hue='Group', multiple='stack', 
                 bins=100, linewidth=0.1)
    
    # Statistics Lines
    l1 = plt.axvline(x=mean_val, color='blue', linestyle='--', linewidth=1, label=f'Mean: {mean_val:.4f}')
    l2 = plt.axvline(x=median_val, color='orange', linestyle=':', linewidth=1.5, label=f'Median: {median_val:.4f}')
    
    plt.title(title, fontsize=15)
    plt.xlabel("IBS with Reference Genome", fontsize=12)
    
    if log_scale:
        plt.yscale('log')
        plt.ylabel("Count (Log Scale)", fontsize=12)
    else:
        plt.ylabel("Count", fontsize=12)
        
    # Combine legends: Seaborn creates a legend for Groups automatically.
    # We want to keep it and add our stats lines.
    # Get handles and labels from seaborn plot
    # Note: ax.get_legend_handles_labels() might return empty if seaborn didn't attach them to ax in standard way
    # or if we removed the legend already.
    # Proper way for seaborn histplot hue legend:
    
    # 1. Get current handles/labels (including hue patches)
    # Seaborn 0.11+ histplot puts legend on the figure or axes.
    # If a legend exists, we can extract its handles.
    handlers = []
    labels = []
    
    if ax.get_legend():
        handlers = ax.get_legend().legend_handles
        labels = [t.get_text() for t in ax.get_legend().get_texts()]
        ax.get_legend().remove()
    
    # 2. Add our stats lines
    handlers.extend([l1, l2])
    labels.extend([l1.get_label(), l2.get_label()])
        
    plt.legend(handles=handlers, labels=labels, title='Groups & Stats', loc='upper right')

    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"Saved {filename}")
    plt.close()

plot_distribution_refined(df_merged, 'IBS_Ref', mean_ibs, median_ibs, 
                  "Distribution of IBS with Reference", f"{OUTPUT_PREFIX}_dist_ibs.png")

plot_distribution_refined(df_merged, 'IBS_Ref', mean_ibs, median_ibs, 
                  "Distribution of IBS with Reference (Log Scale)", f"{OUTPUT_PREFIX}_dist_ibs_log.png", 
                  log_scale=True)

# ==========================================
# Plot 2: Regression
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
    
    # Plotting: JointGrid
    # Filter valid data including Group for plotting
    local_valid_plot = df[[x_col, y_col, 'Group']].replace([np.inf, -np.inf], np.nan).dropna()

    # Use seaborn JointGrid to plot scatter + marginals
    g = sns.JointGrid(data=local_valid_plot, x=x_col, y=y_col, height=10, ratio=5)

    # 1. Main Scatter Plot
    sns.scatterplot(data=local_valid_plot, x=x_col, y=y_col, hue='Group',
                    alpha=0.5, s=15, edgecolor='none', ax=g.ax_joint, legend=True)
    
    # 2. Add Regression Lines to Main Plot
    x_min, x_max = local_valid[x_col].min(), local_valid[x_col].max()
    x_span = x_max - x_min
    if x_span == 0: x_span = 0.01
    x_range = np.linspace(x_min - 0.05*x_span, x_max + 0.05*x_span, 100).reshape(-1, 1)
    
    y_ols = ols.predict(x_range)
    y_huber = huber.predict(x_range)
    
    g.ax_joint.plot(x_range, y_ols, color='blue', linewidth=2, linestyle='--', 
                    label=f'OLS: {ols_eq}, $R^2$={ols_score:.3f}')
    g.ax_joint.plot(x_range, y_huber, color='red', linewidth=2, 
                    label=f'Huber: {huber_eq}, $R^2$={huber_score:.3f}')
    
    # Move legend (Scatter + Lines)
    # The scatter plot legend is automatically handled by sns.scatterplot(legend=True)
    # We need to add lines to it or create a new one.
    handles, labels = g.ax_joint.get_legend_handles_labels()
    # Add regression lines manually since they are not in seaborn hue data
    # (Note: Standard Matplot lines are not automatically picked up by seaborn hue legend unless we merge)
    
    # Simple fix: Let matplotlib legend handle everything
    # But seaborn scatter legend is complex (markers).
    # We can plot lines and just call legend() - it might overwrite or combine.
    g.ax_joint.legend(loc='lower left') # This typically shows regression lines only if seaborn legend is removed or overridden
    
    # 3. Marginals (Stacked Histograms)
    sns.histplot(data=local_valid_plot, x=x_col, hue='Group', multiple='stack', 
                 bins=50, linewidth=0.1, ax=g.ax_marg_x, legend=False)
    sns.histplot(data=local_valid_plot, y=y_col, hue='Group', multiple='stack', 
                 bins=50, linewidth=0.1, ax=g.ax_marg_y, legend=False)

    g.ax_joint.set_xlabel(x_label)
    g.ax_joint.set_ylabel(y_label)
    g.ax_joint.set_xlim(x_min - 0.05*x_span, x_max + 0.05*x_span)
    
    plt.suptitle(f'{y_label} vs {x_label}', y=1.02, fontsize=16)
    
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved {filename}")
    plt.close()

print("Generating Regression Plot (IBS_Ref vs Missing)...")
plot_dual_regression(df_merged, 'Missing_Rate', 'IBS_Ref', 
                     'Missing Rate', 'IBS with Reference Genome', 
                     f"{OUTPUT_PREFIX}_reg_miss_vs_ibs.png")

print("Analysis Complete.")
