import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.linear_model import LinearRegression, HuberRegressor
import sys

# Configuration
GCOUNT_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.variant_stats.gcount"
OUTPUT_PREFIX = "mac_maf_analysis"

def analyze_mac_maf():
    print(f"Reading {GCOUNT_FILE}...")
    try:
        df = pd.read_csv(GCOUNT_FILE, sep=r'\s+')
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # User provided: #CHROM ID REF ALT HOM_REF_CT HET_REF_ALT_CTS TWO_ALT_GENO_CTS HAP_REF_CT HAP_ALT_CTS MISSING_CT
    required = ['HET_REF_ALT_CTS', 'TWO_ALT_GENO_CTS', 'HOM_REF_CT', 'HAP_REF_CT', 'HAP_ALT_CTS']
    if not all(col in df.columns for col in required):
        print(f"Error: Missing required columns. Found: {list(df.columns)}")
        return
    
    # 1. Calculate Allele Counts
    print("Calculating allele counts...")
    # Alt Count = 1*Het + 2*HomAlt + 1*HapAlt
    df['Alt_Count'] = df['HET_REF_ALT_CTS'] + (df['TWO_ALT_GENO_CTS'] * 2) + df['HAP_ALT_CTS']
    # Ref Count = 1*Het + 2*HomRef + 1*HapRef
    df['Ref_Count'] = df['HET_REF_ALT_CTS'] + (df['HOM_REF_CT'] * 2) + df['HAP_REF_CT']
    
    df['Total_Alleles'] = df['Alt_Count'] + df['Ref_Count']
    
    # 2. Calculate MAC (Minor Allele Count) & MAF (Minor Allele Frequency)
    df['MAC'] = df[['Alt_Count', 'Ref_Count']].min(axis=1)
    # MAF calculation
    df['MAF'] = df['MAC'] / df['Total_Alleles']
    
    # --- Task 1: Analyze MAC = 1 ---
    print("\n" + "="*40)
    print("Analysis of Sites with MAC = 1")
    print("="*40)
    
    mac1_df = df[df['MAC'] == 1]
    count_mac1 = len(mac1_df)
    print(f"Number of variants with MAC = 1: {count_mac1}")
    
    if count_mac1 > 0:
        print("\nStatistics of Total Allele Count for these variants:")
        # Display describe() stats
        stats = mac1_df['Total_Alleles'].describe()
        print(stats)
        
        print("\nDetail of Total Alleles for MAC=1 sites (Top 20 examples):")
        # Format output nicely
        print(mac1_df[['ID', 'MAC', 'Total_Alleles', 'MISSING_CT']].head(20).to_string(index=False))
        
        # Check if totals vary significantly
        unique_totals = mac1_df['Total_Alleles'].unique()
        if len(unique_totals) < 20:
            print(f"\nUnique Total Allele counts observed: {sorted(unique_totals)}")
        else:
            print(f"\nNumber of unique Total Allele counts observed: {len(unique_totals)}")
    else:
        print("No variants with MAC=1 found.")

    # --- Task 2: Plot MAC (0-380) vs MAF ---
    print("\n" + "="*40)
    print("Generating MAC (0-380) vs MAF Plot")
    print("="*40)
    
    # Filter Data
    df_plot = df[(df['MAC'] >= 0) & (df['MAC'] <= 380)].copy()
    
    # Remove NaN or Inf if any
    df_plot = df_plot.replace([np.inf, -np.inf], np.nan).dropna(subset=['MAC', 'MAF'])
    
    if len(df_plot) < 10:
        print("Not enough data points in range 0-380 to plot.")
        return
        
    # print(f"Plotting {len(df_plot)} variants...")
    
    # plot_dual_regression_styled(df_plot, 'MAC', 'MAF', 
    #                            'Minor Allele Count (MAC)', 'Minor Allele Frequency (MAF)', 
    #                            f"{OUTPUT_PREFIX}_mac_0_380_vs_maf.png")

    # --- Task 3: Plot MAC vs Total Alleles (0-80 & 0-380) ---
    print("\n" + "="*40)
    print("Generating MAC vs Total Alleles Plots")
    print("="*40)
    
    for limit in [80, 380]:
        range_label = f"0-{limit}"
        print(f"Plotting MAC vs Total Alleles (Range {range_label})...")
        
        # Filter Data
        df_sub = df[(df['MAC'] >= 0) & (df['MAC'] <= limit)].copy()
        
        # Filter NaNs
        df_sub = df_sub.replace([np.inf, -np.inf], np.nan).dropna(subset=['MAC', 'Total_Alleles'])
        
        if len(df_sub) < 10:
            print(f"Not enough data points in range {range_label}.")
            continue
            
        plot_dual_regression_styled(df_sub, 'MAC', 'Total_Alleles', 
                                   'Minor Allele Count (MAC)', 'Total Allele Number', 
                                   f"{OUTPUT_PREFIX}_mac_{range_label}_vs_total_alleles.png")

def plot_dual_regression_styled(df, x_col, y_col, x_label, y_label, filename):
    # Set style
    sns.set_theme(style="ticks")
    
    X = df[x_col].values.reshape(-1, 1)
    y = df[y_col].values
    
    # 1. Regression Models
    # OLS
    ols = LinearRegression()
    ols.fit(X, y)
    ols_slope = ols.coef_[0]
    ols_intercept = ols.intercept_
    ols_score = ols.score(X, y)
    ols_eq = f"y = {ols_slope:.2e}x + {ols_intercept:.2e}"
    
    # Huber
    huber = HuberRegressor(epsilon=1.35)
    huber.fit(X, y)
    huber_slope = huber.coef_[0]
    huber_intercept = huber.intercept_
    huber_score = huber.score(X, y)
    huber_eq = f"y = {huber_slope:.2e}x + {huber_intercept:.2e}"
    
    # 2. JointGrid Plotting
    # We don't have 'Group' for variants usually, so we'll use a single color but keep the marginals
    g = sns.JointGrid(data=df, x=x_col, y=y_col, height=10, ratio=5)
    
    # Main Scatter
    # Using alpha to handle density
    g.plot_joint(sns.scatterplot, alpha=0.3, s=15, edgecolor='none', color='royalblue')
    
    # Regression Lines
    x_min, x_max = df[x_col].min(), df[x_col].max()
    x_range = np.linspace(x_min, x_max, 100).reshape(-1, 1)
    
    y_ols = ols.predict(x_range)
    y_huber = huber.predict(x_range)
    
    g.ax_joint.plot(x_range, y_ols, color='orange', linewidth=2, linestyle='--', 
                    label=f'OLS: {ols_eq}\n$R^2$={ols_score:.3f}')
    g.ax_joint.plot(x_range, y_huber, color='red', linewidth=2, 
                    label=f'Huber: {huber_eq}\n$R^2$={huber_score:.3f}')
    
    g.ax_joint.legend(loc='upper left')
    
    # Marginals (Histograms)
    g.plot_marginals(sns.histplot, kde=False, bins=50, color='royalblue', linewidth=0.1)
    
    # Labels
    g.ax_joint.set_xlabel(x_label, fontsize=12)
    g.ax_joint.set_ylabel(y_label, fontsize=12)
    g.ax_joint.grid(True, linestyle='--', alpha=0.3)
    
    plt.suptitle(f'{y_label} vs {x_label} (Range 0-380)', y=1.02, fontsize=16)
    
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Chart saved to: {filename}")
    plt.close()

if __name__ == "__main__":
    analyze_mac_maf()
