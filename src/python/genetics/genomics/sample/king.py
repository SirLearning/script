from .sample_utils import load_df_from_king, load_df_from_plink2
import numpy as np
import pandas as pd
from infra.utils.io import save_df_to_tsv, save_thresholds
from infra.utils.graph import plot_distribution_with_stats, plot_regression_comparison
import matplotlib.pyplot as plt
import seaborn as sns

def ana_king_kinship(
    input_file,
    output_prefix="king_kinship"
):
    print(f"Processing KING Kinship: {input_file}")
    
    # 1. Read Data
    df = load_df_from_king(input_file)
    if df is None: return

    # Check columns
    if 'KINSHIP' not in df.columns:
        print("Error: KINSHIP column missing.")
        return
    
    # 2. Basic Stats
    stats_dict = {
        'Mean_Kinship': df['KINSHIP'].mean(),
        'Median_Kinship': df['KINSHIP'].median(),
        'Max_Kinship': df['KINSHIP'].max(),
        'Min_Kinship': df['KINSHIP'].min()
    }
    save_thresholds(stats_dict, f"{output_prefix}.stats.tsv")
    print(stats_dict)
    
    # 3. Scatter Plot: IBS0 vs KINSHIP (if IBS0 exists, typical for .kin0)
    if 'IBS0' in df.columns:
        # Base Plot
        plot_regression_comparison(
            df=df, x_col='IBS0', y_col='KINSHIP',
            x_label='Fraction of Zero IBS Sites (IBS0)', y_label='KING Kinship Coefficient',
            filename=f"{output_prefix}.ibs0_vs_kinship.png",
            title="Kinship vs IBS0"
        )
        # Custom annotated plot could be added here similar to tmp script if needed
        # For now, general regression/scatter plot is good.

    # 4. Distribution Plots
    
    # 4.1 Full Range
    plot_distribution_with_stats(
        data=df, col='KINSHIP',
        title="Distribution of KING Kinship (Full Range)",
        filename=f"{output_prefix}.dist.png",
        mean_val=stats_dict['Mean_Kinship'],
        median_val=stats_dict['Median_Kinship'],
        x_label="Kinship Coefficient",
        xlim=(-0.5, 0.5) # Focus on typical range
    )

    # 4.2 Zoomed Range of Interest (-1 to 0.5) if many negatives
    df_zoom = df[(df['KINSHIP'] >= -1) & (df['KINSHIP'] <= 0.5)]
    if not df_zoom.empty:
        plot_distribution_with_stats(
            data=df_zoom, col='KINSHIP',
            title="Distribution of KING Kinship (Zoomed)",
            filename=f"{output_prefix}.dist.zoom.png",
            mean_val=df_zoom['KINSHIP'].mean(),
            median_val=df_zoom['KINSHIP'].median(),
            x_label="Kinship Coefficient",
            xlim=(-1, 0.5)
        )

def ana_derived_het(
    input_file,
    output_prefix="derived_het"
):
    print(f"Processing Derived Heterozygosity (Hi+Hj): {input_file}")
    
    # 1. Read Data (Expecting .kin0 usually)
    df = load_df_from_king(input_file)
    if df is None: return
    
    # Check Cols
    required = ['HETHET', 'IBS0', 'KINSHIP']
    if not all(col in df.columns for col in required):
        print(f"Error: Missing columns for Derived Het calculation: {required}")
        return
        
    # 2. Calc (HETHET - 2*IBS0) / KINSHIP
    epsilon = 1e-9
    df_calc = df[abs(df['KINSHIP']) > epsilon].copy()
    
    df_calc['Hi_Hj'] = (df_calc['HETHET'] - 2 * df_calc['IBS0']) / df_calc['KINSHIP']
    
    # Filter reasonable range?
    df_calc = df_calc[(df_calc['Hi_Hj'] >= -2) & (df_calc['Hi_Hj'] <= 2)]
    
    # Stats
    mean_val = df_calc['Hi_Hj'].mean()
    median_val = df_calc['Hi_Hj'].median()
    
    save_df_to_tsv(
        df_calc[['Sample1', 'Sample2', 'Hi_Hj']], 
        f"{output_prefix}.info.tsv"
    )

    # 3. Plot
    plot_distribution_with_stats(
        data=df_calc, col='Hi_Hj',
        title="Distribution of Derived $H_i + H_j$",
        filename=f"{output_prefix}.dist.png",
        mean_val=mean_val, median_val=median_val,
        x_label="Derived $H_i + H_j$",
        xlim=(-2, 2)
    )

def ana_het_vs_max_kinship(
    scount_file, 
    king_file, 
    output_prefix="het_vs_max_kin"
):
    """
    Compares Heterozygosity Rate (from .scount) vs Max Kinship (from .kin0/.king).
    """
    print(f"Processing Het vs Max Kinship...")
    
    # 1. Load Het
    df_het = load_df_from_plink2(scount_file)
    if df_het is None: return
    
    if 'Het_Rate' not in df_het.columns:
        # Calculate if missing (failsafe)
        if 'HET_SNP_CT' in df_het.columns:
             total = df_het['HOM_REF_CT'] + df_het['HET_SNP_CT'] + df_het['HOM_ALT_SNP_CT']
             df_het['Het_Rate'] = df_het['HET_SNP_CT'] / total
        else:
             print("Error: Could not determine Het Rate")
             return

    # 2. Load King & Find Max
    df_kin = load_df_from_king(king_file)
    if df_kin is None: return
    
    # Stack Sample1 and Sample2 to find max per sample
    k1 = df_kin[['Sample1', 'KINSHIP']].rename(columns={'Sample1': 'Sample'})
    k2 = df_kin[['Sample2', 'KINSHIP']].rename(columns={'Sample2': 'Sample'})
    df_max_kin = pd.concat([k1, k2]).groupby('Sample')['KINSHIP'].max().reset_index()
    df_max_kin.rename(columns={'KINSHIP': 'Max_Kinship'}, inplace=True)
    
    # 3. Merge
    df_merged = pd.merge(df_het, df_max_kin, on='Sample', how='left').fillna(0) # Fill 0 if no kinship pairs
    
    # 4. Plot
    plt.figure(figsize=(10, 7))
    sns.set_style("whitegrid")

    sns.scatterplot(
        data=df_merged, x='Het_Rate', y='Max_Kinship', 
        alpha=0.5, s=25, edgecolor='w', color='teal'
    )
    
    # Reference Lines
    plt.axhline(y=0.354, color='darkred', linestyle='--', label='Duplicate > 0.354')
    
    het_mean = df_merged['Het_Rate'].mean()
    het_std = df_merged['Het_Rate'].std()
    het_limit = het_mean + 3 * het_std
    plt.axvline(x=het_limit, color='orange', linestyle='--', label=f'High Het (> {het_limit:.4f})')
    
    plt.title('Heterozygosity vs. Max Kinship')
    plt.xlabel('Heterozygosity Rate')
    plt.ylabel('Maximum Kinship Coefficient')
    plt.legend()
    
    plt.savefig(f"{output_prefix}.png", dpi=300)
    print(f"Saved {output_prefix}.png")
    plt.close()

if __name__ == "__main__":
    # Simple CLI dispatch
    import sys
    if len(sys.argv) < 2:
        print("Usage: python king.py <command> [args]")
        print("Commands: kinship, derived, het_vs_kin")
        sys.exit(1)
        
    cmd = sys.argv[1]
    
    if cmd == "kinship":
        ana_king_kinship(sys.argv[2], sys.argv[3] if len(sys.argv)>3 else "king_ana")
    elif cmd == "derived":
        ana_derived_het(sys.argv[2], sys.argv[3] if len(sys.argv)>3 else "derived_het")
    elif cmd == "het_vs_kin":
        ana_het_vs_max_kinship(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv)>4 else "het_vs_kin")
    else:
        print(f"Unknown command: {cmd}")
