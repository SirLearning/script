from genetics.genomics.sample.smiss import load_smiss
from infra.utils.io import load_df_from_space_sep, load_df_from_tsv, save_df_to_tsv, save_thresholds
from infra.utils.graph import plot_distribution_with_stats, plot_regression_comparison
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

def load_scount_data(scount_file):
    """
    Reads PLINK .scount file and calculates heterozygosity rate.
    Returns DataFrame with columns ['Sample', 'Het_Rate'] and stats dictionary.
    """
    print(f"Reading .scount data from {scount_file}...")
    df_het = load_df_from_space_sep(scount_file)

    required_cols = ['HOM_REF_CT', 'HET_SNP_CT', 'HOM_ALT_SNP_CT']
    # Check for IID column
    if '#IID' in df_het.columns:
        df_het = df_het.rename(columns={'#IID': 'Sample'})
    elif 'IID' in df_het.columns:
        df_het = df_het.rename(columns={'IID': 'Sample'})
    else:
        print(f"[Warning] Missing IID/Sample column in scount. Found: {df_het.columns}")
        # Try to proceed if strictly formatted or fail
        # return None, None
        pass

    if not all(col in df_het.columns for col in required_cols):
        print(f"[Error] Missing genotype count columns. Found: {df_het.columns}")
        return None, None

    # Calculate Heterozygosity Rate
    df_het['Total_Called'] = df_het['HOM_REF_CT'] + df_het['HET_SNP_CT'] + df_het['HOM_ALT_SNP_CT']
    
    # Avoid division by zero
    df_het = df_het[df_het['Total_Called'] > 0].copy()
    
    df_het['Het_Rate'] = df_het['HET_SNP_CT'] / df_het['Total_Called']
    
    save_df_to_tsv(
        df=df_het[['Sample', 'Het_Rate']],
        output_file=scount_file.replace('.scount', '.het_rates.tsv')
    )
    
    # Stats
    mean_het = df_het['Het_Rate'].mean()
    median_het = df_het['Het_Rate'].median()
    std_het = df_het['Het_Rate'].std()
    
    stats = {
        'Mean_Het': mean_het,
        'Median_Het': median_het,
        'Std_Het': std_het
    }
    
    return df_het[['Sample', 'Het_Rate']], stats


def run_heterozygosity_analysis(
    scount_file,
    het_tsv,
    smiss_file, 
    output_prefix
):
    """
    Main orchestrator for heterozygosity analysis.
    """
    # 1. Load Data
    print("Loading Heterozygosity and Missing Rate data...")
    if het_tsv and os.path.exists(het_tsv):
        print(f"Loading precomputed heterozygosity from {het_tsv}...")
        df_het = load_df_from_tsv(het_tsv, sep='\t')
        het_stats = {
            'Mean_Het': df_het['Het_Rate'].mean(),
            'Median_Het': df_het['Het_Rate'].median(),
            'Std_Het': df_het['Het_Rate'].std()
        }
    else:
        df_het, het_stats = load_scount_data(scount_file)
        if df_het is None: return

    df_miss = load_smiss(smiss_file)
    if df_miss is None: return

    print(f"Heterozygosity - Mean: {het_stats['Mean_Het']:.4f}, Median: {het_stats['Median_Het']:.4f}")

    # 2. Merge
    print("Merging Heterozygosity and Missing Rate data...")
    df_merged = pd.merge(df_het, df_miss, on='Sample', how='inner')
    print(f"Merged samples: {len(df_merged)}")

    # 3. Calculate Thresholds (e.g., Mean + 3SD for Het)
    thresholds = het_stats.copy()
    thresholds['Upper_Threshold_3SD'] = het_stats['Mean_Het'] + 3 * het_stats['Std_Het']
    thresholds['Lower_Threshold_3SD'] = max(0, het_stats['Mean_Het'] - 3 * het_stats['Std_Het'])
    
    save_thresholds(thresholds, f"{output_prefix}.thresholds.tsv")

    # 4. Plots
    # 4.1 Het Distribution
    plot_distribution_with_stats(
        data=df_merged,
        col='Het_Rate',
        title="Distribution of Sample Heterozygosity",
        filename=f"{output_prefix}_dist_het.png",
        mean_val=het_stats['Mean_Het'],
        median_val=het_stats['Median_Het'],
        x_label="Sample Heterozygosity Rate",
        color='forestgreen'
    )
    
    # 4.2 Het Distribution (Log)
    plot_distribution_with_stats(
        data=df_merged,
        col='Het_Rate',
        title="Distribution of Sample Heterozygosity (Log Scale)",
        filename=f"{output_prefix}_dist_het_log.png",
        mean_val=het_stats['Mean_Het'],
        median_val=het_stats['Median_Het'],
        x_label="Sample Heterozygosity Rate",
        log_scale=True,
        color='forestgreen'
    )
    
    # 4.3 Regression: Missing vs Het
    plot_regression_comparison(
        df=df_merged,
        x_col='Missing_Rate',
        y_col='Het_Rate',
        x_label='Missing Rate',
        y_label='Heterozygosity Rate',
        filename=f"{output_prefix}_reg_miss_vs_het.png",
        title="Heterozygosity vs Missing Rate"
    )

    print("Analysis Complete.")

# To allow importing as a module or running as script
if __name__ == "__main__":
    # Example hardcoded path for testing if run directly without args
    # But ideally should use argparse
    pass


