from genetics.germplasm.sample.process import load_df_from_plink2
from infra.utils.io import load_df_from_space_sep, load_df_from_tsv, load_thresholds, save_df_to_tsv, save_thresholds
from infra.utils.graph import plot_distribution_with_stats, plot_regression_comparison
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

def load_scount_data(
    scount_file,
    output_prefix
):
    """
    Reads PLINK .scount file and calculates heterozygosity rate.
    Returns DataFrame with columns ['Sample', 'Het_Rate'] and stats dictionary.
    """
    print(f"Reading .scount data from {scount_file}...")
    df_het = load_df_from_plink2(scount_file)

    required_cols = ['HOM_REF_CT', 'HET_SNP_CT', 'HOM_ALT_SNP_CT']

    if not all(col in df_het.columns for col in required_cols):
        print(f"[Error] Missing genotype count columns. Found: {df_het.columns}")
        return None, None

    # Calculate Heterozygosity Rate
    df_het['Total_Called'] = df_het['HOM_REF_CT'] + df_het['HET_SNP_CT'] + df_het['HOM_ALT_SNP_CT']
    df_het = df_het[df_het['Total_Called'] > 0].copy()
    df_het['Het_Rate'] = df_het['HET_SNP_CT'] / df_het['Total_Called']
    
    save_df_to_tsv(
        df=df_het[['Sample', 'Het_Rate']],
        output_file=f"{output_prefix}.het_rates.tsv"
    )
    
    return df_het[['Sample', 'Het_Rate']]


def ana_heterozygosity(
    scount_file,
    smiss_file, 
    output_prefix,
    tsv_file=None,
    th_tsv_file=None
):
    """
    Main orchestrator for heterozygosity analysis.
    """
    print("Loading Heterozygosity and Missing Rate data...")
    if tsv_file and os.path.exists(tsv_file):
        print(f"Loading precomputed heterozygosity from {tsv_file}...")
        df_het = load_df_from_tsv(tsv_file, sep='\t')
    else:
        # 1. Load Data
        df_het = load_scount_data(scount_file, output_prefix)
        if df_het is None: return
        df_miss = load_df_from_plink2(smiss_file)
        if df_miss is None: return
        # 2. Merge
        print("Merging Heterozygosity and Missing Rate data...")
        df_merged = pd.merge(df_het, df_miss, on='Sample', how='inner')
        print(f"Merged samples: {len(df_merged)}")

    if th_tsv_file and os.path.exists(th_tsv_file):
        print(f"Loading precomputed thresholds from {th_tsv_file}...")
        thresholds = load_thresholds(th_tsv_file)
    else:
        # 3. Calculate Stats & Thresholds
        het_stats = {
                'Mean_Het': df_merged['Het_Rate'].mean(),
                'Median_Het': df_merged['Het_Rate'].median(),
                'Std_Het': df_merged['Het_Rate'].std()
            }
        print(f"Heterozygosity - Mean: {het_stats['Mean_Het']:.4f}, Median: {het_stats['Median_Het']:.4f}")
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
        mean_val=thresholds['Mean_Het'],
        median_val=thresholds['Median_Het'],
        x_label="Sample Heterozygosity Rate",
        color='forestgreen'
    )
    
    # 4.2 Het Distribution (Log)
    plot_distribution_with_stats(
        data=df_merged,
        col='Het_Rate',
        title="Distribution of Sample Heterozygosity (Log Scale)",
        filename=f"{output_prefix}_dist_het_log.png",
        mean_val=thresholds['Mean_Het'],
        median_val=thresholds['Median_Het'],
        x_label="Sample Heterozygosity Rate",
        log_scale=True,
        color='forestgreen'
    )
    
    # 4.3 Regression: Missing vs Het
    plot_regression_comparison(
        df=df_merged,
        x_col='F_MISS',
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


