import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
import os
from infra.plot_utils import plot_distribution_with_stats, plot_regression_comparison
from infra.threshold_utils import save_thresholds

def load_scount_data(scount_file):
    """
    Reads PLINK .scount file and calculates heterozygosity rate.
    Returns DataFrame with columns ['Sample', 'Het_Rate'] and stats dictionary.
    """
    print(f"Reading .scount data from {scount_file}...")
    try:
        df_het = pd.read_csv(scount_file, sep=r'\s+')
    except Exception as e:
        print(f"[Error] Failed to read scount file: {e}")
        return None, None

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

def load_smiss_data(smiss_file):
    """
    Reads PLINK .smiss file.
    Returns DataFrame with columns ['Sample', 'Missing_Rate'].
    """
    print(f"Reading .smiss data from {smiss_file}...")
    try:
        df_mid = pd.read_csv(smiss_file, sep=r'\s+')
    except Exception as e:
        print(f"[Error] Failed to read smiss file: {e}")
        return None

    # Normalizing columns
    if '#IID' in df_mid.columns:
        df_mid = df_mid.rename(columns={'#IID': 'Sample'})
    elif 'IID' in df_mid.columns:
        df_mid = df_mid.rename(columns={'IID': 'Sample'})
    else:
        print(f"[Error] Missing IID column in smiss. Found: {df_mid.columns}")
        return None

    if 'F_MISS' in df_mid.columns:
        df_mid = df_mid.rename(columns={'F_MISS': 'Missing_Rate'})
    else:
        print(f"[Error] Missing F_MISS column in smiss. Found: {df_mid.columns}")
        return None
        
    return df_mid[['Sample', 'Missing_Rate']]


def run_heterozygosity_analysis(
    scount_file, 
    smiss_file, 
    output_prefix
):
    """
    Main orchestrator for heterozygosity analysis.
    """
    # 1. Load Data
    df_het, het_stats = load_scount_data(scount_file)
    if df_het is None: return

    df_miss = load_smiss_data(smiss_file)
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


