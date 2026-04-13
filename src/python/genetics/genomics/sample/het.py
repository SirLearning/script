from .sample_utils import load_df_from_plink2
from infra.utils.io import load_df_from_tsv, load_thresholds, save_df_to_tsv, save_thresholds
from infra.utils.graph import plot_distribution_with_stats, plot_regression_comparison
from infra.utils import plot_stacked_distribution, plot_joint_regression
from genetics.germplasm.sample import anno_group
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
        return None

    # Calculate Heterozygosity Rate
    df_het['Total_Called'] = df_het['HOM_REF_CT'] + df_het['HET_SNP_CT'] + df_het['HOM_ALT_SNP_CT']
    df_het = df_het[df_het['Total_Called'] > 0].copy()
    df_het['Het_Rate'] = df_het['HET_SNP_CT'] / df_het['Total_Called']
    
    save_df_to_tsv(
        df=df_het[['Sample', 'Het_Rate']],
        output_file=f"{output_prefix}.rate.info.tsv"
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
        save_thresholds(thresholds, f"{output_prefix}.th.tsv")

    # 4. Plots
    # 4.1 Het Distribution
    plot_distribution_with_stats(
        data=df_merged,
        col='Het_Rate',
        title="Distribution of Sample Heterozygosity",
        filename=f"{output_prefix}.dist.png",
        mean_val=thresholds['Mean_Het'],
        median_val=thresholds['Median_Het'],
        x_label="Sample Heterozygosity Rate"
    )
    
    # 4.2 Het Distribution (Log)
    plot_distribution_with_stats(
        data=df_merged,
        col='Het_Rate',
        title="Distribution of Sample Heterozygosity (Log Scale)",
        filename=f"{output_prefix}.dist_log.png",
        mean_val=thresholds['Mean_Het'],
        median_val=thresholds['Median_Het'],
        x_label="Sample Heterozygosity Rate",
        log_scale=True
    )
    
    # 4.3 Regression: Missing vs Het
    plot_regression_comparison(
        df=df_merged,
        x_col='F_MISS',
        y_col='Het_Rate',
        x_label='Missing Rate',
        y_label='Heterozygosity Rate',
        filename=f"{output_prefix}.reg_vs_het.png",
        title="Heterozygosity vs Missing Rate"
    )

    print("Analysis Complete.")


def ana_heterozygosity_with_group(
    scount_file,
    smiss_file, 
    group_file,
    output_prefix,
    tsv_file=None,
    th_tsv_file=None
):
    """
    Main orchestrator for heterozygosity analysis incorporating group information.
    Generates stacked distributions and joint regression plots utilizing group labels.
    """
    print("Loading Heterozygosity and Missing Rate data...")
    if tsv_file and os.path.exists(tsv_file):
        print(f"Loading precomputed heterozygosity from {tsv_file}...")
        df_merged = load_df_from_tsv(tsv_file, sep='\t')
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

        # 3. Integrate Group Info
        if group_file:
            df_merged = anno_group(df_merged, group_file)
        
        save_df_to_tsv(df_merged, f"{output_prefix}.group.info.tsv")

    if th_tsv_file and os.path.exists(th_tsv_file):
        print(f"Loading precomputed thresholds from {th_tsv_file}...")
        thresholds = load_thresholds(th_tsv_file)
    else:
        # 4. Calculate Stats & Thresholds
        het_stats = {
                'Mean_Het': df_merged['Het_Rate'].mean(),
                'Median_Het': df_merged['Het_Rate'].median(),
                'Std_Het': df_merged['Het_Rate'].std()
            }
        print(f"Heterozygosity - Mean: {het_stats['Mean_Het']:.4f}, Median: {het_stats['Median_Het']:.4f}")
        thresholds = het_stats.copy()
        thresholds['Upper_Threshold_3SD'] = het_stats['Mean_Het'] + 3 * het_stats['Std_Het']
        thresholds['Lower_Threshold_3SD'] = max(0, het_stats['Mean_Het'] - 3 * het_stats['Std_Het'])
        save_thresholds(thresholds, f"{output_prefix}.group.th.tsv")

    # Set Seaborn Style
    sns.set_theme(style="ticks")

    # 5. Plots
    # 5.1 Het Distribution
    print("Generating Distribution Plots...")
    plot_stacked_distribution(
        df=df_merged,
        col='Het_Rate',
        group_col='Group',
        title="Distribution of Sample Heterozygosity",
        filename=f"{output_prefix}.dist.png",
        mean_val=thresholds['Mean_Het'],
        median_val=thresholds['Median_Het'],
        x_label="Sample Heterozygosity Rate"
    )
    
    # 5.2 Het Distribution (Log)
    plot_stacked_distribution(
        df=df_merged,
        col='Het_Rate',
        group_col='Group',
        title="Distribution of Sample Heterozygosity (Log Scale)",
        filename=f"{output_prefix}.dist_log.png",
        mean_val=thresholds['Mean_Het'],
        median_val=thresholds['Median_Het'],
        x_label="Sample Heterozygosity Rate",
        log_scale=True
    )
    
    # 5.3 Regression: Missing vs Het
    print("Generating Regression Plot (Heterozygosity vs Missing)...")
    plot_joint_regression(
        df=df_merged,
        x_col='Het_Rate',
        y_col='F_MISS',
        group_col='Group',
        x_label='Heterozygosity Rate',
        y_label='Missing Rate',
        filename=f"{output_prefix}.reg_vs_miss.png",
        title="Heterozygosity vs Missing Rate"
    )

    print("Analysis Complete.")

# To allow importing as a module or running as script
if __name__ == "__main__":
    # Example hardcoded path for testing if run directly without args
    # But ideally should use argparse
    pass


