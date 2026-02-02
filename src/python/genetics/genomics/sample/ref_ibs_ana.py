import pandas as pd
import sys
import os
import seaborn as sns

# Add src/python to path if needed to find infra
current_dir = os.path.dirname(os.path.abspath(__file__))
# Navigate up to src/python
# Current: src/python/genetics/genomics/sample/ (4 levels deep from src/python ?)
# src/python is parent of genetics.
# genetics/genomics/sample is 3 levels.
# So ../../../ is correct for src/python.
src_python_dir = os.path.abspath(os.path.join(current_dir, "../../../"))
if src_python_dir not in sys.path:
    sys.path.append(src_python_dir)

from infra.plot_utils import plot_stacked_distribution, plot_joint_regression

def ref_ibs_with_missing(
    scount_file, 
    missing_file, 
    group_file, 
    output_prefix
):
    """
    Analyzes Reference IBS and Missing Rate from PLINK .scount and .smiss files.
    Generates distribution and regression plots.
    """
    print("Processing Reference IBS Analysis...")

    # 1. Read .scount data
    print(f"Reading .scount data from {scount_file}...")
    try:
        df_scount = pd.read_csv(scount_file, sep=r'\s+')
    except Exception as e:
        print(f"Failed to read scount file: {e}")
        return

    # Check columns
    required_cols = ['#IID', 'HOM_REF_CT', 'HET_SNP_CT', 'HOM_ALT_SNP_CT']
    if not all(col in df_scount.columns for col in required_cols):
        print(f"Error: Missing columns in scount. Found: {df_scount.columns}")
        return

    # Calculate IBS_ref
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
    print(f"Reading .smiss data from {missing_file}...")
    try:
        df_missing = pd.read_csv(missing_file, sep=r'\s+')
    except Exception as e:
        print(f"Failed to read smiss file: {e}")
        return

    # Handle IID header variations
    if '#IID' in df_missing.columns:
        iid_col = '#IID'
    elif 'IID' in df_missing.columns:
        iid_col = 'IID'
    else:
        print(f"Error: Missing IID column in smiss. Found: {df_missing.columns}")
        return
        
    missing_col = 'F_MISS'
    if missing_col not in df_missing.columns:
        print(f"Error: Missing '{missing_col}' in smiss. Found: {df_missing.columns}")
        return

    df_missing = df_missing[[iid_col, missing_col]].rename(columns={iid_col: 'Sample', missing_col: 'Missing_Rate'})

    # 3. Merge
    print("Merging IBS and Missing Rate data...")
    df_merged = pd.merge(df_ibs, df_missing, on='Sample', how='inner')
    print(f"Merged samples: {len(df_merged)}")

    # 4. Integrate Group Info
    print("Reading and integrating sample group info...")
    if os.path.exists(group_file):
        try:
            df_group = pd.read_csv(group_file, sep=r'\s+', header=None, names=['Sample', 'Group'])
            # Deduplicate
            df_group = df_group.drop_duplicates(subset=['Sample'])
            
            # Merge
            df_merged = pd.merge(df_merged, df_group, on='Sample', how='left')
            df_merged['Group'] = df_merged['Group'].fillna('Unknown')
        except Exception as e:
            print(f"Warning: Failed to process group file: {e}")
            df_merged['Group'] = 'Unknown'
    else:
        print(f"Warning: Group file not found at {group_file}")
        df_merged['Group'] = 'Unknown'

    print(f"Groups found: {df_merged['Group'].unique()}")

    # Set Seaborn Style
    sns.set_theme(style="ticks")

    # ==========================================
    # Plot 1: Distribution
    # ==========================================
    print("Generating Distribution Plots...")
    
    plot_stacked_distribution(
        data=df_merged, 
        col='IBS_Ref', 
        group_col='Group',
        title="Distribution of IBS with Reference", 
        filename=f"{output_prefix}_dist_ibs.png",
        mean_val=mean_ibs, 
        median_val=median_ibs,
        x_label="IBS with Reference Genome"
    )

    plot_stacked_distribution(
        data=df_merged, 
        col='IBS_Ref', 
        group_col='Group',
        title="Distribution of IBS with Reference (Log Scale)", 
        filename=f"{output_prefix}_dist_ibs_log.png", 
        mean_val=mean_ibs, 
        median_val=median_ibs,
        x_label="IBS with Reference Genome",
        log_scale=True
    )

    # ==========================================
    # Plot 2: Regression
    # ==========================================
    print("Generating Regression Plot (IBS_Ref vs Missing)...")
    
    plot_joint_regression(
        df=df_merged, 
        x_col='Missing_Rate', 
        y_col='IBS_Ref', 
        group_col='Group',
        x_label='Missing Rate', 
        y_label='IBS with Reference Genome', 
        filename=f"{output_prefix}_reg_miss_vs_ibs.png"
    )

    print("Analysis Complete.")

if __name__ == "__main__":
    # Settings
    SCOUNT_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.sample_stats.scount"
    MISSING_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.missing.smiss"
    GROUP_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/sample_group.txt"
    OUTPUT_PREFIX = "ref_ibs_analysis"
    
    ref_ibs_with_missing(SCOUNT_FILE, MISSING_FILE, GROUP_FILE, OUTPUT_PREFIX)
