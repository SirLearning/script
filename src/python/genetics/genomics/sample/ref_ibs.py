from germplasm import integrate_group_info
from genetics.genomics.sample import load_smiss
from infra.utils import plot_stacked_distribution, plot_joint_regression, load_df_from_space_sep, save_sample_df_to_tsv
import pandas as pd
import seaborn as sns

def load_and_calculate_ibs(scount_file):
    """
    Reads .scount file and calculates IBS with reference.
    Returns DataFrame with columns: ['Sample', 'IBS_Ref']
    """
    print(f"[Info] Reading .scount file: {scount_file}...")
    try:
        df_scount = load_df_from_space_sep(scount_file)
    except Exception as e:
        print(f"[Error] Failed to read scount file: {e}")
        return None

    # Check columns
    required_cols = ['#IID', 'HOM_REF_CT', 'HET_SNP_CT', 'HOM_ALT_SNP_CT']
    if not all(col in df_scount.columns for col in required_cols):
        print(f"[Error] Missing columns in scount. Found: {df_scount.columns}")
        return None

    # Calculate IBS_ref
    df_scount['Total_Sites'] = df_scount['HOM_REF_CT'] + df_scount['HOM_ALT_SNP_CT'] + df_scount['HET_SNP_CT']
    # Filter out 0 sites if any
    df_scount = df_scount[df_scount['Total_Sites'] > 0].copy()

    df_scount['Total_Alleles'] = 2 * df_scount['Total_Sites']
    df_scount['Ref_Alleles'] = 2 * df_scount['HOM_REF_CT'] + df_scount['HET_SNP_CT']
    df_scount['IBS_Ref'] = df_scount['Ref_Alleles'] / df_scount['Total_Alleles']

    df_ibs = df_scount[['#IID', 'IBS_Ref']].rename(columns={'#IID': 'Sample'})
    return df_ibs


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

    # 1. Read .scount data & Calculate IBS
    df_ibs = load_and_calculate_ibs(scount_file)
    if df_ibs is None: return

    # Save IBS Data
    ibs_output_file = f"{output_prefix}.ibs.tsv"
    save_sample_df_to_tsv(df_ibs, ibs_output_file)

    # Stats
    mean_ibs = df_ibs['IBS_Ref'].mean()
    median_ibs = df_ibs['IBS_Ref'].median()
    print(f"IBS_Ref - Mean: {mean_ibs:.4f}, Median: {median_ibs:.4f}")

    # 2. Read .smiss data
    print(f"Reading .smiss data from {missing_file}...")
    df_missing = load_smiss(missing_file)
    
    # 3. Merge
    print("Merging IBS and Missing Rate data...")
    df_merged = pd.merge(df_ibs, df_missing, on='Sample', how='inner')
    print(f"Merged samples: {len(df_merged)}")

    # 4. Integrate Group Info
    df_merged = integrate_group_info(group_file, df_merged)
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


def ref_ibs_vs_mapping(
    mapping_file, 
    scount_file, 
    group_file, 
    output_prefix
):
    print("Processing Ref IBS vs Mapping Rate Analysis...")
    
    # 1. Read Mapping Rate Data
    print(f"Reading Mapping Rate file: {mapping_file}...")
    try:
        df_map = pd.read_csv(mapping_file, sep='\t')
    except Exception as e:
        print(f"Error reading mapping file: {e}")
        return
        
    if 'TaxaID' not in df_map.columns or 'Mapping_Rate_Pct' not in df_map.columns:
        print(f"Error: Required columns ('TaxaID', 'Mapping_Rate_Pct') not found in mapping file.")
        print(f"Columns found: {list(df_map.columns)}")
        return
        
    df_map = df_map[['TaxaID', 'Mapping_Rate_Pct']].rename(columns={'TaxaID': 'Sample'})
    
    # 2. Read .scount Data & Calculate IBS
    df_ibs = load_and_calculate_ibs(scount_file)
    if df_ibs is None: return
    
    # 3. Merge Datasets
    print("Merging datasets...")
    df_merged = pd.merge(df_map, df_ibs, on='Sample', how='inner')
    print(f"Matched samples: {len(df_merged)}")
    
    if len(df_merged) == 0:
        print("Error: No intersecting samples found between mapping file and scount file.")
        return
        
    # 4. Integrate Group Info
    df_merged = integrate_group_info(group_file, df_merged)
        
    # Print basic stats
    print(f"Correlation (Pearson): {df_merged['Mapping_Rate_Pct'].corr(df_merged['IBS_Ref']):.4f}")
    
    # Set style
    sns.set_theme(style="ticks")
    
    # 5. Plot
    output_filename = f"{output_prefix}_reg_ref_ibs_vs_map.png"
    plot_joint_regression(
        df=df_merged, 
        x_col='Mapping_Rate_Pct', 
        y_col='IBS_Ref', 
        group_col='Group',
        x_label='Mapping Rate (%)', 
        y_label='IBS with Reference Genome', 
        filename=output_filename
    )
                         
    print("Analysis Complete.")


if __name__ == "__main__":
    # Settings
    SCOUNT_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.sample_stats.scount"
    MISSING_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.missing.smiss"
    GROUP_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/sample_group.txt"
    OUTPUT_PREFIX = "ref_ibs_analysis"
    
    ref_ibs_with_missing(SCOUNT_FILE, MISSING_FILE, GROUP_FILE, OUTPUT_PREFIX)
