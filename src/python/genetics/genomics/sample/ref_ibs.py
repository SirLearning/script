from genetics.genomics.sample.mr import load_df_from_idxstats
from genetics.germplasm.sample.anno import anno_group
from genetics.germplasm.sample.process import load_df_from_plink2
from infra.utils import plot_stacked_distribution, plot_joint_regression
from infra.utils.io import load_df_from_tsv, save_df_to_tsv
import pandas as pd
import seaborn as sns
import os
import argparse

def load_and_calculate_ibs(scount_file):
    """
    Reads .scount file and calculates IBS with reference.
    Returns DataFrame with columns: ['Sample', 'IBS_Ref']
    """
    print(f"[Info] Reading .scount file: {scount_file}...")
    df_scount = load_df_from_plink2(scount_file)

    # Check columns
    required_cols = ['HOM_REF_CT', 'HET_SNP_CT', 'HOM_ALT_SNP_CT']
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

    df_ibs = df_scount[['Sample', 'IBS_Ref']]
    return df_ibs


def ana_ref_ibs(
    scount_file, 
    missing_file, 
    group_file, 
    output_prefix,
    tsv_file=None
):
    """
    Analyzes Reference IBS and Missing Rate from PLINK .scount and .smiss files.
    Generates distribution and regression plots.
    """
    print("Processing Reference IBS Analysis...")

    if tsv_file and os.path.exists(tsv_file):
        print(f"Loading precomputed data from: {tsv_file}...")
        df_merged = load_df_from_tsv(tsv_file)
        if 'IBS_Ref' not in df_merged.columns:
            print("[Error] 'IBS_Ref' column missing in loaded TSV.")
            return
    else:
        # 1. Read .scount data & Calculate IBS
        df_ibs = load_and_calculate_ibs(scount_file)
        if df_ibs is None: return

        # 2. Read .smiss data
        print(f"Reading .smiss data from {missing_file}...")
        df_missing = load_df_from_plink2(missing_file)
        if df_missing is None: return
        
        # 3. Merge
        print("Merging IBS and Missing Rate data...")
        df_merged = pd.merge(df_ibs, df_missing, on='Sample', how='inner')
        print(f"Merged samples: {len(df_merged)}")

        # 4. Integrate Group Info
        if group_file:
             df_merged = anno_group(df_merged, group_file) # Use anno_group logic
        
        # 5. Save
        save_file = f"{output_prefix}.data.tsv"
        print(f"Saving merged data to {save_file}...")
        save_df_to_tsv(df_merged, save_file)

    # Stats
    mean_ibs = df_merged['IBS_Ref'].mean()
    median_ibs = df_merged['IBS_Ref'].median()
    print(f"IBS_Ref - Mean: {mean_ibs:.4f}, Median: {median_ibs:.4f}")
    if 'Group' in df_merged.columns:
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
        x_col='F_MISS', 
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
    output_prefix,
    tsv_file=None
):
    print("Processing Ref IBS vs Mapping Rate Analysis...")
    
    if tsv_file and os.path.exists(tsv_file):
        print(f"Loading precomputed data from: {tsv_file}...")
        df_merged = load_df_from_tsv(tsv_file)
    else:
        # 1. Read Mapping Rate Data
        df_map = load_df_from_idxstats(mapping_file)
        
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
        if group_file:
            df_merged = anno_group(df_merged, group_file)
        
        # 5. Save
        save_file = f"{output_prefix}.ref_ibs.mapping.tsv"
        print(f"Saving merged data to {save_file}...")
        save_df_to_tsv(df_merged, save_file)
        
    # Print basic stats
    if 'Mapping_Rate_Pct' in df_merged.columns and 'IBS_Ref' in df_merged.columns:
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
        filename=output_filename,
        title='Correlation: Reference IBS vs Mapping Rate'
    )
        
        
        
        
def main():
    parser = argparse.ArgumentParser(description="Reference IBS Analysis")
    subparsers = parser.add_subparsers(dest='command', help='Task: basic_ibs (Distribution) or ibs_vs_map (Correlation)')

    # 1. Basic IBS Analysis
    p_ibs = subparsers.add_parser('basic', help='Basic IBS Distribution Analysis')
    p_ibs.add_argument("-c", "--scount", required=True, help="Input .scount file")
    p_ibs.add_argument("-m", "--smiss", required=True, help="Input .smiss file")
    p_ibs.add_argument("-g", "--group", help="Sample group file")
    p_ibs.add_argument("-o", "--output_prefix", default="ref_ibs_analysis", help="Output prefix")
    p_ibs.add_argument("-d", "--data", dest="tsv_file", help="Precomputed data TSV file")

    # 2. IBS vs Mapping
    p_map = subparsers.add_parser('mapping', help='IBS vs Mapping Rate Correlation')
    p_map.add_argument("-idx", "--idxstats", required=True, help="Input mapping rate file (idxstats)")
    p_map.add_argument("-c", "--scount", required=True, help="Input .scount file")
    p_map.add_argument("-g", "--group", help="Sample group file")
    p_map.add_argument("-o", "--output_prefix", default="ref_ibs_vs_map", help="Output prefix")
    p_map.add_argument("-d", "--data", dest="tsv_file", help="Precomputed data TSV file")

    args = parser.parse_args()

    if args.command == 'basic':
        ana_ref_ibs(args.scount, args.smiss, args.group, args.output_prefix, args.tsv_file)
    elif args.command == 'mapping':
        ref_ibs_vs_mapping(args.idxstats, args.scount, args.group, args.output_prefix, args.tsv_file)
    else:
        parser.print_help()
        

if __name__ == "__main__":
    main()
    
