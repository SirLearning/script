import os
import pandas as pd
from infra.utils.io import load_df_from_tsv, save_df_to_tsv
from infra.utils.graph import plot_distribution_with_stats, plot_correlation_with_regression
from genetics.germplasm.sample import load_df_from_tbm, load_df_from_plink2

def coverage_dist(
    input_file="/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap/all.A.taxaBamMap.txt",
    tsv_file=None,
    output_prefix="sample.coverage"
):
    """
    Plots histogram of Individual Mean Coverage from .idepth file.
    Input: .idepth file
    Output: Histogram plot
    """
    if tsv_file and os.path.exists(tsv_file):
        df_coverage = load_df_from_tsv(tsv_file)
        if df_coverage is None: return
    else:
        df = load_df_from_tbm(input_file)
        if df is None: return
        # extract Sample and Coverage columns
        df_coverage = df[['Sample', 'Coverage']].copy()
        # Save processed data
        save_df_to_tsv(df_coverage, f"{output_prefix}.tsv")

    # Use shared plotting
    plot_distribution_with_stats(
        data=df_coverage, col='Coverage',
        title="Individual Mean Coverage",
        filename=f"{output_prefix}.png",
        x_label="Mean Coverage", y_label="Count",
        color="firebrick", bins=50 # Approximate binwidth=1 for typical depth range 0-50
    )
    

def reg_missing_coverage(
    depth_file="/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap/all.A.taxaBamMap.txt",
    miss_file="/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.missing.smiss",
    tsv_file=None,
    output_prefix="sample.coverage_vs_missing"
):
    if tsv_file and os.path.exists(tsv_file):
        merged = load_df_from_tsv(tsv_file)
        if merged is None: return
    else:
        # 1. Read Data
        df_coverage = load_df_from_tbm(depth_file)
        df_miss = load_df_from_plink2(miss_file)
        if df_coverage is None or df_miss is None: 
            print("[Error] Failed to load input files.")
            return
        # 2. Merge Data
        if 'Sample' in df_coverage.columns and 'Sample' in df_miss.columns:
            merged = pd.merge(df_coverage, df_miss, left_on='Sample', right_on='Sample')
        else:
            print("[Error] Columns 'Sample' (depth) or 'Sample' (miss) not found.")
            print(f"Coverage Cols: {df_coverage.columns}")
            print(f"Miss Cols: {df_miss.columns}")
            return
        # 3. Print Stats
        print(f"Coverage samples: {len(df_coverage)}")
        print(f"Missing samples: {len(df_miss)}")
        print(f"Merged samples: {len(merged)}")
        if 'Coverage' in df_coverage.columns:
            mean_coverage = df_coverage['Coverage'].mean()
            print(f"Mean Coverage (All): {mean_coverage:.4f}")
        # Save merged data
        save_df_to_tsv(merged, f"{output_prefix}_data.tsv")

    # 4. Regression Plot
    if 'Coverage' in merged.columns and 'F_MISS' in merged.columns:
        plot_correlation_with_regression(
            data=merged,
            x_col='Coverage', y_col='F_MISS',
            title='Correlation: Sequencing Coverage vs. Missing Rate',
            filename=f'{output_prefix}.png',
            x_label='Coverage (X)', y_label='Missing Rate (F_MISS)',
            color='gray', line_color='red'
        )
        
        # 5. Filter & Print Outliers
        high_depth = merged[merged['Coverage'] > 20]
        print(f"\n========= Samples with Coverage > 20 (Count: {len(high_depth)}) =========")
        print(high_depth[['Sample', 'F_MISS']].to_string(index=False))

        high_miss = merged[merged['F_MISS'] > 0.8]
        print(f"\n========= Samples with Missing Rate > 0.8 (Count: {len(high_miss)}) =========")
        print(high_miss[['Sample', 'F_MISS']].to_string(index=False))
    else:
        print("[Error] Creating plot: Required columns 'Coverage' or 'F_MISS' missing.")



