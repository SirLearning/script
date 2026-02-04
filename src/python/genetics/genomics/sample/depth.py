import pandas as pd
from infra.utils.io import load_df_generic, save_df_to_tsv
from infra.utils.graph import plot_distribution_with_stats, plot_correlation_with_regression

def plot_individual_depth(
    input_file,
    output_prefix="ind_depth"
):
    """
    Plots histogram of Individual Mean Depth from .idepth file.
    Input: .idepth file
    Output: Histogram plot
    """
    print(f"[Info] Processing Individual Depth: {input_file}")
    df = load_df_generic(input_file)
    if df is None: return

    if 'MEAN_DEPTH' not in df.columns:
        print(f"[Error] 'MEAN_DEPTH' column not found in {input_file}")
        return
    
    # Save processed data
    save_df_to_tsv(df, f"{output_prefix}.tsv")

    # Use shared plotting
    plot_distribution_with_stats(
        data=df, col='MEAN_DEPTH',
        title="Individual Mean Depth",
        filename=f"{output_prefix}.png",
        x_label="Mean Depth", y_label="Count",
        color="firebrick", bins=50 # Approximate binwidth=1 for typical depth range 0-50
    )

def missing_vs_depth(
    depth_file="/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap/all.A.taxaBamMap.txt",
    miss_file="/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.missing.smiss",
    output_prefix="depth_vs_missing"
):
    # 1. Read Data
    df_depth = load_df_generic(depth_file)
    df_miss = load_df_generic(miss_file)
    
    if df_depth is None or df_miss is None: 
        print("[Error] Failed to load input files.")
        return

    # 2. Merge Data
    # Automatically handles left_on/right_on if columns exist
    if 'Taxa' in df_depth.columns and '#IID' in df_miss.columns:
        merged = pd.merge(df_depth, df_miss, left_on='Taxa', right_on='#IID')
    elif 'Taxa' in df_depth.columns and 'INDV' in df_miss.columns:
        merged = pd.merge(df_depth, df_miss, left_on='Taxa', right_on='INDV')
    else:
        # Blind merge if columns not standard, or raise error. 
        # But for this specific script, we expect specific columns based on original code.
        # Original code: left_on='Taxa', right_on='#IID'
        if 'Taxa' in df_depth.columns and '#IID' in df_miss.columns:
             merged = pd.merge(df_depth, df_miss, left_on='Taxa', right_on='#IID')
        else:
             print("[Error] Columns 'Taxa' (depth) or '#IID' (miss) not found.")
             print(f"Depth Cols: {df_depth.columns}")
             print(f"Miss Cols: {df_miss.columns}")
             return

    # 3. Print Stats
    print(f"Depth samples: {len(df_depth)}")
    print(f"Missing samples: {len(df_miss)}")
    print(f"Merged samples: {len(merged)}")

    if 'Coverage-Of-All-Bams' in df_depth.columns:
        mean_depth = df_depth['Coverage-Of-All-Bams'].mean()
        print(f"Mean Coverage (All): {mean_depth:.4f}")
    
    # Save merged data
    save_df_to_tsv(merged, f"{output_prefix}_data.tsv")

    # 4. Regression Plot
    if 'Coverage-Of-All-Bams' in merged.columns and 'F_MISS' in merged.columns:
        plot_correlation_with_regression(
            data=merged,
            x_col='Coverage-Of-All-Bams', y_col='F_MISS',
            title='Correlation: Sequencing Depth vs. Missing Rate',
            filename=f'{output_prefix}.png',
            x_label='Coverage (X)', y_label='Missing Rate (F_MISS)',
            color='gray', line_color='red'
        )
        
        # 5. Filter & Print Outliers
        high_depth = merged[merged['Coverage-Of-All-Bams'] > 20]
        print(f"\n========= Samples with Coverage > 20 (Count: {len(high_depth)}) =========")
        print(high_depth[['Taxa', 'F_MISS']].to_string(index=False))

        high_miss = merged[merged['F_MISS'] > 0.8]
        print(f"\n========= Samples with Missing Rate > 0.8 (Count: {len(high_miss)}) =========")
        print(high_miss[['Taxa', 'F_MISS']].to_string(index=False))
    else:
        print("[Error] Creating plot: Required columns 'Coverage-Of-All-Bams' or 'F_MISS' missing.")



