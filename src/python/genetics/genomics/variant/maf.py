import numpy as np
import pandas as pd
from infra.utils.io import load_df_from_space_sep, save_df_to_tsv, save_thresholds
from infra.utils.graph import plot_distribution_with_stats

def maf_dist(
    input_file="/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.maf.rm_germ_dup.afreq", 
    output_prefix="maf_dist_rm_germ_dup"
):
    # 1. Read Data
    df = load_df_from_space_sep(input_file)
    if df is None: return

    # 2. Calculate MAF
    if 'ALT_FREQS' in df.columns:
        df['MAF'] = df['ALT_FREQS'].apply(lambda x: x if x <= 0.5 else 1 - x)
    else:
        col = 'ALT_FREQS' if 'ALT_FREQS' in df.columns else 'AF'
        if col in df.columns:
             df['MAF'] = df[col].apply(lambda x: x if x <= 0.5 else 1 - x)
        else:
             print("Error: Could not find ALT_FREQS or AF column")
             return

    # 3. Save Calculated Data
    save_df_to_tsv(df, f"{output_prefix}_data.tsv")

    # 4. Calculate Stats & Save
    df_sub = df[df['MAF'] <= 0.01].copy()
    mean_maf_all = df['MAF'].mean()
    median_maf_all = df['MAF'].median()
    
    stats_dict = {
        'Total_Variants': len(df),
        'Variants_MAF_le_0.01': len(df_sub),
        'Percent_MAF_le_0.01': (len(df_sub) / len(df)) * 100,
        'Global_Mean_MAF': mean_maf_all,
        'Global_Median_MAF': median_maf_all
    }
    save_thresholds(stats_dict, f"{output_prefix}_stats.tsv")
    print(stats_dict)

    # 5. Plotting
    
    # Common Thresholds
    thresh_lines = [
        {'value': 0.001, 'label': 'Threshold 0.001', 'color': '#55a868', 'linestyle': '--'},
        {'value': 0.005, 'label': 'Threshold 0.005', 'color': '#c44e52', 'linestyle': '--'}
    ]

    # Plot 1: 0-0.01 Linear
    plot_distribution_with_stats(
        data=df_sub, col='MAF',
        title='Distribution of Minor Allele Frequency (0-0.01) - Linear Scale',
        filename=f"{output_prefix}_0.01_linear.png",
        x_label='MAF', y_label='Count of Variants',
        bins=100, color='#4c72b0',
        mean_val=mean_maf_all, median_val=median_maf_all,
        thresholds=thresh_lines,
        xlim=(0, 0.01)
    )

    # Plot 2: 0-0.01 Log
    plot_distribution_with_stats(
        data=df_sub, col='MAF',
        title='Distribution of Minor Allele Frequency (0-0.01) - Log Scale',
        filename=f"{output_prefix}_0.01_log.png",
        x_label='MAF', y_label='Count of Variants',
        bins=100, color='#4c72b0',
        log_scale=True,
        mean_val=mean_maf_all, median_val=median_maf_all,
        thresholds=thresh_lines,
        xlim=(0, 0.01)
    )

    # Plot 3: 0-0.5 Log (Full Range)
    plot_distribution_with_stats(
        data=df, col='MAF',
        title='Distribution of Minor Allele Frequency (0-0.5) - Log Scale',
        filename=f"{output_prefix}_0.5_log.png",
        x_label='MAF', y_label='Count of Variants (Log Scale)',
        bins=500, color='forestgreen',
        log_scale=True,
        mean_val=mean_maf_all, median_val=median_maf_all,
        xlim=(0, 0.5)
    )

    # Plot 4: 0-0.05 Linear
    df_sub_05 = df[df['MAF'] <= 0.05]
    plot_distribution_with_stats(
        data=df_sub_05, col='MAF',
        title='Distribution of Minor Allele Frequency (0-0.05)',
        filename=f"{output_prefix}_0.05_linear.png",
        x_label='MAF', y_label='Count of Variants',
        bins=200, color='forestgreen',
        mean_val=mean_maf_all, median_val=median_maf_all,
        thresholds=[{'value': 0.01, 'label': 'Threshold 0.01', 'color': 'red', 'linestyle': '--'}],
        xlim=(0, 0.05)
    )

def plot_allele_frequency(
    input_file,
    output_prefix="maf_distribution"
):
    """
    Plots histogram of Minor Allele Frequency (MAF) from .frq file.
    Input: .frq file (VCFtools output)
    """
    print(f"[Info] Processing Allele Frequency: {input_file}")
    try:
        # Load Raw using pandas directly for raw reading capability if format is messy, 
        # but standardized space-sep usually works.
        df = load_df_from_space_sep(input_file) # Replaces generic read_csv
        
        # If generic load fail or specialized parsing needed for {ALLELE:FREQ}
        # The logic below relies on columns at index 4 and 5 having 'Allele:Freq' string.
        # This implies standard whitespace separation worked for first few columns.
        if df is None: return

        if df.shape[1] < 6:
            print("[Warning] .frq file has fewer than 6 columns.")
            return

        def parse_maf(val):
            if isinstance(val, str) and ':' in val:
                try:
                    return float(val.split(':')[1])
                except ValueError:
                    return np.nan
            return np.nan

        freq1 = df.iloc[:, 4].apply(parse_maf)
        freq2 = df.iloc[:, 5].apply(parse_maf)
        freqs = pd.concat([freq1, freq2], axis=1)
        maf = freqs.min(axis=1).dropna()
        
        if len(maf) == 0:
            print("[Error] No valid MAF data found.")
            return

        # Create DataFrame for plotting/saving
        df_maf = pd.DataFrame({'MAF': maf})
        save_df_to_tsv(df_maf, f"{output_prefix}.tsv")

        plot_distribution_with_stats(
            data=df_maf, col='MAF',
            title='Minor Allele Frequency (MAF) Distribution',
            filename=f"{output_prefix}.png",
            x_label='MAF', y_label='Count',
            bins=100, color='orange'
        )

    except Exception as e:
        print(f"[Error] Failed to plot allele frequency: {e}")


if __name__ == "__main__":
    maf_dist()
