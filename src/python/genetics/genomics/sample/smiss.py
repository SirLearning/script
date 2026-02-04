from genetics.germplasm.sample.process import load_df_from_plink2
from infra.utils.io import load_df_from_tsv, save_df_to_tsv, load_thresholds, save_thresholds
from infra.utils.graph import plot_distribution_with_stats
import pandas as pd
import os
import argparse

def plot_smiss_dist(
    input_file, 
    output_prefix,
    tsv_file=None,
    th_tsv_file=None,
    use_fixed_thresholds=True
):
    """
    Analyzes sample missing rate from PLINK .smiss file.
    Plots distribution with Mean, Median, and Thresholds.
    Calculates Mean + 3SD threshold and saves it to a TSV file.
    """
    print(f"[Info] Processing Missing Rate Analysis: {input_file}")

    try:
        # Load Data
        if tsv_file and os.path.exists(tsv_file):
            print(f"Loading precomputed missing rates from {tsv_file}...")
            df = load_df_from_tsv(tsv_file)
        else:
            # 1. Read Data
            df = load_df_from_plink2(input_file)
            if df is None: return
            # Save intermediate
            save_df_to_tsv(df[['Sample', 'F_MISS']], f"{output_prefix}.missing_rates.tsv")
            
        col = 'F_MISS'

        # 2. Statistics & Thresholds
        if th_tsv_file and os.path.exists(th_tsv_file):
             print(f"Loading precomputed thresholds from {th_tsv_file}...")
             thresholds_data = load_thresholds(th_tsv_file)
             mean_val = thresholds_data.get('mean_missing', df[col].mean())
             median_val = df[col].median() # Recalculate median usually fast enough, or could save it
             std_val = thresholds_data.get('std_missing', df[col].std())
             calc_threshold = thresholds_data.get('missing_threshold', mean_val + 3 * std_val)
        else:
            mean_val = df[col].mean()
            median_val = df[col].median()
            std_val = df[col].std()
            calc_threshold = mean_val + 3 * std_val
            if pd.isna(calc_threshold):
                calc_threshold = 0.05
            
            # Save Threshold to File
            threshold_file = f"{output_prefix}.threshold.tsv"
            stats_data = {
                'missing_threshold': calc_threshold,
                'mean_missing': mean_val,
                'std_missing': std_val
            }
            save_thresholds(stats_data, threshold_file)
        
        # 3. Determine Thresholds for Plotting
        thresholds_to_plot = []
        
        # Method A: Fixed Thresholds (from original missing_dist)
        if use_fixed_thresholds:
             thresholds_to_plot.append({'value': 0.1, 'label': 'Threshold 0.1', 'color': '#c44e52', 'linestyle': '--'})
             thresholds_to_plot.append({'value': 0.3, 'label': 'Threshold 0.3', 'color': 'orange', 'linestyle': '--'})
        
        # Method B: Statistical Threshold (Mean + 3SD)
        calc_threshold = mean_val + 3 * std_val
        if pd.isna(calc_threshold):
            calc_threshold = 0.05
            
        thresholds_to_plot.append({'value': calc_threshold, 'label': f'Threshold (+3SD): {calc_threshold:.4f}', 'color': 'purple', 'linestyle': '-'})

        # 4. Plot
        plot_distribution_with_stats(
            data=df,
            col=col,
            title='Distribution of Sample Missing Rate',
            filename=f"{output_prefix}.png",
            mean_val=mean_val,
            median_val=median_val,
            std_val=std_val,
            thresholds=thresholds_to_plot,
            x_label='Missing Rate (F_MISS)',
            color='forestgreen',
            xlim=(0, 1.0) if use_fixed_thresholds else None # Use fixed xlim for standard missing plot
        )
        
        # 5. Save Threshold (Already done above if calculating, redundant if loading but safe)
        # threshold_file = f"{output_prefix}.threshold.tsv"
        # ... logic moved up ...

    except Exception as e:
        print(f"[Error] Failed to analyze missing rate: {e}")


def plot_imiss_dist(
    input_file,
    output_prefix="ind_missingness"
):
    """
    Plots histogram of Individual Missingness from .imiss file.
    """
    print(f"[Info] Processing Individual Missingness: {input_file}")
    try:
        df = load_df_from_plink2(input_file)
        col = 'F_MISS'
        if col not in df.columns:
            print(f"[Error] '{col}' column not found in {input_file}")
            return

        # Statistics (Optional but good for plotting)
        mean_val = df[col].mean()
        median_val = df[col].median()

        plot_distribution_with_stats(
            data=df,
            col=col,
            title="Individual Missingness",
            filename=f"{output_prefix}.png",
            mean_val=mean_val,
            median_val=median_val,
            x_label="Fraction of Missing Data",
            color="steelblue",
            bins=100
        )
        
    except Exception as e:
        print(f"[Error] Failed to plot individual missingness: {e}")
        

def main():
    parser = argparse.ArgumentParser(description="Missing Rate Analysis (.smiss, .imiss)")
    subparsers = parser.add_subparsers(dest='command', help='Analysis Type')

    # 1. Sample Missing
    p_smiss = subparsers.add_parser('smiss', help='Analyze Sample Missing Rate')
    p_smiss.add_argument("-i", "--input", required=True, help="Input .smiss file")
    p_smiss.add_argument("-o", "--output_prefix", default="sample_missing_rate_analysis", help="Output prefix")
    p_smiss.add_argument("-d", "--data", dest="tsv_file", help="Precomputed data TSV file")
    p_smiss.add_argument("-th", "--thresholds", dest="th_tsv_file", help="Precomputed thresholds TSV file")
    p_smiss.add_argument("--fixed_thresholds", action='store_true', default=True, help="Use fixed thresholds (0.1, 0.3)")

    # 2. Individual Missing
    p_imiss = subparsers.add_parser('imiss', help='Analyze Individual Missingness')
    p_imiss.add_argument("-i", "--input", required=True, help="Input .imiss file")
    p_imiss.add_argument("-o", "--output_prefix", default="ind_missingness", help="Output prefix")

    args = parser.parse_args()

    if args.command == 'smiss':
        plot_smiss_dist(
            args.input, 
            args.output_prefix, 
            args.tsv_file, 
            args.th_tsv_file, 
            args.fixed_thresholds
        )
    elif args.command == 'imiss':
        plot_imiss_dist(args.input, args.output_prefix)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
