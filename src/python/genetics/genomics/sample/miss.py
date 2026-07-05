from pathlib import Path

from .sample_utils import load_df_from_plink2
from genetics.germplasm.sample import anno_group
from infra.utils.io import load_df_from_tsv, save_df_to_tsv, load_thresholds, save_thresholds
from infra.utils.graph import plot_distribution_with_stats, plot_slope_chart
import pandas as pd
import os
import argparse

THIN_COMMON_SUBGENOMES = ('A', 'B', 'D', 'Others')

def ana_sample_missing(
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
            save_df_to_tsv(df[['Sample', 'F_MISS']], f"{output_prefix}.rate.info.tsv")
            
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
            threshold_file = f"{output_prefix}.th.tsv"
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
            filename=f"{output_prefix}.dist.png",
            mean_val=mean_val,
            median_val=median_val,
            std_val=std_val,
            thresholds=thresholds_to_plot,
            x_label='Missing Rate (F_MISS)',
            # color='forestgreen',
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


def _subgenome_smiss_path(process_dir, subgenome):
    return Path(process_dir) / 'sample' / f'{subgenome}.info.smiss'


def merge_thin_common_sample_missing(
    thin_process_dir,
    compare_process_dir,
    subgenome,
    left_label='test_thin',
    right_label='test_common_thin',
):
    """Inner-join per-sample F_MISS for one subgenome across two process mods."""
    thin_path = _subgenome_smiss_path(thin_process_dir, subgenome)
    compare_path = _subgenome_smiss_path(compare_process_dir, subgenome)
    empty_cols = ['Sample', left_label, right_label]
    if not thin_path.exists():
        print(f'[Warning] Missing smiss: {thin_path}')
        return pd.DataFrame(columns=empty_cols)
    if not compare_path.exists():
        print(f'[Warning] Missing smiss: {compare_path}')
        return pd.DataFrame(columns=empty_cols)

    thin = load_df_from_plink2(str(thin_path))
    compare = load_df_from_plink2(str(compare_path))
    if thin is None or compare is None:
        return pd.DataFrame(columns=empty_cols)

    merged = thin[['Sample', 'F_MISS']].merge(
        compare[['Sample', 'F_MISS']],
        on='Sample',
        how='inner',
        suffixes=('_left', '_right'),
    )
    return merged.rename(columns={
        'F_MISS_left': left_label,
        'F_MISS_right': right_label,
    })


def plot_thin_common_sample_missing_slope(
    thin_process_dir,
    compare_process_dir,
    output_prefix,
    group_file=None,
    subgenomes=THIN_COMMON_SUBGENOMES,
    left_label='test_thin',
    right_label='test_common_thin',
):
    """
    Per-subgenome slope charts: each sample connects left_mod → right_mod F_MISS.

    Writes ``{prefix}.sample.missing_slope.{A|B|D|Others}.png`` and
    ``{prefix}.sample.missing_slope.info.tsv``.
    """
    info_rows = []
    for subgenome in subgenomes:
        merged = merge_thin_common_sample_missing(
            thin_process_dir,
            compare_process_dir,
            subgenome,
            left_label=left_label,
            right_label=right_label,
        )
        if merged.empty:
            print(f'[Warning] No overlapping smiss samples for subgenome {subgenome}; skipping plot.')
            continue

        if group_file:
            merged = anno_group(merged, group_file, save_tsv=False)

        plot_slope_chart(
            merged,
            col1=left_label,
            col2=right_label,
            xlabel1=left_label,
            xlabel2=right_label,
            title=f'Sample missing rate: {left_label} vs {right_label} ({subgenome})',
            filename=f'{output_prefix}.sample.missing_slope.{subgenome}.png',
            group_col='Group' if group_file and 'Group' in merged.columns else None,
            ylabel='Missing Rate (F_MISS)',
            ylim=(0.0, 1.0),
        )

        out = merged.copy()
        out.insert(0, 'subgenome', subgenome)
        info_rows.append(out)

    if info_rows:
        save_df_to_tsv(pd.concat(info_rows, ignore_index=True), f'{output_prefix}.sample.missing_slope.info.tsv')
    return info_rows


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
        ana_sample_missing(
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
