import pandas as pd
import sys
import os

# Add src/python to path if needed to find infra
current_dir = os.path.dirname(os.path.abspath(__file__))
src_python_dir = os.path.abspath(os.path.join(current_dir, "../../../"))
if src_python_dir not in sys.path:
    sys.path.append(src_python_dir)

from infra.plot_utils import plot_distribution_with_stats

def analyze_smiss_distribution(
    input_file, 
    output_prefix,
    use_fixed_thresholds=True
):
    """
    Analyzes sample missing rate from PLINK .smiss file.
    Plots distribution with Mean, Median, and Thresholds.
    Calculates Mean + 3SD threshold and saves it to a TSV file.
    """
    print(f"[Info] Processing Missing Rate Analysis: {input_file}")

    try:
        # 1. Read Data
        df = pd.read_csv(input_file, sep=r'\s+')
        col = 'F_MISS'
        
        if col not in df.columns:
            print(f"[Error] Column '{col}' not found in {input_file}. Columns: {df.columns}")
            return

        # 2. Statistics
        mean_val = df[col].mean()
        median_val = df[col].median()
        std_val = df[col].std()
        
        # 3. Determine Thresholds
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
        
        # 5. Save Threshold to File
        threshold_file = f"{output_prefix}.threshold.tsv"
        out_df = pd.DataFrame({
            'missing_threshold': [calc_threshold],
            'mean_missing': [mean_val],
            'std_missing': [std_val]
        })
        out_df.to_csv(threshold_file, sep='\t', index=False)
        print(f"[Info] Thresholds saved to {threshold_file}")

    except Exception as e:
        print(f"[Error] Failed to analyze missing rate: {e}")


def analyze_imiss_distribution(
    input_file,
    output_prefix="ind_missingness"
):
    """
    Plots histogram of Individual Missingness from .imiss file.
    """
    print(f"[Info] Processing Individual Missingness: {input_file}")
    try:
        df = pd.read_csv(input_file, sep=r'\s+')
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

if __name__ == "__main__":
    # Example usage
    SMISS_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/check_missing.smiss"
    analyze_smiss_distribution(SMISS_FILE, "sample_missing_rate_distribution", use_fixed_thresholds=True)
