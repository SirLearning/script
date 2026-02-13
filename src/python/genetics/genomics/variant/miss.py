import os
from .variant_utils import load_df_from_plink_variant
from infra.utils.io import load_df_from_tsv, save_df_to_tsv
from infra.utils.graph import plot_distribution_with_stats

def ana_variant_missing(
    input_file="/data1/dazheng_tusr1/vmap4.VCF.v1/check_missing.vmiss",
    output_prefix="variant_missingness"
):
    """
    Plots the distribution of variant missing rates from a .vmiss file.
    Input: Suffix .vmiss (from plink --missing)
    """
    # 1. Load Data
    df = load_df_from_plink_variant(input_file)
    if df is None: return

    # 2. Save Intermediate Data
    save_file = f"{output_prefix}.info.tsv"
    save_df_to_tsv(df, save_file)

    # 3. Plot
    # Column 'F_MISS' is standard from load_vmiss
    plot_distribution_with_stats(
        data=df,
        col='F_MISS',
        title='Distribution of Variant Missing Rate',
        filename=f"{output_prefix}.dist.png",
        x_label='Missing Rate (F_MISS)',
        y_label='Count of Variants',
        bins=50,
        thresholds=[{'value': 0.1, 'label': 'Threshold = 0.1', 'color': 'red', 'linestyle': '--'}]
    )


def plot_site_missingness_vcftools(
    input_file,
    output_prefix="site_missingness"
):
    """
    Plots histogram of Site Missingness from .lmiss file.
    Input: .lmiss file (from vcftools --missing-site)
    """
    print(f"[Info] Processing Site Missingness: {input_file}")
    
    # 1. Load Data
    df = load_df_from_tsv(input_file)
    if df is None: return
    
    if 'F_MISS' not in df.columns:
        print(f"[Error] 'F_MISS' column not found in {input_file}")
        return

    # 2. Save Intermediate Data
    save_file = f"{output_prefix}.info.tsv"
    save_df_to_tsv(df, save_file)

    # 3. Plot
    plot_distribution_with_stats(
        data=df,
        col='F_MISS',
        title='Site Missingness',
        filename=f"{output_prefix}.dist.png",
        x_label='Fraction of Missing Data',
        y_label='Count',
        bins=100, # Approximate equivalent to binwidth=0.01 for 0-1 range
        color='forestgreen'
    )



