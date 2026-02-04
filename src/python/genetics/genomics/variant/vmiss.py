import os
from infra.utils.io import load_df_from_space_sep, save_df_to_tsv
from infra.utils.graph import plot_distribution_with_stats

def plot_variant_missingness_distribution_plink(
    input_file="/data1/dazheng_tusr1/vmap4.VCF.v1/check_missing.vmiss",
    output_prefix="variant_missingness"
):
    """
    Plots the distribution of variant missing rates from a .vmiss file.
    Input: Suffix .vmiss (from plink --missing)
    """
    # 1. Load Data
    df = load_vmiss(input_file)
    if df is None: return

    # 2. Save Intermediate Data
    save_file = f"{output_prefix}.tsv"
    save_df_to_tsv(df, save_file)

    # 3. Plot
    # Column 'F_MISS' is standard from load_vmiss
    plot_distribution_with_stats(
        data=df,
        col='F_MISS',
        title='Distribution of Variant Missing Rate',
        filename=f"{output_prefix}.png",
        x_label='Missing Rate (F_MISS)',
        y_label='Count of Variants',
        bins=50,
        color='skyblue',
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
    df = load_df_from_space_sep(input_file)
    if df is None: return
    
    if 'F_MISS' not in df.columns:
        print(f"[Error] 'F_MISS' column not found in {input_file}")
        return

    # 2. Save Intermediate Data
    save_file = f"{output_prefix}.tsv"
    save_df_to_tsv(df, save_file)

    # 3. Plot
    plot_distribution_with_stats(
        data=df,
        col='F_MISS',
        title='Site Missingness',
        filename=f"{output_prefix}.png",
        x_label='Fraction of Missing Data',
        y_label='Count',
        bins=100, # Approximate equivalent to binwidth=0.01 for 0-1 range
        color='forestgreen'
    )


def load_vmiss(filepath):
    """
    Loads variant missingness file (PLINK .vmiss format).
    Extracts Position from ID column if 'Position' column doesn't exist but ID does (e.g. 2-952).
    """
    print(f"[Info] Loading VMISS: {filepath}")
    df = load_df_from_space_sep(filepath)
    if df is None: return None

    # Clean header (remove #)
    df.columns = [c.replace('#', '') if isinstance(c, str) else c for c in df.columns]

    if 'Position' not in df.columns:
        if 'ID' in df.columns:
            # Try to extract position from 'ID' column (e.g., 'Chr-Pos')
            # Assuming ID format: Chr-Pos-...
            try:
                df['Position'] = df['ID'].astype(str).str.split('-').str[1].astype(int)
            except Exception:
                print("[Warning] Could not extract Position from ID column.")
        else:
             pass # Columns might just be CHROM POS etc.

    return df


