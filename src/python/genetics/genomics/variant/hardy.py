from .variant_utils import load_df_from_plink_hardy
import numpy as np
import pandas as pd
from infra.utils.io import save_df_to_tsv, save_thresholds
from infra.utils.graph import plot_stacked_histogram, plot_dual_regression

def ana_hardy_stats(
    input_file, 
    output_prefix="variant_hardy"
):
    print(f"Processing Hardy Stats: {input_file}")
    
    # 1. Read Data
    df = load_df_from_plink_hardy(input_file)
    if df is None: return

    required_cols = ['HOM_A1_CT', 'HET_A1_CT', 'TWO_AX_CT']
    if not all(col in df.columns for col in required_cols):
        print(f"Error: Missing required columns in {input_file}")
        return

    # 2. Calculate Metrics
    df['Total_Samples'] = df['HOM_A1_CT'] + df['HET_A1_CT'] + df['TWO_AX_CT']
    
    # Filter valid
    # Some older plink files might have 0 samples for some reason, filter them out
    df = df[df['Total_Samples'] > 0].copy()

    # Observed Heterozygosity
    het_col = 'O(HET_A1)'
    if het_col not in df.columns:
        df[het_col] = df['HET_A1_CT'] / df['Total_Samples']

    # Freq of A1 (p) & MAF
    df['p'] = (df['HOM_A1_CT'] * 2 + df['HET_A1_CT']) / (df['Total_Samples'] * 2)
    df['MAF'] = df['p'].apply(lambda x: min(x, 1-x))
    
    # Expected Heterozygosity & Inbreeding Coefficient (F)
    df['Hexp'] = 2 * df['p'] * (1 - df['p'])
    
    # Avoid division by zero for F calculation
    df['F'] = np.nan
    mask_hexp = df['Hexp'] > 0
    df.loc[mask_hexp, 'F'] = 1 - (df.loc[mask_hexp, het_col] / df.loc[mask_hexp, 'Hexp'])

    # 3. Save Calculated Data
    save_df_to_tsv(df, f"{output_prefix}.info.tsv")

    # 4. MAF Grouping & Stats
    conditions = [
        (df['MAF'] <= 0.001),
        (df['MAF'] > 0.001) & (df['MAF'] <= 0.05),
        (df['MAF'] > 0.05)
    ]
    choices = ['MAF <= 0.001', '0.001 < MAF <= 0.05', 'MAF > 0.05']
    df['MAF_Group'] = np.select(conditions, choices, default='Unknown')

    mean_het = df[het_col].mean()
    median_het = df[het_col].median()
    
    # Stats Dict
    stats_dict = {
        'Mean_Het': mean_het,
        'Median_Het': median_het,
        'Mean_F': df['F'].mean(),
        'Median_F': df['F'].median(),
        'Valid_Sites_F': int(df['F'].notna().sum())
    }
    save_thresholds(stats_dict, f"{output_prefix}.th.tsv")
    print(stats_dict)

    # 5. Plotting
    PLOT_SAMPLE_SIZE = 1000000
    if len(df) > PLOT_SAMPLE_SIZE:
        df_plot = df.sample(n=PLOT_SAMPLE_SIZE, random_state=42)
    else:
        df_plot = df

    # Plot Settings
    hue_order = ['MAF <= 0.001', '0.001 < MAF <= 0.05', 'MAF > 0.05']
    palette_map = {'MAF <= 0.001': 'lightgray', 
                   '0.001 < MAF <= 0.05': 'teal', 
                   'MAF > 0.05': 'purple'}

    # 5.1 Het Distribution (Log & Linear)
    x_limit_het = df[het_col].quantile(0.99)
    if x_limit_het < 0.01: x_limit_het = 0.01
    x_limit_het = min(x_limit_het * 1.2, 1.0)
    
    plot_stacked_histogram(
        data=df_plot, x_col=het_col, hue_col='MAF_Group',
        hue_order=hue_order, palette_map=palette_map,
        title='Distribution of Variant Heterozygosity - Stacked by MAF (Linear)',
        filename=f"{output_prefix}.het.dist.linear.png",
        xlabel='Observed Heterozygosity', ylabel='Number of Variants',
        mean_val=mean_het, median_val=median_het,
        xlim=(0, x_limit_het), log_scale=False
    )

    plot_stacked_histogram(
        data=df_plot, x_col=het_col, hue_col='MAF_Group',
        hue_order=hue_order, palette_map=palette_map,
        title='Distribution of Variant Heterozygosity - Stacked by MAF (Log)',
        filename=f"{output_prefix}.het.dist.log.log.png",
        xlabel='Observed Heterozygosity', ylabel='Number of Variants',
        mean_val=mean_het, median_val=median_het,
        xlim=(0, x_limit_het), log_scale=True
    )

    # 5.2 F Distribution
    df_plot_f = df_plot.dropna(subset=['F'])
    f_high_maf = df[(df['F'].notna()) & (df['MAF'] > 0.05)]['F']
    median_f_high = f_high_maf.median() if not f_high_maf.empty else None

    plot_stacked_histogram(
        data=df_plot_f, x_col='F', hue_col='MAF_Group',
        hue_order=hue_order, palette_map=palette_map,
        title='Distribution of Inbreeding Coefficient (F)',
        filename=f"{output_prefix}.f.dist.png",
        xlabel='Inbreeding Coefficient (F)', ylabel='Number of Variants',
        mean_val=df['F'].mean(), median_val=median_f_high,
        xlim=(-1.0, 1.0), log_scale=True 
    )

    # 5.3 Regression MAF vs Het
    plot_dual_regression(
        df_valid=df.dropna(subset=['MAF', het_col]), 
        df_plot=df_plot.dropna(subset=['MAF', het_col]),
        x_col='MAF', y_col=het_col,
        x_label='Minor Allele Frequency (MAF)',
        y_label='Observed Heterozygosity',
        filename=f"{output_prefix}.reg_maf_het.png"
    )

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        ana_hardy_stats(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "hardy_out")
    else:
        print("Usage: python hardy.py <input.hardy> <output_prefix>")
