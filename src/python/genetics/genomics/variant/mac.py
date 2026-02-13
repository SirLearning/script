from .variant_utils import load_df_from_plink_gcount
import numpy as np
import pandas as pd
from infra.utils.io import save_df_to_tsv
from infra.utils.graph import plot_dual_regression, plot_heatmap_custom

def ana_mac_stats(
    input_file,
    output_prefix="variant_mac"
):
    print(f"Processing MAC Stats: {input_file}")
    
    # 1. Load Data
    df = load_df_from_plink_gcount(input_file)
    if df is None: return

    required = ['HET_REF_ALT_CTS', 'TWO_ALT_GENO_CTS', 'HOM_REF_CT', 'HAP_REF_CT', 'HAP_ALT_CTS']
    if not all(col in df.columns for col in required):
        print(f"Error: Missing required columns in {input_file}")
        return

    # 2. Calculate Counts
    # Alt Count = 1*Het + 2*HomAlt + 1*HapAlt
    df['Alt_Count'] = df['HET_REF_ALT_CTS'] + (df['TWO_ALT_GENO_CTS'] * 2) + df['HAP_ALT_CTS']
    # Ref Count = 1*Het + 2*HomRef + 1*HapRef
    df['Ref_Count'] = df['HET_REF_ALT_CTS'] + (df['HOM_REF_CT'] * 2) + df['HAP_REF_CT']
    
    df['Total_Alleles'] = df['Alt_Count'] + df['Ref_Count']
    df['MAC'] = df[['Alt_Count', 'Ref_Count']].min(axis=1)
    df['MAF'] = df['MAC'] / df['Total_Alleles']
    
    # Hobs
    df['Total_Samples'] = df['HOM_REF_CT'] + df['HET_REF_ALT_CTS'] + df['TWO_ALT_GENO_CTS'] + df['HAP_REF_CT'] + df['HAP_ALT_CTS']
    df['Hobs'] = df['HET_REF_ALT_CTS'] / df['Total_Samples']

    # Het Fraction (Prop of MAC that is heterozygous)
    df['Het_Fraction'] = np.nan
    mask_mac = df['MAC'] > 0
    df.loc[mask_mac, 'Het_Fraction'] = df.loc[mask_mac, 'HET_REF_ALT_CTS'] / df.loc[mask_mac, 'MAC']

    save_df_to_tsv(df, f"{output_prefix}.info.tsv")

    # 3. MAC=1 Analysis
    mac1_df = df[df['MAC'] == 1]
    with open(f"{output_prefix}.mac1.stats.txt", "w") as f:
        f.write(f"MAC=1 Count: {len(mac1_df)}\n")
        f.write(str(mac1_df['Total_Alleles'].describe()))
    
    # 4. Plots
    
    # 4.1 MAC vs MAF (Regression) - Range 0-380
    df_sub = df[(df['MAC'] >= 0) & (df['MAC'] <= 380)].dropna(subset=['MAC', 'MAF'])
    if len(df_sub) > 10:
        plot_dual_regression(
            df_valid=df_sub, df_plot=df_sub,
            x_col='MAC', y_col='MAF',
            x_label='Minor Allele Count (MAC)', y_label='Minor Allele Frequency (MAF)',
            filename=f"{output_prefix}.reg_mac_maf.380.png"
        )
    
    # 4.2 MAC vs Hobs (Regression) - Range 1-100
    df_sub_het = df[(df['MAC'] >= 1) & (df['MAC'] <= 100) & (df['Hobs'] > 0)].dropna()
    if len(df_sub_het) > 10:
        plot_dual_regression(
            df_valid=df_sub_het, df_plot=df_sub_het,
            x_col='MAC', y_col='Hobs',
            x_label='Minor Allele Count (1-100)', y_label='Observed Heterozygosity',
            filename=f"{output_prefix}.reg_mac_hobs.100.png"
        )

    # 4.3 Heatmap: MAC vs Het Fraction
    mac_num = 200
    df_frac = df[(df['MAC'] >= 1) & (df['MAC'] <= mac_num)].dropna(subset=['Het_Fraction'])
    
    if len(df_frac) > 10:
        # Construct Matrix
        # X: MAC (1..mac_num), Y: Fraction (0..1, bin 0.05)
        bin_width = 0.05
        y_bins = np.arange(0, 1.0 + bin_width, bin_width) # 21 bins
        n_y_bins = len(y_bins) - 1 # 20 bins
        
        # We need a matrix of shape (n_y_bins, mac_num)
        matrix = np.zeros((n_y_bins, mac_num))
        
        for mac_val in range(1, mac_num + 1):
            subset = df_frac[df_frac['MAC'] == mac_val]
            if not subset.empty:
                # Histogram of fractions for this MAC
                counts, _ = np.histogram(subset['Het_Fraction'], bins=y_bins)
                # Normalize
                total = counts.sum()
                if total > 0:
                    matrix[:, mac_val-1] = counts / total
        
        # Flip for plotting (top is 1.0)
        matrix = np.flipud(matrix)
        
        # Labels
        x_labels = [str(i) for i in range(1, mac_num + 1)]
        y_labels = [f"{(1.0 - i*bin_width):.1f}" for i in range(n_y_bins)] # Approx labels

        plot_heatmap_custom(
            data_matrix=matrix, x_labels=x_labels, y_labels=y_labels,
            title=f'Het Fraction Distribution per MAC (Proportion)',
            filename=f"{output_prefix}.heatmap_mac_hetfrac.png",
            xlabel='MAC', ylabel='Het Fraction',
            cbar_label='Proportion'
        )

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        ana_mac_stats(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "mac_out")
    else:
        print("Usage: python mac.py <input.gcount> <output_prefix>")
