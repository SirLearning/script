import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
import os

# Set global style
sns.set_theme(style="whitegrid")






def plot_site_depth(
    input_file,
    output_prefix="site_depth"
):
    """
    Plots Site Depth stats from .ldepth.mean file.
    1. Histogram of Mean Depth
    2. Scatter of Mean vs SD (if VAR_DEPTH exists)
    3. Histogram of SD Depth (if VAR_DEPTH exists)
    Input: .ldepth.mean file
    """
    print(f"[Info] Processing Site Depth: {input_file}")
    try:
        df = pd.read_csv(input_file, sep=r'\s+')
        if 'MEAN_DEPTH' not in df.columns:
            print(f"[Error] 'MEAN_DEPTH' column not found in {input_file}")
            return

        # Plot 1: Mean Depth Histogram
        plt.figure(figsize=(10, 8))
        sns.histplot(data=df, x='MEAN_DEPTH', binwidth=1, color="darkcyan", edgecolor="black", kde=False)
        plt.title("Site Mean Depth")
        plt.xlabel("Mean Depth")
        plt.ylabel("Count")
        
        outfile_mean = f"{output_prefix}_mean.png"
        plt.tight_layout()
        plt.savefig(outfile_mean, dpi=300)
        plt.close()
        print(f"[Success] Mean Depth plot saved to {outfile_mean}")
        
        # Check for VAR_DEPTH
        if 'VAR_DEPTH' in df.columns:
            # Calculate SD
            df['SD_DEPTH'] = np.sqrt(df['VAR_DEPTH'])
            
            # Plot 2: Mean vs SD Scatter
            plt.figure(figsize=(10, 8))
            sns.scatterplot(data=df, x='MEAN_DEPTH', y='SD_DEPTH', alpha=0.1, color="darkblue", edgecolor=None)
            plt.title("Site Depth: Mean vs SD")
            plt.xlabel("Mean Depth")
            plt.ylabel("Standard Deviation")
            
            outfile_scatter = f"{output_prefix}_mean_vs_sd.png"
            plt.tight_layout()
            plt.savefig(outfile_scatter, dpi=300)
            plt.close()
            print(f"[Success] Mean vs SD plot saved to {outfile_scatter}")
            
            # Plot 3: SD Histogram
            plt.figure(figsize=(10, 8))
            sns.histplot(data=df, x='SD_DEPTH', binwidth=1, color="cyan", edgecolor="black", kde=False)
            plt.title("Site Depth Standard Deviation")
            plt.xlabel("SD Depth")
            plt.ylabel("Count")
            
            outfile_sd = f"{output_prefix}_sd.png"
            plt.tight_layout()
            plt.savefig(outfile_sd, dpi=300)
            plt.close()
            print(f"[Success] SD Depth plot saved to {outfile_sd}")

    except Exception as e:
        print(f"[Error] Failed to plot site depth: {e}")
