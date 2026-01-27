import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
import fastcluster

# Settings
FILE_PATHS = {
    "Normal": {
        "matrix": "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.mibs",
        "id": "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.mibs.id",
        "output_prefix": "normal_ibs"
    },
    "RawSample": {
        "matrix": "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.raw_sample.mibs",
        "id": "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.raw_sample.mibs.id",
        "output_prefix": "raw_sample_ibs"
    },
    "RawMissingFilteredSample": {
        "matrix": "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.raw_missing.mibs",
        "id": "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.raw_missing.mibs.id",
        "output_prefix": "raw_missing_filtered_ibs"
    }
}

GROUP_FILE = "/data1/dazheng_tusr1/vmap4.VCF.v1/sample_group.txt"

def process_and_plot_heatmap(name, config):
    print(f"\nProcessing {name} IBS Matrix for Heatmap...")
    mat_file = config["matrix"]
    id_file = config["id"]
    out_prefix = config["output_prefix"]
    
    # Group Color Palette
    GROUP_COLORS = {
        'A': 'red',
        'AB': 'orange',
        'ABD': 'green',
        'Watkins': 'blue',
        'HZNU': 'purple',
        'WAP': 'brown',
        'Nature': 'cyan',
        'Others': 'grey',
        'S': 'pink',
        'D': 'yellow'
    }
    
    # 1. Read IDs
    print(f"Reading IDs from {id_file}...")
    try:
        df_id = pd.read_csv(id_file, sep=r'\s+', header=None)
        if df_id.shape[1] >= 2:
            sample_ids = df_id[1].values
        else:
            sample_ids = df_id[0].values
    except Exception as e:
        print(f"Failed to read ID file: {e}")
        return

    n_samples = len(sample_ids)
    print(f"Number of samples: {n_samples}")

    # 2. Read Groups & Assign Colors
    print(f"Reading Groups from {GROUP_FILE}...")
    sample_to_group = {}
    if os.path.exists(GROUP_FILE):
        try:
            df_group_file = pd.read_csv(GROUP_FILE, sep=r'\s+', header=None, names=['Sample', 'Group'])
            df_group_file = df_group_file.drop_duplicates(subset=['Sample'])
            sample_to_group = dict(zip(df_group_file['Sample'], df_group_file['Group']))
        except Exception as e:
            print(f"Warning: Failed to read group file: {e}")
    else:
        print("Warning: Group file not found.")

    # Create Color Mapping Vector
    sample_colors = []
    for sid in sample_ids:
        grp = sample_to_group.get(sid, 'Others')
        color = GROUP_COLORS.get(grp, 'grey')
        sample_colors.append(color)
    
    row_colors = pd.Series(sample_colors, index=sample_ids, name='Group')

    # 3. Read Matrix
    print(f"Reading Matrix from {mat_file} (this may take time)...")
    try:
        df_mat = pd.read_csv(mat_file, sep=r'\s+', header=None)
        # Set index and columns for clustermap
        df_mat.index = sample_ids
        df_mat.columns = sample_ids
    except Exception as e:
        print(f"Failed to read Matrix file: {e}")
        return

    if df_mat.shape != (n_samples, n_samples):
        print(f"Error: Matrix shape {df_mat.shape} does not match ID count {n_samples}")
        return
    
    # Check for NaNs
    if df_mat.isnull().values.any():
        nan_count = df_mat.isnull().sum().sum()
        print(f"Warning: Found {nan_count} NaN values in matrix. Filling with 0 to proceed with clustering.")
        df_mat = df_mat.fillna(0)

    # 4. Plot Heatmap

    print(f"Plotting Clustered Heatmap for {name}...")
    # Since N=7675 is large, we adjust figsize and disable xticklabels/yticklabels
    # 'row_colors' adds the side bar. col_colors for top bar.
    
    # Limit max size or performance might suffer?
    # seaborn clustermap can be slow for N=7000+. 
    # Use method='average' (UPGMA) or 'ward'. 'average' is standard for UPGMA.
    # Metric: euclidean or correlation? 
    # Usually IBS -> Distance = 1-IBS.
    # standard_scale=None since IBS is normalized.
    
    try:
        # Increase recursion limit just in case for linkage
        sys.setrecursionlimit(20000)
        
        g = sns.clustermap(
            df_mat,
            method='average',      # UPGMA
            metric='euclidean',    # metric for linkage (on IBS rows? or convert to dist?)
                                   # Ideally we cluster on Distance = 1-IBS.
                                   # But clustermap metric='euclidean' computes dist between IBS vectors.
                                   # Standard approach: provide pre-calculated linkage or let sns do it.
                                   # "correlation" is also common.
            row_colors=row_colors, # Add group colors
            col_colors=row_colors, # Add group colors
            cmap="viridis",        # Color map
            figsize=(15, 15),
            xticklabels=False,     # Hide labels
            yticklabels=False,     # Hide labels
            cbar_kws={'label': 'IBS Score'},
            dendrogram_ratio=(0.15, 0.15), # Ratio of dendrogram size
            colors_ratio=0.03      # Ratio of color bar size
        )
        
        # Add Legend for Groups
        # Create custom legend entries
        from matplotlib.patches import Patch
        
        # Get present groups
        present_groups = set(sample_to_group.values())
        if 'Others' in sample_colors: # Check if 'Others' was actually used
             if any(c == GROUP_COLORS['Others'] for c in sample_colors):
                 present_groups.add('Others')
                 
        legend_elements = [Patch(facecolor=color, edgecolor='none', label=grp) 
                           for grp, color in GROUP_COLORS.items() 
                           if grp in present_groups] # Only show groups present in data
                           
        # Add legend to the figure
        # Place it to the right of the figure content
        # Using figure coordinates, (1.0, 0.5) is the right edge center
        legend = g.fig.legend(handles=legend_elements, title='Sample Groups', 
                   loc='center left', bbox_to_anchor=(1.0, 0.5), 
                   borderaxespad=0.)

        out_file = f"{out_prefix}_heatmap.png"
        
        # Ensure the legend is saved by using bbox_extra_artists (though tight usually handles fig.legend)
        g.savefig(out_file, dpi=300, bbox_inches='tight')
        print(f"Saved {out_file}")
        plt.close() # Close the current figure explicitly if needed, but g.savefig manages its own.
        
    except Exception as e:
        print(f"Error during plotting: {e}")
        import traceback
        traceback.print_exc()

# Run for both
# process_and_plot_heatmap("Normal IBS", FILE_PATHS["Normal"])
# process_and_plot_heatmap("Raw Sample IBS", FILE_PATHS["RawSample"])
process_and_plot_heatmap("Raw Missing Filtered Sample IBS", FILE_PATHS["RawMissingFilteredSample"])

print("\nAnalysis Complete.")
