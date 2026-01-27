import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

def show_ibs(
    tag="ibs_analysis",
    matrix_file="/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.raw_missing.mibs",
    id_file="/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.raw_missing.mibs.id",
    output_prefix="ibs_analysis",
    group_file="/data1/dazheng_tusr1/vmap4.VCF.v1/sample_group.txt"
):
    file_config = {
        "matrix": matrix_file,
        "id": id_file,
        "output_prefix": output_prefix
    }

    # Run for both
    process_and_plot(tag, file_config, group_file)

    print("\nAnalysis Complete.")

def process_and_plot(name, config, group_file):
    print(f"\nProcessing {name} IBS Matrix...")
    mat_file = config["matrix"]
    id_file = config["id"]
    out_prefix = config["output_prefix"]
    
    # 1. Read IDs
    print(f"Reading IDs from {id_file}...")
    try:
        # PLINK .id files usually have FID IID or just IID. 
        # mibs.id usually acts as FID IID matching the matrix rows
        df_id = pd.read_csv(id_file, sep=r'\s+', header=None)
        # Assuming 2nd column is IID if 2 cols, or 1st if 1.
        if df_id.shape[1] >= 2:
            sample_ids = df_id[1].values
        else:
            sample_ids = df_id[0].values
    except Exception as e:
        print(f"Failed to read ID file: {e}")
        return

    n_samples = len(sample_ids)
    print(f"Number of samples: {n_samples}")

    # 2. Read Groups
    print(f"Reading Groups from {group_file}...")
    sample_to_group = {}
    if os.path.exists(group_file):
        try:
            df_group_file = pd.read_csv(group_file, sep=r'\s+', header=None, names=['Sample', 'Group'])
            # Creating a dictionary map
            # Deduplicate just in case
            df_group_file = df_group_file.drop_duplicates(subset=['Sample'])
            sample_to_group = dict(zip(df_group_file['Sample'], df_group_file['Group']))
        except Exception as e:
            print(f"Warning: Failed to read group file: {e}")
    else:
        print("Warning: Group file not found.")

    # Convert sample_ids to groups (vectorized later, but listing them now)
    # Default group 'Others' (actually map unlisted to 'Others' or specific logic?)
    # User said: "A-others". So if not A or AB, treat as Others.
    # But first we need the raw group name (e.g. ABD, Watkins).
    # If a sample is NOT in sample_group.txt, we probably treat it as Others too?
    # Let's map strict names first.
    
    sample_groups_list = [sample_to_group.get(sid, 'Others') for sid in sample_ids]
    # Verify what groups we have
    # We only care if they are 'A', 'AB' or something else.
    # So let's normalize the list for faster processing?
    # No, keep original for debugging maybe. But for classification we need A/AB check.
    
    # 3. Read Matrix
    print(f"Reading Matrix from {mat_file} (this may take time)...")
    try:
        # matrix is N x N, space separated, no header?
        # Usually .mibs is just numbers.
        # Use numpy loadtxt or pandas read_csv. Pandas is usually faster for big csv.
        df_mat = pd.read_csv(mat_file, sep=r'\s+', header=None)
        matrix = df_mat.values
    except Exception as e:
        print(f"Failed to read Matrix file: {e}")
        return

    if matrix.shape != (n_samples, n_samples):
        print(f"Error: Matrix shape {matrix.shape} does not match ID count {n_samples}")
        return

    # 4. Extract Upper Triangle (excluding diagonal)
    print("Extracting upper triangle pairs...")
    iu1 = np.triu_indices(n_samples, k=1)
    ibs_values = matrix[iu1]
    
    # Get indices to map back to groups
    idx1 = iu1[0]
    idx2 = iu1[1]
    
    # This might be large. 7675^2 / 2 ~ 29 million.
    # Vectorized group lookup
    groups_arr = np.array(sample_groups_list)
    g1_list = groups_arr[idx1]
    g2_list = groups_arr[idx2]
    
    print("Classifying pairs...")
    # Vectorized classification is hard with custom logic strings. 
    # Use pandas for easier mapping
    df_pairs = pd.DataFrame({
        'IBS': ibs_values,
        'G1': g1_list,
        'G2': g2_list
    })
    
    # Vectorized categorization
    # Create booleans
    is_a1 = (df_pairs['G1'] == 'A')
    is_ab1 = (df_pairs['G1'] == 'AB')
    is_a2 = (df_pairs['G2'] == 'A')
    is_ab2 = (df_pairs['G2'] == 'AB')
    
    # Initialize as Others-Others
    df_pairs['Category'] = 'Others-Others'
    
    # A-A
    mask_AA = is_a1 & is_a2
    df_pairs.loc[mask_AA, 'Category'] = 'A-A'
    
    # AB-AB
    mask_ABAB = is_ab1 & is_ab2
    df_pairs.loc[mask_ABAB, 'Category'] = 'AB-AB'
    
    # A-AB (mixed order A-AB or AB-A)
    mask_A_AB = (is_a1 & is_ab2) | (is_ab1 & is_a2)
    df_pairs.loc[mask_A_AB, 'Category'] = 'A-AB'
    
    # A-Others (A-!A!AB or !A!AB-A)
    # Others means NOT A and NOT AB
    mask_others1 = ~(is_a1 | is_ab1)
    mask_others2 = ~(is_a2 | is_ab2)
    
    mask_A_Others = (is_a1 & mask_others2) | (mask_others1 & is_a2)
    df_pairs.loc[mask_A_Others, 'Category'] = 'A-Others'
    
    # AB-Others
    mask_AB_Others = (is_ab1 & mask_others2) | (mask_others1 & is_ab2)
    df_pairs.loc[mask_AB_Others, 'Category'] = 'AB-Others'
    
    # Others-Others (already default, but can check)
    # mask_OO = mask_others1 & mask_others2
    # df_pairs.loc[mask_OO, 'Category'] = 'Others-Others'
    
    # 5. Stats
    print("Calculating Statistics...")
    mean_val = df_pairs['IBS'].mean()
    median_val = df_pairs['IBS'].median()
    min_val = df_pairs['IBS'].min()
    max_val = df_pairs['IBS'].max()
    
    print(f"IBS Stats ({name}):")
    print(f"  Mean  : {mean_val:.6f}")
    print(f"  Median: {median_val:.6f}")
    print(f"  Min   : {min_val:.6f}")
    print(f"  Max   : {max_val:.6f}")
    
    # 6. Plotting
    print("Plotting...")
    
    # Define hue order for stacking (Bottom to Top)
    # User requested: A-A, A-AB, A-Others, AB-AB, AB-Others, Others-Others
    # Seaborn/Matplotlib stack order: First in list is usually at the bottom.
    # hue_order = ['A-A', 'A-AB', 'A-Others', 'AB-AB', 'AB-Others', 'Others-Others']
    hue_order = ['Others-Others', 'AB-Others', 'AB-AB', 'A-Others', 'A-AB', 'A-A']
    
    # Helper to plot one distribution
    def plot_dist_variant(log_scale, filename_suffix):
        plt.figure(figsize=(12, 8))
        
        ax = sns.histplot(data=df_pairs, x='IBS', hue='Category', hue_order=hue_order,
                     multiple='stack', bins=100, linewidth=0.1, palette='tab10')
        
        # Add stats lines
        plt.axvline(x=mean_val, color='blue', linestyle='--', linewidth=1.5, label=f'Mean: {mean_val:.4f}')
        plt.axvline(x=median_val, color='orange', linestyle=':', linewidth=2, label=f'Median: {median_val:.4f}')
        plt.plot([], [], ' ', label=f'Min: {min_val:.4f}')
        plt.plot([], [], ' ', label=f'Max: {max_val:.4f}')
        
        log_txt = " (Log Scale)" if log_scale else ""
        plt.title(f'Pairwise IBS Distribution ({name}){log_txt}', fontsize=16)
        plt.xlabel('IBS Score', fontsize=14)
        
        if log_scale:
            plt.yscale('log')
            plt.ylabel('Count of Pairs (Log Scale)', fontsize=14)
        else:
            plt.ylabel('Count of Pairs', fontsize=14)
        
        # Refine Legend
        if ax.get_legend():
            h_sns = ax.get_legend().legend_handles
            l_sns = [t.get_text() for t in ax.get_legend().get_texts()]
            ax.get_legend().remove()
        else:
            h_sns, l_sns = [], []

        from matplotlib.lines import Line2D
        h_stats = [
            Line2D([0], [0], color='blue', linestyle='--', linewidth=1.5),
            Line2D([0], [0], color='orange', linestyle=':', linewidth=2),
            Line2D([0], [0], linestyle='none'),
            Line2D([0], [0], linestyle='none')
        ]
        l_stats = [
            f'Mean: {mean_val:.4f}',
            f'Median: {median_val:.4f}',
            f'Min: {min_val:.4f}',
            f'Max: {max_val:.4f}'
        ]
        
        final_handles = h_sns + h_stats
        final_labels = l_sns + l_stats
        
        plt.legend(handles=final_handles, labels=final_labels, title='Category & Stats', loc='upper left', bbox_to_anchor=(1, 1))
        
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        out_file = f"{out_prefix}_{filename_suffix}.png"
        plt.savefig(out_file, dpi=300)
        print(f"Saved {out_file}")
        plt.close()

    # Plot Linear
    plot_dist_variant(log_scale=False, filename_suffix="distribution_linear")
    
    # Plot Log
    plot_dist_variant(log_scale=True, filename_suffix="distribution_log")

    # New Plotting Function: Vertical Subplots (Category-wise)
    def plot_vertical_categories(filename_suffix):
        # Define order: Logical ordering for inspection
        cats_to_plot = ['A-A', 'A-AB', 'A-Others', 'AB-AB', 'AB-Others', 'Others-Others']
        
        # Filter only existing categories
        existing_cats = set(df_pairs['Category'].unique())
        cats_to_plot = [c for c in cats_to_plot if c in existing_cats]
        
        num_cats = len(cats_to_plot)
        if num_cats == 0: return

        # Create subplots
        fig, axes = plt.subplots(nrows=num_cats, ncols=1, figsize=(12, 3 * num_cats), sharex=True)
        if num_cats == 1: axes = [axes] # Handle single case

        # Color mapping to match the stack plot
        palette = sns.color_palette('tab10', len(hue_order))
        cat_color_map = dict(zip(hue_order, palette))
        
        for i, cat in enumerate(cats_to_plot):
            ax = axes[i]
            subset = df_pairs[df_pairs['Category'] == cat]
            
            if len(subset) == 0: continue

            # Plot Density (KDE) + Histogram
            # stat='density' makes the area sum to 1, removing count magnitude differences
            sns.histplot(data=subset, x='IBS', color=cat_color_map.get(cat, 'blue'), 
                         bins=100, linewidth=0.1, ax=ax, stat='density', alpha=0.4, label='Hist')
            sns.kdeplot(data=subset, x='IBS', color=cat_color_map.get(cat, 'blue'), 
                        ax=ax, linewidth=2, label='KDE')
            
            # Add stats for this category
            l_mean = subset['IBS'].mean()
            l_median = subset['IBS'].median()
            
            ax.axvline(x=l_mean, color='blue', linestyle='--', linewidth=1.5)
            ax.axvline(x=l_median, color='orange', linestyle=':', linewidth=2)
            
            # Labeling
            ax.set_ylabel('Density')
            ax.set_title(f"Category: {cat} (Count: {len(subset)})", loc='left', fontsize=12, fontweight='bold')
            
            # Stats text
            stats_text = f"Mean: {l_mean:.4f}\nMedian: {l_median:.4f}"
            ax.text(0.98, 0.85, stats_text, transform=ax.transAxes, ha='right', va='top', 
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))
            
            ax.grid(axis='y', alpha=0.3)
            
            # Only legend on first plot to avoid clutter? Or no legend needed as title explains it.
        
        plt.xlabel('IBS Score', fontsize=14)
        plt.tight_layout()
        
        out_file = f"{out_prefix}_{filename_suffix}.png"
        plt.savefig(out_file, dpi=300)
        print(f"Saved {out_file}")
        plt.close()

    plot_vertical_categories("distribution_split_categories")


