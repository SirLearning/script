import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


def ab_ibs_trend(
    matrix_file="/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.raw_missing.mibs",
    id_file="/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.raw_missing.mibs.id", 
    group_file="/data1/dazheng_tusr1/vmap4.VCF.v1/sample_group.txt",
    output_prefix="ab_ibs_trend_plot"
):
    # Map arguments to legacy variable names
    MATRIX_FILE = matrix_file
    ID_FILE = id_file
    GROUP_FILE = group_file
    OUTPUT_FILE = f"{output_prefix}.png"

    print("Starting AB IBS Trend Analysis...")
    
    # 1. Read IDs
    print(f"Reading IDs from {ID_FILE}...")
    try:
        df_id = pd.read_csv(ID_FILE, sep=r'\s+', header=None)
        if df_id.shape[1] >= 2:
            sample_ids = df_id[1].values
        else:
            sample_ids = df_id[0].values
    except Exception as e:
        print(f"Failed to read ID file: {e}")
        return

    n_samples = len(sample_ids)
    print(f"Total samples: {n_samples}")

    # 2. Read Groups & Define Masks
    print(f"Reading Groups from {GROUP_FILE}...")
    sample_to_group = {}
    if os.path.exists(GROUP_FILE):
        df_group_file = pd.read_csv(GROUP_FILE, sep=r'\s+', header=None, names=['Sample', 'Group'])
        df_group_file = df_group_file.drop_duplicates(subset=['Sample'])
        sample_to_group = dict(zip(df_group_file['Sample'], df_group_file['Group']))
    
    # Define groups based on previous logic (A, AB, Others)
    # Target: AB Population
    # Comparison: Others Population (Not A, Not AB)
    
    # Groups vector aligned with matrix rows/cols
    groups_vector = np.array([sample_to_group.get(sid, 'Others') for sid in sample_ids])
    
    # Create Boolean Masks
    mask_ab = (groups_vector == 'AB')
    mask_a = (groups_vector == 'A')
    # Others is everything that is NOT A and NOT AB
    mask_others = ~(mask_a | mask_ab)
    
    indices_ab = np.where(mask_ab)[0]
    indices_a = np.where(mask_a)[0]
    indices_others = np.where(mask_others)[0]
    
    count_ab = len(indices_ab)
    count_a = len(indices_a)
    count_others = len(indices_others)
    
    print(f"Count of 'AB' samples: {count_ab}")
    print(f"Count of 'A' samples: {count_a}")
    print(f"Count of 'Others' samples: {count_others}")
    
    if count_ab == 0:
        print("Error: No AB samples found.")
        return
    # We continue even if Others or A is 0, but check before plotting

    # 3. Read Matrix
    print(f"Reading Matrix from {MATRIX_FILE} (this is large, please wait)...")
    try:
        # Optimization: We only need rows corresponding to AB samples.
        # But standard read_csv reads all.
        # Since we need to calculate means against columns (Others and AB), we strictly need:
        # Rows: indices_ab
        # Cols: indices_ab AND indices_others
        
        # Loading full matrix is simplest given likely pandas implementation efficiency vs chunking text
        df_mat = pd.read_csv(MATRIX_FILE, sep=r'\s+', header=None)
        matrix = df_mat.values
    except Exception as e:
        print(f"Failed to read Matrix file: {e}")
        return

    if matrix.shape != (n_samples, n_samples):
        print(f"Error: Matrix shape {matrix.shape} does not match ID count {n_samples}")
        return

    # 4. Calculate Means
    print("Calculating Mean IBS Scores...")
    
    # Store results: [ (SampleName, Mean_AB, Mean_Others, Mean_A), ... ]
    plotting_data = []
    
    # Iterate only through AB samples
    for i in indices_ab:
        sample_name = sample_ids[i]
        
        # Row for this sample
        row_values = matrix[i, :]
        
        # 1. Mean with AB (excluding self)
        current_mask_ab = mask_ab.copy()
        current_mask_ab[i] = False # exclude self
        
        vals_ab_clean = row_values[current_mask_ab]
        
        if len(vals_ab_clean) > 0:
            mean_ab = np.mean(vals_ab_clean)
        else:
            mean_ab = np.nan
            
        # 2. Mean with Others
        vals_others = row_values[mask_others]
        
        if len(vals_others) > 0:
            mean_others = np.mean(vals_others)
        else:
            mean_others = np.nan

        # 3. Mean with A
        vals_a = row_values[mask_a]
        
        if len(vals_a) > 0:
            mean_a = np.mean(vals_a)
        else:
            mean_a = np.nan
            
        plotting_data.append({
            'Sample': sample_name,
            'Mean_AB': mean_ab,
            'Mean_Others': mean_others,
            'Mean_A': mean_a
        })
        
    df_res = pd.DataFrame(plotting_data)
    
    # Remove NaNs
    # But only drop if the columns we NEED are NaN. 
    # For now, let's keep all and drop inside plotting function depending on what we plot.
    print(f"Calculated means for {len(df_res)} AB samples.")

    # 5. Plotting Function
    def plot_slope_chart(data, col_start, col_end, label_start, label_end, title, filename):
        # Drop NaNs for these specific columns
        df_plot = data.dropna(subset=[col_start, col_end])
        
        if len(df_plot) == 0:
            print(f"No valid data to plot for {filename}")
            return

        print(f"Generating Slope Chart: {filename}...")
        
        plt.figure(figsize=(8, 10))
        
        # X coordinates
        x_start = 0
        x_end = 1
        
        color_decrease = '#1f77b4' # Muted Blue
        color_increase = '#d62728' # Muted Red
        alpha_line = 0.4
        
        n_inc = 0
        n_dec = 0

        # Plot lines
        for _, row in df_plot.iterrows():
            val_start = row[col_start]
            val_end = row[col_end]
            
            if val_end > val_start:
                c = color_increase # Rising
                z = 2
                n_inc += 1
            else:
                c = color_decrease # Falling
                z = 1
                n_dec += 1
                
            plt.plot([x_start, x_end], [val_start, val_end], 
                     color=c, alpha=alpha_line, linewidth=1, marker='o', 
                     markersize=4, zorder=z)

        # Decorate Plot
        plt.xticks([x_start, x_end], [label_start, label_end], fontsize=12)
        plt.ylabel('Mean IBS Score', fontsize=12)
        plt.title(title, fontsize=14)
        
        # Add manual legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color=color_increase, lw=2, label='Higher End Value (Increasing)'),
            Line2D([0], [0], color=color_decrease, lw=2, label='Lower End Value (Decreasing)')
        ]
        plt.legend(handles=legend_elements, loc='best')
        
        # Add stats annotation
        stats_txt = (f"Total Samples: {len(df_plot)}\n"
                     f"Increasing Trend: {n_inc}\n"
                     f"Decreasing Trend: {n_dec}")
        
        plt.text(0.05, 0.95, stats_txt, transform=plt.gca().transAxes, 
                 fontsize=10, verticalalignment='top', bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

        plt.grid(axis='y', alpha=0.2, linestyle='--')
        plt.xlim(-0.2, 1.2)
        plt.tight_layout()
        plt.savefig(filename, dpi=300)
        print(f"Plot saved to {filename}")

    # Plot 1: AB -> Others
    if count_others > 0:
        plot_slope_chart(df_res, 'Mean_AB', 'Mean_Others', 
                         'Mean IBS (vs AB)', 'Mean IBS (vs Others)',
                         'IBS Trend: AB Samples (vs AB -> vs Others)',
                         'ab_ibs_trend_ab_vs_others.png')
    
    # Plot 2: A -> AB
    if count_a > 0:
        plot_slope_chart(df_res, 'Mean_A', 'Mean_AB', 
                         'Mean IBS (vs A)', 'Mean IBS (vs AB)',
                         'IBS Trend: AB Samples (vs A -> vs AB)',
                         'ab_ibs_trend_a_vs_ab.png')

    # 5b. New Function: 3-Axis Plot (A -> AB -> Others)
    def plot_three_axis_chart(data, col1, col2, col3, labels, title, filename):
        # Drop NaNs
        df_plot = data.dropna(subset=[col1, col2, col3])
        if len(df_plot) == 0:
            print(f"No valid data for {filename}")
            return
            
        print(f"Generating 3-Axis Chart: {filename}...")
        plt.figure(figsize=(10, 8))
        
        # Coordinates
        x_coords = [0, 1, 2]
        
        # Colors Strategy
        # 1. Up -> Up (Red)
        # 2. Up -> Down (Orange)
        # 3. Down -> Up (Green)
        # 4. Down -> Down (Blue)
        
        c_up_up = '#d62728'      # Red
        c_up_down = '#ff7f0e'    # Orange
        c_down_up = '#2ca02c'    # Green
        c_down_down = '#1f77b4'  # Blue
        
        alpha_line = 0.5
        
        counts = {'Up-Up': 0, 'Up-Down': 0, 'Down-Up': 0, 'Down-Down': 0}
        
        for _, row in df_plot.iterrows():
            v1, v2, v3 = row[col1], row[col2], row[col3]
            
            # Determine Pattern
            # Stage 1: col1 -> col2
            increase1 = (v2 > v1)
            # Stage 2: col2 -> col3
            increase2 = (v3 > v2)
            
            if increase1 and increase2:
                color = c_up_up
                counts['Up-Up'] += 1
                z = 3
            elif increase1 and not increase2:
                color = c_up_down
                counts['Up-Down'] += 1
                z = 2
            elif not increase1 and increase2:
                color = c_down_up
                counts['Down-Up'] += 1
                z = 2
            else: # Down Down
                color = c_down_down
                counts['Down-Down'] += 1
                z = 1
                
            plt.plot(x_coords, [v1, v2, v3], color=color, alpha=alpha_line, 
                     linewidth=1, marker='o', markersize=4, zorder=z)
                     
        # Decoration
        plt.xticks(x_coords, labels, fontsize=12)
        plt.ylabel('Mean IBS Score', fontsize=12)
        plt.title(title, fontsize=14)
        
        # Legend (Outside to right)
        from matplotlib.lines import Line2D
        legend_lines = [
            Line2D([0], [0], color=c_up_up, lw=2, label='Idx: A < AB < Others (Up-Up)'),
            Line2D([0], [0], color=c_up_down, lw=2, label='Idx: A < AB > Others (Up-Down)'),
            Line2D([0], [0], color=c_down_up, lw=2, label='Idx: A > AB < Others (Down-Up)'),
            Line2D([0], [0], color=c_down_down, lw=2, label='Idx: A > AB > Others (Down-Down)')
        ]
        plt.legend(handles=legend_lines, loc='upper left', bbox_to_anchor=(1.02, 1), title="Trend Patterns")
        
        # Stats box
        stats_txt = (f"Total: {len(df_plot)}\n"
                     f"Up-Up: {counts['Up-Up']}\n"
                     f"Up-Down: {counts['Up-Down']}\n"
                     f"Down-Up: {counts['Down-Up']}\n"
                     f"Down-Down: {counts['Down-Down']}")
        plt.text(0.02, 0.98, stats_txt, transform=plt.gca().transAxes,
                 fontsize=10, va='top', bbox=dict(boxstyle='round', fc='white', alpha=0.9))
                 
        plt.grid(axis='y', alpha=0.3, linestyle='--')
        plt.tight_layout()
        plt.savefig(filename, dpi=300)
        print(f"Saved {filename}")

    # Plot 3: 3 Stages
    if count_a > 0 and count_others > 0:
        plot_three_axis_chart(df_res, 
                              col1='Mean_A', col2='Mean_AB', col3='Mean_Others',
                              labels=['Mean IBS (vs A)', 'Mean IBS (vs AB)', 'Mean IBS (vs Others)'],
                              title='IBS Trend (A -> AB -> Others)',
                              filename=f'{output_prefix}_three_axis.png')

if __name__ == "__main__":
    ab_ibs_trend()
