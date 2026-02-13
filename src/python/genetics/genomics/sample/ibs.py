from .sample_utils import load_df_from_plink2, load_matrix_file, load_ids_file
from genetics.germplasm.sample import anno_group
import numpy as np
import pandas as pd
from infra.utils.io import save_df_to_tsv
from infra.utils.graph import plot_clustermap, plot_regression_comparison, plot_scatter_with_thresholds, plot_slope_chart
import seaborn as sns
import sys

def ana_ibs_matrix(
    matrix_file,
    id_file,
    group_file=None,
    output_prefix="ibs_ana",
    plot_heatmap=True
):
    print(f"Processing IBS Matrix: {matrix_file}")
    
    # 1. Load IDs
    ids = load_ids_file(id_file)
    if not ids:
        print("Error: No IDs found.")
        return
    n_samples = len(ids)
    
    # 2. Load Matrix
    df_mat = load_matrix_file(matrix_file)
    if df_mat is None: return
    
    # Validate shape
    if df_mat.shape != (n_samples, n_samples):
        print(f"Error: Matrix shape {df_mat.shape} != {n_samples}x{n_samples}")
        return
        
    df_mat.index = ids
    df_mat.columns = ids
    
    # Fill NaNs
    if df_mat.isnull().values.any():
        df_mat = df_mat.fillna(0)
    
    # 3. Heatmap
    if plot_heatmap:
        print("Generating Heatmap...")
        # Prepare Colors  if group file offered
        row_colors = None
        if group_file:
            # Use anno_group to get groups
            df_ids = pd.DataFrame({'Sample': ids})
            df_grouped = anno_group(df_ids, group_file)
            
            if df_grouped is not None and 'Group' in df_grouped.columns:
                # Map Groups to Colors
                unique_groups = df_grouped['Group'].unique()
                palette = sns.color_palette("hsv", len(unique_groups) + 1)
                color_map = dict(zip(unique_groups, palette))
                color_map['Unknown'] = (0.5, 0.5, 0.5)
                
                # Make sure order matches ids (which is index of matrix)
                # reindex df_grouped to match ids order if merge shuffled it?
                # anno_group does left join on input df, so order should be preserved if we pass df_ids
                # but pandas replace/merge might change index. 
                # Let's verify by setting index
                df_grouped = df_grouped.set_index('Sample').reindex(ids)
                
                sample_colors = [color_map.get(grp, (0.5,0.5,0.5)) for grp in df_grouped['Group']]
                row_colors = pd.Series(sample_colors, index=ids, name='Group')
            
        plot_clustermap(
            data_matrix=df_mat,
            row_colors=row_colors,
            title=f"IBS Heatmap ({n_samples} samples)",
            filename=f"{output_prefix}.heatmap.png",
            xticklabels=False, yticklabels=False
        )

def ana_het_vs_max_ibs(
    scount_file, 
    matrix_file,
    id_file,
    output_prefix="het_vs_max_ibs"
):
    print(f"Processing Het vs Max IBS...")
    
    # 1. Load Het
    df_het = load_df_from_plink2(scount_file)
    if df_het is None: return
    
    if 'Het_Rate' not in df_het.columns:
        if 'HET_SNP_CT' in df_het.columns:
             total = df_het['HOM_REF_CT'] + df_het['HET_SNP_CT'] + df_het['HOM_ALT_SNP_CT']
             df_het['Het_Rate'] = df_het['HET_SNP_CT'] / total
    
    # 2. Load IBS Matrix & IDs
    ids = load_ids_file(id_file)
    df_mat = load_matrix_file(matrix_file)
    if df_mat is None: return
    
    mat = df_mat.values
    
    # 3. Calculate Max IBS (excluding diagonal)
    np.fill_diagonal(mat, 0)
    max_ibs = np.max(mat, axis=1)
    # partner_idx = np.argmax(mat, axis=1)
    
    df_max_ibs = pd.DataFrame({
        'Sample': ids,
        'Max_IBS': max_ibs
    })
    
    # 4. Merge
    df_merged = pd.merge(df_het, df_max_ibs, on='Sample', how='inner')
    
    # 5. Plot
    het_mean = df_merged['Het_Rate'].mean()
    het_std = df_merged['Het_Rate'].std()
    het_limit = het_mean + 3 * het_std
    
    thresholds_h = [{'value': 0.95, 'label': 'Duplicate > 0.95', 'color': 'red', 'linestyle': '--'}]
    thresholds_v = [{'value': het_limit, 'label': f'High Het (> {het_limit:.4f})', 'color': 'orange', 'linestyle': '--'}]

    plot_scatter_with_thresholds(
        data=df_merged, 
        x_col='Het_Rate', 
        y_col='Max_IBS',
        title='Heterozygosity vs. Max IBS',
        filename=f"{output_prefix}.png",
        thresholds_h=thresholds_h,
        thresholds_v=thresholds_v,
        xlabel='Heterozygosity Rate',
        ylabel='Maximum IBS Similarity'
    )

def ana_ibs_trend(
    matrix_file,
    id_file,
    group_file,
    output_prefix="ibs_trend"
):
    """
    Analyzes IBS trends between groups (A, AB, Others).
    Ideally needs customized plotting for slope charts or multi-axis charts.
    """
    print("Processing IBS Trend Analysis...")
    
    # 1. Load Data
    ids = load_ids_file(id_file)
    df_mat = load_matrix_file(matrix_file)
    if df_mat is None: return
    
    if len(ids) != df_mat.shape[0]:
        print("Size mismatch.")
        return
    
    mat = df_mat.values
    
    # 2. Define Groups with anno_group
    df_ids = pd.DataFrame({'Sample': ids})
    df_grouped = anno_group(df_ids, group_file)
    
    if df_grouped is None or 'Group' not in df_grouped.columns:
        print("Group annotation failed.")
        return
        
    # Re-align with ids list ordering
    df_grouped = df_grouped.set_index('Sample').reindex(ids)
    
    # Map each sample to A, AB, or Others
    labels = []
    # Fillna just in case
    df_grouped['Group'] = df_grouped['Group'].fillna('Others')
    
    for grp in df_grouped['Group']:
        if grp == 'A' or grp == 'AB': 
            labels.append(grp)
        else:
            labels.append('Others')
            
    labels = np.array(labels)
    
    idx_a = np.where(labels == 'A')[0]
    idx_ab = np.where(labels == 'AB')[0]
    idx_others = np.where(labels == 'Others')[0]
    
    if len(idx_ab) == 0:
        print("No AB samples found.")
        return
        
    # 3. Calculate Means for AB samples
    # For each AB sample, calculate mean IBS with A, with Others, and with other ABs
    
    results = []
    for i in idx_ab:
        # Avoid self if needed, but mean usually okay with self for population stats unless size is small
        # Usually exclude self for "within group" but here we compare to populations.
        
        # Mean to A
        mean_a = np.nan
        if len(idx_a) > 0:
            mean_a = np.mean(mat[i, idx_a])
            
        # Mean to Others
        mean_others = np.nan
        if len(idx_others) > 0:
            mean_others = np.mean(mat[i, idx_others])
            
        # Mean to AB (exclude self)
        # Construct mask for self exclusion
        local_ab_idx = idx_ab[idx_ab != i]
        mean_ab = np.nan
        if len(local_ab_idx) > 0:
            mean_ab = np.mean(mat[i, local_ab_idx])
            
        results.append({
            'Sample': ids[i],
            'Mean_A': mean_a,
            'Mean_Others': mean_others,
            'Mean_AB': mean_ab
        })
        
    df_res = pd.DataFrame(results)
    save_df_to_tsv(df_res, f"{output_prefix}.info.tsv")
    
    # 4. Simple Slope Plot (A -> Others) similar to tmp script
    if 'Mean_A' in df_res.columns and 'Mean_Others' in df_res.columns:
        df_plot = df_res.dropna(subset=['Mean_A', 'Mean_Others']).copy()
        
        plot_slope_chart(
             df=df_plot,
             col1='Mean_A',
             col2='Mean_Others',
             xlabel1='Mean IBS with A',
             xlabel2='Mean IBS with Others',
             title='IBS Trend: AB Samples affinity to A vs Others',
             filename=f"{output_prefix}.slope.png"
        )

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python ibs.py <command> [args]")
        print("Commands: matrix, het_vs_ibs, trend")
        sys.exit(1)
        
    cmd = sys.argv[1]
    
    if cmd == "matrix":
        # python ibs.py matrix mat.mibs mat.mibs.id [group.txt]
        grp = sys.argv[4] if len(sys.argv) > 4 else None
        ana_ibs_matrix(sys.argv[2], sys.argv[3], grp, "ibs_ana")
    elif cmd == "het_vs_ibs":
        ana_het_vs_max_ibs(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5] if len(sys.argv)>5 else "het_vs_ibs")
    elif cmd == "trend":
        ana_ibs_trend(sys.argv[2], sys.argv[3], sys.argv[4], "ibs_trend")
    else:
        print(f"Unknown command: {cmd}")
