
from turtle import pd
import pandas as pd
import matplotlib.pyplot as plt

def plot_individual_depth(
    input_file,
    output_prefix="ind_depth"
):
    """
    Plots histogram of Individual Mean Depth from .idepth file.
    Input: .idepth file
    Output: Histogram plot
    """
    print(f"[Info] Processing Individual Depth: {input_file}")
    try:
        df = pd.read_csv(input_file, sep=r'\s+')
        if 'MEAN_DEPTH' not in df.columns:
            print(f"[Error] 'MEAN_DEPTH' column not found in {input_file}")
            return

        plt.figure(figsize=(10, 8))
        sns.histplot(data=df, x='MEAN_DEPTH', binwidth=1, color="firebrick", edgecolor="black", kde=False)
        plt.title("Individual Mean Depth")
        plt.xlabel("Mean Depth")
        plt.ylabel("Count")
        
        outfile = f"{output_prefix}.png"
        plt.tight_layout()
        plt.savefig(outfile, dpi=300)
        plt.close()
        print(f"[Success] Plot saved to {outfile}")
    except Exception as e:
        print(f"[Error] Failed to plot individual depth: {e}")


def missing_vs_depth(
    depth_file="/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap/all.A.taxaBamMap.txt",
    miss_file="/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.missing.smiss",
    output_prefix="depth_vs_missing"
):
    # 1. 读取数据
    # 深度数据：注意第一列叫 Taxa
    df_depth = pd.read_csv(depth_file, sep='\t')
    # 缺失率数据：注意第一列叫 #IID
    df_miss = pd.read_csv(miss_file, sep='\s+')

    # 2. 合并数据
    # left_on 指第一个表的 ID 列，right_on 指第二个表的 ID 列
    # 这会自动处理顺序不一致的问题，只有 ID 完全相同的行才会被合并
    merged = pd.merge(df_depth, df_miss, left_on='Taxa', right_on='#IID')

    # 3. 打印检查
    print(f"深度文件样本数: {len(df_depth)}")
    print(f"缺失率文件样本数: {len(df_miss)}")
    print(f"成功匹配的样本数: {len(merged)}")

    # 计算并打印平均深度
    mean_depth = df_depth['Coverage-Of-All-Bams'].mean()
    print(f"所有样本的平均测序深度: {mean_depth:.4f}")

    # 4. 回归分析与绘图
    # 计算回归参数
    clean_data = merged.dropna(subset=['Coverage-Of-All-Bams', 'F_MISS'])
    x = clean_data['Coverage-Of-All-Bams']
    y = clean_data['F_MISS']
    slope, intercept = np.polyfit(x, y, 1)
    r_squared = np.corrcoef(x, y)[0, 1] ** 2

    plt.figure(figsize=(10, 6))
    # 使用 regplot 自动计算回归线 $y = ax + b$
    sns.regplot(x='Coverage-Of-All-Bams', y='F_MISS', data=merged, 
                scatter_kws={'alpha':0.4, 's':10, 'color':'gray'}, 
                line_kws={'color':'red', 'label':'Linear Regression'})

    # 标注方程和R方 (坐标 0.5, 0.9 表示在图上方中间位置)
    plt.text(0.5, 0.9, f'$y = {slope:.4f}x + {intercept:.4f}$\n$R^2 = {r_squared:.4f}$', 
            transform=plt.gca().transAxes, fontsize=12, color='darkred',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray', boxstyle='round'))

    plt.title('Correlation: Sequencing Depth vs. Missing Rate')
    plt.xlabel('Coverage (X)')
    plt.ylabel('Missing Rate (F_MISS)')
    plt.ylim(0, 1)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.savefig(f'{output_prefix}.png')
    print(f"Plot saved to {output_prefix}.png")

    # 5. 打印深度大于 20 的样本名
    high_depth_samples = merged[merged['Coverage-Of-All-Bams'] > 20]
    print(f"\n========= Samples with Coverage > 20 (Count: {len(high_depth_samples)}) =========")
    # 打印 Taxa 和 F_MISS 两列
    print(high_depth_samples[['Taxa', 'F_MISS']].to_string(index=False))

    # 6. 打印缺失率高于 0.8 的样本名
    high_miss_samples = merged[merged['F_MISS'] > 0.8]
    print(f"\n========= Samples with Missing Rate > 0.8 (Count: {len(high_miss_samples)}) =========")
    # 打印 Taxa 和 F_MISS 两列
    print(high_miss_samples[['Taxa', 'F_MISS']].to_string(index=False))


def ref_ibs_vs_mapping(
    mapping_file="/data1/dazheng_tusr1/vmap4.VCF.v1/vmap4_v1_idxstat_summary.txt", 
    scount_file="/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.sample_stats.scount", 
    group_file="/data1/dazheng_tusr1/vmap4.VCF.v1/sample_group.txt", 
    output_prefix="ref_ibs_vs_mapping_zoomed"
):
    print("Processing Ref IBS vs Mapping Rate Analysis...")
    
    # 1. Read Mapping Rate Data
    print(f"Reading Mapping Rate file: {mapping_file}...")
    try:
        df_map = pd.read_csv(mapping_file, sep='\t')
    except Exception as e:
        print(f"Error reading mapping file: {e}")
        return
        
    if 'TaxaID' not in df_map.columns or 'Mapping_Rate_Pct' not in df_map.columns:
        print(f"Error: Required columns ('TaxaID', 'Mapping_Rate_Pct') not found in mapping file.")
        print(f"Columns found: {list(df_map.columns)}")
        return
        
    df_map = df_map[['TaxaID', 'Mapping_Rate_Pct']].rename(columns={'TaxaID': 'Sample'})
    
    # 2. Read .scount Data & Calculate IBS
    print(f"Reading .scount file: {scount_file}...")
    try:
        df_scount = pd.read_csv(scount_file, sep=r'\s+')
    except Exception as e:
        print(f"Failed to read scount file: {e}")
        return

    required_cols = ['#IID', 'HOM_REF_CT', 'HET_SNP_CT', 'HOM_ALT_SNP_CT']
    if not all(col in df_scount.columns for col in required_cols):
        print(f"Error: Missing columns in scount. Found: {df_scount.columns}")
        return

    # Calculate IBS_ref
    df_scount['Total_Sites'] = df_scount['HOM_REF_CT'] + df_scount['HOM_ALT_SNP_CT'] + df_scount['HET_SNP_CT']
    df_scount = df_scount[df_scount['Total_Sites'] > 0].copy()

    df_scount['Total_Alleles'] = 2 * df_scount['Total_Sites']
    df_scount['Ref_Alleles'] = 2 * df_scount['HOM_REF_CT'] + df_scount['HET_SNP_CT']
    df_scount['IBS_Ref'] = df_scount['Ref_Alleles'] / df_scount['Total_Alleles']

    df_ibs = df_scount[['#IID', 'IBS_Ref']].rename(columns={'#IID': 'Sample'})
    
    # 3. Merge Datasets
    print("Merging datasets...")
    df_merged = pd.merge(df_map, df_ibs, on='Sample', how='inner')
    print(f"Matched samples: {len(df_merged)}")
    
    if len(df_merged) == 0:
        print("Error: No intersecting samples found between mapping file and scount file.")
        return
        
    # 4. Integrate Group Info
    print(f"Reading Group info: {group_file}...")
    if os.path.exists(group_file):
        try:
            df_group = pd.read_csv(group_file, sep=r'\s+', header=None, names=['Sample', 'Group'])
            df_group = df_group.drop_duplicates(subset=['Sample'])
            df_merged = pd.merge(df_merged, df_group, on='Sample', how='left')
            df_merged['Group'] = df_merged['Group'].fillna('Unknown')
        except Exception as e:
            print(f"Warning: Failed to process group file: {e}")
            df_merged['Group'] = 'Unknown'
    else:
        print(f"Warning: Group file not found at {group_file}")
        df_merged['Group'] = 'Unknown'
        
    # Print basic stats
    print(f"Correlation (Pearson): {df_merged['Mapping_Rate_Pct'].corr(df_merged['IBS_Ref']):.4f}")
    
    # Set style
    sns.set_theme(style="ticks")
    
    # 5. Plot
    output_filename = f"{output_prefix}_reg_ref_ibs_vs_map.png"
    plot_dual_regression(df_merged, 'Mapping_Rate_Pct', 'IBS_Ref',
                         'Mapping Rate (%)', 'IBS with Reference Genome',
                         output_filename)
                         
    print("Analysis Complete.")


