import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def missing_dist(
    input_file="/data1/dazheng_tusr1/vmap4.VCF.v1/check_missing.smiss", 
    output_prefix="sample_missing_rate_distribution"
):
    # 1. 读取数据 (PLINK2生成的 .smiss 文件)
    # smiss 文件通常也是空格分隔
    # 请确保文件名指向正确的 .smiss 文件路径
    # Use raw string for regex separator to avoid warnings
    df = pd.read_csv(input_file, sep=r'\s+')

    # 2. 检查列名
    # 样本缺失率文件 (.smiss) 中，缺失率列通常也是 'F_MISS'
    missing_col = 'F_MISS'

    # Calculate Statistics
    mean_val = df[missing_col].mean()
    median_val = df[missing_col].median()
    print(f"Sample Missing Rate - Mean: {mean_val:.5f}, Median: {median_val:.5f}")

    # 3. 绘图设置
    plt.figure(figsize=(10, 6))

    # Style based on maf_dist.py (forestgreen for general distribution)
    sns.histplot(df[missing_col], bins=100, kde=False, color='forestgreen', edgecolor='none', alpha=0.8)

    # 4. Add Vertical Lines (Thresholds & Stats)
    # Thresholds 0.1 and 0.3
    plt.axvline(x=0.1, color='#c44e52', linestyle='--', linewidth=1.5, label='Threshold 0.1') # Red-ish
    plt.axvline(x=0.3, color='orange', linestyle='--', linewidth=1.5, label='Threshold 0.3')

    # Stats
    plt.axvline(x=mean_val, color='black', linestyle='-', linewidth=1.5, label=f'Mean: {mean_val:.5f}')
    plt.axvline(x=median_val, color='purple', linestyle='-.', linewidth=1.5, label=f'Median: {median_val:.5f}')

    # 5. 添加标题和标签
    plt.title('Distribution of Sample Missing Rate', fontsize=15)
    plt.xlabel('Missing Rate (F_MISS)', fontsize=12)
    plt.ylabel('Count of Samples', fontsize=12)

    # Set X-axis range 0-1
    plt.xlim(0, 1.0)

    # 6. Grid, Legend, Save
    plt.legend(loc='upper right')
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    plt.tight_layout()

    plt.savefig(f'{output_prefix}.png', dpi=300)
    print(f"Plot saved to {output_prefix}.png")

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


def calculate_missing_threshold(input_file, output_prefix):
    """
    Calculates missing rate threshold (Mean + 3SD), generates a distribution plot,
    and prints the threshold to stdout for Nextflow capture.
    Other logs are sent to stderr to avoid polluting stdout.
    """
    import sys
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    try:
        df = pd.read_csv(input_file, sep=r'\s+')
    except Exception as e:
        sys.stderr.write(f"[Error] Failed to read {input_file}: {e}\n")
        print("0.05", end="")
        return

    col = 'F_MISS'
    if col not in df.columns:
        sys.stderr.write(f"[Error] Column {col} not found in {input_file}\n")
        print("0.05", end="")
        return

    # Stats
    mean_val = df[col].mean()
    std_val = df[col].std()
    threshold = mean_val + 3 * std_val
    
    # Fallback if calculation fails (e.g. NaN)
    if pd.isna(threshold):
        threshold = 0.05
    
    # Plotting
    try:
        plt.figure(figsize=(10, 6))
        sns.histplot(df[col], bins=50, kde=True, color='skyblue', edgecolor='black', alpha=0.7)
        
        plt.axvline(mean_val, color='blue', linestyle='--', label=f'Mean: {mean_val:.4f}')
        plt.axvline(threshold, color='red', linestyle='-', linewidth=2, label=f'Threshold (+3SD): {threshold:.4f}')
        
        plt.title(f'Sample Missing Rate Distribution\nThreshold: {threshold:.4f}')
        plt.xlabel('Missing Rate (F_MISS)')
        plt.ylabel('Count')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.3)
        plt.tight_layout()
        
        plot_path = f"{output_prefix}.png"
        plt.savefig(plot_path, dpi=300)
        sys.stderr.write(f"[Info] Plot saved to {plot_path}\n")
    except Exception as e:
        sys.stderr.write(f"[Warning] Plotting failed: {e}\n")

    # Critical: Output ONLY the numeric value to stdout without newline
    print(f"{threshold:.5f}", end="")



