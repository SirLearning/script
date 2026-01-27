import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.linear_model import LinearRegression, HuberRegressor

def analyze_mac_vs_heterozygosity(gcount_file):
    print(f"正在读取 {gcount_file} ...")
    # 读取 .gcount 文件
    try:
        df = pd.read_csv(gcount_file, sep=r'\s+')
    except Exception as e:
        print(f"读取文件失败: {e}")
        return
    
    # --- 列名匹配 ---
    # User provided: #CHROM ID REF ALT HOM_REF_CT HET_REF_ALT_CTS TWO_ALT_GENO_CTS HAP_REF_CT HAP_ALT_CTS MISSING_CT
    required = ['HET_REF_ALT_CTS', 'TWO_ALT_GENO_CTS', 'HOM_REF_CT', 'HAP_REF_CT', 'HAP_ALT_CTS']
    if not all(col in df.columns for col in required):
        print(f"错误: 文件缺少必要列。现有列: {list(df.columns)}")
        return
    
    # 1. 计算 Counts
    df['Alt_Count'] = df['HET_REF_ALT_CTS'] + (df['TWO_ALT_GENO_CTS'] * 2) + df['HAP_ALT_CTS']
    df['Ref_Count'] = df['HET_REF_ALT_CTS'] + (df['HOM_REF_CT'] * 2) + df['HAP_REF_CT']
    df['Total_Alleles'] = df['Alt_Count'] + df['Ref_Count']
    
    # 2. 计算 MAC (Minor Allele Count)
    df['MAC'] = df[['Alt_Count', 'Ref_Count']].min(axis=1)
    
    # 3. 计算 观察杂合度 (Observed Heterozygosity)
    # Hobs = Het_Count / Total_Samples (Total Samples = HomRef + Het + HomAlt + HapRef + HapAlt)
    # 这里注意： plink2 输出的 genotype counts 是样本数，不是 allele 数
    # Total Samples = HOM_REF_CT + HET_REF_ALT_CTS + TWO_ALT_GENO_CTS + HAP_REF_CT + HAP_ALT_CTS
    df['Total_Samples'] = df['HOM_REF_CT'] + df['HET_REF_ALT_CTS'] + df['TWO_ALT_GENO_CTS'] + df['HAP_REF_CT'] + df['HAP_ALT_CTS']
    df['Hobs'] = df['HET_REF_ALT_CTS'] / df['Total_Samples']
    
    # --- 统计部分 ---
    total_variants = len(df)
    print("-" * 45)
    print(f"总位点数: {total_variants}")
    print("-" * 45)
    print(f"{'MAC':<10} {'Count':<15} {'Percentage':<15}")
    print("-" * 45)
    
    for i in range(1, 11):
        count = (df['MAC'] == i).sum()
        pct = count / total_variants
        print(f"{i:<10} {count:<15} {pct:.4%}")
        
    print("-" * 45)

    # --- 绘图部分: MAC 1-100 vs Heterozygosity ---
    print("正在生成 MAC (1-100) vs Heterozygosity 回归图...")
    
    # 筛选数据: MAC 在 1 到 100 之间
    df_plot = df[(df['MAC'] >= 1) & (df['MAC'] <= 100) & (df['Hobs'] > 0)].copy()
    
    if len(df_plot) < 10:
        print("数据过少，无法绘图。")
        return

    # Regression
    X = df_plot['MAC'].values.reshape(-1, 1)
    y = df_plot['Hobs'].values
    
    # 1. OLS
    ols = LinearRegression()
    ols.fit(X, y)
    ols_slope = ols.coef_[0]
    ols_intercept = ols.intercept_
    ols_score = ols.score(X, y)
    ols_eq = f"y = {ols_slope:.4f}x + {ols_intercept:.4f}"
    
    # 2. Huber
    huber = HuberRegressor(epsilon=1.35)
    huber.fit(X, y)
    huber_slope = huber.coef_[0]
    huber_intercept = huber.intercept_
    huber_score = huber.score(X, y)
    huber_eq = f"y = {huber_slope:.4f}x + {huber_intercept:.4f}"
    
    # Plotting
    plt.figure(figsize=(10, 8))
    
    # Scatter (带透明度)
    sns.scatterplot(data=df_plot, x='MAC', y='Hobs', 
                    alpha=0.3, s=15, color='royalblue', edgecolor=None, label='Data')
    
    # Regression Lines
    x_range = np.linspace(1, 100, 100).reshape(-1, 1)
    y_ols = ols.predict(x_range)
    y_huber = huber.predict(x_range)
    
    plt.plot(x_range, y_ols, color='orange', linewidth=2, linestyle='--', 
             label=f'OLS: {ols_eq}, $R^2$={ols_score:.3f}')
    plt.plot(x_range, y_huber, color='red', linewidth=2, 
             label=f'Huber: {huber_eq}, $R^2$={huber_score:.3f}')
    
    plt.title('Minor Allele Count (MAC 1-100) vs Heterozygosity')
    plt.xlabel('MAC (Minor Allele Count)')
    plt.ylabel('Observed Heterozygosity')
    plt.xlim(0, 105)
    plt.ylim(0, df_plot['Hobs'].max() * 1.1)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    output_file = "mac_vs_het_regression.png"
    plt.savefig(output_file, dpi=300)
    print(f"图表已保存至: {output_file}")

    mac_num = 200
    # --- 绘图部分 2: MAC vs Het Fraction (杂合成分占比) ---
    print(f"正在生成 MAC (1-{mac_num}) vs Het Fraction 热图...")

    # 计算 Het Fraction: 多少比例的 MAC 是由 Heterozygotes 贡献的
    # Het Fraction = HET_REF_ALT_CTS / MAC
    df['Het_Fraction'] = df['HET_REF_ALT_CTS'] / df['MAC']
    
    # 筛选数据: MAC 在 1 到 mac_num 之间
    # 注意：这里我们保留 Het_Fraction = 0 的点
    df_frac = df[(df['MAC'] >= 1) & (df['MAC'] <= mac_num)].copy()
    # 移除空值 (如果有)
    df_frac = df_frac.dropna(subset=['Het_Fraction'])
    
    if len(df_frac) < 10:
        print("数据过少，无法绘制 Het Fraction 图。")
        return

    # --- 构造热图数据 ---
    from matplotlib.colors import LogNorm
    
    # X轴: MAC (1-200)
    # Y轴: Het Fraction (0-1), 步长 0.005 -> 200个bins
    
    mac_min, mac_max = 1, mac_num
    bin_width = 0.05
    y_bins = np.arange(0, 1.0 + bin_width, bin_width) # 0, 0.005, ..., 1.0
    
    # 结果矩阵列表
    heatmap_norm_cols = []
    heatmap_count_cols = []
    
    # 遍历每个 MAC 值 (从 1 到 mac_num)
    for mac in range(mac_min, mac_max + 1):
        subset = df_frac[df_frac['MAC'] == mac]
        if len(subset) == 0:
            # 如果该 MAC 没有变异，全填0
            counts = np.zeros(len(y_bins) - 1)
            norm_counts = counts
        else:
            # 计算直方图
            counts, _ = np.histogram(subset['Het_Fraction'], bins=y_bins)
            # 归一化：计算该 MAC 下，各 Fraction bin 占总数的比例
            norm_counts = counts / len(subset)
            
        heatmap_count_cols.append(counts)
        heatmap_norm_cols.append(norm_counts)
    
    # === 1. 绘制比例热图 (Proportion) ===
    # 转换为 numpy 数组并转置
    data_matrix_norm = np.array(heatmap_norm_cols).T
    # 翻转矩阵: 1.0 在顶部
    data_matrix_norm = np.flipud(data_matrix_norm)
    
    plt.figure(figsize=(14, 10))
    ax1 = sns.heatmap(data_matrix_norm, cmap='viridis', 
                     cbar_kws={'label': 'Proportion of Variants within MAC'})
    
    # 设置 X 轴标签
    x_ticks = np.arange(0, mac_num, mac_num // 10) # 0, 20, ... 180
    x_labels = [str(x + 1) for x in x_ticks]
    
    # 设置 Y 轴标签
    n_bins = data_matrix_norm.shape[0]
    label_step = 0.1
    tick_step = int(label_step / bin_width)
    if tick_step < 1: tick_step = 1
    y_ticks = np.arange(0, n_bins + 1, tick_step)
    y_ticks = y_ticks[y_ticks <= n_bins]
    y_labels = [f"{1.0 - (i * bin_width):.1f}" for i in y_ticks]
    
    def set_axis_labels(ax):
        ax.set_xticks(x_ticks + 0.5)
        ax.set_xticklabels(x_labels, rotation=0)
        ax.set_xlabel('MAC (Minor Allele Count)')
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_labels, rotation=0)
        ax.set_ylabel('Het Fraction (HET_Count / MAC)')

    set_axis_labels(ax1)
    
    plt.title('Heatmap: Het Fraction Distribution per MAC (Proportion)')
    
    output_file_frac = f"mac_{mac_num}_vs_het_fraction_heatmap_0.05.png"
    plt.savefig(output_file_frac, dpi=300)
    print(f"图表已保存至: {output_file_frac}")

    # === 2. 绘制计数热图 (Count) ===
    print(f"正在生成 MAC (1-{mac_num}) vs Het Fraction (Count) 热图...")
    
    data_matrix_count = np.array(heatmap_count_cols).T
    data_matrix_count = np.flipud(data_matrix_count)
    
    plt.figure(figsize=(14, 10))
    
    # 使用 LogNorm 来展示数量差异
    # vmin=1 确保 log(1)=0, 并且忽略 0 值 (sns heatmap 会将其视为无数据或根据cmap处理)
    ax2 = sns.heatmap(data_matrix_count, cmap='viridis', 
                      norm=LogNorm(vmin=1),
                      cbar_kws={'label': 'Count of Variants (Log Scale)'})
    
    set_axis_labels(ax2)
    plt.title('Heatmap: Het Fraction Distribution per MAC (Count)')
    
    output_file_count = f"mac_{mac_num}_vs_het_fraction_count_heatmap_0.05.png"
    plt.savefig(output_file_count, dpi=300)
    print(f"图表已保存至: {output_file_count}")

# 运行
if __name__ == "__main__":
    analyze_mac_vs_heterozygosity("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.variant_stats.gcount")