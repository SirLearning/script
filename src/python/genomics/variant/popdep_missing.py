import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import sys

# 全局设置：用于绘图的降采样数量
PLOT_SAMPLE_SIZE = 50000 

print("正在读取缺失率数据 (vmiss)...")
# 1. 读取 vmiss 文件
try:
    df_vmiss = pd.read_csv("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.missing.vmiss", sep='\s+')
    df_vmiss.columns = [c.replace('#', '') for c in df_vmiss.columns]
    
    # 提取 ID 中的位置信息 (例如 '2-952' -> 952)
    # 假设格式始终为 '染色体-位置-...'，取第二个字段
    df_vmiss['Pos_key'] = df_vmiss['ID'].str.split('-').str[1].astype(int)
except Exception as e:
    print(f"读取 vmiss 文件失败: {e}")
    sys.exit(1)

print("正在读取深度数据 (popdep)...")
# 2. 读取深度数据文件
try:
    df_popdep = pd.read_csv("/data/home/tusr1/01projects/vmap4/05reliable.lib/01test/chr002_popdep.txt", 
                            sep='\s+', 
                            usecols=['Position', 'Depth_SD', 'Depth_Mean'])
except Exception as e:
    print(f"读取 popdep 文件失败: {e}")
    sys.exit(1)

print("正在匹配数据位置...")
# 3. 合并数据 (内连接)
merged = pd.merge(df_vmiss[['Pos_key', 'F_MISS']], 
                    df_popdep, 
                    left_on='Pos_key', 
                    right_on='Position', 
                    how='inner')

if merged.empty:
    print("错误：两个文件之间没有匹配的位置，请检查 ID 格式。")
    sys.exit(1)

print(f"匹配完成，共有 {len(merged)} 个共享位点。开始分析...")

# 定义一个绘图函数
def plot_regression(x_col, x_label, output_filename):
    print(f"-->正在分析 {x_col} vs F_MISS ...")
    
    # 1. 准备数据 (Drop NA)
    # 仅提取需要的列，减少内存拷贝
    df_clean = merged[[x_col, 'F_MISS']].dropna()
    x_full = df_clean[x_col]
    y_full = df_clean['F_MISS']
    
    # 2. 统计分析 (使用全量数据确保准确)
    print("   计算全量数据回归参数...")
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_full, y_full)
    r_squared = r_value**2
    line_eq = f"y = {slope:.4f}x + {intercept:.4f}"

    # 3. 绘图 (使用降采样数据加速)
    plt.figure(figsize=(12, 8))
    
    if len(df_clean) > PLOT_SAMPLE_SIZE:
        print(f"   位点过多 ({len(df_clean)})，随机降采样至 {PLOT_SAMPLE_SIZE} 个用于绘图...")
        df_plot = df_clean.sample(n=PLOT_SAMPLE_SIZE, random_state=42)
    else:
        df_plot = df_clean

    # 绘制散点 (Scatter)
    sns.scatterplot(x=df_plot[x_col], y=df_plot['F_MISS'], 
                    alpha=0.2, s=3, color='royalblue', edgecolor=None, label='Data Points (Sampled)')

    # 绘制回归线 (基于全量数据的方程)
    print("   绘制回归线...")
    x_range = np.linspace(df_clean[x_col].min(), df_clean[x_col].max(), 100)
    y_pred = slope * x_range + intercept
    plt.plot(x_range, y_pred, color='red', linewidth=2, label='Regression Line (All Data)')

    plt.text(0.05, 0.95, f"Equation: {line_eq}\n$R^2 = {r_squared:.4f}$\nP-value = {p_value:.2e}\n(Stats based on {len(df_clean):,} sites)", 
                transform=plt.gca().transAxes, 
                verticalalignment='top', 
                fontsize=12,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    plt.title(f'Correlation between {x_label} and Variant Missing Rate (Chr 002)', fontsize=15)
    plt.xlabel(f'{x_label} ({x_col})', fontsize=12)
    plt.ylabel('Missing Rate (F_MISS)', fontsize=12)
    plt.xlim(0, 10)
    plt.ylim(0, 1)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    print(f"   图片已保存: {output_filename}")
    plt.close()

# 4. 导入numpy (如果之前没导入) 并绘图
import numpy as np

# 计算变异系数 (CV = Depth_SD / Depth_Mean)
# 注意处理 Depth_Mean 为 0 的情况，避免除以零错误
merged['Depth_CV'] = merged['Depth_SD'] / merged['Depth_Mean']
merged.loc[merged['Depth_Mean'] == 0, 'Depth_CV'] = np.nan # 或者设为 0，取决于业务逻辑，这里跳过这些点

# 绘制 Depth_SD vs F_MISS
plot_regression('Depth_SD', 'Depth Standard Deviation', "missing_vs_depth_sd_regression_0.2.png")

# 绘制 Depth_Mean vs F_MISS
plot_regression('Depth_Mean', 'Mean Depth', "missing_vs_depth_mean_regression_0.2.png")

# 绘制 Depth_CV (变异系数) vs F_MISS
plot_regression('Depth_CV', 'Depth Coefficient of Variation (SD/Mean)', "missing_vs_depth_cv_regression_0.2.png")

print("所有分析完成！")
