import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import sys
import numpy as np

# 全局设置：用于绘图的降采样数量
PLOT_SAMPLE_SIZE = 1000000 
OUTPUT_PREFIX = "variant_depth_mq_"

print("正在读取深度数据 (popdep)...")
# 1. 读取深度数据文件
try:
    # 假设深度文件格式：Position ... Depth_Mean, Depth_SD
    df_popdep = pd.read_csv("/data/home/tusr1/01projects/vmap4/05reliable.lib/01test/chr002_popdep.txt", 
                            sep='\s+', 
                            usecols=['Position', 'Depth_Mean', 'Depth_SD'])
except Exception as e:
    print(f"读取 popdep 文件失败: {e}")
    sys.exit(1)

# 计算 Variance 和 CV
print("计算 Variance 和 CV...")
df_popdep['Depth_Var'] = df_popdep['Depth_SD'] ** 2
# 处理 Depth_Mean 为 0 导致的 CV 无穷大或 NaN
df_popdep['Depth_CV'] = np.where(df_popdep['Depth_Mean'] > 0, 
                                 df_popdep['Depth_SD'] / df_popdep['Depth_Mean'], 
                                 np.nan)
# 移除 CV 计算失败的行 (可选)
# df_popdep = df_popdep.dropna(subset=['Depth_CV'])

print("正在读取 MQ 数据...")
# 2. 读取 MQ 数据文件
try:
    # 修正：MQ 文件无表头，列顺序：Chrom (0), Pos (1), MQ (2)
    mq_file = "/data/home/tusr1/01projects/vmap4/05reliable.lib/01test/chr002.50.site_mq.txt"
    print(f"Reading MQ file: {mq_file}")
    
    # 显式指定列名
    # 增加 dtype 转换 (强制 MQ 转 numeric, 失败转 NaN)
    df_mq = pd.read_csv(mq_file, sep=r'\s+', header=None, names=['Chrom', 'Position', 'MQ'], low_memory=False)
    
    # 强制将 MQ 列转换为数值型 (将非数字转为 NaN)
    df_mq['MQ'] = pd.to_numeric(df_mq['MQ'], errors='coerce')
    
    # 删除 MQ 为 NaN 的行
    df_mq = df_mq.dropna(subset=['MQ'])
    
    mq_col = 'MQ' # 目标 MQ 列名
    print(f"MQ 数据已读取，样本行: \n{df_mq.head(3)}")

except Exception as e:
    print(f"读取 MQ 文件失败: {e}")
    sys.exit(1)

print("正在匹配数据位置...")
# 3. 合并数据 (内连接)
# 假设两个 DataFrame 都有 Position 列
merged = pd.merge(df_popdep, 
                    df_mq, 
                    on='Position', 
                    how='inner')

if merged.empty:
    print("错误：两个文件之间没有匹配的位置 (Position)，请检查。")
    sys.exit(1)

print(f"匹配完成，共有 {len(merged)} 个共享位点。开始分析...")

# --- 功能 1: 绘制 MQ 分布图 ---
def plot_distribution(df, col, label, filename):
    print(f"--> 正在绘制 {label} 分布图...")
    plt.figure(figsize=(10, 6))
    
    # 降采样
    if len(df) > PLOT_SAMPLE_SIZE:
        df_plot = df.sample(n=PLOT_SAMPLE_SIZE, random_state=42)
    else:
        df_plot = df

    # 直方图
    # 修改：对于 MQ 这种整数型数据 (0-60)，使用 discrete=True 是最准确的绘图方式。
    # 它不会进行 bin 聚合，而是让每个整数值独立成为一个柱子 (相当于 bin 宽度为 1 且中心对齐)。
    sns.histplot(df_plot[col], discrete=True, kde=True, color='teal', stat="density")
    
    # 统计值 (基于全量数据)
    mean_val = df[col].mean()
    median_val = df[col].median()
    min_val = df[col].min()
    max_val = df[col].max()
    
    plt.axvline(mean_val, color='red', linestyle='--', label=f'Mean: {mean_val:.2f}')
    plt.axvline(median_val, color='orange', linestyle='--', label=f'Median: {median_val:.2f}')
    
    plt.plot([], [], ' ', label=f'Min: {min_val:.2f}')
    plt.plot([], [], ' ', label=f'Max: {max_val:.2f}')
    
    plt.title(f'Distribution of {label}')
    plt.xlabel(label)
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"   图片已保存: {filename}")
    plt.close()

# --- 功能 2: 绘制 Regression 散点图 ---
def plot_regression(df, x_col, y_col, x_label, y_label, filename):
    print(f"--> 正在分析 {y_label} vs {x_label} ...")
    
    # 清洗数据
    df_clean = df[[x_col, y_col]].dropna()
    x_full = df_clean[x_col]
    y_full = df_clean[y_col]
    
    # 统计 (全量)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_full, y_full)
    r_squared = r_value**2
    line_eq = f"y = {slope:.4f}x + {intercept:.4f}"
    
    # 降采样
    if len(df_clean) > PLOT_SAMPLE_SIZE:
        df_plot = df_clean.sample(n=PLOT_SAMPLE_SIZE, random_state=42)
    else:
        df_plot = df_clean
        
    plt.figure(figsize=(12, 8))
    
    # 散点
    sns.scatterplot(x=df_plot[x_col], y=df_plot[y_col], 
                    alpha=0.2, s=5, color='royalblue', edgecolor=None, label='Data (Sampled)')
    
    # 回归线
    x_range = np.linspace(df_clean[x_col].min(), df_clean[x_col].max(), 100)
    y_pred = slope * x_range + intercept
    plt.plot(x_range, y_pred, color='red', linewidth=2, label='Regression (All Data)')
    
    # 标注
    plt.text(0.05, 0.95, f"{line_eq}\n$R^2 = {r_squared:.4f}$", 
             transform=plt.gca().transAxes, 
             fontsize=12, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
             
    plt.title(f'{y_label} vs {x_label}')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    # 限制范围：Depth Mean 0-10, MQ 0-70 (通常 MQ 最大 60, 有时 255)
    # 根据数据自动调整，但给个上限防止极值
    x_limit = df_clean[x_col].quantile(0.999)
    y_limit = df_clean[y_col].max() * 1.05
    if y_limit > 70: y_limit = 70 # 强制 MQ 上限为 70 (因为有效范围通常是0-60)
    
    plt.xlim(0, max(10, x_limit)) # 至少显示到 10，或者看数据 99.9%
    plt.ylim(0, y_limit)
    
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"   图片已保存: {filename}")
    plt.close()

# 执行
# 1. MQ 分布
plot_distribution(merged, mq_col, 'Mapping Quality (MQ)', f"{OUTPUT_PREFIX}dist_mq.png")

# 2. Depth Mean vs MQ
plot_regression(merged, 'Depth_Mean', mq_col, 'Mean Depth', 'Mapping Quality (MQ)', 
                f"{OUTPUT_PREFIX}reg_depth_mean_vs_mq.png")

# 3. Depth SD vs MQ
plot_regression(merged, 'Depth_SD', mq_col, 'Depth Standard Deviation', 'Mapping Quality (MQ)', 
                f"{OUTPUT_PREFIX}reg_depth_sd_vs_mq.png")

# 4. Depth Variance vs MQ
plot_regression(merged, 'Depth_Var', mq_col, 'Depth Variance', 'Mapping Quality (MQ)', 
                f"{OUTPUT_PREFIX}reg_depth_var_vs_mq.png")

# 5. Depth CV vs MQ
plot_regression(merged, 'Depth_CV', mq_col, 'Depth CV', 'Mapping Quality (MQ)', 
                f"{OUTPUT_PREFIX}reg_depth_cv_vs_mq.png")

print("所有分析完成！")
