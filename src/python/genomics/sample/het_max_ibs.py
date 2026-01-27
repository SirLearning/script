import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

# 设置文件路径
scount_file = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.sample_stats.scount"
ibs_matrix_file = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.mibs"
ibs_id_file = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.ibs.mibs.id"
output_file = "het_vs_max_ibs_final.png"

# 1. 读取杂合度数据 (.scount)
print("正在读取 .scount 数据...")
try:
    df_het = pd.read_csv(scount_file, sep=r'\s+')
except Exception as e:
    print(f"读取 scount 文件失败: {e}")
    sys.exit(1)

# column check
required_cols = ['#IID', 'HOM_REF_CT', 'HET_SNP_CT', 'HOM_ALT_SNP_CT']
missing = [c for c in required_cols if c not in df_het.columns]
if missing:
    print(f"错误: .scount 文件缺少列: {missing}")
    print(f"现有列: {df_het.columns.tolist()}")
    sys.exit(1)

# 计算已观测到的位点总数 (不含缺失)
df_het['Total_Called'] = df_het['HOM_REF_CT'] + df_het['HET_SNP_CT'] + df_het['HOM_ALT_SNP_CT']

# 计算杂合率
df_het['Het_Rate'] = df_het['HET_SNP_CT'] / df_het['Total_Called']

# 准备合并用的 dataframe
df_het_clean = df_het[['#IID', 'Het_Rate']].rename(columns={'#IID': 'Sample'})

# 2. 读取 IBS 数据
print("正在读取 IBS 矩阵数据 (可能需要一些时间)...")
try:
    ibs_matrix = pd.read_csv(ibs_matrix_file, sep=r'\s+', header=None)
    ibs_ids = pd.read_csv(ibs_id_file, sep=r'\s+', header=None)
except Exception as e:
    print(f"读取 IBS 文件失败: {e}")
    sys.exit(1)

# 确保矩阵维度与ID数量匹配
if len(ibs_matrix) != len(ibs_ids):
    print(f"错误: IBS 矩阵行数 ({len(ibs_matrix)}) 与 ID 文件行数 ({len(ibs_ids)}) 不匹配")
    sys.exit(1)

# 3. 提取每个样本的最大亲缘关系值 (Max Kinship)
print("正在计算每行最大 IBS (Diagonal excluded)...")
matrix_values = ibs_matrix.values
# 将对角线（自交相似度）设为 0，避免最大值总是自己 (通常为1)
np.fill_diagonal(matrix_values, 0)
max_ibs = np.max(matrix_values, axis=1)
max_ibs_indices = np.argmax(matrix_values, axis=1)

# 构造绘图数据
# ibs_ids[1] 是 IID 列
ids_vec = ibs_ids[1].values
df_max_ibs = pd.DataFrame({
    'Sample': ids_vec,
    'Max_IBS': max_ibs,
    'Partner_Sample': ids_vec[max_ibs_indices]
})

# 4. 合并数据
# 注意：IBS ID 和 scount ID 需要对应
print("合并数据...")
plot_data = pd.merge(df_het_clean, df_max_ibs, on='Sample', how='inner')

if plot_data.empty:
    print("警告: 合并后数据为空，请检查 Sample ID 是否一致")
    sys.exit(1)

print(f"绘图点数: {len(plot_data)}")

# 5. 绘图
plt.figure(figsize=(10, 8))
sns.set_style("white") # Remove background grid/shading

sns.scatterplot(
    data=plot_data, 
    x='Het_Rate', 
    y='Max_IBS', 
    alpha=0.5, 
    s=25,
    edgecolor='w',
    color='royalblue'
)

# 绘制互为最大IBS的连线 (Reciprocal Best Hits)
print("正在绘制互为最大IBS的连线...")
coord_map = {row['Sample']: (row['Het_Rate'], row['Max_IBS']) for _, row in plot_data.iterrows()}
partner_map = {row['Sample']: row['Partner_Sample'] for _, row in plot_data.iterrows()}
processed_pairs = set()

for sample_a, pair_b in partner_map.items():
    # 检查互为最佳匹配: A的最爱是B，且B的最爱是A
    if pair_b in partner_map and partner_map[pair_b] == sample_a:
        #防止重复绘制
        pair_key = tuple(sorted((sample_a, pair_b)))
        if pair_key not in processed_pairs:
            processed_pairs.add(pair_key)
            
            # 获取坐标
            xa, ya = coord_map[sample_a]
            xb, yb = coord_map[pair_b]
            
            plt.plot([xa, xb], [ya, yb], color='gray', alpha=0.5, linewidth=0.5, zorder=1)

# print statistics about RBH
num_rbh_samples = len(processed_pairs) * 2
total_samples = len(plot_data)
rbh_ratio = num_rbh_samples / total_samples if total_samples > 0 else 0
print(f"Reciprocal Best Hit (RBH) Statistics:")
print(f"  Total Samples: {total_samples}")
print(f"  Samples in RBH pairs: {num_rbh_samples} ({len(processed_pairs)} pairs)")
print(f"  RBH Ratio: {rbh_ratio:.2%}")

# 阈值线
plt.axhline(y=0.95, color='red', linestyle='--', label='Duplicate Threshold (0.95)')

# 自动计算小麦的杂合率异常线 (均值 + 3倍标准差)
het_mean = plot_data['Het_Rate'].mean()
het_std = plot_data['Het_Rate'].std()
het_limit = het_mean + 3 * het_std

plt.axvline(x=het_limit, color='orange', linestyle='--', label=f'High Het (> {het_limit:.4f})')

plt.title('Wheat Vmap4 QC: Heterozygosity vs. Max IBS')
plt.xlabel('Individual Heterozygosity Rate')
plt.ylabel('Maximum IBS Similarity (with any other sample)')
plt.legend()
plt.tight_layout()

plt.savefig(output_file, dpi=300)
print(f"绘图成功！已保存到 {output_file}")
print(f"异常杂合率阈值设定为: {het_limit:.4f} (Mean: {het_mean:.4f}, Std: {het_std:.4f})")
