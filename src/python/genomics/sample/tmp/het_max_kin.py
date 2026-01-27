import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 1. 读取杂合度数据 (.scount)
print("正在读取 .scount 数据...")
# Use raw string for regex separator
df_het = pd.read_csv("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.sample_stats.scount", sep=r'\s+')

# 计算已观测到的位点总数 (不含缺失)
# Column names corrected: HET_CT -> HET_SNP_CT, HOM_ALT_CT -> HOM_ALT_SNP_CT
df_het['Total_Called'] = df_het['HOM_REF_CT'] + df_het['HET_SNP_CT'] + df_het['HOM_ALT_SNP_CT']

# 计算杂合率
df_het['Het_Rate'] = df_het['HET_SNP_CT'] / df_het['Total_Called']

# 准备合并用的 dataframe
df_het_clean = df_het[['#IID', 'Het_Rate']].rename(columns={'#IID': 'Sample'})

# 2. 读取亲缘关系数据 (.kin0)
print("正在计算每个样本的最大 Kinship...")
# Use raw string for regex separator
df_kin = pd.read_csv("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.king.kin0", sep=r'\s+')

# 提取所有对子中的最大值
# Column names corrected: ID1 -> #IID1, ID2 -> IID2
k1 = df_kin[['#IID1', 'KINSHIP']].rename(columns={'#IID1': 'Sample'})
k2 = df_kin[['IID2', 'KINSHIP']].rename(columns={'IID2': 'Sample'})
df_max_kin = pd.concat([k1, k2]).groupby('Sample')['KINSHIP'].max().reset_index()
df_max_kin.rename(columns={'KINSHIP': 'Max_Kinship'}, inplace=True)

# 3. 合并数据
plot_data = pd.merge(df_het_clean, df_max_kin, on='Sample', how='left').fillna(0)

# 4. 绘图
plt.figure(figsize=(10, 7))
sns.set_style("whitegrid")

sns.scatterplot(
    data=plot_data, 
    x='Het_Rate', 
    y='Max_Kinship', 
    alpha=0.5, 
    s=25, 
    edgecolor='w',
    color='teal'
)

# 阈值线
plt.axhline(y=0.354, color='darkred', linestyle='--', label='Duplicate (Kinship > 0.354)')

# 自动计算小麦的杂合率异常线 (均值 + 3倍标准差)
het_mean = plot_data['Het_Rate'].mean()
het_std = plot_data['Het_Rate'].std()
het_limit = het_mean + 3 * het_std

plt.axvline(x=het_limit, color='orange', linestyle='--', label=f'High Het (> {het_limit:.4f})')

plt.title('Wheat Vmap4 QC: Heterozygosity vs. Max Kinship', fontsize=14)
plt.xlabel('Heterozygosity Rate [HET_SNP_CT / Total_Called]', fontsize=12)
plt.ylabel('Maximum Kinship Coefficient', fontsize=12)
plt.legend()

plt.tight_layout()
output_file = "het_vs_max_kinship_final.png"
plt.savefig(output_file, dpi=300)
print(f"绘图成功！已保存到 {output_file}")
print(f"异常杂合率阈值设定为: {het_limit:.4f} (Mean: {het_mean:.4f}, Std: {het_std:.4f})")
