import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. 读取数据 (PLINK2生成的 .smiss 文件)
# smiss 文件通常也是空格分隔
# 请确保文件名指向正确的 .smiss 文件路径
file_path = "/data1/dazheng_tusr1/vmap4.VCF.v1/check_missing.smiss" 
# Use raw string for regex separator to avoid warnings
df = pd.read_csv(file_path, sep=r'\s+')

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

plt.savefig('sample_missing_rate_distribution.png', dpi=300)
print("Plot saved to sample_missing_rate_distribution.png")