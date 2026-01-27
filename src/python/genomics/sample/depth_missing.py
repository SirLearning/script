import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# 1. 读取数据
# 深度数据：注意第一列叫 Taxa
df_depth = pd.read_csv("/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap/all.A.taxaBamMap.txt", sep='\t')
# 缺失率数据：注意第一列叫 #IID
df_miss = pd.read_csv("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.missing.smiss", sep='\s+')

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
plt.savefig('depth_vs_missing.png')
print("Plot saved to depth_vs_missing.png")

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

