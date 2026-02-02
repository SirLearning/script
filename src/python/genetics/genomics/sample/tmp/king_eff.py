import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# 1. 加载数据
print("Loading data...")
# 如果是方阵格式，需要先将其展平为一维向量
king_df = pd.read_csv("/data1/dazheng_tusr1/vmap4.VCF.v1/chr001_king.king", sep='\s+', header=None)
kinship_values = king_df.values[np.triu_indices_from(king_df.values, k=1)] # 只取上三角，避开对角线

# 统计原始数据范围
print(f"Original data range: min={kinship_values.min()}, max={kinship_values.max()}")

# 2. 数据清洗与绘图
plt.figure(figsize=(10, 6))

# 过滤掉极端的负值 (例如 < -1)，这些通常是噪音，会严重压缩X轴
# 正常的Kinship系数范围通常在 -0.5 到 0.5 之间
valid_range_mask = kinship_values > -1
filtered_data = kinship_values[valid_range_mask]
n_removed = len(kinship_values) - len(filtered_data)
print(f"Removed {n_removed} extreme negative values (< -1) to improve visualization.")

sns.histplot(filtered_data, bins=100, color='teal', kde=False)

# 3. 装饰
plt.axvline(x=0, color='black', linestyle='-', label='Unrelated (~0)') # 0点线
plt.axvline(x=0.0442, color='red', linestyle='--', label='3rd Degree')
plt.axvline(x=0.177, color='orange', linestyle='--', label='1st Degree')

plt.title('Distribution of KING Kinship Coefficients')
plt.xlabel('Kinship Coefficient')
plt.ylabel('Frequency (Pairs)')
plt.yscale('log') # 建议使用对数纵坐标，因为无关对子数远多于亲属对子

# 限制X轴范围，聚焦于有意义的区域
plt.xlim(-0.5, 0.6) 

plt.legend()
plt.savefig('king_coefficient_distribution.png')
print("Plot saved to king_coefficient_distribution.png")
