import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. 读取数据 (注意：PLINK2的文件通常以 # 开头，read_csv会自动处理或需指定 sep)
# 如果是空格分隔，使用 sep='\s+'
file_path = "/data1/dazheng_tusr1/vmap4.VCF.v1/check_missing.vmiss" 
df = pd.read_csv(file_path, sep='\s+')

# 2. 检查列名，PLINK2 经常在第一列前加 '#'
# 我们需要的是 'F_MISS' 这一列
missing_col = 'F_MISS'

# 3. 绘图设置
plt.figure(figsize=(10, 6))
sns.histplot(df[missing_col], bins=50, kde=False, color='skyblue', edgecolor='black')

# 4. 添加标题和标签
plt.title('Distribution of Variant Missing Rate', fontsize=15)
plt.xlabel('Missing Rate (F_MISS)', fontsize=12)
plt.ylabel('Count of Variants', fontsize=12)

# 5. 可选：添加一条红线表示常用的阈值（如 0.1）
plt.axvline(x=0.1, color='red', linestyle='--', label='Threshold = 0.1')
plt.legend()

plt.grid(axis='y', alpha=0.3)
plt.savefig('missing_rate_distribution.png')
print("Plot saved to missing_rate_distribution.png")