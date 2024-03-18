import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

chr1A_allTE_summ = pd.read_table('data/01chr1A.anno.stats', sep='\s+', header=0, index_col=0)
chr1A_allTE_summ.index.name = 'Classification'

print(chr1A_allTE_summ.head())

# 1. plot length
plt.figure(figsize=(10, 10))
sns.barplot(data=chr1A_allTE_summ, x='Classification', y='size')
plt.xticks(rotation=45)
plt.show()

# 2. plot proportion
plt.figure(figsize=(10, 10))
sns.barplot(data=chr1A_allTE_summ, x='Classification', y='percent')
plt.xticks(rotation=45)
plt.show()

# 3. plot number
plt.figure(figsize=(10, 10))
sns.barplot(data=chr1A_allTE_summ, x='Classification', y='count')
plt.xticks(rotation=45)
plt.show()
