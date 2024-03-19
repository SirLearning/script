import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

plt.style.use('bmh')

chr1A_summ = pd.read_table('data/chr1A/stats.nl.txt', sep='\s+', header=0, index_col=0)
chr1ANM_summ = pd.read_table('data/chr1ANM/stats.nl.txt', sep='\s+', header=0, index_col=0)
chr1A_summ.index.name = 'Classification'
chr1ANM_summ.index.name = 'Classification'

chr1A_summ.index = chr1A_summ.index.fillna('NULL')
chr1ANM_summ.index = chr1ANM_summ.index.fillna('NULL')

mpl.rcParams['font.size'] = 15

# 1. plot length
plt.figure(figsize=(10, 10))
sns.barplot(data=chr1A_summ, x='Classification', y='size', label='chr1A')
sns.barplot(data=chr1ANM_summ, x='Classification', y='size', label='chr1ANM', alpha=0.75)
plt.xticks(rotation=45)
plt.legend()
plt.show()

# 2. plot proportion
plt.figure(figsize=(10, 10))
sns.barplot(data=chr1A_summ, x='Classification', y='percent', label='chr1A')
sns.barplot(data=chr1ANM_summ, x='Classification', y='percent', label='chr1ANM', alpha=0.75)
plt.xticks(rotation=45)
plt.legend()
plt.show()

# 3. plot number
plt.figure(figsize=(10, 10))
sns.barplot(data=chr1A_summ, x='Classification', y='count', label='chr1A')
sns.barplot(data=chr1ANM_summ, x='Classification', y='count', label='chr1ANM', alpha=0.75)
plt.xticks(rotation=45)
plt.legend()
plt.show()