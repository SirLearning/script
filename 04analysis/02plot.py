import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

chr1A_allTE_summ = pd.read_table('data/chr1A.anno.stats.all', sep='\s+', header=0, index_col=0)
chr1ANM_allTE_summ = pd.read_table('data/chr1ANM.anno.stats', sep='\s+', header=0, index_col=0)
chr1A_allTE_summ.index.name = 'Classification'
chr1ANM_allTE_summ.index.name = 'Classification'

chr1A_allTE_summ.index = chr1A_allTE_summ.index.fillna('NULL')
chr1ANM_allTE_summ.index = chr1ANM_allTE_summ.index.fillna('NULL')

mpl.rcParams['font.size'] = 15

# 1. plot length
plt.figure(figsize=(10, 10))
sns.barplot(data=chr1A_allTE_summ, x='Classification', y='size', label='chr1A', color='#6780B6')
sns.barplot(data=chr1ANM_allTE_summ, x='Classification', y='size', label='chr1ANM', color='#A3A1FB', alpha=0.5)
plt.xticks(rotation=45)
plt.legend()
plt.show()

# 2. plot proportion
plt.figure(figsize=(10, 10))
sns.barplot(data=chr1A_allTE_summ, x='Classification', y='percent', label='chr1A', color='#6780B6')
sns.barplot(data=chr1ANM_allTE_summ, x='Classification', y='percent', label='chr1ANM', color='#A3A1FB', alpha=0.5)
plt.xticks(rotation=45)
plt.legend()
plt.show()

# 3. plot number
plt.figure(figsize=(10, 10))
sns.barplot(data=chr1A_allTE_summ, x='Classification', y='count', label='chr1A', color='#6780B6')
sns.barplot(data=chr1ANM_allTE_summ, x='Classification', y='count', label='chr1ANM', color='#A3A1FB', alpha=0.5)
plt.xticks(rotation=45)
plt.legend()
plt.show()