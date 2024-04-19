import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

plt.style.use('seaborn-v0_8-deep')

intact_summ = pd.read_table('data/np.stats', sep='\s+', header=0, index_col=0)
intact_summ.index.name = 'Classification'
intact_summ['size (Mb)'] = intact_summ['size'] / 1000000
intact_summ.index = intact_summ.index.fillna('NULL')

TElib_summ = pd.read_table('data/TElib.stats', sep='\s+', header=0, index_col=0)
TElib_summ.index.name = 'Classification'
TElib_summ['size (Mb)'] = TElib_summ['size'] / 1000000
TElib_summ.index = TElib_summ.index.fillna('NULL')

mpl.rcParams['font.size'] = 24
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20


# 1. plot proportion
plt.figure(figsize=(16, 12))
sns.barplot(data=intact_summ, x='Classification', y='percent', label='intact chr1A')
sns.barplot(data=TElib_summ, x='Classification', y='percent', label='TElib chr1A', alpha=0.6)
plt.xticks(rotation=45)
plt.title('Proportion of TEs by different methods')
plt.legend()
plt.show()

# 2. plot number
plt.figure(figsize=(16, 12))
sns.barplot(data=intact_summ, x='Classification', y='count', label='intact chr1A')
sns.barplot(data=TElib_summ, x='Classification', y='count', label='TElib chr1A', alpha=0.6)
plt.xticks(rotation=45)
plt.title('Counts of TEs by different methods')
plt.legend()
plt.show()