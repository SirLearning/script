import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

plt.style.use('fast')

chr1A_summ = pd.read_table('data/chr1A/stats.nl.txt', sep='\s+', header=0, index_col=0)
chr1ANM_summ = pd.read_table('data/CS/stats.nl.txt', sep='\s+', header=0, index_col=0)
chr1A_summ.index.name = 'Classification'
chr1ANM_summ.index.name = 'Classification'

chr1A_summ['size (Mb)'] = chr1A_summ['size'] / 1000000
chr1ANM_summ['size (Mb)'] = chr1ANM_summ['size'] / 1000000

chr1A_summ.index = chr1A_summ.index.fillna('NULL')
chr1ANM_summ.index = chr1ANM_summ.index.fillna('NULL')

mpl.rcParams['font.size'] = 24
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20

length_summ = sum(chr1ANM_summ['size (Mb)'][:])
print(length_summ/1000000)

# # 1. plot length
# plt.figure(figsize=(16, 12))
# sns.scatterplot(data=chr1A_summ, x='Classification', y='size (Mb)', label='intact chr1A', s=300)
# sns.scatterplot(data=chr1ANM_summ, x='Classification', y='size (Mb)', label='N-sliced chr1A', alpha=0.6, s=300)
# plt.xticks(rotation=45)
# plt.title('Size of TEs by different methods')
# plt.legend()
# plt.show()
#
# # 2. plot proportion
# plt.figure(figsize=(16, 12))
# sns.scatterplot(data=chr1A_summ, x='Classification', y='percent', label='intact chr1A', s=300)
# sns.scatterplot(data=chr1ANM_summ, x='Classification', y='percent', label='N-sliced chr1A', alpha=0.6, s=300)
# plt.xticks(rotation=45)
# plt.title('Proportion of TEs by different methods')
# plt.legend()
# plt.show()
#
# # 3. plot number
# plt.figure(figsize=(16, 12))
# sns.scatterplot(data=chr1A_summ, x='Classification', y='count', label='intact chr1A', s=300)
# sns.scatterplot(data=chr1ANM_summ, x='Classification', y='count', label='N-sliced chr1A', alpha=0.6, s=300)
# plt.xticks(rotation=45)
# plt.title('Counts of TEs by different methods')
# plt.legend()
# plt.show()


# # 4. plot proportion CS
# plt.figure(figsize=(16, 12))
# sns.scatterplot(data=chr1A_summ, x='Classification', y='percent', label='EDTA', s=300)
# sns.scatterplot(data=chr1ANM_summ, x='Classification', y='percent', label='CS v1.0', alpha=0.6, s=300)
# plt.xticks(rotation=45)
# plt.title('Proportion of TEs by different methods')
# plt.legend()
# plt.show()
#
# # 5. plot number CS
# plt.figure(figsize=(16, 12))
# sns.scatterplot(data=chr1A_summ, x='Classification', y='count', label='EDTA', s=300)
# sns.scatterplot(data=chr1ANM_summ, x='Classification', y='count', label='CS v1.0', alpha=0.6, s=300)
# plt.xticks(rotation=45)
# plt.title('Counts of TEs by different methods')
# plt.legend()
# plt.show()

# 4. plot proportion CS
plt.figure(figsize=(16, 12))
sns.barplot(data=chr1A_summ, x='Classification', y='percent', label='EDTA')
sns.barplot(data=chr1ANM_summ, x='Classification', y='percent', label='CS v1.0', alpha=0.6)
plt.xticks(rotation=45)
plt.title('Proportion of TEs by different methods')
plt.legend()
plt.show()

# 5. plot number CS
plt.figure(figsize=(16, 12))
sns.barplot(data=chr1A_summ, x='Classification', y='count', label='EDTA')
sns.barplot(data=chr1ANM_summ, x='Classification', y='count', label='CS v1.0', alpha=0.6)
plt.xticks(rotation=45)
plt.title('Counts of TEs by different methods')
plt.legend()
plt.show()