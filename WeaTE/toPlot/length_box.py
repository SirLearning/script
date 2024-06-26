import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# def plot_length(triticeae):
#     length = pd.read_table('data/' + triticeae + '/stats.length.old.txt', sep='\t')
#     length.columns = ['Classification', 'width']
#     length['width (kb)'] = length['width'] / 1000
#     length['triticeae'] = triticeae  # 添加一个新的列 'triticeae'
#     sns.boxplot(x='Classification', y='width (kb)', hue='triticeae', data=length, showfliers=False,
#                 palette='Set3', linewidth=2.5, ax=ax)
#     plt.xticks(rotation=45)
#
#
# fig, ax = plt.subplots()
# ax.figure.set_size_inches(24, 12)
#
# plot_length('CS')
# plot_length('EDTA')
# plt.title('TE Length Distribution in different annotation')
# plt.legend(title='Name')  # 添加图例
# plt.show()

# output_name = 'data/stats.length.txt'

def plot_length(name):
    length = pd.read_table('../data/old/' + name + '/stats.length.txt', sep='\t')
    length.columns = ['Classification', 'width']
    length['length (kb)'] = length['width'] / 1000
    length['triticeae'] = name  # 添加一个新的列 'triticeae'
    return length


fig, ax = plt.subplots()
plt.rcParams['font.size'] = 24
ax.figure.set_size_inches(26, 10)

# Concatenate the dataframes
length = pd.concat([plot_length('CS'), plot_length('EDTA'), plot_length('curated_lib')])
grouped = length.groupby(['Classification', 'triticeae']).mean()
print(grouped.to_markdown())
# grouped.to_csv(output_name, sep='\t', header=False)


length.reset_index(drop=True, inplace=True)

sns.boxplot(x='Classification', y='length (kb)', hue='triticeae', data=length, showfliers=False,
            palette='Set3', linewidth=2, dodge=True, ax=ax, width=0.8)
plt.xticks(rotation=45)
plt.title('TE Length Distribution in different annotation')
plt.legend(title='Name')  # 添加图例
plt.tick_params(axis='both', labelsize=18)
plt.show()