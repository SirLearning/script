import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# def plot_length(triticeae):
#     length = pd.read_table('data/' + triticeae + '/stats.length.txt', sep='\t')
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


def plot_length(name):
    length = pd.read_table('data/' + name + '/stats.length.txt', sep='\t')
    length.columns = ['Classification', 'width']
    length['width (kb)'] = length['width'] / 1000
    length['triticeae'] = name  # 添加一个新的列 'triticeae'
    return length


fig, ax = plt.subplots()
ax.figure.set_size_inches(24, 12)

# Concatenate the dataframes
length = pd.concat([plot_length('CS'), plot_length('EDTA'), plot_length('trep')])

length.reset_index(drop=True, inplace=True)

sns.boxplot(x='Classification', y='width (kb)', hue='triticeae', data=length, showfliers=False,
            palette='Set3', linewidth=2, dodge=True, ax=ax, width=0.8)
plt.xticks(rotation=45)
plt.title('TE Length Distribution in different annotation')
plt.legend(title='Name')  # 添加图例
plt.show()