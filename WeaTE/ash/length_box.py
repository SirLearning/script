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




