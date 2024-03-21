import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def plot_length(name):
    length = pd.read_table('data/' + name + '/stats.length.txt', sep='\t')
    length.columns = ['Classification', 'width']
    length['width (kb)'] = length['width'] / 1000

    plt.figure(figsize=(16, 12))
    sns.boxplot(x='Classification', y='width (kb)', data=length, showfliers=False, palette='Set3', linewidth=2.5)
    plt.xticks(rotation=45)
    plt.title('TE Length Distribution in' + name + 'annotation')
    plt.show()


plot_length('CS')
plot_length('EDTA')
