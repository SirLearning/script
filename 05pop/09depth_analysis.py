import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

summ_name = 'data/1B.mosdepth.summary.txt'
TEcode = 'data/TEcode'


def main():
    depth = pd.concat([plot_depth('1A'), plot_depth('1B'), plot_depth('1D')])
    depth.reset_index(drop=True, inplace=True)

    # 1. depth difference
    fig, ax = plt.subplots()
    plt.rcParams['font.size'] = 24
    ax.figure.set_size_inches(18, 9)

    sns.boxplot(x='type', y='depth', hue='chr', data=depth, showfliers=False,
                palette='Set3', linewidth=2, dodge=True, ax=ax, width=0.8)
    plt.title('TE Depth Distribution in different chromosome')
    plt.legend(title='Name')
    plt.show()

    # 2. proportion of TEs
    grouped = depth.groupby(['type', 'chr'])
    proportion = grouped.agg(
        sum_depth=('depth', 'sum'),
    )
    proportion['percent'] = proportion['sum_depth'] / proportion.groupby('chr')['sum_depth'].transform('sum') * 100
    proportion.reset_index(inplace=True)

    plt.figure(figsize=(18, 9))
    sns.barplot(x='type', y='percent', hue='chr', data=proportion, palette='Set3')
    plt.title('Proportion of TEs in CS chromosomes')
    plt.legend()
    plt.show()

    # 3. retry
    proportion_wide = proportion.pivot(index='chr', columns='type', values='percent')
    bottom = np.zeros(len(proportion_wide))

    plt.figure(figsize=(14, 16))
    plt.style.use('seaborn-v0_8-deep')
    for column in proportion_wide.columns:
        plt.bar(proportion_wide.index, proportion_wide[column], bottom=bottom, label=column, alpha=0.8)
        bottom += proportion_wide[column]
    plt.title('Proportion of TEs in CS chromosomes')
    plt.xlabel('Chromosome')
    plt.ylabel('Proportion (%)')
    plt.legend()
    plt.show()

def plot_depth(chr):
    summ = pd.read_table("data/" + chr + ".mosdepth.summary.txt", sep="\t", header=0)
    summ['type'] = summ['chrom'].str.split('#').str[1]
    tecode = pd.read_table(TEcode, sep=",", header=None)
    tecode.columns = ['cls', 'new_cls']
    for i in range(0, len(tecode)):
        summ.loc[summ['type'] == tecode['cls'][i], 'type'] = tecode['new_cls'][i]
    summ['depth'] = summ['mean']
    summ['chr'] = chr
    summ.drop(['chrom', 'mean', 'min', 'max'], axis=1, inplace=True)
    # summ now: chr, type, depth, length, bases
    return summ

def nl_stats(chr):
    allte = pd.read_table("data/stats."+ chr + ".length.txt", sep='\t', header=None)

    return te_length, te_summ

if __name__ == '__main__':
    main()