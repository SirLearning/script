import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

TEcode = 'data/TEcode'


def main():
    depth = pd.concat([plot_depth('1A'), plot_depth('1B'), plot_depth('1D')])
    depth.reset_index(drop=True, inplace=True)
    groupe_reads = depth.groupby(['type', 'chr'])
    reads_pct = to_group(groupe_reads)
    ref_pct = pd.concat([nl_stats('chr1A'), nl_stats('chr1B'), nl_stats('chr1D')])

    # 1. ref_depth difference
    fig, ax = plt.subplots()
    plt.rcParams['font.size'] = 24
    ax.figure.set_size_inches(18, 9)

    sns.boxplot(x='type', y='ref_depth', hue='chr', data=depth, showfliers=False,
                palette='Set3', linewidth=2, dodge=True, ax=ax, width=0.8)
    plt.title('TE Depth Distribution in different chromosome')
    plt.legend(title='Name')
    plt.show()

    # 2. reads_pct of TEs
    fig, ax = plt.subplots()
    ax.figure.set_size_inches(18, 9)

    sns.barplot(x='type', y='percent', hue='chr', data=reads_pct, palette='Set3', ax=ax)
    sns.barplot(x='type', y='percent', hue='chr', data=ref_pct, ax=ax, alpha=0.5)
    plt.title('Proportion of TEs in CS chromosomes')
    plt.legend()
    plt.show()

    # 3. retry
    reads_wide = reads_pct.pivot(index='chr', columns='type', values='percent')
    bottom = np.zeros(len(reads_wide))

    plt.figure(figsize=(14, 16))
    plt.style.use('seaborn-v0_8-deep')
    for column in reads_wide.columns:
        plt.bar(reads_wide.index, reads_wide[column], bottom=bottom, label=column, alpha=0.8)
        bottom += reads_wide[column]
    plt.title('Proportion of TEs in CS chromosomes')
    plt.xlabel('Chromosome')
    plt.ylabel('Proportion (%)')
    plt.legend()
    plt.show()

def plot_depth(chr):
    summ = pd.read_table("data/old_reads_depth/" + chr + ".mosdepth.summary.txt", sep="\t", header=0)
    # summ = pd.read_table("data/vu_reads_depth/" + chr + ".mosdepth.summary.txt", sep="\t", header=0)
    summ['type'] = summ['chrom'].str.split('#').str[1]
    tecode = pd.read_table(TEcode, sep=",", header=None)
    tecode.columns = ['cls', 'new_cls']
    for i in range(0, len(tecode)):
        summ.loc[summ['type'] == tecode['cls'][i], 'type'] = tecode['new_cls'][i]
    summ['ref_depth'] = summ['mean']
    summ['chr'] = chr
    summ.drop(['chrom', 'mean', 'min', 'max'], axis=1, inplace=True)
    # summ now: chr, type, ref_depth, length, bases
    return summ

def to_group(grouped):
    pct = grouped.agg(
        sum_depth=('ref_depth', 'sum'),
    )
    pct['percent'] = pct['sum_depth'] / pct.groupby('chr')['sum_depth'].transform('sum') * 100
    pct.reset_index(inplace=True)
    return pct

def nl_stats(chr):
    te = pd.read_table("data/ref_depth/stats." + chr + ".length.txt", sep='\t', header=None)
    te.columns = ['type', 'length']
    te['length'] = te['length'].astype(int)
    te = te[~te['type'].str.contains('X')]
    te_length = te.groupby('type').agg(
        sum_length=('length', 'sum'),
    )
    te_length['percent'] = te_length['sum_length'] / te_length['sum_length'].sum() * 100
    te_length['chr'] = chr
    te_length.reset_index(inplace=True)
    return te_length

if __name__ == '__main__':
    main()