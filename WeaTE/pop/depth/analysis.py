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
    ref_pct = pd.concat([nl_stats('1A'), nl_stats('1B'), nl_stats('1D')])

    depth = depth.sort_values(by='type')
    reads_pct = reads_pct.sort_values(by='type')
    ref_pct = ref_pct.sort_values(by='type')

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

    sns.barplot(x='type', y='percent', hue='chr', data=reads_pct, palette='Set3', ax=ax, legend=['reads_1A', 'reads_1B', 'reads_1D'])
    sns.barplot(x='type', y='percent', hue='chr', data=ref_pct, ax=ax, alpha=0.5, legend=['ref_1A', 'ref_1B', 'ref_1D'])
    plt.title('Proportion of TEs in CS chromosomes')
    plt.legend()
    plt.show()

    # 3. MSE of TEs
    fig, ax = plt.subplots()
    ax.figure.set_size_inches(18, 9)
    merged = calculate_cv(reads_pct, ref_pct)
    sns.barplot(x='type', y='coefficient of variation', hue='chr', data=merged, palette='Set3', ax=ax)
    plt.title('Deviation of TEs proportion')
    plt.ylabel('Deviation (%)')
    plt.legend()
    plt.show()

    # 4. accumulated proportion
    reads_wide = reads_pct.pivot(index='chr', columns='type', values='percent')
    bottom = np.zeros(len(reads_wide))

    plt.figure(figsize=(18, 9))
    plt.style.use('seaborn-v0_8-deep')
    al_tag = 0
    for column in reads_wide.columns:
        plt.barh(reads_wide.index, reads_wide[column], left=bottom, label=column, alpha=(0.8 - al_tag), height=0.6)
        bottom += reads_wide[column]
        al_tag += 0.05
    plt.title('Proportion of TEs in CS chromosomes')
    plt.ylabel('Chromosome')
    plt.xlabel('Proportion (%)')
    plt.legend()
    plt.show()

    # 5. reads_pct of TEs after threshold
    threshold = pd.read_table('data/vu_reads_depth/threshold.txt', header=None)
    threshold.columns = ['chrom', 'chr']
    th = depth
    for i in range(len(threshold)):
        th = th[~((th['chrom'] == threshold.loc[i, 'chrom']) & (th['chr'] == threshold.loc[i, 'chr']))]
    th_groupe = th.groupby(['type', 'chr'])
    th_reads = to_group(th_groupe)

    fig, ax = plt.subplots()
    ax.figure.set_size_inches(18, 9)
    sns.barplot(x='type', y='percent', hue='chr', data=th_reads, palette='Set3', ax=ax,
                legend=['reads_1A', 'reads_1B', 'reads_1D'])
    sns.barplot(x='type', y='percent', hue='chr', data=ref_pct, ax=ax, alpha=0.5, legend=['ref_1A', 'ref_1B', 'ref_1D'])
    plt.title('Proportion of TEs in CS chromosomes')
    plt.legend()
    plt.show()

    fig, ax = plt.subplots()
    ax.figure.set_size_inches(18, 9)
    th_merged = calculate_cv(th_reads, ref_pct)
    sns.barplot(x='type', y='coefficient of variation', hue='chr', data=th_merged, palette='Set3', ax=ax)
    plt.title('Deviation of TEs proportion after threshold')
    plt.ylabel('Deviation (%)')
    plt.legend()
    plt.show()


def calculate_cv(reads_pct, ref_pct):
    merged = pd.merge(reads_pct, ref_pct, on=['type', 'chr'], suffixes=('_reads', '_ref'))
    merged['coefficient of variation'] = (merged['percent_reads'] - merged['percent_ref']) / merged['percent_ref'] * 100
    return merged


def plot_depth(chr):
    summ = pd.read_table("data/vu_reads_depth/" + chr + ".mosdepth.summary.txt", sep="\t", header=0)
    # summ = pd.read_table("data/vu_reads_depth/" + chr + ".mosdepth.summary.txt", sep="\t", header=0)
    summ['type'] = summ['chrom'].str.split('#').str[1]
    tecode = pd.read_table(TEcode, sep=",", header=None)
    tecode.columns = ['cls', 'new_cls']
    for i in range(0, len(tecode)):
        summ.loc[summ['type'] == tecode['cls'][i], 'type'] = tecode['new_cls'][i]
    summ['ref_depth'] = summ['mean']
    summ['chr'] = chr
    summ.drop(['mean', 'min', 'max'], axis=1, inplace=True)
    # summ now: chr, type, ref_depth, length, bases
    return summ

def to_group(grouped):
    pct = grouped.agg(
        sum_depth=('ref_depth', 'sum'),
    )
    pct['percent'] = pct['sum_depth'] / pct.groupby('chr')['sum_depth'].transform('sum') * 100
    pct.drop('sum_depth', axis=1, inplace=True)
    pct.reset_index(inplace=True)
    return pct

def nl_stats(chr):
    te = pd.read_table("data/ref_depth/stats.chr" + chr + ".length.txt", sep='\t', header=None)
    te.columns = ['type', 'length']
    te['length'] = te['length'].astype(int)
    te = te[~te['type'].str.contains('X')]
    te_length = te.groupby('type').agg(
        sum_length=('length', 'sum'),
    )
    te_length['percent'] = te_length['sum_length'] / te_length['sum_length'].sum() * 100
    te_length['chr'] = chr
    te_length.drop('sum_length', axis=1, inplace=True)
    te_length.reset_index(inplace=True)
    return te_length

if __name__ == '__main__':
    main()