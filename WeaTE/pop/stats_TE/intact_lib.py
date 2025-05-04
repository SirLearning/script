# stats of TE number and length
import sys
import pandas as pd

allTE_name = 'transposon/np.gff3'
fai_name = 'transposon/chr1A.fa.fai'
output_name = 'transposon/np.stats'
TEcode_name = 'transposon/TEcode'

TE_code = pd.read_table(TEcode_name, sep=',', header=None)
TE_code.columns = ['cls', 'new_cls']


def nl_stats(allte, fai):
    # 1. TE number into dataframe
    genome_size = pd.read_table(fai, header=None)
    genome_size.columns = ['chr', 'size', 'start', 'line', 'width']
    size_sum = genome_size['size'].sum()

    # 2. TE annotation into dataframe
    # 2.1 all TE annotation
    allte = pd.read_table(allte, sep='\t', header=None)
    allte.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    allte['width'] = allte['end'] - allte['start'] + 1
    # 2.2 TE attributes
    attributes = allte['attributes'].str.split(';', expand=True)
    attributes.columns = ['ID', 'Name', 'Classification', 'Sequence_ontology', 'Identity', 'Method', 'others1',
                                'others2']
    # 2.3 merge all TE annotation and TE attributes
    allte = pd.concat([allte, attributes], axis=1)
    return allte, size_sum


def sum_allte(te, size):
    # 1.1 show classification
    te['Classification'] = te['Classification'].str.split('=').str[1]

    # 1.2 no Parent
    # chr1A_allTE['Parent'] = chr1A_allTE['Name'].str.split('=').str[0]
    index = ~te['Name'].str.contains('Parent')
    te = te[index]

    # 1.3 Unspecified annotation
    te.loc[te['Classification'] == 'Unspecified', 'Classification'] = all_TE['Name'].str.split('=').str[1]
    te.loc[:, 'Classification'] = te['Classification'].str.split('_').str[0]
    for i in range(0, len(TE_code)):
        te.loc[te['Classification'] == TE_code['cls'][i], 'Classification'] = TE_code['new_cls'][i]

    # 1.4 get the TE number in each type
    grouped = te.query(
        "`type` not in ['target_site_duplication', 'long_terminal_repeat']").groupby('Classification')
    te_summ = grouped.agg(
        count=('type', 'count'),
        mean_size=('width', 'mean'),
        size=('width', 'sum'),
        percent=('width', lambda x: x.sum() / size * 100)
    )
    return te_summ

def main():
    # main
    [all_TE, sum_size] = nl_stats(allTE_name, fai_name)  # perform stats
    TE_summ = sum_allte(all_TE, sum_size)  # get sum
    # show sum
    pd.set_option('display.max_rows', None)
    with open(output_name, 'w') as f:
        f.write(TE_summ.to_string(index_names=False))
    print(TE_summ)
    print(TE_summ.index)
