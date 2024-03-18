# stats of TE number and length
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# input file 1, 2: chr1A.fa.mod.EDTA.TEanno.gff3, 01chr1A.fa.fai
# allTE_name = 'data/01chr1A.anno.gff3'
# fai_name = 'data/01chr1A.fa.fai'
# output_name = 'data/01chr1A.anno.stats'
allTE_name = sys.argv[1]
fai_name = sys.argv[2]
output_name = sys.argv[3]
TE_code = pd.read_table('data/TEcode', sep=',', header=None)
TE_code.columns = ['cls', 'new_cls']

def nl_stats(fai, allte):
    # 1. TE number into dataframe
    genome_size = pd.read_table(fai, header=None)
    genome_size.columns = ['chr', 'size', 'start', 'line', 'width']
    chr1a_size = genome_size['size'].sum()

    # 2. TE annotation into dataframe
    # 2.1 all TE annotation
    chr1a_allte = pd.read_table(allte, sep='\t', header=None)
    chr1a_allte.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    chr1a_allte['width'] = chr1a_allte['end'] - chr1a_allte['start'] + 1
    # 2.2 TE attributes
    chr1A_attributes = chr1a_allte['attributes'].str.split(';', expand=True)
    chr1A_attributes.columns = ['ID', 'Name', 'Classification', 'Sequence_ontology', 'Identity', 'Method', 'others1',
                                'others2', 'others3']
    # 2.3 merge all TE annotation and TE attributes
    chr1a_allte = pd.concat([chr1a_allte, chr1A_attributes], axis=1)
    return chr1a_allte, chr1a_size


def sum_allte(allte, size):
    # 1.1 show classification
    allte['Classification'] = allte['Classification'].str.split('=').str[1]

    # 1.2 no Parent
    # chr1A_allTE['Parent'] = chr1A_allTE['Name'].str.split('=').str[0]
    index = ~allte['Name'].str.contains('Parent')
    allte = allte[index]

    # 1.3 Unspecified annotation
    allte.loc[allte['Classification'] == 'Unspecified', 'Classification'] = chr1A_allTE['Name'].str.split('=').str[1]
    allte.loc[:, 'Classification'] = allte['Classification'].str.split('_').str[0]
    for i in range(0, len(TE_code)):
        allte.loc[allte['Classification'] == TE_code['cls'][i], 'Classification'] = TE_code['new_cls'][i]
    # Create a dictionary from TE_code DataFrame
    # cls_dict = TE_code.set_index('cls')['new_cls'].to_dict()
    # # Use map function to replace 'Classification' column values
    # allte['Classification'] = allte['Classification'].map(cls_dict)

    # 1.3 get the TE number in each type
    grouped = allte.query(
        "`type` not in ['target_site_duplication', 'long_terminal_repeat']").groupby('Classification')
    allte_summ = grouped.agg(
        count=('type', 'count'),
        mean_size=('width', 'mean'),
        size=('width', 'sum'),
        percent=('width', lambda x: x.sum() / size * 100)
    )
    # 删除名为'Classification'的行
    # allte_summ = allte_summ.drop('Classification')
    return allte_summ


# main
[chr1A_allTE, sum_size] = nl_stats(fai_name, allTE_name)  # perform stats
chr1A_allTE_summ = sum_allte(chr1A_allTE, sum_size)  # get sum
# show sum
pd.set_option('display.max_rows', None)
with open(output_name, 'w') as f:
    f.write(chr1A_allTE_summ.to_string(index_names=False))
print(chr1A_allTE_summ)
print(chr1A_allTE_summ.index)