# stats of TE number and length
import sys
import pandas as pd

# allTE_name = sys.argv[1]
# fai_name = sys.argv[2]
# output_name = sys.argv[3]
# TE_length_name = sys.argv[4]

allTE_name = 'data/np.gff3'
fai_name = 'data/chr1A.fa.fai'
output_name = 'data/EDTA/stats.length.txt'
# TE_length_name = 'data/'
TEcode_name = 'data/TEcode'

TE_code = pd.read_table(TEcode_name, sep=',', header=None)
TE_code.columns = ['cls', 'new_cls']

# def nl_stats(allte, fai):
#     # 1. TE number into dataframe
#     genome_size = pd.read_table(fai, header=None)
#     genome_size.columns = ['chr', 'size', 'start', 'line', 'width']
#     size_sum = genome_size['size'].sum()
#
#     # 2. TE annotation into dataframe
#     # 2.1 all TE annotation
#     allte = pd.read_table(allte, sep='\t', header=None)
#     allte.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes', 'Classification']
#     allte['width'] = allte['end'] - allte['start'] + 1
#
#     # 1.4 get the TE number in each type
#     grouped = allte.groupby('Classification')
#     te_summ = grouped.agg(
#         count=('type', 'count'),
#         mean_size=('width', 'mean'),
#         size=('width', 'sum'),
#         percent=('width', lambda x: x.sum() / size_sum * 100)
#     )
#
#     te_length = pd.concat([allte['Classification'], allte['width']], axis=1)
#     return te_length, te_summ
#

# EDTA length
def nl_stats(allte_name, fai):
    # 1. TE number into dataframe
    genome_size = pd.read_table(fai, header=None)
    genome_size.columns = ['chr', 'size', 'start', 'line', 'width']
    size_sum = genome_size['size'].sum()

    # 2. TE annotation into dataframe
    # 2.1 all TE annotation
    allte = pd.read_table(allte_name, sep='\t', header=None)
    allte.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    allte['width'] = allte['end'] - allte['start'] + 1
    # 2.2 TE attributes
    attributes = allte['attributes'].str.split(';', expand=True)
    # attributes.columns = ['ID', 'Name', 'Classification', 'Sequence_ontology', 'Identity', 'Method', 'others1',
    #                             'others2', 'others3']
    attributes.columns = ['ID', 'Name', 'Classification', 'Sequence_ontology', 'Identity', 'Method', 'others1',
                          'others2']
    # 2.3 merge all TE annotation and TE attributes
    allte = pd.concat([allte, attributes], axis=1)
    # 1.1 show classification
    allte['Classification'] = allte['Classification'].str.split('=').str[1]

    # 1.3 Unspecified annotation
    allte.loc[allte['Classification'] == 'Unspecified', 'Classification'] = allte['Name'].str.split('=').str[1]
    allte.loc[:, 'Classification'] = allte['Classification'].str.split('_').str[0]
    for i in range(0, len(TE_code)):
        allte.loc[allte['Classification'] == TE_code['cls'][i], 'Classification'] = TE_code['new_cls'][i]

    # 1.4 get the TE number in each type
    grouped = allte.groupby('Classification')
    te_summ = grouped.agg(
        count=('type', 'count'),
        mean_size=('width', 'mean'),
        size=('width', 'sum'),
        percent=('width', lambda x: x.sum() / size_sum * 100)
    )

    te_length = pd.concat([allte['Classification'], allte['width']], axis=1)
    return te_length, te_summ


# main
[TE_length, TE_summ] = nl_stats(allTE_name, fai_name)  # perform stats
# show sum
pd.set_option('display.max_rows', None)
with open(output_name, 'w') as f:
    f.write(TE_summ.to_string(index_names=False))
TE_length.to_csv(output_name, sep='\t', header=False, index=False)