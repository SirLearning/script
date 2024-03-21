# stats of TE number and length
import sys
import pandas as pd

allTE_name = sys.argv[1]
fai_name = sys.argv[2]
output_name = sys.argv[3]

# allTE_name = 'data/test.ln.gff3'
# fai_name = 'data/chr1A.fa.fai'
# output_name = 'data/test.gff3'


def nl_stats(allte, fai):
    # 1. TE number into dataframe
    genome_size = pd.read_table(fai, header=None)
    genome_size.columns = ['chr', 'size', 'start', 'line', 'width']
    size_sum = genome_size['size'].sum()

    # 2. TE annotation into dataframe
    # 2.1 all TE annotation
    allte = pd.read_table(allte, sep='\t', header=None)
    allte.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes', 'Classification']
    allte['width'] = allte['end'] - allte['start'] + 1

    # 1.4 get the TE number in each type
    grouped = allte.groupby('Classification')
    te_summ = grouped.agg(
        count=('type', 'count'),
        mean_size=('width', 'mean'),
        size=('width', 'sum'),
        percent=('width', lambda x: x.sum() / size_sum * 100)
    )
    return te_summ


# main
TE_summ = nl_stats(allTE_name, fai_name)  # perform stats
# show sum
pd.set_option('display.max_rows', None)
with open(output_name, 'w') as f:
    f.write(TE_summ.to_string(index_names=False))
print(TE_summ)
print(TE_summ.index)
