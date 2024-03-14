# stats of TE number and length
import sys
import pandas as pd

# 1. TE number in each genome
# 1.1 get the genome annotation
allTE_name = sys.argv[1]
fai_name = sys.argv[2]

genome_size = pd.read_table(fai_name, header=None)
genome_size.columns = ['chr', 'size', 'start', 'line', 'width']
sum_size = genome_size['size'].sum()

chr1A_allTE = pd.read_table(allTE_name, sep='\t', header=None)
chr1A_allTE.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
chr1A_allTE['width'] = chr1A_allTE['end'] - chr1A_allTE['start'] + 1
chr1A_attributes = chr1A_allTE['attributes'].str.split(';', expand=True)
print(chr1A_attributes.shape[1])
chr1A_attributes.columns = ['ID', 'Name', 'Classification', 'Sequence_ontology', 'Identity', 'Method', 'others', 'others', 'others']
chr1A_allTE = pd.concat([chr1A_allTE, chr1A_attributes], axis=1)
chr1A_classification = chr1A_allTE['Classification'].str.split('=', expand=True)

# 1.2 get the TE number in each type
grouped = chr1A_allTE.query("`type` not in ['repeat_region', 'target_site_duplication', 'long_terminal_repeat']").groupby('Classification')
chr1A_allTE_summ = grouped.agg(
    count=('type', 'count'),
    mean_size=('width', 'mean'),
    size=('width', 'sum'),
    percent=('width', lambda x: x.sum() / sum_size * 100)
)
pd.set_option('display.max_rows', None)
print(chr1A_allTE_summ)