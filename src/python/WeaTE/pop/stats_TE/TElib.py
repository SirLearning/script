import pandas as pd

fai_name = 'transposon/CS.1A.TElib.fa.fai'
chr_name = 'transposon/chr1A.fa.fai'
output_name = 'transposon/TElib.stats'
TEcode_name = 'transposon/TEcode'

chr = pd.read_table(chr_name, sep='\t', header=None)
chr.columns = ['chr', 'size', 'start', 'line', 'width']

TE_code = pd.read_table(TEcode_name, sep=',', header=None)
TE_code.columns = ['cls', 'new_cls']

fai = pd.read_table(fai_name, sep='\t', header=None)
fai.columns = ['chr', 'size', 'start', 'line', 'width']
fai['classification'] = fai['chr'].str.split('#').str[1]

for i in range(0, len(TE_code)):
    fai.loc[fai['classification'] == TE_code['cls'][i], 'classification'] = TE_code['new_cls'][i]

grouped = fai.groupby('classification')
fai_summ = grouped.agg(
    count=('chr', 'count'),
    size=('size', 'sum'),
    mean_size=('size', 'mean'),
    percent=('size', lambda x: x.sum() / chr['size'].sum() * 100)
)
with open(output_name, 'w') as f:
    f.write(fai_summ.to_string(index_names=False))
print(fai_summ)
print(fai_summ['count'].sum())