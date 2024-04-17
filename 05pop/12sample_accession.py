import pandas as pd

sample_name = 'data/sample.txt'
out_name = 'data/accession.txt'

sample = pd.read_table(sample_name, sep='::', engine='python', header=None)
sample.columns = ['accession0', 'name0']
sample['accession'] = sample['accession0'].str.split(',').str[0]
sample['name'] = sample['name0'].str.split(',').str[0]
sample.drop(['accession0', 'name0'], axis=1, inplace=True)
sample.to_csv(out_name, sep='\t', index=False)