import sys

import pandas as pd

anno_name = sys.argv[1]
tecode_name = sys.argv[2]
output_name = sys.argv[3]

anno = pd.read_table(anno_name, sep='\t', header=None)
anno.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
attributes = anno['attributes'].str.split(';', expand=True)
attributes.columns = ['ID', 'Name', 'Classification', 'Sequence_ontology', 'Identity', 'Method', 'others1',
                      'others2', 'others3']
anno = pd.concat([anno, attributes], axis=1)
anno['length'] = anno['end'] - anno['start'] + 1

# 1. rearrange classification
anno['Classification'] = anno['Classification'].str.split('=').str[1]
# 1.2 no Parent
index = ~anno['Name'].str.contains('Parent')
anno = anno[index]
# 1.3 Unspecified annotation
TEcode = pd.read_table(tecode_name, sep=",", header=None)
TEcode.columns = ['cls', 'new_cls']
anno.loc[anno['Classification'] == 'Unspecified', 'Classification'] = anno['Name'].str.split('=').str[1]
anno.loc[:, 'Classification'] = anno['Classification'].str.split('_').str[0]
anno['lib_correct'] = 2
for i in range(0, len(TEcode)):
    anno.loc[anno['Classification'] == TEcode['cls'][i], 'Classification'] = TEcode['new_cls'][i]
    anno.loc[anno['Classification'] == TEcode['new_cls'][i], 'lib_correct'] = 4
anno.drop(anno['lib_correct'] == 2, axis=0)

te_length = pd.concat([anno['Classification'], anno['length']], axis=1)
with open(output_name, 'w') as f:
    te_length.to_csv(f, sep='\t', index=False)
