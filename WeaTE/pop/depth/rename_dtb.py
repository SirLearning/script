import sys

import pandas as pd

depth_name = sys.argv[1]
tecode_name = sys.argv[2]
output_name = sys.argv[3]

# anno_name = 'transposon/test.gff3'
# tecode_name = 'transposon/TEcode'
# output_name = 'transposon/length.txt'

depth = pd.read_table(depth_name, sep='\t', header=None)
depth.columns = ['TE', 'depth', 'prob']

# 1. rearrange classification
depth['Classification'] = depth['Classification'].str.split('=').str[1]
# 1.2 no Parent
index = ~depth['Name'].str.contains('Parent')
depth = depth[index]
# 1.3 Unspecified annotation
TEcode = pd.read_table(tecode_name, sep=",", header=None)
TEcode.columns = ['cls', 'new_cls']
depth.loc[depth['Classification'] == 'Unspecified', 'Classification'] = depth['Name'].str.split('=').str[1]
depth.loc[:, 'Classification'] = depth['Classification'].str.split('_').str[0]
depth['lib_correct'] = False
for i in range(0, len(TEcode)):
    depth.loc[depth['Classification'] == TEcode['cls'][i], 'Classification'] = TEcode['new_cls'][i]
    depth.loc[depth['Classification'] == TEcode['new_cls'][i], 'lib_correct'] = True
depth = depth[depth['lib_correct'] == True]

te_length = pd.concat([depth['Classification'], depth['length']], axis=1)
with open(output_name, 'w') as f:
    te_length.to_csv(f, sep='\t', index=False)
