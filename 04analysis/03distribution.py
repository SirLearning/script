import sys
import pandas as pd

# allte_name = 'data/01chr1A.anno.gff3'
# TEcode_name = 'data/TEcode'
# output_name = 'data/chr1A.anno.mc.gff3'

allte_name = sys.argv[1]
TEcode_name = sys.argv[2]
output_name = sys.argv[3]

# 1. TE annotation into dataframe
allte = pd.read_table(allte_name, sep='\t', header=None)
allte.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
# rearrange the attributes
attributes = allte['attributes'].str.split(';', expand=True)
attributes.columns = ['ID', 'Name', 'Classification', 'Sequence_ontology', 'Identity', 'Method', 'others1', 'others2']
allte = pd.concat([allte, attributes], axis=1)

# 2. change the classification
allte['Classification'] = allte['Classification'].str.split('=').str[1]
# 2.2 no Parent
index = ~allte['Name'].str.contains('Parent')
allte = allte[index]
# 2.3 Unspecified annotation
allte.loc[allte['Classification'] == 'Unspecified', 'Classification'] = allte['Name'].str.split('=').str[1]
allte.loc[:, 'Classification'] = allte['Classification'].str.split('_').str[0]
TE_code = pd.read_table(TEcode_name, sep=',', header=None)
TE_code.columns = ['cls', 'new_cls']
for i in range(0, len(TE_code)):
    allte.loc[allte['Classification'] == TE_code['cls'][i], 'Classification'] = TE_code['new_cls'][i]

# # 3. merge the attributes
# allte = allte.fillna('')
# allte['attributes'] = allte['ID'] + ';' + allte['Name'] + ';' + 'Classification=' + allte['Classification'] + ';' + allte['Sequence_ontology'] + ';' + allte['Identity'] + ';' + allte['Method'] + ';' + allte['others1'] + ';' + allte['others2'] + ';' + allte['others3']
allte = allte.drop(['ID', 'Name', 'Sequence_ontology', 'Identity', 'Method', 'others1', 'others2'], axis=1)
# print(allte['attributes'])

# 4. output the gff3 file
allte.to_csv(output_name, sep='\t', header=False, index=False)