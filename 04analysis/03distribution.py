import sys
import pandas as pd

# allte_name = 'data/01chr1A.anno.gff3'
# TEcode_name = 'data/TEcode'
# output_name = 'data/chr1A.anno.mc.gff3'

all_TE_name = sys.argv[1]
TEcode_name = sys.argv[2]
output_name = sys.argv[3]

# 1. TE annotation into dataframe
all_TE = pd.read_table(all_TE_name, sep='\t', header=None)
all_TE.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
# rearrange the attributes
attributes = all_TE['attributes'].str.split(';', expand=True)
attributes.columns = ['ID', 'Name', 'Classification', 'Sequence_ontology', 'Identity', 'Method', 'others1', 'others2', 'others3']
all_TE = pd.concat([all_TE, attributes], axis=1)

# 2. change the classification
all_TE['Classification'] = all_TE['Classification'].str.split('=').str[1]
# 2.2 no Parent
index = ~all_TE['Name'].str.contains('Parent')
all_TE = all_TE[index]
# 2.3 Unspecified annotation
all_TE.loc[all_TE['Classification'] == 'Unspecified', 'Classification'] = all_TE['Name'].str.split('=').str[1]
all_TE.loc[:, 'Classification'] = all_TE['Classification'].str.split('_').str[0]
TE_code = pd.read_table(TEcode_name, sep=',', header=None)
TE_code.columns = ['cls', 'new_cls']
for i in range(0, len(TE_code)):
    all_TE.loc[all_TE['Classification'] == TE_code['cls'][i], 'Classification'] = TE_code['new_cls'][i]

# # 3. merge the attributes
# allte = allte.fillna('')
# allte['attributes'] = allte['ID'] + ';' + allte['Name'] + ';' + 'Classification=' + allte['Classification'] + ';' + allte['Sequence_ontology'] + ';' + allte['Identity'] + ';' + allte['Method'] + ';' + allte['others1'] + ';' + allte['others2'] + ';' + allte['others3']
all_TE = all_TE.drop(['ID', 'Name', 'Sequence_ontology', 'Identity', 'Method', 'others1', 'others2', 'others3'], axis=1)
# print(allte['attributes'])

# 4. output the gff3 file
empty_cols = all_TE.columns[all_TE.isnull().all()]
all_TE = all_TE.drop(empty_cols, axis=1)
all_TE.to_csv(output_name, sep='\t', header=False, index=False)