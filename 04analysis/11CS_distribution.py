import sys
import pandas as pd
import re

# all_TE_name = 'data/test.cs.anno.gff3'
# output_name = 'data/test.gff3'

all_TE_name = sys.argv[1]
output_name = sys.argv[2]

# 1. TE annotation into dataframe
all_TE = pd.read_table(all_TE_name, sep='\t', header=None)
all_TE.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
# rearrange the attributes
attributes = all_TE['attributes'].str.split(';', expand=True)
attributes.columns = ['ID', 'Name', 'Ontology_term', 'compo', 'soloLTR', 'status']
all_TE = pd.concat([all_TE, attributes], axis=1)

# 2. change the classification
all_TE['Classification'] = all_TE['compo'].str.split('=').str[1]
for j in range(0, len(all_TE['Classification'])):
    line = all_TE.loc[j, 'Classification']
    line = re.split('\s', line)
    max_value = 0
    classification = ''
    # print(line)
    for i in range(0, len(line) -1):
        if re.match('(\w+_\w+|no_match)', line[i]):
            if max_value < float(line[i + 1]):
                max_value = float(line[i + 1])
                classification = line[i]
        else:

    all_TE.loc[j, 'Classification'] = classification
all_TE['Classification'] = all_TE['Classification'].str.replace('no_match', 'XXX')
all_TE['Classification'] = all_TE['Classification'].str.split('_').str[0]
# 3. merge the attributes
all_TE = all_TE.drop(['ID', 'Name', 'Ontology_term', 'compo', 'soloLTR', 'status'], axis=1)
# print(allte['attributes'])

# 4. output the gff3 file
empty_cols = all_TE.columns[all_TE.isnull().all()]
all_TE = all_TE.drop(empty_cols, axis=1)
all_TE.to_csv(output_name, sep='\t', header=False, index=False)
