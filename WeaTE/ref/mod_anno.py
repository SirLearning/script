import sys
import pandas as pd
import re

"""
Module to modify the EDTA TE annotation file
"""

""":parameter
anno = '{seq}.mod.EDTA.TEanno.gff3' # EDTA
anno = 'iwgsc_refseqv1.0_TransposableElements_2017Mar13.gff3' # CS
nn = 'N1.bed'
fai = '{seq}.fai'
output = 'mod.edta.gff3'
output = 'mod.cs.gff3'
output = 'mod.chrom.site.gff3'
output = 'mod.n.site.gff3'
"""


def edta(anno_name):
    # 1. read anno
    anno = pd.read_table(anno_name, sep='\t', header=None, comment='#')
    anno.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    anno['width'] = anno['end'] - anno['start'] + 1
    attributes = anno['attributes'].str.split(';', expand=True)
    if len(attributes.columns) == 9:
        attributes.columns = ['ID', 'Name', 'Classification', 'Sequence_ontology', 'Identity', 'Method',
                              'others1', 'others2', 'others3']
        attributes.drop(['others1', 'others2', 'others3'], axis=1, inplace=True)
    else:
        attributes.columns = ['ID', 'Name', 'Classification', 'Sequence_ontology', 'Identity', 'Method',
                              'others1', 'others2']
        attributes.drop(['others1', 'others2'], axis=1, inplace=True)
    anno = pd.concat([anno, attributes], axis=1)
    # 2. delete lines with 'Parent'
    anno = no_parent(anno)
    # 3. change the classification according to TEcode
    anno = te_code(anno)
    anno = anno.drop(['ID', 'Name', 'Sequence_ontology', 'Identity', 'Method'], axis=1)
    anno.reset_index(drop=True, inplace=True)
    return anno


def cs(anno_name):
    # 1. read anno
    anno = pd.read_table(anno_name, sep='\t', header=None, comment='#')
    anno.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    attributes = anno['attributes'].str.split(';', expand=True)
    attributes.columns = ['ID', 'Name', 'Ontology_term', 'compo', 'soloLTR', 'status']
    anno = pd.concat([anno, attributes], axis=1)
    # 2. most identity classification
    anno['Classification'] = anno['compo'].str.split('=').str[1]
    for j in range(0, len(anno['Classification'])):
        line = anno.loc[j, 'Classification']
        line = re.split('\s', line)
        max_value = 0
        classification = ''
        # print(line)
        for i in range(0, len(line) - 1):
            if re.match('(\w+_\w+|no_match)', line[i]):
                if max_value < float(line[i + 1]):
                    max_value = float(line[i + 1])
                    classification = line[i]
        anno.loc[j, 'Classification'] = classification
    anno['Classification'] = anno['Classification'].str.replace('no_match', 'XXX')
    anno['Classification'] = anno['Classification'].str.split('_').str[0]
    # 3. merge the attributes
    anno = anno.drop(['ID', 'Name', 'Ontology_term', 'compo', 'soloLTR', 'status'], axis=1)
    empty_cols = anno.columns[anno.isnull().all()]
    anno = anno.drop(empty_cols, axis=1)
    anno.reset_index(drop=True, inplace=True)
    return anno


def cs_v2(anno_name):
    # 1. read anno
    anno = pd.read_table(anno_name, sep='\t', header=None, comment='#')
    anno.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    anno = no_parent(anno)
    attributes = anno['attributes'].str.split(';', expand=True)
    attributes.columns = ['ID', 'compo', 'copie', 'post', 'range', 'status']
    anno = pd.concat([anno, attributes], axis=1)
    # 2. most identity classification
    anno['Classification'] = anno['ID'].str.split('_').str[1]
    anno['Classification'] = anno['Classification'].str.replace('no', 'XXX')
    # 3. merge the attributes
    anno = anno.drop(['ID', 'compo', 'copie', 'post', 'range', 'status'], axis=1)
    empty_cols = anno.columns[anno.isnull().all()]
    anno = anno.drop(empty_cols, axis=1)
    anno.reset_index(drop=True, inplace=True)
    return anno


def no_parent(anno):
    index = ~anno['attributes'].str.contains('Parent')
    anno = anno[index]
    return anno


def te_code(anno):
    data = {
        'cls': ['DNA/DTA', 'DNA/DTC', 'DNA/DTH', 'DNA/DTM', 'DNA/DTT', 'DNA/Helitron', 'LINE/unknown', 'LTR/Copia',
                'LTR/Gypsy', 'LTR/unknown', 'MITE/DTA', 'MITE/DTC', 'MITE/DTH', 'MITE/DTM', 'MITE/DTT', 'Unknown'],
        'new_cls': ['DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DHH', 'RIX', 'RLC', 'RLG', 'RLX', 'DTA', 'DTC', 'DTH', 'DTM',
                    'DTT', 'XXX']
    }
    TEcode = pd.DataFrame(data)
    anno['Classification'] = anno['Classification'].str.split('=').str[1]
    anno.loc[anno['Classification'] == 'Unspecified', 'Classification'] = anno['Name'].str.split('=').str[1]
    anno.loc[:, 'Classification'] = anno['Classification'].str.split('_').str[0]
    for i in range(0, len(TEcode)):
        anno.loc[anno['Classification'] == TEcode['cls'][i], 'Classification'] = TEcode['new_cls'][i]
    return anno


def chrom_site(fai_name, anno_name, output_name):
    fai = pd.read_table(fai_name, sep='\t', header=None)
    fai.columns = ['seq', 'length', 'start', 'lb', 'lw']
    with open(anno_name, 'r') as anno_file, open(output_name, 'w') as output_file:
        for line in anno_file:
            element = line.split('\t')
            chrom = element[0]
            start = int(element[3])
            end = int(element[4])
            if chrom == str(fai['seq'][1]):
                element[3] = str(start + int(fai['length'][0]))
                element[4] = str(end + int(fai['length'][0]))
            output_file.write('\t'.join(element))


def n_site(nn_name, anno_name, output_name):
    nn = pd.read_table(nn_name, sep='\t', header=None)
    nn.columns = ['chrom', 'start', 'end']
    with open(anno_name, 'r') as anno_file, open(output_name, 'w') as output_file:
        for line in anno_file:
            element = line.split('\t')
            chrom = element[0]
            start = int(element[3])
            end = int(element[4])
            for i in range(0, len(nn)):
                if chrom == 'chr1A_' + str(i + 2):
                    element[3] = str(start + int(nn['end'][i]))
                    element[4] = str(end + int(nn['end'][i]))
                    output_file.write('\t'.join(element))
                    break
            else:
                output_file.write(line)


def seq_region(anno_name, region):
    anno = edta(anno_name)
    anno = anno.drop(['Classification'], axis=1)
    region = region.split(':')
    anno = anno.loc[(int(anno['start']) >= int(region[0])) & (int(anno['end']) <= int(region[1]))]
    return anno


def mod_type(mod, type, anno):
    if mod == 'edta':
        if re.match(r'^\d+:\d:$', type):
            anno = seq_region(anno, type)
        elif type == 'TE':
            anno = edta(anno)
    elif mod == 'cs':
        anno = cs(anno)
    elif mod == 'cs_v2':
        anno = cs_v2(anno)
    return anno


def main():
    anno = sys.argv[1]
    output = sys.argv[2]
    mod = sys.argv[3]
    type = sys.argv[4]
    mod_type(mod, type, anno).to_csv(output, sep='\t', header=False, index=False)
    # anno.to_csv(output, sep='\t', header=False, index=False)
    # columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes', 'width', 'classification']


if __name__ == '__main__':
    main()