import sys
import pandas as pd

"""
Module to modify the EDTA TE annotation file
"""

""":parameter
anno = '{seq}.mod.EDTA.TEanno.gff3'
fai = '{seq}.fai'
output = 'mod.edta.gff3'
"""


def mod_anno(anno):
    # 1. read anno
    anno = pd.read_table(anno, sep='\t', header=None)
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
    anno.drop(['ID', 'Name', 'Sequence_ontology', 'Identity', 'Method'], axis=1, inplace=True)
    anno.reset_index(drop=True, inplace=True)
    return anno


def no_parent(anno):
    index = ~anno['Name'].str.contains('Parent')
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
    anno.loc[anno['Classification'] == 'Unspecified', 'Classification'] = anno['Name'].str.split('=').str[1]
    anno.loc[:, 'Classification'] = anno['Classification'].str.split('_').str[0]
    for i in range(0, len(TEcode)):
        anno.loc[anno['Classification'] == TEcode['cls'][i], 'Classification'] = TEcode['new_cls'][i]
    return anno


def main():
    anno = sys.argv[1]
    output = sys.argv[2]
    anno = mod_anno(anno)
    anno.to_csv(output, sep='\t', header=False, index=False)
    # columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes', 'width', 'classification']


if __name__ == '__main__':
    main()