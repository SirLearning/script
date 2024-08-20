# stats of TE number and length
import sys
import pandas as pd

"""
Module to statistically measure the annotation files.
"""

""":parameter
fai = '{seq}.fai'
anno = 'mod.anno.gff3'
output1 = 'stats.summ.txt'
output2 = 'stats.length.txt'
output3 = 'stats.overlap.gff3'
"""


def summary(fai_name, anno_name):
    # 1. read data
    fai = pd.read_table(fai_name, header=None)
    fai.columns = ['chr', 'size', 'start', 'line', 'width']
    size = fai['size'].sum()
    if isinstance(anno_name, str):
        anno = pd.read_table(anno_name, sep='\t', header=None, comment='#')
        anno.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes',
                        'width', 'classification']
    else:
        anno = anno_name
    # 2. get the TE number in each type
    grouped = anno.query(
        "`type` not in ['target_site_duplication', 'long_terminal_repeat']").groupby('classification')
    summ = grouped.agg(
        count=('type', 'count'),
        mean_size=('width', 'mean'),
        size=('width', 'sum'),
        percent=('width', lambda x: x.sum() / size * 100)
    )
    return summ


def length(anno_name):
    if isinstance(anno_name, str):
        anno = pd.read_table(anno_name, sep='\t', header=None)
        anno.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes',
                        'width', 'classification']
    else:
        anno = anno_name
    TE_length = pd.concat([anno['classification'], anno['width']], axis=1)
    return TE_length


def self_overlap(anno_name, mode):
    if isinstance(anno_name, str):
        anno = pd.read_table(anno_name, sep='\t', header=None)
        anno.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes',
                        'width', 'classification']
    else:
        anno = anno_name
    # check overlap
    anno['overlap'] = False
    end0 = 0
    classification0 = ''
    width0 = 0
    if mode == 1:  # mode 1: 0.8 mean length overlap with before in same classification
        for i in range(0, len(anno)):
            if (end0 - anno['start'][i]) > (width0 + anno['width'][i]) * 0.4 and anno['classification'][i] == classification0:
                anno.loc[i, 'overlap'] = True
            end0 = anno['end'][i]
            classification0 = anno['classification'][i]
            width0 = anno['width'][i]
    elif mode == 2:  # mode 2: 0.8 overlap with before
        for i in range(0, len(anno)):
            if (end0 - anno['start'][i]) > (width0 + anno['width'][i]) * 0.4:
                anno.loc[i, 'overlap'] = True
            end0 = anno['end'][i]
            width0 = anno['width'][i]
    elif mode == 3:  # mode 3: 0.8 overlap with before in same classification
        for i in range(0, len(anno)):
            if (end0 - anno['start'][i]) > anno['width'][i] * 0.8 and anno['classification'][i] == classification0:
                anno.loc[i, 'overlap'] = True
            end0 = anno['end'][i]
            classification0 = anno['classification'][i]
    elif mode == 4:  # mode 4: 0.8 overlap with longest before
        for i in range(0, len(anno)):
            if (end0 - anno['start'][i]) > anno['width'][i] * 0.8:
                anno.loc[i, 'overlap'] = True
            if anno['end'][i] > end0:
                end0 = anno['end'][i]
    return anno


def main():
    fai = sys.argv[1]
    anno = sys.argv[2]
    output1 = sys.argv[3]
    output2 = sys.argv[4]
    output3 = sys.argv[5]

    summ = summary(fai, anno)   # get summary
    pd.set_option('display.max_rows', None)
    with open(output1, 'w') as f:
        f.write(summ.to_string(index_names=False))

    TE_length = length(anno)  # get length
    TE_length.to_csv(output2, sep='\t', header=False,  index=False)
    # columns = ['Classification', 'length']

    TE_overlap = self_overlap(anno, 4)  # get overlap
    TE_overlap.to_csv(output3, sep='\t', header=False,  index=False)
    # columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes', 'width', 'classification', 'overlap']


if __name__ == '__main__':
    main()
