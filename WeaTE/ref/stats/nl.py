# stats of TE number and length
import sys
import pandas as pd

"""
Module to summarize the number and length of TEs in a genome.
"""

""":parameter
fai = '{seq}.fai'
anno = 'stats.edta.gff3'
output = 'stats.nl.txt'
"""


def stats_anno(fai_name, anno_name):
    # 1. read data
    fai = pd.read_table(fai_name, header=None)
    fai.columns = ['chr', 'size', 'start', 'line', 'width']
    size = fai['size'].sum()
    anno = pd.read_table(anno_name, sep='\t', header=None)
    # 2. get the TE number in each type
    grouped = anno.query(
        "`type` not in ['target_site_duplication', 'long_terminal_repeat']").groupby('Classification')
    summ = grouped.agg(
        count=('type', 'count'),
        mean_size=('width', 'mean'),
        size=('width', 'sum'),
        percent=('width', lambda x: x.sum() / size * 100)
    )
    return summ


def main():
    fai = sys.argv[1]
    anno = sys.argv[2]
    output = sys.argv[3]
    summ = stats_anno(fai, anno)  # get sum
    # show sum
    pd.set_option('display.max_rows', None)
    with open(output, 'w') as f:
        f.write(summ.to_string(index_names=False))


if __name__ == '__main__':
    main()
