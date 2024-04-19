# stats of TE number and length
import sys
import pandas as pd

"""
Module to summarize the number and length of TEs in a genome.
"""

""": parameter
fai = 'lib.fa.fai'
cu_code = 'curatedLib.txt'
output_name = 'data/curated_lib/stats.length.txt'
"""


def mod_clib(fai_name, curated_name):   # fai to classified length
    # 1. read lib.fai
    fai = pd.read_table(fai_name, sep="\t", header=None)
    fai.columns = ['TE', 'length', 'start', 'line', 'width']
    cu_code = pd.read_table(curated_name, sep=":", header=None)
    cu_code.columns = ['Classification', 'TE', 'species']
    # 2. classify lib
    fai['Classification'] = ''
    for i in range(len(fai)):
        for j in range(len(cu_code)):
            if fai['TE'][i] == cu_code['TE'][j]:
                fai.loc[i, 'Classification'] = cu_code.loc[j, 'Classification']
                break
    to_stats = pd.concat([fai['Classification'], fai['length']], axis=1)
    return to_stats


def mod_trep(fai_name):     # fai to classified length
    # 1. read lib.fai
    fai = pd.read_table(fai_name, header=None)
    fai.columns = ['TE', 'length', 'start', 'line', 'width']
    # 2. classify lib
    fai['Classification'] = fai['TE'].str.split('_').str[0]
    to_stats = pd.concat([fai['Classification'], fai['length']], axis=1)
    return to_stats


def main():
    fai = sys.argv[1]
    cu_code = sys.argv[2]
    output_name = sys.argv[3]
    TE_length = mod_clib(fai, cu_code)  # fai to classified length
    # show sum
    pd.set_option('display.max_rows', None)
    TE_length.to_csv(output_name, sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()
