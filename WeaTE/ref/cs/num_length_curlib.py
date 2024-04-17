# stats of TE number and length
import sys
import pandas as pd

fai_name = 'data/abd.fa.fai'
curated_name = 'data/curatedLib.txt'
output_name = 'data/curated_lib/stats.length.txt'


def nl_stats(fai, curated):
    # 1. TE number into dataframe
    te = pd.read_table(fai, sep="\t", header=None)
    te.columns = ['TE', 'length', 'start', 'line', 'width']

    cu = pd.read_table(curated, sep=":", header=None)
    cu.columns = ['Classification', 'TE', 'species']

    for i in range(len(te)):
        for j in range(len(cu)):
            if te['TE'][i] == cu['TE'][j]:
                te.loc[i, 'TE'] = cu.loc[j, 'Classification']
                break
    te['Classification'] = te['TE'].str.split('_').str[0]

    te_length = pd.concat([te['Classification'], te['length']], axis=1)
    print(te_length)
    return te_length


# main
TE_length = nl_stats(fai_name, curated_name)  # perform stats
# show sum
pd.set_option('display.max_rows', None)
TE_length.to_csv(output_name, sep='\t', header=False, index=False)