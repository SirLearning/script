# stats of TE number and length
import sys
import pandas as pd

fai_name = 'data/tr2.fa.fai'
output_name = 'data/stats.length.old.txt'


def nl_stats(fai):
    # 1. TE number into dataframe
    te = pd.read_table(fai, header=None)
    te.columns = ['TE', 'length', 'start', 'line', 'width']
    te['Classification'] = te['TE'].str.split('_').str[0]

    te_length = pd.concat([te['Classification'], te['length']], axis=1)
    return te_length


# main
TE_length = nl_stats(fai_name)  # perform stats
# show sum
pd.set_option('display.max_rows', None)
TE_length.to_csv(output_name, sep='\t', header=False, index=False)