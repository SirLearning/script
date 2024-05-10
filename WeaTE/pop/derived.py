import sys

import matplotlib.pyplot as plt
import pandas as pd
import upsetplot

"""
Module to determine the derived TEs in species.
"""

""": parameter
fai = 'lib.fa.fai'
output_name = 'lib.bed'
"""


def fai_to_bed(fai_name):
    fai = pd.read_table(fai_name, header=None)
    fai.columns = ['chr', 'size', 'start', 'line', 'width']
    fai['start'] = 1
    fai['end'] = fai['size']
    bed = fai[['chr', 'start', 'end']]
    return bed


def to_plot(spc):
    cor_A = {
        'CS_A': [True, True, True, True, True, True, True, True],
        'barley': [False, True, False, False, False, False, False, False],
        'rye': [False, False, True, False, False, False, False, False],
        'CS_D': [False, False, False, True, False, False, False, False],
        'CS_B': [False, False, False, False, True, False, False, False],
        'durum': [False, False, False, False, False, True, False, False],
        'wild_emmer': [False, False, False, False, False, False, True, False],
        'Lancer': [False, False, False, False, False, False, False, True],
        'count': [3996, 3721, 3298, 3105, 3080, 2269, 2015, 1966]
    }
    cor_B = {
        'CS_B': [True, True, True, True, True, True, True, True],
        'barley': [False, True, False, False, False, False, False, False],
        'rye': [False, False, True, False, False, False, False, False],
        'CS_D': [False, False, False, True, False, False, False, False],
        'CS_A': [False, False, False, False, True, False, False, False],
        'durum': [False, False, False, False, False, True, False, False],
        'wild_emmer': [False, False, False, False, False, False, True, False],
        'Lancer': [False, False, False, False, False, False, False, True],
        'count': [4638, 4328, 3908, 3718, 3673, 2604, 2441, 2212]
    }
    cor_D = {
        'CS_D': [True, True, True, True, True, True, True],
        'barley': [False, True, False, False, False, False, False],
        'rye': [False, False, True, False, False, False, False],
        'CS_B': [False, False, False, True, False, False, False],
        'CS_A': [False, False, False, False, True, False, False],
        'tauchii': [False, False, False, False, False, True, False],
        'Lancer': [False, False, False, False, False, False, True],
        'count': [3623, 3282, 2952, 2812, 2776, 1997, 1455]
    }

    se_A = cor_to_se(cor_A)
    se_B = cor_to_se(cor_B)
    se_D = cor_to_se(cor_D)
    upsetplot.plot(se_A, show_counts=True, facecolor="#348ABD")
    plt.show()
    upsetplot.plot(se_B, show_counts=True, facecolor="#7A68A6")
    plt.show()
    upsetplot.plot(se_D, show_counts=True, facecolor="#467821")
    plt.show()


def cor_to_se(cor):
    df = pd.DataFrame(cor)
    se = df.set_index([col for col in df.columns if col != 'count']).squeeze()
    se[:] = df['count']
    return se


def main():
    # bed = fai_to_bed(sys.argv[1])
    # bed.to_csv(sys.argv[2], sep='\t', header=False, index=False)
    to_plot('cs')


if __name__ == '__main__':
    main()
