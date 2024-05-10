import sys
import seaborn as sns
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


def upset_plot(spc):
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


def te_code(sample):
    data = {
        'cls': ['DNA/DTA', 'MITE/DTA', 'DNA/DTC', 'MITE/DTC', 'DNA/DTH', 'MITE/DTH', 'DNA/DTM', 'MITE/DTM',
                'DNA/DTT', 'MITE/DTT', 'TIR/Tc1_Mariner', 'TIR/Tc1', 'DNA/Helitron', 'LTR/Copia', 'LTR/Gypsy',
                'LTR/unknown', 'LINE/unknown', 'Unknown', 'pararetrovirus'],
        'new_cls': ['DTA', 'DTA', 'DTC', 'DTC', 'DTH', 'DTH', 'DTM', 'DTM',
                    'DTT', 'DTT', 'DTT', 'DTT', 'DHH', 'RLC', 'RLG',
                    'RLX', 'RIX', 'XXX', 'RLR']
    }
    TEcode = pd.DataFrame(data)
    for i in range(0, len(TEcode)):
        sample.loc[sample['Classification'] == TEcode['cls'][i], 'Classification'] = TEcode['new_cls'][i]
    return sample


def choose_sample(spc):
    samples = []
    if spc == 'tauchii':
        samples = ['CRR072405', 'CRR072406', 'CRR072407', 'CRR072408', 'CRR072409']
    if spc == 'landrace':
        samples = ['CRR072401', 'SRR7164576', 'SRR7164580', 'SRR7164572', 'SRR7164606']
    if spc == 'cultivar':
        samples = ['SRR7164620', 'SRR7164628', 'SRR7164670', 'SRR7164669', 'SRR7164604']
    if spc == 'durum':
        samples = ['CRR072337', 'CRR072338', 'CRR072340', 'CRR072341', 'CRR072347']
    if spc == 'wild_emmer':
        samples = ['CRR072247', 'CRR072248', 'CRR072249', 'CRR072250', 'CRR072251']
    return samples


def extract_derived(genome):
    species = []
    if genome == 'd':
        species = ['tauchii', 'landrace', 'cultivar']
    if genome == 'ab':
        species = ['landrace', 'cultivar', 'durum', 'wild_emmer']
    summary = pd.DataFrame()
    for spc in species:
        sample_name = choose_sample(spc)
        for name in sample_name:
            sample = pd.read_table('data/' + genome + '/' + genome + '.' + name + '.sum', sep='\t', header=None)
            sample.columns = ['seq', 'length', 'bases', 'mean', 'min', 'max']
            cu_code = pd.read_table('data/curatedLib.txt', sep=":", header=None)
            cu_code.columns = ['Classification', 'TE', 'species']
            # 1. classify
            sample['Classification'] = ''
            for i in range(len(sample)):
                for j in range(len(cu_code)):
                    if sample['seq'][i] == cu_code['TE'][j]:
                        sample.loc[i, 'Classification'] = cu_code.loc[j, 'Classification']
                        break
                if sample['Classification'][i] == '':
                    sample.loc[i, 'Classification'] = sample['seq'][i].split('#')[1]
                    sample = te_code(sample)
            # 2. stats
            grouped = sample.groupby('Classification')
            summ = grouped.agg(
                count=('Classification', 'count'),
                size=('bases', 'sum'),
            )
            summ['species'] = spc
            summary = pd.concat([summary, summ])
    return summary


def to_boxplot():
    summ = extract_derived('ab')
    summ['size'] = summ['size'] / 10**9
    fig, ax = plt.subplots()
    plt.rcParams['font.size'] = 24
    ax.figure.set_size_inches(18, 9)
    sns.boxplot(x='Classification', y='size', hue='species', data=summ, showfliers=False,
                palette='Set3', linewidth=2, dodge=True, ax=ax, width=0.8)
    plt.title('TE content in different species')
    plt.ylabel('size (Gb)')
    plt.legend(title='Name')
    plt.show()


def main():
    # bed = fai_to_bed(sys.argv[1])
    # bed.to_csv(sys.argv[2], sep='\t', header=False, index=False)
    # upset_plot('cs')

    # 1. run stats
    to_boxplot()


if __name__ == '__main__':
    main()
