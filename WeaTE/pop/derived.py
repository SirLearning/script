import sys
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import upsetplot
from upsetplot import from_contents

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


def cs_age(genome):
    spc = []
    age = []
    sign = []
    seq = []
    if genome == 'a':
        spc = ['durum', 'wild_emmer', 'D_lineage', 'B_lineage']
        age = ['durum(8,441_ya)', 'wild_emmer(10,041_ya)', 'D_lineage(5.69_Mya)', 'B_lineage(6.67_Mya)']
        sign = ['a-dw', 'a-we', 'a-dcs', 'a-bcs']
        seq = ['cs-dw', 'cs-we', 'cs-dcs', 'cs-bcs']
    if genome == 'b':
        spc = ['durum', 'wild_emmer', 'A_lineage', 'D_lineage']
        age = ['durum(8,441_ya)', 'wild_emmer(10,041_ya)', 'A_lineage(6.67_Mya)', 'D_lineage(6.67_Mya)']
        sign = ['b-dw', 'b-we', 'b-acs', 'b-dcs']
        seq = ['cs-dw', 'cs-we', 'cs-acs', 'cs-dcs']
    if genome == 'd':
        spc = ['tauchii', 'A_lineage', 'B_lineage']
        age = ['tauchii(175,987_ya)', 'A_lineage(5.69_Mya)', 'B_lineage(6.67_Mya)']
        sign = ['d-at', 'd-acs', 'd-bcs']
        seq = ['cs-at', 'cs-acs', 'cs-bcs']
    return age, spc, sign, seq


def mod_sample(sample):
    cu_code = pd.read_table('data/curatedLib.txt', sep=":", header=None)
    cu_code.columns = ['Classification', 'TE', 'species']
    sample['Classification'] = ''
    for i in range(len(sample)):
        for j in range(len(cu_code)):
            if sample['seq'][i] == cu_code['TE'][j]:
                sample.loc[i, 'Classification'] = cu_code.loc[j, 'Classification']
                break
        if sample['Classification'][i] == '':
            sample.loc[i, 'Classification'] = sample['seq'][i].split('#')[1]
            sample = te_code(sample)
    return sample


def read_temp():
    summ = pd.read_table('data/temp', sep='\t', header=0)
    return summ


def read_record(sign):
    summ = pd.read_table('data/' + sign, sep='\t', header=0)
    return summ


def write_temp(df):
    df.to_csv('data/temp', sep='\t', header=True, index=True)


def extract_derived(genome):
    species = []
    if genome == 'd':
        species = ['landrace', 'cultivar', 'tauchii']
    if genome == 'ab':
        species = ['landrace', 'cultivar', 'durum', 'wild_emmer']
    summary = pd.DataFrame()
    for spc in species:
        sample_name = choose_sample(spc)
        for name in sample_name:
            sample = pd.read_table('data/' + genome + '/' + genome + '.' + name + '.sum', sep='\t', header=None)
            sample.columns = ['seq', 'length', 'bases', 'mean', 'min', 'max']
            sample = mod_sample(sample)
            # 2. stats
            grouped = sample.groupby('Classification')
            summ = grouped.agg(
                count=('Classification', 'count'),
                size=('bases', 'sum'),
            )
            summ['species'] = spc
            summ['sample'] = name
            summary = pd.concat([summary, summ])
    return summary


def cs_div(genome):
    [ages, species, signals] = cs_age(genome)
    sample_name = choose_sample('landrace')
    summary = pd.DataFrame()
    for i in range(len(species)):
        for name in sample_name:
            sample = pd.read_table('data/cs/' + signals[i] + '.' + name + '.sum', sep='\t', header=None)
            sample.columns = ['seq', 'length', 'bases', 'mean', 'min', 'max']
            sample = mod_sample(sample)
            grouped = sample.groupby('Classification')
            summ = grouped.agg(
                count=('Classification', 'count'),
                size=('bases', 'sum'),
            )
            summ['species'] = species[i]
            summ['sample'] = name
            summ['age'] = ages[i]
            summary = pd.concat([summary, summ])
    return summary


def cs_cor():
    cor = pd.DataFrame()
    genome = 'd'
    [ages, species, signals, seqs] = cs_age(genome)
    cor = {}
    for i in range(len(species)):
        data_list = pd.read_table('data/cs/' + genome + '/' + seqs[i], sep='\t', header=None)[0]
        cor['not in ' + species[i]] = data_list
    # print(cor)
    df = from_contents(cor)
    # df = df.dropna()
    if genome == 'a': upsetplot.plot(df, show_counts=True, facecolor="#348ABD")
    if genome == 'b': upsetplot.plot(df, show_counts=True, facecolor="#7A68A6")
    if genome == 'd': upsetplot.plot(df, show_counts=True, facecolor="#467821")
    plt.show()
    plt.show()


def lib_boxplot():
    # summ = extract_derived('ab')
    # write_temp(summ)
    # summ = read_temp()
    summ = read_record('lib_d')
    summ['size'] = summ['size'] / 10**9
    fig, ax = plt.subplots()
    plt.rcParams['font.size'] = 20
    ax.figure.set_size_inches(16, 7.5)
    sns.boxplot(x='Classification', y='size', hue='species', data=summ, showfliers=True,
                palette='Set3', linewidth=1.5, dodge=True, ax=ax, width=0.75)
    plt.title('TE content in different species (genome D)')
    plt.ylabel('size (Gb)')
    plt.legend(title='Species', framealpha=0.5, fontsize=20, title_fontsize=20)
    plt.show()


def cs_boxplot():
    # summ = cs_div('d')
    # write_temp(summ)
    # summ = read_temp()
    summ = read_record('cs_d')
    summ['size'] = summ['size'] / 10**9
    fig, ax = plt.subplots()
    plt.rcParams['font.size'] = 24
    ax.figure.set_size_inches(16, 7.5)
    sns.boxplot(x='Classification', y='size', hue='age', data=summ, showfliers=True,
                palette='Set3', linewidth=1.5, dodge=True, ax=ax, width=0.75)
    plt.title('TE mutation load & Diverge time (genome D)')
    # plt.xlabel('TE type')
    plt.ylabel('TE Mutation Load (Gb)')
    plt.legend( title='BW diverge time with:', framealpha=0.5, fontsize=18, title_fontsize=22)
    plt.show()


def main():
    # bed = fai_to_bed(sys.argv[1])
    # bed.to_csv(sys.argv[2], sep='\t', header=False, index=False)

    # upset_plot('cs')

    # 1. run stats
    # lib_boxplot()
    # cs_boxplot()
    cs_cor()


if __name__ == '__main__':
    main()
