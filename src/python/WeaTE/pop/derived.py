import sys
import string
import scipy
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
    order = []
    region = []
    if genome == 'a':
        spc = ['durum', 'wild_emmer', 'D_lineage', 'B_lineage']
        age = ['durum(8,441_ya)', 'wild_emmer(10,041_ya)', 'D_lineage(5.69_Mya)', 'B_lineage(6.67_Mya)']
        sign = ['a-dw', 'a-we', 'a-dcs', 'a-bcs']
        seq = ['cs-dw', 'cs-we', 'cs-dcs', 'cs-bcs']
        order = ['a-0', 'a-1', 'a-2', 'a-3']
        region = ['ancestor-WE', 'WE-DW', 'DW-BW', 'BW-now']
    if genome == 'b':
        spc = ['durum', 'wild_emmer', 'A_lineage', 'D_lineage']
        age = ['durum(8,441_ya)', 'wild_emmer(10,041_ya)', 'A_lineage(6.67_Mya)', 'D_lineage(6.67_Mya)']
        sign = ['b-dw', 'b-we', 'b-acs', 'b-dcs']
        seq = ['cs-dw', 'cs-we', 'cs-acs', 'cs-dcs']
        order = ['b-0', 'b-1', 'b-2', 'b-3']
        region = ['ancestor-WE', 'WE-DW', 'DW-BW', 'BW-now']
    if genome == 'd':
        spc = ['tauchii', 'A_lineage', 'B_lineage']
        age = ['tauchii(175,987_ya)', 'A_lineage(5.69_Mya)', 'B_lineage(6.67_Mya)']
        sign = ['d-at', 'd-acs', 'd-bcs']
        seq = ['cs-at', 'cs-acs', 'cs-bcs']
        order = ['d-0', 'd-1', 'd-2']
        region = ['ancestor-AT', 'AT-BW', 'BW-now']
    if genome == 'we':
        order = ['we-0', 'we-1', 'we-2']
        region = ['ancestor-WE', 'WE-DW', 'WE-now']
    if genome == 'dw':
        order = ['dw-0', 'dw-1', 'dw-2']
        region = ['ancestor-WE', 'WE-DW', 'DW-now']
    if genome == 'ctv':
        order = ['ctv-0', 'ctv-1', 'ctv-2']
        region = ['ancestor-landrace', 'landrace-cultivar', 'cultivar-now']
    if genome == 'ldr':
        order = ['ldr-0', 'ldr-1', 'ldr-2']
        region = ['ancestor-landrace', 'landrace-cultivar', 'landrace-now']
    return age, spc, sign, seq, order, region


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
    summ = pd.read_table('transposon/' + sign, sep='\t', header=0)
    return summ


def write_temp(df):
    df.to_csv('transposon/temp', sep='\t', header=True, index=True)


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
            sample = pd.read_table('transposon/' + genome + '/' + genome + '.' + name + '.sum', sep='\t', header=None)
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
    [ages, species, signals, region1, region2] = cs_age(genome)
    sample_name = choose_sample('landrace')
    summary = pd.DataFrame()
    for i in range(len(species)):
        for name in sample_name:
            sample = pd.read_table('transposon/cs/' + signals[i] + '.' + name + '.sum', sep='\t', header=None)
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


def cs_region(genome):
    [ages, species, signals, seqs, orders, regions] = cs_age(genome)
    sample_name = choose_sample('landrace')
    summary = pd.DataFrame()
    for i in range(len(regions)):
        for name in sample_name:
            sample = pd.read_table('transposon/cs/' + orders[i] + '.' + name + '.sum', sep='\t', header=None)
            sample.columns = ['seq', 'length', 'bases', 'mean', 'min', 'max']
            sample = mod_sample(sample)
            grouped = sample.groupby('Classification')
            summ = grouped.agg(
                count=('Classification', 'count'),
                size=('bases', 'sum'),
            )
            summ['region'] = regions[i]
            summ['sample'] = name
            summary = pd.concat([summary, summ])
    return summary


def dm_region(genome):
    [ages, species, signals, seqs, orders, regions] = cs_age(genome)
    sample_name = []
    if genome == 'we':
        sample_name = choose_sample('wild_emmer')
    if genome == 'dw':
        sample_name = choose_sample('durum')
    if genome == 'ctv':
        sample_name = choose_sample('cultivar')
    if genome == 'ldr':
        sample_name = choose_sample('landrace')
    summary = pd.DataFrame()
    for i in range(len(regions)):
        for name in sample_name:
            sample = pd.read_table('transposon/lc/' + orders[i] + '.' + name + '.sum', sep='\t', header=None)
            sample.columns = ['seq', 'length', 'bases', 'mean', 'min', 'max']
            sample = mod_sample(sample)
            grouped = sample.groupby('Classification')
            summ = grouped.agg(
                count=('Classification', 'count'),
                size=('bases', 'sum'),
            )
            summ['region'] = regions[i]
            summ['sample'] = name
            summary = pd.concat([summary, summ])
    return summary


def cs_cor():
    cor = pd.DataFrame()
    genome = 'd'
    [ages, species, signals, seqs] = cs_age(genome)
    cor = {}
    for i in range(len(species)):
        data_list = pd.read_table('transposon/cs/' + genome + '/' + seqs[i], sep='\t', header=None)[0]
        cor['not in ' + species[i]] = data_list
    # print(cor)
    df = from_contents(cor)
    # df = df.dropna()
    if genome == 'a': upsetplot.plot(df, show_counts=True, facecolor="#348ABD")
    if genome == 'b': upsetplot.plot(df, show_counts=True, facecolor="#7A68A6")
    if genome == 'd': upsetplot.plot(df, show_counts=True, facecolor="#467821")
    plt.show()
    plt.show()


def count_p(file1, file2, region1, region2, Classification):
    df1 = read_record(file1)
    df2 = read_record(file2)
    stat, p_value = scipy.stats.ttest_ind(df1[(df1["region"] == region1) & (df1["Classification"] == Classification)]["size"],
                                          df2[(df2["region"] == region2) & (df2["Classification"] == Classification)]["size"],
                                          equal_var=False)
    print(stat, p_value)
    return p_value


def sm_cp(spc, regions):
    sample_name = []
    if spc == 'we':
        sample_name = choose_sample('wild_emmer')
    if spc == 'dw':
        sample_name = choose_sample('durum')
    RLC = pd.DataFrame()
    for i in range(len(regions)):
        for j in range(len(sample_name)):
            sample = pd.read_table('transposon/dm/' + spc + '-' + str(i) + '.' + sample_name[j] + '.sum', sep='\t', header=None)
            sample.columns = ['seq', 'length', 'bases', 'mean', 'min', 'max']
            sample = mod_sample(sample)
            sample_RLC = sample[sample['Classification'] == 'RLC']
            sample_RLC['region'] = regions[i]
            sample_RLC['sample'] = spc + '-' + str(j+1)
            sample_RLC['species'] = spc
            RLC = pd.concat([RLC, sample_RLC])
    return RLC


def order_region(genome,spc):
    order = []
    region = []
    if genome == 'd':
        if spc == 'at':
            order = ['at-0', 'at-1', 'at-2']
            region = ['ancestor-AT', 'AT-BW', 'AT-now']
        if spc == 'cs':
            order = ['cs-0', 'cs-1', 'cs-2']
            region = ['ancestor-AT', 'AT-BW', 'BW-now']
    if genome == 'ab':
        if spc == 'cs':
            order = ['cs-0', 'cs-1', 'cs-2']
            region = ['ancestor-DW', 'DW-BW', 'BW-now']
        if spc == 'dw':
            order = ['dw-0', 'dw-1', 'dw-2']
            region = ['ancestor-DW', 'DW-BW', 'DW-now']
    return order, region


def poly_region(genome, spc):
    [orders, regions] = order_region(genome, spc)
    sample_name = []
    if spc == 'cs':
        sample_name = choose_sample('landrace')
    if spc == 'at':
        sample_name = choose_sample('tauchii')
    if spc == 'we':
        sample_name = choose_sample('wild_emmer')
    if spc == 'dw':
        sample_name = choose_sample('durum')
    summary = pd.DataFrame()
    for i in range(len(regions)):
        for name in sample_name:
            sample = pd.read_table('transposon/poly/' + genome + '/' + orders[i] + '.' + name + '.sum', sep='\t', header=None)
            sample.columns = ['seq', 'length', 'bases', 'mean', 'min', 'max']
            sample = mod_sample(sample)
            grouped = sample.groupby('Classification')
            summ = grouped.agg(
                count=('Classification', 'count'),
                size=('bases', 'sum'),
            )
            summ['region'] = regions[i]
            summ['sample'] = name
            summary = pd.concat([summary, summ])
    return summary


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
    summ = cs_region('d')
    write_temp(summ)
    # summ = read_temp()
    # summ = read_record('cs_ra')
    summ['size'] = summ['size'] / 10**9
    fig, ax = plt.subplots()
    plt.rcParams['font.size'] = 24
    ax.figure.set_size_inches(16, 7.5)
    # sns.boxplot(x='Classification', y='size', hue='age', transposon=summ, showfliers=True,
    sns.boxplot(x='Classification', y='size', hue='region', data=summ, showfliers=True,
                palette='Set3', linewidth=1.5, dodge=True, ax=ax, width=0.75)
    plt.title('TE mutation load during wheat evolution (genome D)')
    # plt.xlabel('TE type')
    plt.ylabel('TE Mutation Load (Gb)')
    plt.legend( title='BW evolution period', framealpha=0.5, fontsize=18, title_fontsize=22)
    plt.show()


def dm_boxplot():
    # summ = dm_region('dw')
    # summ = poly_region('d', 'at')
    # summ = dm_region('ldr')
    # write_temp(summ)
    # summ = read_temp()
    summ = read_record('poly_d-at')
    summ['size'] = summ['size'] / 10**9
    summ['region'] = string.capwords(summ['region'])
    fig, ax = plt.subplots()
    plt.rcParams['font.size'] = 24
    ax.figure.set_size_inches(16, 7.5)
    sns.boxplot(x='Classification', y='size', hue='region', data=summ, showfliers=True,
                palette='Set3', linewidth=1.5, dodge=True, ax=ax, width=0.75)
    plt.title('TE mutation load during landrace evolution')
    # plt.title('TE mutation load during durum wheat (DW) evolution')
    # plt.xlabel('TE type')
    plt.ylabel('TE Mutation Load (Gb)')
    plt.legend( title='Evolution period', framealpha=0.5, fontsize=18, title_fontsize=22)
    plt.show()


def cp_num_boxplot():
    regions = ['ancestor-WE', 'WE-DW']
    we_RLC = sm_cp('we', regions)
    dw_RLC = sm_cp('dw', regions)
    RLC = pd.concat([we_RLC, dw_RLC])
    write_temp(RLC)
    # RLC = read_temp()
    # RLC = read_record('dm_RLC')
    # print(RLC[RLC['mean'] > 20000])
    fig, ax = plt.subplots()
    plt.rcParams['font.size'] = 24
    ax.figure.set_size_inches(16, 7.5)
    sns.boxplot(x='region', y='mean', hue='sample', data=RLC, showfliers=False,
                palette='Set3', linewidth=1.5, dodge=True, ax=ax, width=0.75)
    plt.title('RLC copy number in different samples before domestication')
    plt.ylabel('TE copy number')
    plt.legend( title='Samples', framealpha=0.5, fontsize=18, title_fontsize=22)
    plt.show()


def main():
    # bed = fai_to_bed(sys.argv[1])
    # bed.to_csv(sys.argv[2], sep='\t', header=False, index=False)

    # upset_plot('cs')

    # 1. run stats
    # lib_boxplot()
    # cs_boxplot()
    # cs_cor()

    # 2. analysis
    # dm_boxplot()
    # cp_num_boxplot()

    # 3. count p-value
    count_p('dm_we', 'dm_dw', 'ancestor-WE', 'RLC')
    # count_p('lc_ctv', 'lc_ldr', 'ancestor-landrace', 'ancestor-landrace', 'DHH')
    # count_p('lc_ctv', 'lc_ldr', 'ancestor-landrace', 'ancestor-landrace', 'RLC')
    # count_p('lc_ctv', 'lc_ldr', 'landrace-cultivar', 'landrace-cultivar', 'RLC')
    # count_p('lc_ctv', 'lc_ldr', 'landrace-cultivar', 'landrace-cultivar', 'RLG')
    # count_p('lc_ctv', 'lc_ldr', 'cultivar-now', 'landrace-now', 'DHH')
    # count_p('lc_ctv', 'lc_ldr', 'cultivar-now', 'landrace-now', 'RLG')
    # df = read_record('dm_RLC')
    # stat1, p_value1 = scipy.stats.ttest_ind(
    #     # df[(df["region"] == 'ancestor-WE') & (df["species"] == 'we')]["mean"],
    #     # df[(df["region"] == 'ancestor-WE') & (df["species"] == 'dw')]["mean"],
    #     df[(df["region"] == 'ancestor-WE') & (df["sample"] == 'we-1')]["mean"],
    #     df[(df["region"] == 'ancestor-WE') & (df["sample"] == 'dw-1')]["mean"],
    #     equal_var=False)
    # stat, p_value = scipy.stats.ttest_ind(
    #     # df[(df["region"] == 'WE-DW') & (df["species"] == 'we')]["mean"],
    #     # df[(df["region"] == 'WE-DW') & (df["species"] == 'dw')]["mean"],
    #     df[(df["region"] == 'WE-DW') & (df["sample"] == 'we-1')]["mean"],
    #     df[(df["region"] == 'WE-DW') & (df["sample"] == 'dw-1')]["mean"],
    #     equal_var=False)
    # print(p_value1, p_value)


if __name__ == '__main__':
    main()
