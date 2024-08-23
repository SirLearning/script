import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd


def dtb(file_name):
    fig, ax = plt.subplots()
    ax.figure.set_size_inches(12, 8)
    for i in range(0, len(file_name)):
        te_plot(ax, file_name.iloc[i, 0], 'chr1A', i)
    ax.set_title("chr1A")
    ax.set_xlabel('chromosome (Mb)')
    ax.set_ylabel('density')
    ax.legend(fontsize=14, framealpha=0.5)
    plt.show()


def plot_te_hist(axs, file, alpha, method):
    file = str(file)
    density = pd.read_table('data/distribution/' + method + '.' + file + '.density.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    axs.hist(density['win_start']/10**6, weights=density['density'], bins=density['win_num'].iloc[-1], alpha=1 - alpha / 36,
             label=file)


def te_plot(axs, file, chrom, i):
    file = str(file)
    density = pd.read_table('data/' + chrom + '/stats.' + file + '.dtb.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    axs.plot(density['win_start'] / 10**6, density['density'], label=file, alpha=1-i/36)


def compare_te_hist(axs, file, method):
    file = str(file)
    density = pd.read_table('data/distribution/' + method + '/' + method + '.' + file + '.density.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    axs.hist(density['win_start']/1000000, weights=density['density'], bins=density['win_num'].iloc[-1], alpha=0.8, label=method)


def compare_te(axs, file, name, method, color, fig_type):
    file = str(file)
    density = pd.read_table('data/old/' + method + '/stats.' + file + '.dtb.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    if fig_type == 'hist':
        axs.hist(density['win_start']/1000000, weights=density['density'], bins=density['win_num'].iloc[-1], alpha=0.7,
                 label=name, color=color)
    elif fig_type == 'plot':
        axs.plot(density['win_start']/1000000, density['density'], label=method, alpha=0.8, color=color)
        # axs.plot(density['win_start']/1000000, density['density'], label=triticeae, alpha=0.8, color=color)


def difference_te_hist(axs, file, method):
    file = str(file)
    density = pd.read_table('data/' + method + '/chr1ANM.' + file + '.density.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    axs.hist(density['win_start']/1000000, weights=density['density'], bins=density['win_num'].iloc[-1], alpha=0.4, label=method)


def difference_te(axs, file, name, method, color, fig_type):
    file = str(file)
    density = pd.read_table('data/old/' + method + '/stats.' + file + '.dtb.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    if fig_type == 'hist':
        axs.hist(density['win_start']/1000000, weights=density['density'], bins=density['win_num'].iloc[-1], alpha=0.8,
                 label=name, color=color)
    elif fig_type == 'plot':
        axs.plot(density['win_start']/1000000, density['density'], label=name, alpha=0.8, color=color)


def centromere(axs, spc, method):
    if spc == 'cs_n' or spc == 'ldr':
        if method == 'chr1A':
            axs.axvspan(209.9, 216.2, color='green', alpha=0.01)
        if method == 'chr1B':
            axs.axvspan(237.7, 244.5, color='green', alpha=0.01)
        if method == 'chr1D':
            axs.axvspan(166.0, 174.6, color='green', alpha=0.01)
    if spc == 'we':
        if method == 'chr1A':
            axs.axvspan(116.0, 116.7, color='green', alpha=0.01)
            axs.axvspan(212.1, 215.7, color='green', alpha=0.01)
            axs.axvspan(217.0, 218.7, color='green', alpha=0.01)
        if method == 'chr1B':
            axs.axvspan(234.6, 235.2, color='green', alpha=0.01)
            axs.axvspan(237.2, 238.3, color='green', alpha=0.01)
            axs.axvspan(277.9, 278.7, color='green', alpha=0.01)
    if spc == 'dw':
        if method == 'chr1A':
            axs.axvspan(208.0, 211.4, color='green', alpha=0.01)
        if method == 'chr1B':
            axs.axvspan(235.4, 237.2, color='green', alpha=0.01)
    if spc == 'at':
        axs.axvspan(175.9, 178.5, color='green', alpha=0.01)
        axs.axvspan(179.1, 180.9, color='green', alpha=0.01)


def genome_dtb(axs, spc, name, method, i):
    cover = pd.read_table('data/' + spc + '/' + method + '/coverage.' + name + '.dtb.txt', sep='\s+', header=None)
    cover.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'coverage']
    axs.plot(cover['win_start']/1000000, cover['coverage']*100, label=name, alpha=1-i/20, linewidth=4)
    centromere(axs, spc, method)
    # axs.hist(cover['win_start'] / 1000000, weights=cover['coverage'], bins=cover['win_num'].iloc[-1], alpha=0.4, label=name)


def do_plot(axs, file_name, chr):
    axs.figure.set_size_inches(16, 7)
    for i in range(0, 12):
        genome_dtb(axs, 'cs_n', file_name[i], chr, i)
    axs.set_xlabel('chromosome (Mb)')
    axs.set_ylabel('coverage (%)')
    # axs.set_title("TE distribution along " + chr)
    # axs.legend(fontsize=14, framealpha=0.5)
    plt.show()


def cs_dtb(axs, spc, name, chr, i):
    cover = pd.read_table('data/' + spc + '/coverage.' + name + '.dtb.txt', sep='\t', header=None)
    cover.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'coverage']
    cover = cover[cover['chrom'] == chr]
    axs.plot(cover['win_start']/1000000, cover['coverage']*100, label=name, alpha=1-i/20, linewidth=4)
    centromere(axs, spc, chr)
    # axs.hist(cover['win_start'] / 1000000, weights=cover['coverage'], bins=cover['win_num'].iloc[-1], alpha=0.4, label=name


def new_plot(axs, file_name, chr):
    axs.figure.set_size_inches(16, 7)
    for i in range(0, 11):
        cs_dtb(axs, 'cs_v2', file_name[i], chr, i)
    axs.set_xlabel('chromosome (Mb)')
    axs.set_ylabel('coverage (%)')
    # axs.set_title("TE distribution along " + chr)
    # axs.legend(fontsize=14, framealpha=0.5)
    plt.show()


def plot_length(name):
    length = pd.read_table('data/old/' + name + '/stats.length.txt', sep='\t')
    length.columns = ['Classification', 'width']
    length['length (kb)'] = length['width'] / 1000
    length['triticeae'] = name  # 添加一个新的列 'triticeae'
    return length


def main():
    file_name = ['DHH', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'RIJ', 'RIL', 'RIR', 'RLC', 'RLG']
    plt.style.use('seaborn-v0_8-deep')
    # plt.style.use('fast')
    mpl.rcParams['font.size'] = 24
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20
    fig, axA = plt.subplots()
    new_plot(axA, file_name, 'Chr1A')
    fig, axB = plt.subplots()
    new_plot(axB, file_name, 'Chr1B')
    fig, axD = plt.subplots()
    new_plot(axD, file_name, 'Chr1D')
    # fig, axR = plt.subplots()
    # do_plot(axR, file_name, 'chr1R')
    # fig, axH = plt.subplots()
    # do_plot(axH, file_name, 'chr1H')




if __name__ == '__main__':
    main()

# # plot all TEs in chr1A
# fig, ax = plt.subplots()
# ax.figure.set_size_inches(12, 8)
# for i in range(0, len(file_name)):
#     te_plot(ax, file_name.iloc[i, 0], 'chr1A', i)
# ax.set_title("TE distribution along chr1A")
# ax.set_xlabel('chromosome (Mb)')
# ax.set_ylabel('density')
# ax.legend(fontsize=11.5, framealpha=0.5)
# plt.show()
# # plot all TEs in chr1B
# fig, ax = plt.subplots()
# ax.figure.set_size_inches(12, 8)
# for i in range(0, len(file_name)):
#     te_plot(ax, file_name.iloc[i, 0], 'chr1B', i)
# ax.set_title("TE distribution along chr1B")
# ax.set_xlabel('chromosome (Mb)')
# ax.set_ylabel('density')
# ax.legend(fontsize=11.5, framealpha=0.5)
# plt.show()
# # plot all TEs in chr1D
# fig, ax = plt.subplots()
# ax.figure.set_size_inches(12, 8)
# for i in range(0, len(file_name)):
#     te_plot(ax, file_name.iloc[i, 0], 'chr1D', i)
# ax.set_title("TE distribution along chr1D")
# ax.set_xlabel('chromosome (Mb)')
# ax.set_ylabel('density')
# ax.legend(fontsize=11.5, framealpha=0.5)
# plt.show()
#
#
# fig, ax = plt.subplots()
# ax.figure.set_size_inches(16, 8)
# for i in range(0, len(file_name)):
#     te_plot(ax, file_name.iloc[i, 0], 'chr1D1', i)
# ax.set_title("TE distribution along chr1D")
# ax.set_xlabel('chromosome (Mb)')
# ax.set_ylabel('density')
# ax.legend(fontsize=11.5, framealpha=0.5)
# plt.show()

# # plot strange region in chr1D
# fig, ax = plt.subplots()
# ax.figure.set_size_inches(12, 8)
# for i in range(0, len(file_name)):
#     te_plot(ax, file_name.iloc[i, 0], '01stg', i)
# ax.set_title("TE distribution along chr1D (93-96 Mb)")
# ax.set_xlabel('chromosome (Mb)')
# ax.set_ylabel('density')
# ax.legend(fontsize=14, framealpha=0.5)
# plt.show()
# fig, ax = plt.subplots()
# ax.figure.set_size_inches(12, 8)
# for i in range(0, len(file_name)):
#     te_plot(ax, file_name.iloc[i, 0], '02stg', i)
# ax.set_title("TE distribution along chr1D (358-360 Mb)")
# ax.set_xlabel('chromosome (Mb)')
# ax.set_ylabel('density')
# ax.legend(fontsize=14, framealpha=0.5)
# plt.show()




# fig, ax = plt.subplots()
# ax.figure.set_size_inches(10, 5)
# for i in range(0, len(file_name)):
#     plot_te_hist(ax, file_name.iloc[i, 0], i, 'chr1ANM')
# ax.set_title("TE distribution along the chromosome")
# ax.set_xlabel('chromosome (Mb)')
# ax.set_ylabel('density')
# ax.legend(fontsize=8, framealpha=0.5)
# plt.show()

# # plot every TE in different method
# for i in range(0, 18):
#     fig, ax = plt.subplots()
#     ax.figure.set_size_inches(10, 5)
#     compare_te_hist(ax, file_name.iloc[i, 0], 'chr1A')
#     compare_te_hist(ax, file_name.iloc[i, 0], 'chr1ANM')
#     ax.set_title(file_name.iloc[i, 0] + " distribution along the chromosome")
#     ax.set_xlabel('chromosome (Mb)')
#     ax.set_ylabel('density')
#     ax.legend(fontsize=12, framealpha=0.5)
#     plt.show()

# # difference in different method
# for i in range(0, 18):
#     fig, ax = plt.subplots()
#     ax.figure.set_size_inches(10, 5)
#     compare_te_hist(ax, file_name.iloc[i, 0], 'chr1A')
#     compare_te_hist(ax, file_name.iloc[i, 0], 'chr1ANM')
#     difference_te_hist(ax, file_name.iloc[i, 0], 'ANv')
#     difference_te_hist(ax, file_name.iloc[i, 0], 'NAv')
#     ax.set_title(file_name.iloc[i, 0] + " distribution along the chromosome")
#     ax.set_xlabel('chromosome (Mb)')
#     ax.set_ylabel('density')
#     ax.legend(fontsize=12, framealpha=0.5)
#     plt.show()
