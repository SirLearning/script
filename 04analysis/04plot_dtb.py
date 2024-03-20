import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

file_name = pd.read_table('data/dtb_files.name', header=None)
file_name = file_name.fillna('NULL')


mpl.rcParams['font.size'] = 20
plt.style.use('seaborn-v0_8-deep')

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15


def plot_te_hist(axs, file, alpha, method):
    file = str(file)
    density = pd.read_table('data/distribution/' + method + '.' + file + '.density.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    axs.hist(density['win_start']/1000000, weights=density['density'], bins=density['win_num'].iloc[-1], alpha=1 - alpha / 36,
             label=file)


def te_plot(axs, file, chrom, i):
    file = str(file)
    density = pd.read_table('data/' + chrom + '/stats.' + file + '.dtb.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    axs.plot(density['win_start'] / 1000000, density['density'], label=file, alpha=1-i/36)


def compare_te_hist(axs, file, method):
    file = str(file)
    density = pd.read_table('data/distribution/' + method + '/' + method + '.' + file + '.density.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    axs.hist(density['win_start']/1000000, weights=density['density'], bins=density['win_num'].iloc[-1], alpha=0.8, label=method)


def difference_te_hist(axs, file, method):
    file = str(file)
    density = pd.read_table('data/' + method + '/chr1ANM.' + file + '.density.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    axs.hist(density['win_start']/1000000, weights=density['density'], bins=density['win_num'].iloc[-1], alpha=0.4, label=method)


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
# ax.figure.set_size_inches(12, 8)
# for i in range(0, len(file_name)):
#     te_plot(ax, file_name.iloc[i, 0], 'chr1D1', i)
# ax.set_title("TE distribution along chr1D")
# ax.set_xlabel('chromosome (Mb)')
# ax.set_ylabel('density')
# ax.legend(fontsize=11.5, framealpha=0.5)
# plt.show()

fig, ax = plt.subplots()
ax.figure.set_size_inches(12, 8)
for i in range(0, len(file_name)):
    te_plot(ax, file_name.iloc[i, 0], '01stg', i)
ax.set_title("TE distribution along chr1D (93-96 Mb)")
ax.set_xlabel('chromosome (Mb)')
ax.set_ylabel('density')
ax.legend(fontsize=11.5, framealpha=0.5)
plt.show()
fig, ax = plt.subplots()
ax.figure.set_size_inches(12, 8)
for i in range(0, len(file_name)):
    te_plot(ax, file_name.iloc[i, 0], '02stg', i)
ax.set_title("TE distribution along chr1D (358-360 Mb)")
ax.set_xlabel('chromosome (Mb)')
ax.set_ylabel('density')
ax.legend(fontsize=11.5, framealpha=0.5)
plt.show()


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
