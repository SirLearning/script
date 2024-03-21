import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl

file_name = pd.read_table('data/dtb_files.name', header=None)
file_name = file_name.fillna('NULL')
plt.style.use('seaborn-v0_8-deep')


def compare_te_hist(axs, file, name, method, color):
    file = str(file)
    density = pd.read_table('data/' + method + '/stats.' + file + '.dtb.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    # axs.hist(density['win_start']/1000000, weights=density['density'], bins=density['win_num'].iloc[-1], alpha=0.7,
    #          label=name, color=color)
    # axs.plot(density['win_start'] / 1000000, density['density'], label=method, alpha=0.8, color=color)
    axs.plot(density['win_start']/1000000, density['density'], label=name, alpha=0.8, color=color)

def difference_te_hist(axs, file, name, method, color):
    file = str(file)
    density = pd.read_table('data/distribution/' + method + '/stats.' + file + '.dtb.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    axs.hist(density['win_start'] / 1000000, weights=density['density'], bins=density['win_num'].iloc[-1], alpha=0.8,
             label=name, color=color)
    # axs.plot(density['win_start']/1000000, density['density'], label=name, alpha=0.8, color=color)


mpl.rcParams['font.size'] = 24
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20

# # difference in different method
# for i in range(0, 18):
#     fig, ax = plt.subplots()
#     ax.figure.set_size_inches(10, 8)
#     compare_te_hist(ax, file_name.iloc[i, 0], 'intact chr1A', 'chr1A', '#348ABD')
#     compare_te_hist(ax, file_name.iloc[i, 0], 'N-sliced chr1A', 'chr1ANM', '#7A68A6')
#     difference_te_hist(ax, file_name.iloc[i, 0], 'Only in intact chr1A', 'ANv', '#A60628')
#     difference_te_hist(ax, file_name.iloc[i, 0], 'Only in N-sliced chr1A', 'NAv', '#467821')
#     ax.set_xlabel('chromosome (Mb)')
#     ax.set_ylabel('density')
#     ax.set_title(file_name.iloc[i, 0] + " distribution along the chromosome")
#     ax.legend(fontsize=14, framealpha=0.5)
#     plt.show()

# correlation between EDTA and CS
for i in range(0, 18):
    fig, ax = plt.subplots()
    ax.figure.set_size_inches(10, 8)
    compare_te_hist(ax, file_name.iloc[i, 0], 'intact chr1A', 'EDTA', '#348ABD')
    compare_te_hist(ax, file_name.iloc[i, 0], 'N-sliced chr1A', 'CS', '#7A68A6')
    difference_te_hist(ax, file_name.iloc[i, 0], 'Only in intact chr1A', 'ANv', '#A60628')
    difference_te_hist(ax, file_name.iloc[i, 0], 'Only in N-sliced chr1A', 'NAv', '#467821')
    ax.set_xlabel('chromosome (Mb)')
    ax.set_ylabel('density')
    ax.set_title(file_name.iloc[i, 0] + " distribution along the chromosome")
    ax.legend(fontsize=14, framealpha=0.5)
    plt.show()