import matplotlib.pyplot as plt
import pandas as pd

file_name = pd.read_table('data/dtb_files.name', header=None)
file_name = file_name.fillna('NULL')
plt.style.use('bmh')


def plot_te_hist(axs, file, alpha, method):
    file = str(file)
    density = pd.read_table('data/distribution/' + method + '.' + file + '.density.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    axs.hist(density['win_start'], weights=density['density'], bins=density['win_num'].iloc[-1], alpha=1 - alpha / 36,
             label=file)
    axs.set_xlabel('chromosome')
    axs.set_ylabel('density')


def compare_te_hist(axs, file, method):
    file = str(file)
    density = pd.read_table('data/distribution/' + method + '.' + file + '.density.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    axs.hist(density['win_start'], weights=density['density'], bins=density['win_num'].iloc[-1], alpha=0.6, label=method)
    axs.set_xlabel('chromosome')
    axs.set_ylabel('density')


# plot all TEs in one method
fig, ax = plt.subplots()
ax.figure.set_size_inches(10, 5)
for i in range(0, len(file_name) - 1):
    plot_te_hist(ax, file_name.iloc[i, 0], i, 'chr1A')
ax.set_title("TE distribution along the chromosome")
ax.legend(fontsize=8, framealpha=0.5)
plt.show()

fig, ax = plt.subplots()
ax.figure.set_size_inches(10, 5)
for i in range(0, len(file_name) - 1):
    plot_te_hist(ax, file_name.iloc[i, 0], i, 'chr1ANM')
ax.set_title("TE distribution along the chromosome")
ax.legend(fontsize=8, framealpha=0.5)
plt.show()

# plot every TE in different method
for i in range(0, 17):
    fig, ax = plt.subplots()
    ax.figure.set_size_inches(10, 5)
    compare_te_hist(ax, file_name.iloc[i, 0], 'chr1A')
    compare_te_hist(ax, file_name.iloc[i, 0], 'chr1ANM')
    ax.set_title(file_name.iloc[i, 0] + " distribution along the chromosome")
    ax.legend(fontsize=12, framealpha=0.5)
    plt.show()
