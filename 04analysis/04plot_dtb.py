import matplotlib.pyplot as plt
import pandas as pd

plt.style.use('bmh')


def plot_te_hist(axs, file, alpha):
    file = str(file)
    density = pd.read_table('data/distribution/chr1A.' + file + '.density.txt', sep='\s+', header=None)
    density.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    axs.hist(density['win_start'], weights=density['density'], bins=density['win_num'].iloc[-1], alpha=1-alpha/36, label=file)
    axs.set_xlabel('chromosome')
    axs.set_ylabel('density')


fig, ax = plt.subplots()
ax.figure.set_size_inches(10, 5)
file_name = pd.read_table('data/dtb_files.name', header=None)
file_name = file_name.fillna('NULL')
for i in range(0, len(file_name)-1):
    plot_te_hist(ax, file_name.iloc[i, 0], i)

ax.set_title("TE distribution along the chromosome")
ax.legend(fontsize=8, framealpha = 0.5)
plt.show()

fig, ax = plt.subplots()
ax.figure.set_size_inches(10, 5)
