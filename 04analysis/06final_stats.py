import pandas as pd

file_name = pd.read_table('data/dtb_files.name', header=None)
for i in [0, len(file_name) - 1]:
    chr1A = pd.read_table('data/distribution/chr1A.' + file_name.iloc[i, 0] + '.density.txt', sep='\s+', header=None)
    chr1A.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
    chr1ANM = pd.read_table('data/distribution/chr1ANM.' + file_name.iloc[i, 0] + '.density.txt', sep='\s+', header=None)
    chr1ANM.columns = ['chrom', 'win_start', 'win_end', 'win_num', 'none', 'strand', 'density']
