import pandas as pd

lib_name = 'data/intact.fa.fai'

# Read the library file
lib = pd.read_table(lib_name, sep="\t", header=None)
lib.columns = ['seq', 'length', 'start', 'line_length', 'line_bytes']
lib['name'] = lib['seq'].str.split('|').str[0]
lib['locus'] = lib['seq'].str.split('|').str[1]

for i in range(0, len(lib)):
    if lib.loc[i, 'name'].startswith('TE'):
        lib.loc[i, 'type'] = 'TElib'
    else:
        lib.loc[i, 'type'] = 'intacLib'

grouped = lib.groupby('type')
summary = grouped.agg(
    count=('seq', 'count'),
)
print(summary)
