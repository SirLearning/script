import pandas as pd


def main():
    lib_name = 'transposon/intact.fa.fai'
    TElib_name = 'transposon/TElib.fai'
    intactTE_name = 'transposon/intactTE.fai'

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
        mean_length=('length', 'mean'),
        sum_length=('length', 'sum'),
        max_length=('length', 'max'),
        min_length=('length', 'min')
    )
    print(summary)

    lib[lib['type'] == 'TElib'].to_csv(TElib_name, sep='\t', index=False)
    lib[lib['type'] == 'intacLib'].to_csv(intactTE_name, sep='\t', index=False)


if __name__ == '__main__':
    main()
