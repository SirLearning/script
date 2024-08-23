# stats of TE number and length
import sys
import pandas as pd

"""
Module to summarize the number and length of TEs in a genome.
"""

""": parameter
fai = 'lib.fa.fai'
cu_code = 'curatedLib.txt'
output_name = 'data/curated_lib/stats.length.txt'
"""


def mod_clib(fai_name, curated_name):   # fai to classified length
    # 1. read lib.fai
    fai = pd.read_table(fai_name, sep="\t", header=None)
    fai.columns = ['TE', 'length', 'start', 'line', 'width']
    cu_code = pd.read_table(curated_name, sep=":", header=None)
    cu_code.columns = ['Classification', 'TE', 'species']
    # 2. classify lib
    fai['Classification'] = ''
    for i in range(len(fai)):
        for j in range(len(cu_code)):
            if fai['TE'][i] == cu_code['TE'][j]:
                fai.loc[i, 'Classification'] = cu_code.loc[j, 'Classification']
                break
    to_stats = pd.concat([fai['Classification'], fai['length']], axis=1)
    return to_stats


def mod_trep(fai_name):     # fai to classified length
    # 1. read lib.fai
    fai = pd.read_table(fai_name, header=None)
    fai.columns = ['TE', 'length', 'start', 'line', 'width']
    # 2. classify lib
    fai['Classification'] = fai['TE'].str.split('_').str[0]
    to_stats = pd.concat([fai['Classification'], fai['length']], axis=1)
    return to_stats


def RMout_to_bed(input_file):
    with open(input_file, 'r') as te, open(f"{input_file}.bed", 'w') as bed:
        for line in te:
            line = line.strip()
            if not any(keyword in line for keyword in ["Score", "score", "SW"]) and not line.startswith("#"):
                if line:
                    line = line.lstrip()
                    fields = line.split()
                    fields[8] = fields[8].replace("C", "-")
                    chr, start, end, strand = fields[4], fields[5], fields[6], fields[8]
                    ID = ";".join(fields)
                    start = str(int(start) + 1)
                    bed.write(f"{chr}\t{start}\t{end}\t{ID}\t.\t{strand}\n")


def fai_to_bed(fai_name):
    fai = pd.read_table(fai_name, header=None)
    fai.columns = ['chr', 'size', 'start', 'line', 'width']
    fai['start'] = 1
    fai['end'] = fai['size']
    bed = fai[['chr', 'start', 'end']]
    return bed


def tree_anno(fai_name):
    anno = pd.read_table(fai_name, header=None)
    anno.columns = ['ID', 'length', 'start', 'line', 'width']
    anno['Species'] = anno['ID'].str.split('_').str[0]
    anno.drop(['length', 'start', 'line', 'width'], axis=1, inplace=True)
    return anno


def main():
    if sys.argv[1] == 'RM':
        RMout_to_bed(sys.argv[2])
        bed = fai_to_bed(sys.argv[3])
        bed.to_csv(sys.argv[4], sep='\t', header=False, index=False)
    elif sys.argv[1] == 'tree':
        fai = sys.argv[2]
        anno = tree_anno(fai)
        anno.to_csv(sys.argv[3], sep='\t', header=True, index=False)
    else:
        fai = sys.argv[1]
        cu_code = sys.argv[2]
        output_name = sys.argv[3]
        TE_length = mod_clib(fai, cu_code)  # fai to classified length
        # show sum
        pd.set_option('display.max_rows', None)
        TE_length.to_csv(output_name, sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()
