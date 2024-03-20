import sys
import pandas as pd

# fai_name = 'data/cp.chr1B.fa.fai'
# anno_name = 'data/test.gff3'
# output_name = 'data/stats.anno.gff3'

fai_name = sys.argv[1]
anno_name = sys.argv[2]
output_name = sys.argv[3]

fai = pd.read_table(fai_name, sep='\t', header=None)
fai.columns = ['seq', 'length', 'start', 'lb', 'lw']

with open(anno_name, 'r') as anno_file, open(output_name, 'w') as output_file:
    for line in anno_file:
        element = line.split('\t')
        chrom = element[0]
        start = int(element[3])
        end = int(element[4])
        if chrom == fai['seq'][1]:
            element[3] = str(start + int(fai['end'][1]))
            element[4] = str(end + int(fai['end'][1]))
        output_file.write('\t'.join(element))
