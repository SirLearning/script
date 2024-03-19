import sys
import pandas as pd

# nn_name = 'data/N1.bed'
# anno_name = 'data/04chr1ANM.anno.gff3'
# output_name = 'data/04chr1ANM.anno.N.gff3'

nn_name = sys.argv[1]
anno_name = sys.argv[2]
output_name = sys.argv[3]

nn = pd.read_table(nn_name, sep='\t', header=None)
nn.columns = ['chrom', 'start', 'end']

with open(anno_name, 'r') as anno_file, open(output_name, 'w') as output_file:
    for line in anno_file:
        element = line.split('\t')
        chrom = element[0]
        start = int(element[3])
        end = int(element[4])
        for i in range(0, len(nn) - 1):
            if chrom == 'chr1A_' + str(i + 2):
                element[3] = str(start + int(nn['end'][i]))
                element[4] = str(end + int(nn['end'][i]))
                output_file.write('\t'.join(element))
                break
        else:
            output_file.write(line)
