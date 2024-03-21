import sys

import pandas as pd

# anno_name = 'data/test.ln.gff3'
# output_name = 'data/test.gff3'

anno_name = sys.argv[1]
output_name = sys.argv[2]

anno = pd.read_table(anno_name, sep='\t', header=None)
anno.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes', 'Classification']
anno['width'] = anno['end'] - anno['start'] + 1

anno['overlap'] = 0

fm_end = 0
fm_class = ''
fm_width = 0
# flexible in 80% of the width
for i in range(0, len(anno)):
    # if (fm_end - anno['start'][i]) > (fm_width + anno['width'][i]) / 10 and anno['Classification'][i] == fm_class:
    if (fm_end - anno['start'][i]) > (fm_width + anno['width'][i]) / 10:
        anno.loc[i, 'overlap'] = 1
    fm_end = anno['end'][i]
    fm_class = anno['Classification'][i]
    fm_width = anno['width'][i]
anno.to_csv(output_name, sep='\t', index=False)
