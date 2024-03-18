import sys

import pandas as pd

intactLTR_name = 'data/02intact.LTR.gff3'  # sys.argv[1]
intactLTR_age = pd.read_table('data/chr1A.fa.mod.pass.list', header=None)  # sys.argv[2]

