# retrieve RT domain sequences
# ../../31.intactLTR/Aly/intact.LTR.fa.rexdb.dom.tsv
# ../../31.intactLTR/Aly/intact.LTR.fa.rexdb.dom.faa
grep "Ty1_copia" $1 \
  | grep -P "\-RT\t" \
  | cut -f 1 \
  | seqtk subseq $2 \
  - > dom.RT.raw.fa

# simplify the ID
awk -F':' '{print $1}' Aly.RT.raw.fa \
  | sed 's/\Class_I\/LTR\///' \
  | sed 'ss/>/>Aly-/' > Aly.RT.fa
awk -F':' '{print $1}' Ath.RT.raw.fa \
  | sed 's/\Class_I\/LTR\///' \
  | sed 'ss/>/>Ath-/' > Ath.RT.fa
# combine the RT sequences
cat Aly.RT.fa Ath.RT.fa > RT.fa

# align the RT sequences
mafft --auto RT.fa > RT.aln

# build phylogenetic tree using IQ-TREE
iqtree -s RT.aln -bb 1000 \
    -pre Ty1_copia.iqtree 1>iqtree.log 2>&1 &
# polish the tree using iTOL: https://itol.embl.de/