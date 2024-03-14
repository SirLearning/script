# retrieve RT domain sequences
grep "Ty1_copia" ../../31.intactLTR/Aly/intact.LTR.fa.rexdb.dom.tsv \
  | grep -P "\-RT\t" \
  | cut -f 1 \
  | seqtk subseq ../../31.intactLTR/Aly/intact.LTR.fa.rexdb.dom.faa \
  - > Aly.RT.raw.fa

grep "Ty1_copia" ../../31.intactLTR/Ath/intact.LTR.fa.rexdb.dom.tsv \
  | grep -P "\-RT\t" \
  | cut -f 1 \
  | seqtk subseq ../../31.intactLTR/Ath/intact.LTR.fa.rexdb.dom.faa \
  - > Ath.RT.raw.fa

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