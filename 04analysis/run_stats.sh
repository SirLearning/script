# step 1 number and length
python 01num_length.py data/01chr1A.fa.mod.EDTA.TEanno.gff3 data/01chr1A.fa.fai data/01chr1A.anno.stats

# step 2 intact LTR
#awk '$3~"LTR_retrotransposon"' \
#  data/01chr1A.fa.mod.EDTA.TEanno.gff3 \
#  | grep "Method=structural" > data/02intact.LTR.gff3
#awk -F';' '{print $1}' data/02intact.LTR.gff3\
#  | sed 's/ID=//' | awk '$7!="?"' \
#  | awk '{print $1"\t"$4-1"\t"$5"\t"$9"\t.\t"$7}' >> data/02intact.LTR.bed6
#awk -F';' '{print $1}' data/02intact.LTR.gff3\
#  | sed 's/ID=//' | awk '$7=="?"' \
#  | awk '{print $1"\t"$4-1"\t"$5"\t"$9"\t.\t+"}' >> data/02intact.LTR.bed6
#bedtools getfasta -fi data/01chr1A.fa \
#  -bed intact.LTR.bed6 -fo - -name -s \
#  | awk -F'::' '{print $1}' | seqtk seq -l 60 - > data/02intact.LTR.fa
#TEsorter -db rexdb data/02intact.LTR.fa

# step 3 LTR insertion time
python 03data_orgnz.py data/01chr1A.fa.mod.EDTA.TEanno.gff3 data/01chr1A.fa.fai