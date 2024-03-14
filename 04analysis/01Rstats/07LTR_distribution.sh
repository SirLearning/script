# fai index to extract genome size
# generate windows
samtools faidx ../../12.ref/Aly.genome.fa.
awk '{print $1"\t"$2}' ../../12.ref/Aly.genome.fa.fai > genome.txt
#   generate bed3 and bed6
bedtools makewindows -g genome.txt -w 1000000 -s 200000 > chr_sw.bed3
awk '{print $1"\t"$2"\t"$3"\t"NR"\t.\t+"}' chr_sw.bed3 > chr_sw.bed6
# calculate the number of genes in each window
awk '$3=="gene"' ../../12.ref/Aly.genome.gtf | bedtools coverage -a chr_sw.bed6 \
  -b - -counts -F 0.5 > chr_sw.gene.density.txt
# calculate the number of TE in each window
awk '$3=="Copia_LTR_retrotransposon"' ../../31.intactLTR/Aly/intact.LTR.gff3 \
  | bedtools coverage -a chr_sw.bed6 -b - -counts -F 0.5 > chr_sw.Copia.density.txt
awk '$3=="Gypsy_LTR_retrotransposon"' ../../31.intactLTR/Aly/intact.LTR.gff3 \
  | bedtools coverage -a chr_sw.bed6 -b - -counts -F 0.5 > chr_sw.Gypsy.density.txt