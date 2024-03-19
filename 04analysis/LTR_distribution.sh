# process the N_cut annotation
#python ~/script/04analysis/04n_site_redirect.py N1.bed
# generate windows 1. ../../12.ref/Aly.genome.fa.fai
awk '{print $1"\t"$2}' $1 > genome.txt
#   generate bed3 and bed6
bedtools makewindows -g genome.txt -w 5000000 -s 5000000 > chr1A.bed3
awk '{print $1"\t"$2"\t"$3"\t"NR"\t.\t+"}' chr1A.bed3 > chr1A.bed6
# calculate the number of genes in each window
#awk '$3=="gene"' ../../12.ref/Aly.genome.gtf | bedtools coverage -a chr_sw.bed6 \
#  -b - -counts -F 0.5 > chr_sw.gene.density.txt
python ~/script/04analysis/03distribution.py $2 ~/script/04analysis/data/TEcode mc.gff3
# calculate the number of TE in each window 2. ../../31.intactLTR/Aly/intact.LTR.gff3
declare -a arr=("DHH" "DTA" "DTC" "DTH" "DTM" "DTT" "DTX" "DXX" "NULL" "RIJ" "RIL" "RIR" "RIX" "RLC" "RLG" "RLX" "RSX" "XXX")
for i in "${arr[@]}"
do
   awk -v var="$i" '$10==var' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > "chr1ANM.${i}.density.txt"
done