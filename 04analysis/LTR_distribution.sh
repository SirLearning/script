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
awk '{$1="chr1A"}1' OFS='\t' mc.gff3 > chr1A.mc.gff3
# calculate the number of TE in each window 2. ../../31.intactLTR/Aly/intact.LTR.gff3
awk '$10=="DHH"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.DHH.density.txt
awk '$10=="DTA"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.DTA.density.txt
awk '$10=="DTC"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.DTC.density.txt
awk '$10=="DTH"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.DTH.density.txt
awk '$10=="DTM"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.DTM.density.txt
awk '$10=="DTT"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.DTT.density.txt
awk '$10=="DTX"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.DTX.density.txt
awk '$10=="DXX"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.DXX.density.txt
awk '$10=="NULL"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.NULL.density.txt
awk '$10=="RIJ"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.RIJ.density.txt
awk '$10=="RIL"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.RIL.density.txt
awk '$10=="RIR"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.RIR.density.txt
awk '$10=="RIX"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.RIX.density.txt
awk '$10=="RLC"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.RLC.density.txt
awk '$10=="RLG"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.RLG.density.txt
awk '$10=="RLX"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.RLX.density.txt
awk '$10=="RSX"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.RSX.density.txt
awk '$10=="XXX"' chr1A.mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1ANM.XXX.density.txt

#declare -a arr=("DHH" "DTA" "DTC" "DTH" "DTM" "DTT" "DTX" "DXX" "NULL" "RIJ" "RIL" "RIR" "RIX" "RLC" "RLG" "RLX" "RSX" "XXX")
#
#for i in "${arr[@]}"
#do
#   awk -v var="$i" '$10==var' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > "chr1A.${i}.density.txt"
#done