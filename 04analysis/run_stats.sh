# step 1 number and length
# output file
touch stats.nl.txt
## $1: {seq}.mod.EDTA.TEanno.gff3
## $2: {seq}.fai
python ~/script/04analysis/01num_length.py $1 $2 stats.nl.txt ~/script/04analysis/data/TEcode
printf "1\n"
# step 2 distribution
#python ~/script/04analysis/04n_site_redirect.py N1.bed # N_cut annotation
# generate windows
awk '{print $1"\t"$2}' $2 > stats.genome
bedtools makewindows -g stats.genome -w 5000000 -s 5000000 > stats.dtb.bed3
awk '{print $1"\t"$2"\t"$3"\t"NR"\t.\t+"}' stats.dtb.bed3 > stats.dtb.bed6
# output file
touch stats.dtb.gff3
python ~/script/04analysis/03distribution.py $1 ~/script/04analysis/data/TEcode stats.dtb.gff3
#awk '{$1="chr1A"}1' OFS='\t' stats.dtb.gff3 > stats.dtb.gff3 # when N1.bed is used, need to reconsider the file name
# calculate the number of TE in each window 2. ../../31.intactLTR/Aly/intact.LTR.gff3
#awk '$10=="DHH"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DHH.dtb.txt
#awk '$10=="DTA"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTA.dtb.txt
#awk '$10=="DTC"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTC.dtb.txt
#awk '$10=="DTH"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTH.dtb.txt
#awk '$10=="DTM"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTM.dtb.txt
#awk '$10=="DTT"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTT.dtb.txt
#awk '$10=="DTX"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTX.dtb.txt
#awk '$10=="DXX"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DXX.dtb.txt
#awk '$10=="NULL"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.NULL.dtb.txt
#awk '$10=="RIJ"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RIJ.dtb.txt
#awk '$10=="RIL"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RIL.dtb.txt
#awk '$10=="RIR"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RIR.dtb.txt
#awk '$10=="RIX"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RIX.dtb.txt
#awk '$10=="RLC"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RLC.dtb.txt
#awk '$10=="RLG"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RLG.dtb.txt
#awk '$10=="RLX"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RLX.dtb.txt
#awk '$10=="RSX"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RSX.dtb.txt
#awk '$10=="XXX"' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.XXX.dtb.txt
# another way
declare -a arr=("DHH" "DTA" "DTC" "DTH" "DTM" "DTT" "DTX" "DXX" "NULL" "RIJ" "RIL" "RIR" "RIX" "RLC" "RLG" "RLX" "RSX" "XXX")
for i in "${arr[@]}"
do
   awk -v var="$i" '$10==var' stats.dtb.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > "stats.${i}.dtb.txt"
done
# step 3 test N_cut
# redirect N_cut annotation
#python ~/script/04analysis/04n_site_redirect.py N1.bed ../chr1ANM.fa.mod.EDTA.TEanno.gff3 chr1ANM.anno.gff3