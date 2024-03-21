## $1: stats.dtb.gff3
## $2: {seq}.fai
## $3: window size
awk '{print $1"\t"$2}' $2 > stats.genome
bedtools makewindows -g stats.genome -w $3 -s $3 > stats.dtb.bed3
awk '{print $1"\t"$2"\t"$3"\t"NR"\t.\t+"}' stats.dtb.bed3 > stats.dtb.bed6
declare -a arr=("DHH" "DTA" "DTC" "DTH" "DTM" "DTT" "DTX" "DXX" "NULL" "RIJ" "RIL" "RIR" "RIX" "RLC" "RLG" "RLX" "RSX" "XXX")
for i in "${arr[@]}"
do
   awk -v var="$i" '$10==var' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > "stats.${i}.dtb.txt"
done