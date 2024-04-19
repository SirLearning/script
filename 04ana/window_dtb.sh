## $1: stats.dtb.gff3
## $2: {seq}.fai
## $3: window size
awk '{print $1"\t"$2}' $2 > stats.genome
bedtools makewindows -g stats.genome -w $3 -s $3 > stats.dtb.bed3
awk '{print $1"\t"$2"\t"$3"\t"NR"\t.\t+"}' stats.dtb.bed3 > stats.dtb.bed6

awk -F'\t' '$10=="DHH"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DHH.dtb.txt
awk -F'\t' '$10=="DTA"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTA.dtb.txt
awk -F'\t' '$10=="DTC"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTC.dtb.txt
awk -F'\t' '$10=="DTH"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTH.dtb.txt
awk -F'\t' '$10=="DTM"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTM.dtb.txt
awk -F'\t' '$10=="DTT"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTT.dtb.txt
awk -F'\t' '$10=="DTX"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTX.dtb.txt
awk -F'\t' '$10=="DXX"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DXX.dtb.txt
awk -F'\t' '$10=="RIJ"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RIJ.dtb.txt
awk -F'\t' '$10=="RIL"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RIL.dtb.txt
awk -F'\t' '$10=="RIR"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RIR.dtb.txt
awk -F'\t' '$10=="RIX"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RIX.dtb.txt
awk -F'\t' '$10=="RLC"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RLC.dtb.txt
awk -F'\t' '$10=="RLG"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RLG.dtb.txt
awk -F'\t' '$10=="RLX"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RLX.dtb.txt
awk -F'\t' '$10=="RSX"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RSX.dtb.txt
awk -F'\t' '$10=="XXX"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.XXX.dtb.txt
awk -F'\t' '$10=="NULL"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.NULL.dtb.txt