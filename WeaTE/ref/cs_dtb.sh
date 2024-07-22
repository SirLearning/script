## $1: {seq}.fai
## $2: {seq}.mod.EDTA.TEanno.gff3
## $3: window size 5000000
touch mod.cs.gff3
python ~/script/WeaTE/ref/mod_anno.py $2 mod.cs.gff3 cs_v2
awk '{print $1"\t"$2}' $1 > stats.genome
bedtools makewindows -g stats.genome -w $3 -s $3 > stats.dtb.bed3
awk '{print $1"\t"$2"\t"$3"\t"NR"\t.\t+"}' stats.dtb.bed3 > stats.dtb.bed6

awk -F'\t' '$10=="DHH"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.DHH.dtb.txt
awk -F'\t' '$10=="DTA"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.DTA.dtb.txt
awk -F'\t' '$10=="DTC"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.DTC.dtb.txt
awk -F'\t' '$10=="DTH"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.DTH.dtb.txt
awk -F'\t' '$10=="DTM"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.DTM.dtb.txt
awk -F'\t' '$10=="DTT"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.DTT.dtb.txt
awk -F'\t' '$10=="RIJ"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.RIJ.dtb.txt
awk -F'\t' '$10=="RIL"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.RIL.dtb.txt
awk -F'\t' '$10=="RIR"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.RIR.dtb.txt
awk -F'\t' '$10=="RLC"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.RLC.dtb.txt
awk -F'\t' '$10=="RLG"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.RLG.dtb.txt
awk -F'\t' '$10=="RSX"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.RSX.dtb.txt

#awk -F'\t' '$16=="DHH"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > coverage.DHH.dtb.txt
#awk -F'\t' '$16=="DTA"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > coverage.DTA.dtb.txt
#awk -F'\t' '$16=="DTC"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > coverage.DTC.dtb.txt
#awk -F'\t' '$16=="DTH"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > coverage.DTH.dtb.txt
#awk -F'\t' '$16=="DTM"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > coverage.DTM.dtb.txt
#awk -F'\t' '$16=="DTT"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > coverage.DTT.dtb.txt
#awk -F'\t' '$16=="RIJ"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > coverage.RIJ.dtb.txt
#awk -F'\t' '$16=="RIL"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > coverage.RIL.dtb.txt
#awk -F'\t' '$16=="RIR"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > coverage.RIR.dtb.txt
#awk -F'\t' '$16=="RLC"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > coverage.RLC.dtb.txt
#awk -F'\t' '$16=="RLG"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > coverage.RLG.dtb.txt
#awk -F'\t' '$16=="RSX"' mod.cs.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > coverage.RSX.dtb.txt


