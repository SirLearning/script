## $1: {seq}.fai
## $2: {seq}.mod.EDTA.TEanno.gff3
## $3: window size 5000000
touch mod.edta.gff3
python ~/script/WeaTE/ref/mod_anno.py $2 mod.edta.gff3
awk '{print $1"\t"$2}' $1 > stats.genome
bedtools makewindows -g stats.genome -w $3 -s $3 > stats.dtb.bed3
awk '{print $1"\t"$2"\t"$3"\t"NR"\t.\t+"}' stats.dtb.bed3 > stats.dtb.bed6
#
#awk -F'\t' '$10=="DHH"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DHH.dtb.txt
#awk -F'\t' '$10=="DTA"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTA.dtb.txt
#awk -F'\t' '$10=="DTC"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTC.dtb.txt
#awk -F'\t' '$10=="DTH"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTH.dtb.txt
#awk -F'\t' '$10=="DTM"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTM.dtb.txt
#awk -F'\t' '$10=="DTT"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTT.dtb.txt
#awk -F'\t' '$10=="DTX"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DTX.dtb.txt
#awk -F'\t' '$10=="DXX"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.DXX.dtb.txt
#awk -F'\t' '$10=="RIJ"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RIJ.dtb.txt
#awk -F'\t' '$10=="RIL"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RIL.dtb.txt
#awk -F'\t' '$10=="RIR"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RIR.dtb.txt
#awk -F'\t' '$10=="RIX"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RIX.dtb.txt
#awk -F'\t' '$10=="RLC"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RLC.dtb.txt
#awk -F'\t' '$10=="RLG"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RLG.dtb.txt
#awk -F'\t' '$10=="RLX"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RLX.dtb.txt
#awk -F'\t' '$10=="RSX"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.RSX.dtb.txt
#awk -F'\t' '$10=="XXX"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.XXX.dtb.txt
#awk -F'\t' '$10=="NULL"' $1 | bedtools coverage -a stats.dtb.bed6 -b - -counts -F 0.5 > stats.NULL.dtb.txt

awk -F'\t' '$11=="DHH"' mod.edta.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.DHH.dtb.txt
awk -F'\t' '$11=="DTA"' mod.edta.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.DTA.dtb.txt
awk -F'\t' '$11=="DTC"' mod.edta.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.DTC.dtb.txt
awk -F'\t' '$11=="DTH"' mod.edta.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.DTH.dtb.txt
awk -F'\t' '$11=="DTM"' mod.edta.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.DTM.dtb.txt
awk -F'\t' '$11=="DTT"' mod.edta.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.DTT.dtb.txt
awk -F'\t' '$11=="RIJ"' mod.edta.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.RIJ.dtb.txt
awk -F'\t' '$11=="RIL"' mod.edta.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.RIL.dtb.txt
awk -F'\t' '$11=="RIR"' mod.edta.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.RIR.dtb.txt
awk -F'\t' '$11=="RLC"' mod.edta.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.RLC.dtb.txt
awk -F'\t' '$11=="RLG"' mod.edta.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.RLG.dtb.tx
awk -F'\t' '$11=="RSX"' mod.edta.gff3 | bedtools coverage -a stats.dtb.bed6 -b - -mean -F 0.5 > coverage.RSX.dtb.txt
