# process the N_cut annotation
python ~/script/04analysis/03n_site_redirect.py N1.bed
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
awk '$10=="DHH"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.DHH.density.txt
awk '$10=="DTA"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.DTA.density.txt
awk '$10=="DTC"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.DTC.density.txt
awk '$10=="DTH"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.DTH.density.txt
awk '$10=="DTM"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.DTM.density.txt
awk '$10=="DTT"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.DTT.density.txt
awk '$10=="DTX"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.DTX.density.txt
awk '$10=="DXX"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.DXX.density.txt
awk '$10=="NULL"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.NULL.density.txt
awk '$10=="RIJ"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.RIJ.density.txt
awk '$10=="RIL"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.RIL.density.txt
awk '$10=="RIR"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.RIR.density.txt
awk '$10=="RIX"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.RIX.density.txt
awk '$10=="RLC"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.RLC.density.txt
awk '$10=="RLG"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.RLG.density.txt
awk '$10=="RLX"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.RLX.density.txt
awk '$10=="RSX"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.RSX.density.txt
awk '$10=="XXX"' mc.gff3 | bedtools coverage -a chr1A.bed6 -b - -counts -F 0.5 > chr1A.XXX.density.txt