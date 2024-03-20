# 1. split chrom version
# step 1
touch stats.anno.gff3
## $1: {seq}.fai
## $2: {seq}.mod.EDTA.TEanno.gff3
python ~/script/04analysis/09site_redirect_chr.py $1 $2 stats.anno.gff3
awk '{$1="chr1D"}1' OFS='\t' stats.anno.gff3 > stats.anno.mod.gff3

## 2. N version
## step 1 redirect N_cut annotation
## output file
#touch stats.anno.gff3
### $1: {seq}.mod.EDTA.TEanno.gff3
#python ~/script/04analysis/05site_redirect_N.py ~/script/04analysis/data/N1.bed $1 stats.anno.gff3
#awk '{$1="chr1A"}1' OFS='\t' stats.anno.gff3 > stats.anno.mod.gff3