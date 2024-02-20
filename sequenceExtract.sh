# move chr1A to certain file 
cat iwgsc_refseqv1.0_TransposableElements_2017Mar13.gff3 | grep -w 'chr1A' > chr1A/chr1A.gff3
seqkit grep -p 1 abd_iwgscV1.fa -o chr1A/chr1A1.fa

seqkit stats 

# change .gff file to .bed file by awk
awk -F'\t' '$3 == "repeat_region" {split($9, a, ";"); for (i in a) {split(a[i], b, "="); if (b[1] == "ID") TE_id = b[2];}  print $1, $4-1, $5, TE_id, ".", $7}' chr1A.gff3 > annotations.bed

gff2bed <JM44.repeat.masked.gff > annotation.bed
bedtools getfasta -s -fi JM44.repeat.masked.fasta -bed annotation.bed -fo teseq.fasta
seqkit stats teseq.fasta
python DataProcessing.py > call.txt
bedtools merge -i annotation.bed > merge.bed

# seqkit to parent
cat annotation.bed | grep -w 'Parent' > parent.bed
bedtools getfasta -s -fi 01data/JM44.repeat.masked.fasta -bed new.bed -fo parent.fasta
cat annotation.bed | grep -v 'Parent' > vParent.bed
bedtools getfasta -s -fi 01data/JM44.repeat.masked.fasta -bed new.bed -fo vParent.fasta
seqkit stats 02result/teSeq.fasta parent.fasta vParent.fasta
cat 01data/
seqkit fx2tab -nliH parent.fasta -o new.txt
python DataProcessing.py > table.txt