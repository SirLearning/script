# gff to bed
gff2bed <JM44.repeat.masked.gff > annotation.bed
awk -F'\t' '$3 == "repeat_region" {split($9, a, ";"); for (i in a) {split(a[i], b, "="); if (b[1] == "ID") TE_id = b[2];}  print $1, $4-1, $5, TE_id, ".", $7}' chr1A.gff3 > annotations.bed
# bed to fasta
bedtools getfasta -s -fi JM44.repeat.masked.fasta -bed annotation.bed -fo teseq.fasta
# embl to fasta
