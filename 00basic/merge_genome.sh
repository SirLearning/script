# chr6A for example
line_num=$(grep -n '>32' chr6A.fa | cut -d : -f 1)
sed -i "${line_num}d" chr6A.fa
sed -i 's/>31/>chr6A/g' chr6A.fa
awk '/^>/&&NR>1{print "";}{printf "%s",/^>/?$0"\n":$0}' chr6A.fa > chr6A1.fa
seqkit seq chr6A1.fa > chr6A.fa
rm chr6A1.fa
samtools faidx chr6A.fa