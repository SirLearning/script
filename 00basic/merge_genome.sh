# chr1D for example
line_num=$(grep -n '>6' chr1D.fa | cut -d : -f 1)
sed -i "${line_num}d" chr1D.fa
sed -i 's/>5/>chr1D/g' chr1D.fa
awk '/^>/&&NR>1{print "";}{printf "%s",/^>/?$0"\n":$0}' chr1D.fa > chr1D1.fa
seqkit seq chr1D1.fa > chr1D.fa
rm chr1D1.fa
samtools faidx chr1D.fa