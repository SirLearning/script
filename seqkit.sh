# move chr1A to certain file 
cat iwgsc_refseqv1.0_TransposableElements_2017Mar13.gff3 | grep -w 'chr1A' > chr1A/chr1A.gff3
seqkit grep -p 1 abd_iwgscV1.fa -o chr1A/chr1A1.fa

seqkit stats 
