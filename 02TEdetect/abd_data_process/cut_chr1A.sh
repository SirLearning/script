# test 2
# cat chr1A.fa | seqkit subseq -r 1:1000000 > chr1A1M.fa

# cut fasta into different length
faToTwoBit $1 genome.2bit
echo "1"
twoBitInfo -nBed genome.2bit N.bed
echo "2"
awk '{print $1 "\t" $2}' $1.fai > $1.genome
bedtools complement -i N.bed -g $1.genome | bedtools getfasta -fo chr1AN.fa -fi $1 -bed -