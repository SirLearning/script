# generate chr1A1M.fa
cat $1 | seqkit subseq -r 1:1000000 > chr1A1M.fa

# cut fasta into different length, generate chr1AN.fa
faToTwoBit $1 genome.2bit
echo "1"
twoBitInfo -nBed genome.2bit N.bed
echo "2"
samtools faidx $1
# cut chr1A.fa into 1Mbp by different length of N, generate chr1ANM.fa
python ~/script/02TEdetect/02ABD/01cut_by_N/02cut_Nbed.py N.bed N1.bed $2
awk '{print $1 "\t" $2}' $1.fai > $1.genome
bedtools complement -i N1.bed -g $1.genome | bedtools getfasta -fo chr1AN.fa -fi $1 -bed -