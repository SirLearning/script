export 1LANGUAGE = "en_US.UTF-8"
export LANGUAGE="en_US.UTF-8"
export LC_ALL="en_US.UTF-8"
export LC_CTYPE="UTF-8"
export LANG="en_US.UTF-8"
# test 1
perl ~/transposon/tools/EDTA/EDTA.pl --genome chr1A.fa --curatedlib ../00data/trep-db_complete_Rel-19.fasta --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 25 > test.log 2>&1 &
# test 2
cat chr1A.fa | seqkit subseq -r 1:1000000 > chr1A1M.fa
## cut fasta into different length
faToTwoBit $workPath/genome.fa $workPath/genome.2bit
twoBitInfo -nBed $workPath/genome.2bit $workPath/N.bed

bedtools complement -i $workPath/N.bed -g $workPath/genome.genome | bedtools getfasta -fo $workPath/genome.temp -fi $workPath/genome.fa -bed -

# embl to fasta
