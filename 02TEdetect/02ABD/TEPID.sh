# $1 test number
# $2 cs region
# $3 kg region

# 0. 20M test
seqkit subseq -r $2 ../00map/cs.1A.fa -o cs.$1.fa
seqkit subseq -r $3 ../00map/kg.1A.fa -o kg.$1.fa
iss generate -p 10 -g kg.$1.fa -n 2666667  --mode basic -o ./kg.$1

# 1. tepid-map
## preparing
bowtie2-build --threads 20 cs.$1.fa cs.$1
yaha -t 10 -g cs.$1.fa
## mapping
#nohup tepid-map -x /data1/home/dazheng/transposon/load/01test/abd.cs.fa  -p 20 -s 2000 -n CRR072401 -1 f1.CRR072401.fq.gz -2 r2.CRR072401.fq.gz &
tepid-map -x ./cs.$1 -y ./cs.$1.X*  -p 20 -s 2000 -n kg.$1 -1 kg.$1_R1.fastq -2 kg.$1_R2.fastq
#samtools index -c CRR072401.bam # csi index

# 2. tepid-discover
#nohup tepid-discover -k -p 20 -n CRR072401 -c CRR072401.bam -s CRR072401.split.bam -t subD.cs.gff3 &
gff2bed < cs.$1.gff3 > cs.$1.bed
tepid-discover -p 20 -n kg.$1  -c kg.$1.bam -s kg.$1.split.bam -t cs.$1.bed
