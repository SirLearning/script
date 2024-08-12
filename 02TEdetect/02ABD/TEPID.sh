# 1. tepid-map
## preparing
bowtie2-build --threads 20 ../00data/abd.cs.fa abd.cs
yaha -g ../00data/abd.cs.fa
## mapping
nohup tepid-map -x /data1/home/dazheng/transposon/load/01test/abd.cs.fa  -p 20 -s 2000 -n CRR072401 -1 f1.CRR072401.fq.gz -2 r2.CRR072401.fq.gz &
samtools index -c CRR072401.bam # csi index
# 2. tepid-discover
nohup tepid-discover -k -p 20 -n CRR072401 -c CRR072401.bam -s CRR072401.split.bam -t subD.cs.gff3 &