bwa-mem2 mem -t 15 $2.known.fa f1.$1.fq.gz r2.$1.fq.gz | samtools view -S -b - > $1.raw.bam

samtools sort -n -m 4G -@ 10 -o $1.sort.bam -O bam $1.raw.bam && rm $1.raw.bam
samtools fixmate -@ 10 -m $1.sort.bam $1.fixmate.bam && rm $1.sort.bam
samtools sort -m 4G -@ 10 -o $1.fixmate.pos.bam -O bam $1.fixmate.bam && rm $1.fixmate.bam
samtools flagstat $1.fixmate.pos.bam > stats.$1.txt
samtools index $1.fixmate.pos.bam && mosdepth -t 10 -n -Q 20 $1 $1.fixmate.pos.bam