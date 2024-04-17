
bwa-mem2 mem -t 15 $2.known.fa ~/transposon/ABD/00data/03reads/f1.$1.fq ~/transposon/ABD/00data/03reads/r2.$1.fq | samtools view -S -b - > $1.raw.bam

samtools sort -n -m 4G -@ 10 -o $1.sort.bam -O bam $1.raw.bam && rm $1.raw.bam
samtools fixmate -@ 10 -m $1.sort.bam $1.fixmate.bam && rm $1.sort.bam
samtools sort -m 4G -@ 10 -o $1.fixmate.pos.bam -O bam $1.fixmate.bam && rm $1.fixmate.bam
samtools markdup -@ 10 -r $1.fixmate.pos.bam $1.rmdup.bam && $1.fixmate.pos.bam
samtools index $1.rmdup.bam && mosdepth -t 10 -n -Q 20 $1 $1.rmdup.bam
samtools flagstat $1.rmdup.bam > $1.stats.txt