## bwa-mem2 mem -> samtools view
#name=trboe;
#echo "bwa-mem2 mem -t 72 -R '@RG\tID:${name}_simulated\tPL:illumina\tSM:${name}' abd_iwgscV1.fa.gz ${name}_f1.fq.gz ${name}_r2.fq.gz | samtools view -S -b - > ${name}_raw.bam" > bwa.sh
## samtools sort -> fixmate -> sort -> markdup
#i=psjun;
#nohup samtools sort -n -m 4G -@ 10 -o ${i}.namesort.bam -O bam ${i}_raw.bam && samtools fixmate -@ 10 -m ${i}.namesort.bam ${i}.fixmate.bam && samtools sort -m 4G -@ 10 -o ${i}.fixmate.pos.bam -O bam ${i}.fixmate.bam && rm -f ${i}.namesort.bam && samtools markdup -@ 10 -r ${i}.fixmate.pos.bam ${i}.rmdup.bam && rm -f ${i}.fixmate.bam && rm -f ${i}.fixmate.pos.bam &

# split version
bwa-mem2 index $1.fa
bwa-mem2 mem -t 5 $1.fa $1.sr_R1.fastq $1.sr_R2.fastq | samtools view -S -b - > $1.raw.bam

samtools sort -n -m 4G -@ 5 -o $1.sort.bam -O bam $1.raw.bam
samtools fixmate -@ 5 -m $1.sort.bam $1.fixmate.bam
samtools sort -m 4G -@ 5 -o $1.fixmate.pos.bam -O bam $1.fixmate.bam
samtools markdup -@ 5 -r $1.fixmate.pos.bam $1.rmdup.bam

samtools index $1.rmdup.bam && mosdepth -t 5 -n -Q 20 $1 $1.rmdup.bam
rm $1.raw.bam $1.sort.bam $1.fixmate.bam $1.fixmate.pos.bam
