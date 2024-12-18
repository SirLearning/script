for i in `cat fqlist.txt`
do
bwa mem -t 20 -R '@RG\tID:Aegilops\tPL:illumina\tSM:Aegilops' /data1/publicData/wheat/reference/v1.0/D/bwaLib/d_iwgscV1.fa.gz ${i}_1.fq.gz ${i}_2.fq.gz | samtools view -S -b -> ${i}.bam
samtools sort -n -m 4G -@ 20 -o ${i}.namesort.bam -O bam ${i}.bam \
 && samtools fixmate -@ 20 -m ${i}.namesort.bam ${i}.fixmate.bam \
 && samtools sort -m 4G -@ 20 -o ${i}.fixmate.pos.bam -O bam ${i}.fixmate.bam \
 && rm -f ${i}.namesort.bam \
 && samtools markdup -@ 20 -r ${i}.fixmate.pos.bam ${i}.rmdup.bam \
 && rm -f ${i}.fixmate.bam \
 && rm -f ${i}.fixmate.pos.bam
done

# can use 30 threads and run 5 samples at the same time