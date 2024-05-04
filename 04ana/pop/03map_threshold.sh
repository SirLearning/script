# nohup sh ~/script/04ana/pop/03map_threshold.sh 20 1A &
bwa-mem2 mem -t 15 $2.known.fa ../f1.$1.fq.gz ../r2.$1.fq.gz | samtools view -S -b - > $1.raw.bam

#samtools sort -n -m 4G -@ 10 -o $1.sort.bam -O bam $1.raw.bam && rm $1.raw.bam
#samtools fixmate -@ 10 -m $1.sort.bam $1.fixmate.bam && rm $1.sort.bam
#samtools sort -m 4G -@ 10 -o $1.fixmate.pos.bam -O bam $1.fixmate.bam && rm $1.fixmate.bam
#samtools markdup -@ 10 -r $1.fixmate.pos.bam $1.rmdup.bam && rm $1.fixmate.pos.bam
#samtools index $1.rmdup.bam && mosdepth -t 10 -n -Q 20 $1 $1.rmdup.bam
#samtools flagstat $1.rmdup.bam >

samtools sort -n -m 4G -@ 10 -o $1.sort.bam -O bam $1.raw.bam && samtools flagstat $1.raw.bam > t1.$1.stats.txt
samtools fixmate -@ 10 -m $1.sort.bam $1.fixmate.bam && samtools flagstat $1.sort.bam > t2.$1.stats.txt
samtools sort -m 4G -@ 10 -o $1.fixmate.pos.bam -O bam $1.fixmate.bam && samtools flagstat $1.fixmate.bam > t3.$1.stats.txt
samtools markdup -@ 10 -r $1.fixmate.pos.bam $1.rmdup.bam && samtools flagstat $1.rmdup.bam > t4.$1.stats.txt
samtools index $1.raw.bam && mosdepth -t 10 -n -Q 20 $1.raw $1.raw.bam
samtools index $1.sort.bam && mosdepth -t 10 -n -Q 20 $1.sort $1.sort.bam
samtools index $1.fixmate.bam && mosdepth -t 10 -n -Q 20 $1.fixmate $1.fixmate.bam
samtools index $1.rmdup.bam && mosdepth -t 10 -n -Q 20 $1.rmdup $1.rmdup.bam

# get mapping rate
#find . -name "*stats.txt" -exec awk 'FNR==8{print FILENAME ": " $0}' {} \;

#mv CRR072247 CRR072248 CRR072249 CRR072250 CRR072251 /data1/home/dazheng/transposon/pop/00data/01WEW
#mv CRR072337 CRR072338 CRR072340 CRR072341 CRR072347 /data1/home/dazheng/transposon/pop/00data/02DW
#mv CRR072405 CRR072406 CRR072407 CRR072408 CRR072409 /data1/home/dazheng/transposon/pop/00data/03AT