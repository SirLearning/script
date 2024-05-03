# reads download
fastq-dump --gzip --split-3 --defline-qual '+' --defline-seq '@\$ac-\$si/\$ri'   SRR_ID

samtools view filename.bam | less -S # view the bam file

prefetch -X 200G -O ./ SRR_ID # download the SRR_ID file
fasterq-dump --split-3 --defline-seq '@$sn[_$rn]/$ri' --defline-qual '+' SRR_ID # split the SRR_ID file
nohup fasterq-dump -p -e 20 --split-3 -O ./ SRR7164628/SRR7164628.sra

#nohup fasterq-dump -e 20 --split-3 -O ./ SRR7164604/SRR7164604.sra && pigz -3 SRR7164604_1.fastq.gz SRR7164604_1.fastq.gz && rm SRR7164604/SRR7164604.sra &
#nohup fasterq-dump -e 20 --split-3 -O ./ SRR7164620/SRR7164620.sra && pigz -3 SRR7164620_1.fastq.gz SRR7164620_1.fastq.gz && rm SRR7164620/SRR7164620.sra &
#nohup fasterq-dump -e 20 --split-3 -O ./ SRR7164628/SRR7164628.sra && pigz -3 SRR7164628_1.fastq SRR7164628_2.fastq && rm SRR7164628/SRR7164628.sra &
#nohup fasterq-dump -e 20 --split-3 -O ./ SRR7164669/SRR7164669.sra && pigz -3 SRR7164669_1.fastq SRR7164669_2.fastq && rm SRR7164669/SRR7164669.sra &
#nohup fasterq-dump -e 20 --split-3 -O ./ SRR7164670/SRR7164670.sra && pigz -3 SRR7164670_1.fastq SRR7164670_2.fastq && rm SRR7164670/SRR7164670.sra &

java -cp ~/js/pop/src/ ssReads 1 SRR7164604/SRR7164604_1.fastq.gz f1.SRR7164604.fq SRR7164604/SRR7164604_2.fastq.gz r2.SRR7164604.fq 906 && pigz -3 f1.SRR7164604.fq r2.SRR7164604.fq
java -cp ~/js/pop/src/ ssReads 1 SRR7164620/SRR7164620_1.fastq.gz f1.SRR7164620.fq SRR7164620/SRR7164620_2.fastq.gz r2.SRR7164620.fq 878 && pigz -3 f1.SRR7164620.fq r2.SRR7164620.fq
java -cp ~/js/pop/src/ ssReads 1 SRR7164628/SRR7164628_1.fastq.gz f1.SRR7164628.fq SRR7164628/SRR7164628_2.fastq.gz r2.SRR7164628.fq 905 && pigz -3 f1.SRR7164628.fq r2.SRR7164628.fq
java -cp ~/js/pop/src/ ssReads 1 SRR7164669/SRR7164669_1.fastq.gz f1.SRR7164669.fq SRR7164669/SRR7164669_2.fastq.gz r2.SRR7164669.fq 902 && pigz -3 f1.SRR7164669.fq r2.SRR7164669.fq
java -cp ~/js/pop/src/ ssReads 1 SRR7164670/SRR7164670_1.fastq.gz f1.SRR7164670.fq SRR7164670/SRR7164670_2.fastq.gz r2.SRR7164670.fq 917 && pigz -3 f1.SRR7164670.fq r2.SRR7164670.fq

10000/11.05259133 =
10000/11.39205023 =
10000/11.0573459 =
10000/11.09383727 =
10000/10.90905583 =
