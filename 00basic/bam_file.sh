# reads download
fastq-dump --gzip --split-3 --defline-qual '+' --defline-seq '@\$ac-\$si/\$ri'   SRR_ID

samtools view filename.bam | less -S # view the bam file

prefetch -X 200G -O ./ SRR_ID # download the SRR_ID file
fasterq-dump --split-3 --defline-seq '@$sn[_$rn]/$ri' --defline-qual '+' SRR_ID # split the SRR_ID file
pigz