awk '$3~"LTR_retrotransposon"' \
  $1 \
  | grep "Method=structural" > intact.LTR.gff3

# 有strand信息的LTR 23.EDTA/Aly/genome.fa.mod.EDTA.TEanno.gff3
awk -F';' '{print $1' intact.LTR.gff3\
  | sed 's/ID=//' | awk '$7="?"' \
  | awk '{print $1"\t"$4-1"\t"$5"\t"$9"\t.\t"$7}' >> intact.LTR.bed6

# 无strand信息的LTR
awk -F';' '{print $1' intact.LTR.gff3\
  | sed 's/ID=//' | awk '$7="?"' \
  | awk '{print $1"\t"$4-1"\t"$5"\t"$9"\t.\t+"}' >> intact.LTR.bed6

# 根据gff提取序列 12.ref/Aly.genome.fa
bedtools getfasta -fi $2 \
  -bed intact.LTR.bed6 -fo - -name -s \
  | awk -F'::' '{print $1}' | seqtk seq -l 60 - > intact.LTR.fa

# TEsorter进一步分类
TEsorter -db rexdb intact.LTR.fa