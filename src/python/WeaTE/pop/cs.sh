cat /data1/home/dazheng/transposon/pop/03mapping/$1/CRR072401.mosdepth.summary.txt | grep -f $3-$2 > $1-$2.CRR072401.sum
cat /data1/home/dazheng/transposon/pop/03mapping/$1/SRR7164576.mosdepth.summary.txt | grep -f $3-$2 > $1-$2.SRR7164576.sum
cat /data1/home/dazheng/transposon/pop/03mapping/$1/SRR7164580.mosdepth.summary.txt | grep -f $3-$2 > $1-$2.SRR7164580.sum
cat /data1/home/dazheng/transposon/pop/03mapping/$1/SRR7164572.mosdepth.summary.txt | grep -f $3-$2 > $1-$2.SRR7164572.sum
cat /data1/home/dazheng/transposon/pop/03mapping/$1/SRR7164606.mosdepth.summary.txt | grep -f $3-$2 > $1-$2.SRR7164606.sum

