#cat $3/$1.mosdepth.summary.txt | grep -f $3.$2 > $3.$1.sum
cat /data1/home/dazheng/transposon/pop/03mapping/$3/$1.mosdepth.summary.txt | grep -f $2-0 > $2-0.$1.sum
cat /data1/home/dazheng/transposon/pop/03mapping/$3/$1.mosdepth.summary.txt | grep -f $2-1 > $2-1.$1.sum
cat /data1/home/dazheng/transposon/pop/03mapping/$3/$1.mosdepth.summary.txt | grep -f $2-2 > $2-2.$1.sum