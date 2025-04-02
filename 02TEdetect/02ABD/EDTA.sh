export LANGUAGE="en_US.UTF-8"
export LC_ALL="en_US.UTF-8"
export LC_CTYPE="UTF-8"
export LANG="en_US.UTF-8"

## $1: {seq}
## $2: thread


# run EDTA
## 2>&1 & 的意思是：将标准错误重定向到标准输出，并且将命令放到后台执行
# perl ~/transposon/tools/EDTA/EDTA.pl --genome $1 --curatedlib ~/transposon/ABD/00data/00library/curated_lib-db_complete_Rel-19.fasta --overwrite 1 --sensitive 1 --stats 1 --evaluate 1 --threads $2 > test.log 2>&1 &

# 继续运行已经断开的程序，从上次断开的地方开始运行
perl ~/transposon/tools/EDTA/EDTA.pl \
--genome $1 \
--curatedlib ~/transposon/ABD/00data/00library/trep-db_complete_Rel-19.fasta \
--overwrite 0 \
--sensitive 1 \
--anno 1 \
--threads $2 > test.log 2>&1 &

nohup perl /data/home/dazheng/software/EDTA/EDTA.pl --genome $1 --overwrite 0 --sensitive 1 --anno 1 --threads $2 --curatedlib /data/home/dazheng/transposon/00library/00curatedLib/$3.fa > test.log 2>&1 &

