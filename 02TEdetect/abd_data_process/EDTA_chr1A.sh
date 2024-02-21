export 1LANGUAGE = "en_US.UTF-8"
export LANGUAGE="en_US.UTF-8"
export LC_ALL="en_US.UTF-8"
export LC_CTYPE="UTF-8"
export LANG="en_US.UTF-8"

# test 1
perl ~/transposon/tools/EDTA/EDTA.pl --genome $1 --curatedlib ../00data/trep-db_complete_Rel-19.fasta --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 25 > test.log 2>&1 &
## 2>&1 & 的意思是：将标准错误重定向到标准输出，并且将命令放到后台执行