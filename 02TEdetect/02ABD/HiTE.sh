#python main.py  --genome ${genome}  --thread ${thread}  --outdir ${output_dir}
# --genome /data1/home/dazheng/transposon/ABD/00data/05refer/a.cs.fa \
# --thread 20 \
# --outdir /data1/home/dazheng/transposon/ABD/09HiTE/01A \

cd /data1/home/dazheng/transposon/tools/HiTE
# $1: genome file
# $2: output directory
# $3: thread
python main.py \
 --genome $1 \
 --outdir $2 \
 --thread $3 \
 --domain 1 \
 --recover 1 \
 --annotate 1 \
 --intact_anno 1 \
 --BM_EDTA 1 \
 --EDTA_home /data1/home/dazheng/transposon/tools/EDTA \
 --coverage_threshold 0.8 \
 --is_wicker 1
