# to h
RepeatMasker -e ncbi -pa 10 -q -no_is -norna -nolow -div 20 -lib $2.fa $1.fa # 80% same
perl ~/transposon/tools/parseRM/RMout_to_bed.pl $1.fa.out base0
rm *.out *.cat *.masked *.tbl
touch $1.bed
samtools faidx $1.fa
python ~/script/WeaTE/pop/derived.py $1.fa.fai $1.bed
bedtools intersect -f 0.8 -r -v -a $1.bed -b $1.fa.out.bed > $1-$2.bed
rm $1.bed $1.fa.out.bed
# to r

