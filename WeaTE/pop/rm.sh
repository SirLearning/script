# to h
RepeatMasker -e ncbi -pa 10 -q -no_is -norna -nolow -div 20 -lib cs.fa la.fa # 80% same
perl ~/transposon/tools/parseRM/RMout_to_bed.pl la.fa.out base0
rm *.out *.cat *.masked *.tbl
touch la.bed
samtools faidx la.fa
python ~/script/WeaTE/pop/derived.py la.fa.fai la.bed
bedtools intersect -f 0.8 -r -v -a la.bed -b la.fa.out.bed > la.cs.bed
# to r

