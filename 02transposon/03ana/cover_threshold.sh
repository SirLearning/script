bedtools intersect -f 1 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.95 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.9 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.85 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.8 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.75 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.7 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.65 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.6 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.55 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.5 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.45 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.4 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.35 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.3 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.25 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.2 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.15 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.1 -r -v -a $1 -b $2 | wc -l >> $3
bedtools intersect -f 0.05 -r -v -a $1 -b $2 | wc -l >> $3

bedtools intersect -f 1 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.95 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.9 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.85 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.8 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.75 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.7 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.65 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.6 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.55 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.5 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.45 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.4 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.35 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.3 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.25 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.2 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.15 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.1 -r -v -a $2 -b $1 | wc -l >> $4
bedtools intersect -f 0.05 -r -v -a $2 -b $1 | wc -l >> $4