# scan the number of '=' char in the file
awk '{ print gsub(/=/,"") }' 01cs.anno.gff3 > s.txt
# list the line number of specific char
grep -n '6' s.txt