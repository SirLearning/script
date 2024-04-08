awk -F'\t' '{print $1}' vu.fai > name
seqkit grep -f name 1A.lib.fa > 1A.known.fa