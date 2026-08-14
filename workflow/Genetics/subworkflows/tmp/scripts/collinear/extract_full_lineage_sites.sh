#!/usr/bin/env bash
# Emit full-lineage syntenic site coordinates from one matrix (gzip TSV).
# Usage: extract_full_lineage_sites.sh <chr_label> <matrix.tsv.gz> <out.tsv.gz>
set -euo pipefail

chr_label=$1
matrix=$2
out=$3
sg=${chr_label: -1}

case "$sg" in
  A) lin='traesA trboe trmon trura trdicA trdurA' ;;
  B) lin='traesB trdicB trdurB aespe' ;;
  D) lin='traesD aetau aecom aeuni aelon aebic' ;;
  *) echo "bad chr $chr_label" >&2; exit 1 ;;
esac

tmp=$(mktemp)
{
  echo -e "chr_label\tref_pos_1based"
  pigz -dc -- "$matrix" | mawk -F'\t' -v chr="$chr_label" -v lincols="$lin" '
  BEGIN { nl = split(lincols, L, " ") }
  NR == 1 {
    for (i = 1; i <= NF; i++) h[$i] = i
    for (i = 1; i <= nl; i++) {
      if (!(L[i] in h)) { print "Missing column " L[i] > "/dev/stderr"; exit 2 }
      li[i] = h[L[i]]
    }
    next
  }
  {
    ok = 1
    for (i = 1; i <= nl; i++) if ($(li[i]) == "-") { ok = 0; break }
    if (ok) print chr "\t" ($(h["ref_pos"]) + 1)
  }
  '
} | pigz -c >"$tmp"
mv "$tmp" "$out"
echo "DONE $chr_label -> $out ($(pigz -dc "$out" | tail -n +2 | wc -l) sites)"
