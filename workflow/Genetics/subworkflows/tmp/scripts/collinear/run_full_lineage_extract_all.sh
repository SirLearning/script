#!/usr/bin/env bash
# Parallel extract of full-lineage site lists for all syntenic matrices.
set -euo pipefail

MATRIX_DIR=${MATRIX_DIR:-/data/home/tusr1/01projects/vmap4/00data/06maf/syntenic_sites_k2}
OUT_DIR=${OUT_DIR:?OUT_DIR required}
WORKERS=${WORKERS:-7}
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

mkdir -p "$OUT_DIR"

echo "matrix_dir=$MATRIX_DIR out_dir=$OUT_DIR workers=$WORKERS" >&2
mapfile -t MATS < <(ls -1 "$MATRIX_DIR"/*.syntenic.matrix.tsv.gz | sort)
echo "n_matrices=${#MATS[@]}" >&2

printf '%s\n' "${MATS[@]}" | xargs -P "$WORKERS" -I{} bash -c '
  m="$1"
  script="$2"
  outdir="$3"
  base=$(basename "$m")
  chr=${base%%.*}
  out="$outdir/${chr}.tsv.gz"
  if [[ -s "$out" ]]; then
    echo "SKIP $chr (exists)"
    exit 0
  fi
  "$script" "$chr" "$m" "$out"
' _ {} "$SCRIPT_DIR/extract_full_lineage_sites.sh" "$OUT_DIR"

echo "ALL_DONE full_lineage_sites"
