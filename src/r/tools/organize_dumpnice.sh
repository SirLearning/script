#!/usr/bin/env bash
set -euo pipefail

# Organize and copy R scripts from note/04DumpNice into src/r/dumpnice by category
# Usage: src/r/tools/organize_dumpnice.sh [source_dir] [dest_root]

REPO_ROOT="$(cd "$(dirname "$0")"/../../.. && pwd)"
SOURCE_DIR="${1:-$REPO_ROOT/note/04DumpNice}"
DEST_ROOT="${2:-$REPO_ROOT/src/r/dumpnice}"

mkdir -p "$DEST_ROOT"
MANIFEST="$DEST_ROOT/_manifest.csv"
SUMMARY="$DEST_ROOT/_summary.txt"

echo "source,target,category,md5_source,md5_target" > "$MANIFEST"

categorize() {
  local path="$1"
  local base="$(basename "$path")"
  local parent="$(basename "$(dirname "$path")")"
  local s="${base,,}"
  local p="${parent,,}"
  if [[ "$s" =~ vcf ]] || [[ "$p" =~ vcf ]]; then echo vcf; return; fi
  if [[ "$s" =~ haplo ]] || [[ "$p" =~ haplo ]]; then echo haplotype; return; fi
  if [[ "$s" =~ gwas|gapit|rmvp|mvp ]] || [[ "$p" =~ gwas ]]; then echo gwas; return; fi
  if [[ "$s" =~ plot|pheatmap|heatmap|venn|density|qmatrix|barcode|(^|[^a-z])map([^a-z]|$)|genestructure|chrbin ]] ; then echo plot; return; fi
  if [[ "$s" =~ (^|[^a-z])fst([^a-z]|$) ]] || [[ "$p" =~ fst ]]; then echo fst; return; fi
  if [[ "$s" =~ (^|[^a-z])ibs([^a-z]|$) ]] || [[ "$p" =~ ibs ]]; then echo ibs; return; fi
  if [[ "$s" =~ (^|[^a-z])maf([^a-z]|$) ]] || [[ "$p" =~ maf ]]; then echo maf; return; fi
  if [[ "$s" =~ (^|[^a-z])pi([^a-z]|$) ]] || [[ "$p" =~ (^|[^a-z])pi([^a-z]|$) ]]; then echo pi; return; fi
  if [[ "$s" =~ (^|[^a-z])ld|ldjump ]] || [[ "$p" =~ ld ]]; then echo ld; return; fi
  if [[ "$s" =~ admixture|structure|qmatrix ]] || [[ "$p" =~ admixture ]]; then echo admixture; return; fi
  if [[ "$s" =~ enrich ]] || [[ "$p" =~ enrich ]]; then echo enrichment; return; fi
  if [[ "$s" =~ (^|[^a-z])pca([^a-z]|$) ]] || [[ "$p" =~ pca ]]; then echo pca; return; fi
  if [[ "$s" =~ (^|[^a-z])rda([^a-z]|$) ]] || [[ "$p" =~ rda ]]; then echo rda; return; fi
  if [[ "$s" =~ migration ]] || [[ "$p" =~ migration ]]; then echo migration; return; fi
  if [[ "$s" =~ permut ]] || [[ "$p" =~ permut ]]; then echo permutation; return; fi
  if [[ "$s" =~ depth ]] || [[ "$p" =~ depth ]]; then echo depth; return; fi
  if [[ "$s" =~ (^|[^a-z])bwa([^a-z]|$) ]] || [[ "$p" =~ bwa ]]; then echo bwa; return; fi
  if [[ "$s" =~ (^|[^a-z])tree([^a-z]|$)|nucdiff ]] || [[ "$p" =~ tree ]]; then echo tree; return; fi
  if [[ "$s" =~ statistic|basic_statistic|seperationsite|tajimasd|nesize|maf-dis ]] || [[ "$p" =~ statistic ]]; then echo statistic; return; fi
  if [[ "$s" =~ qc ]] || [[ "$p" =~ qc ]]; then echo qc; return; fi
  if [[ "$s" =~ 流程|vf_pipe|vmap3|storage|mcmc|bayenv|lfmm|variants_age|method|testgf|39_42-21chr ]] ; then echo pipeline; return; fi
  echo utils
}

copy_one() {
  local src="$1"
  local cat="$2"
  local dst_dir="$DEST_ROOT/$cat"
  mkdir -p "$dst_dir"
  local base="$(basename "$src")"
  local target="$dst_dir/$base"
  local md5_src="$(md5sum "$src" | awk '{print $1}')"
  if [[ -f "$target" ]]; then
    local md5_tgt="$(md5sum "$target" | awk '{print $1}')"
    if [[ "$md5_src" == "$md5_tgt" ]]; then
      echo "$src,$target,$cat,$md5_src,$md5_tgt" >> "$MANIFEST"
      return
    fi
    local stem="${base%.*}"
    local ext=".${base##*.}"
    local i=1
    while :; do
      local cand="$dst_dir/${stem}__$(printf '%02d' "$i")$ext"
      if [[ ! -f "$cand" ]]; then target="$cand"; break; fi
      i=$((i+1))
      if [[ $i -gt 99 ]]; then echo "Too many name collisions for $base" >&2; exit 1; fi
    done
  fi
  cp -p "$src" "$target"
  local md5_tgt="$(md5sum "$target" | awk '{print $1}')"
  echo "$src,$target,$cat,$md5_src,$md5_tgt" >> "$MANIFEST"
}

export -f categorize
export -f copy_one

mapfile -t FILES < <(find "$SOURCE_DIR" -type f \( -iname "*.r" -o -iname "*.R" \))

echo "Found ${#FILES[@]} R files under $SOURCE_DIR"
for f in "${FILES[@]}"; do
  catg="$(categorize "$f")"
  copy_one "$f" "$catg"
done

declare -A COUNTS
while IFS=, read -r src tgt cat md5s md5t; do
  if [[ "$src" == "source" ]]; then continue; fi
  COUNTS[$cat]=$(( ${COUNTS[$cat]:-0} + 1 ))
done < "$MANIFEST"

{
  echo "DumpNice organization summary"
  echo "Source: $SOURCE_DIR"
  echo "Dest:   $DEST_ROOT"
  echo
  echo "Counts by category:"
  tmpfile="$(mktemp)"
  for k in "${!COUNTS[@]}"; do
    printf "%s,%d\n" "$k" "${COUNTS[$k]}" >> "$tmpfile"
  done
  sort "$tmpfile" | while IFS=, read -r k v; do
    echo "- $k: $v"
  done
  rm -f "$tmpfile"
} > "$SUMMARY"

echo "Done. Manifest: $MANIFEST"
