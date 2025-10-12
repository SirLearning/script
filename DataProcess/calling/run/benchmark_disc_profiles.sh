#!/usr/bin/env bash
# Benchmark different fastcall3_disc resource profiles.
# Usage: bash benchmark_disc_profiles.sh <base_cmd_args>
# Example:
#   bash benchmark_disc_profiles.sh "--home_dir /data/dazheng/01projects/vmap4 --java_lib /data/dazheng/lib/jvm --pop chr1 --job test_ABD --workflow_mode disc_only --tiger_jar TIGER_F3_20250915.jar"

set -euo pipefail

BASE_ARGS="$*"
if [[ -z "$BASE_ARGS" ]]; then
  echo "Provide base Nextflow arguments (everything after the .nf path)" >&2
  exit 1
fi

NF_SCRIPT="/data/dazheng/git/script/DataProcess/calling/run/runFastCall3.nf"
CONFIG_FILE="/data/dazheng/git/script/DataProcess/calling/run/nextflow.config"
PROFILES=(disc_c8_f12 disc_c12_f10 disc_c16_f8 disc_c24_f5 disc_c32_f4)
RESULTS_DIR="benchmark_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS_DIR"

cat >"$RESULTS_DIR"/README.txt <<EOT
Benchmark of fastcall3_disc profiles
Date: $(date)
Host: $(hostname)
Profiles: ${PROFILES[*]}
Command base args: $BASE_ARGS
EOT

summary_csv="$RESULTS_DIR/summary.csv"
echo "profile,total_tasks,avg_cpu,max_cpu,avg_rss_gb,max_rss_gb,avg_read_mb,avg_write_mb,wall_minutes" > "$summary_csv"

for profile in "${PROFILES[@]}"; do
  echo "=== Running profile: $profile ==="
  run_dir="$RESULTS_DIR/run_$profile"
  mkdir -p "$run_dir"
  nextflow run "$NF_SCRIPT" -c "$CONFIG_FILE" $BASE_ARGS -profile $profile -with-trace "$run_dir/trace.txt" -ansi-log false > "$run_dir/pipeline.log" 2>&1 || echo "Run failed for $profile (continuing)" >&2
  if [[ -f "$run_dir/trace.txt" ]]; then
    # Extract metrics (assuming tab-delimited trace)
    # Skip header for averages
    awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++){h[$i]=i}} NR>1 {cpu+=$h["%cpu"]; rss+=$h["rss"]; r+=$h["read_bytes"]; w+=$h["write_bytes"]; if($h["%cpu"]>maxcpu)maxcpu=$h["%cpu"]; if($h["rss"]>maxrss)maxrss=$h["rss"]; n++} END {if(n>0) printf("%s,%d,%.2f,%.2f,%.3f,%.3f,%.2f,%.2f,NA\n", ENVIRON["profile"], n, cpu/n, maxcpu, rss/n/1024/1024/1024, maxrss/1024/1024/1024, r/n/1024/1024, w/n/1024/1024);}' profile="$profile" "$run_dir/trace.txt" >> "$summary_csv"
  else
    echo "$profile,0,NA,NA,NA,NA,NA,NA,NA" >> "$summary_csv"
  fi
 done

echo "Benchmark complete. Summary: $summary_csv"
