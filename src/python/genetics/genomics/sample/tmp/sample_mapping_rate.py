import os
import re
import sys
import subprocess
import argparse
import concurrent.futures
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(description="Parallel samtools flagstat parser for Vmap4 project.")
    parser.add_argument("-i", "--input", required=True, help="Path to all.ALL.taxaBamMap.txt")
    parser.add_argument("-o", "--output", default="flagstat_summary.txt", help="Output tab-separated txt file")
    parser.add_argument("-t", "--threads", type=int, default=32, help="Number of parallel threads (default: 32)")
    return parser.parse_args()

def run_flagstat(task):
    """
    执行 samtools flagstat 并解析输出
    task = (taxa_id, bam_path)
    """
    taxa_id, bam_path = task
    if not os.path.exists(bam_path):
        return [taxa_id, bam_path, "FileMissing", "0.00", "0", "0", "0"]

    try:
        # 调用 samtools flagstat
        proc = subprocess.run(
            ["samtools", "flagstat", bam_path],
            capture_output=True, text=True, check=True
        )
        stdout = proc.stdout

        # 正则解析
        # 1. Mapped % -> mapped (100.00% : N/A)
        m_pct = re.search(r"mapped \((\d+\.\d+)%", stdout)
        # 2. Properly paired % -> properly paired (79.43% : N/A)
        pp_pct = re.search(r"properly paired \((\d+\.\d+)%", stdout)
        # 3. Primary mapped count -> 94300301 + 0 primary mapped
        pm_cnt = re.search(r"(\d+) \+ \d+ primary mapped", stdout)
        # 4. Singletons count -> 370837 + 0 singletons
        s_cnt = re.search(r"(\d+) \+ \d+ singletons", stdout)
        # 5. Different chr count -> with mate mapped to a different chr (非 mapQ>=5 的那行)
        dc_cnt = re.search(r"(\d+) \+ \d+ with mate mapped to a different chr\n", stdout)

        return [
            taxa_id,
            m_pct.group(1) if m_pct else "0.00",
            pp_pct.group(1) if pp_pct else "0.00",
            pm_cnt.group(1) if pm_cnt else "0",
            s_cnt.group(1) if s_cnt else "0",
            dc_cnt.group(1) if dc_cnt else "0"
        ]
    except Exception as e:
        return [taxa_id, "Error", "0.00", "0", "0", "0"]

def main():
    args = parse_args()
    
    tasks = []
    print(f"[{sys.argv[0]}] Reading input: {args.input}")

    # 解析输入文件
    try:
        with open(args.input, 'r') as f:
            header = f.readline() # 跳过表头
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) < 3: continue
                
                taxa_id = cols[0]
                # 提取第三列及其后所有BAM路径（防止一个Taxa有多个BAM）
                for bam_path in cols[2:]:
                    if bam_path.strip():
                        tasks.append((taxa_id, bam_path.strip()))
    except FileNotFoundError:
        print(f"Error: Input file {args.input} not found.")
        return

    print(f"Total BAMs to process: {len(tasks)}")
    print(f"Running with {args.threads} parallel processes...")

    # 输出表头
    out_columns = [
        "TaxaID", "Mapped_Pct", "Properly_Paired_Pct", 
        "Primary_Mapped_Count", "Singletons_Count", "Diff_Chr_Count"
    ]

    with open(args.output, 'w') as f_out:
        f_out.write("\t".join(out_columns) + "\n")
        f_out.flush()

        # 使用进程池执行并行任务
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
            future_to_task = {executor.submit(run_flagstat, task): task for task in tasks}
            
            # as_completed yields futures as they complete
            for future in tqdm(concurrent.futures.as_completed(future_to_task), total=len(tasks)):
                try:
                    res = future.result()
                    f_out.write("\t".join(res) + "\n")
                    f_out.flush()
                except Exception as exc:
                    print(f"Task generated an exception: {exc}")

    print(f"\nSuccess! Results saved to: {args.output}")

if __name__ == "__main__":
    main()