import os
import sys
import subprocess
import argparse
import concurrent.futures
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(description="Parallel samtools idxstats parser for Vmap4 project.")
    parser.add_argument("-i", "--input", required=True, help="Path to all.ALL.taxaBamMap.txt")
    parser.add_argument("-o", "--output", default="idxstats_mapping_rate.txt", help="Output tab-separated txt file")
    parser.add_argument("-t", "--threads", type=int, default=32, help="Number of parallel threads (default: 32)")
    return parser.parse_args()

def run_idxstats(task):
    """
    Execute samtools idxstats and calculate mapping rate
    task = (taxa_id, bam_path)
    """
    taxa_id, bam_path = task
    
    # Check if BAM exists
    if not os.path.exists(bam_path):
        return [taxa_id, bam_path, "FileMissing", "0", "0", "0.00"]

    try:
        # Run samtools idxstats
        # idxstats requires the .bai/.csi index file to be present.
        # It is extremely fast because it reads the index header.
        proc = subprocess.run(
            ["samtools", "idxstats", bam_path],
            capture_output=True, text=True, check=True
        )
        
        total_mapped = 0
        total_unmapped = 0
        
        # Parse output
        # Format: refName refLength mappedReads unmappedReads
        for line in proc.stdout.splitlines():
            parts = line.split('\t')
            if len(parts) >= 4:
                # ref_name = parts[0]
                # ref_len = int(parts[1])
                mapped = int(parts[2])
                unmapped = int(parts[3])
                
                total_mapped += mapped
                total_unmapped += unmapped
        
        total_reads = total_mapped + total_unmapped
        
        if total_reads > 0:
            mapping_rate = (total_mapped / total_reads) * 100
        else:
            mapping_rate = 0.0
            
        return [
            taxa_id,
            bam_path,
            str(total_mapped),
            str(total_unmapped),
            f"{mapping_rate:.4f}"
        ]

    except subprocess.CalledProcessError:
        return [taxa_id, bam_path, "Error_Samtools", "0", "0", "0.00"]
    except Exception as e:
        return [taxa_id, bam_path, f"Error: {str(e)}", "0", "0", "0.00"]

def main():
    args = parse_args()
    
    tasks = []
    print(f"[{sys.argv[0]}] Reading input: {args.input}")

    # Parse input file
    try:
        with open(args.input, 'r') as f:
            header = f.readline() 
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) < 3: continue
                
                taxa_id = cols[0]
                # Extract 3rd column and onwards for BAM paths
                for bam_path in cols[2:]:
                    if bam_path.strip():
                        tasks.append((taxa_id, bam_path.strip()))
    except FileNotFoundError:
        print(f"Error: Input file {args.input} not found.")
        return

    print(f"Total BAMs to process: {len(tasks)}")
    print(f"Running with {args.threads} parallel processes...")

    # Output columns
    out_columns = ["TaxaID", "BamPath", "Mapped_Reads", "Unmapped_Reads", "Mapping_Rate_Pct"]

    with open(args.output, 'w') as f_out:
        f_out.write("\t".join(out_columns) + "\n")
        f_out.flush() 

        # Use ProcessPoolExecutor for parallelism
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
            task_futures = {executor.submit(run_idxstats, task): task for task in tasks}
            
            # as_completed yields futures as they complete
            for future in tqdm(concurrent.futures.as_completed(task_futures), total=len(tasks)):
                try:
                    res = future.result()
                    f_out.write("\t".join(res) + "\n")
                    f_out.flush() # Flush immediately
                except Exception as exc:
                    print(f"Task generated an exception: {exc}")

    print(f"\nSuccess! Results saved to: {args.output}")

if __name__ == "__main__":
    main()