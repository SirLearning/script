from genetics.germplasm import integrate_group_info
from genetics.germplasm.sample.anno import anno_group
from genetics.germplasm.sample.process import load_df_from_plink2, load_df_from_tbm
from infra.utils import plot_joint_regression
from infra.utils.io import load_df_from_tsv, save_df_to_tsv
import os
import re
import argparse
import subprocess
import concurrent.futures
import pandas as pd
import seaborn as sns
from tqdm import tqdm


def load_df_from_idxstats(idxstats_file):
    print(f"Reading Mapping Rate file: {idxstats_file}...")
    try:
        df = load_df_from_tsv(idxstats_file)
    except Exception as e:
        print(f"Error reading mapping file: {e}")
        return
    if 'TaxaID' not in df.columns or 'Mapping_Rate_Pct' not in df.columns:
        print(f"Error: Required columns ('TaxaID', 'Mapping_Rate_Pct') not found in mapping file.")
        print(f"Columns found: {list(df.columns)}")
        return
    return df

def parse_bam_map_file(input_file):
    """
    Parses the taxaBamMap file using shared utility.
    Returns a list of tasks: (taxa_id, bam_path)
    """
    tasks = []
    # Load using centralized logic
    df = load_df_from_tbm(input_file)
    
    if df is None or df.empty:
        return []

    # Iterate rows
    # Expectation: df has 'Sample' and 'Bam_Path' (list)
    if 'Sample' not in df.columns or 'Bam_Path' not in df.columns:
        print(f"[Error] Loaded dataframe missing required columns (Sample, Bam_Path). Columns: {df.columns}")
        return []

    for _, row in df.iterrows():
        taxa_id = str(row['Sample']) # Ensure string
        bam_list = row['Bam_Path']
        
        if isinstance(bam_list, list):
            for p in bam_list:
                if p: tasks.append((taxa_id, str(p).strip()))
        
    return tasks

# ---------------------------------------------------------
# Idxstats Logic
# ---------------------------------------------------------

def run_idxstats_task(task):
    """
    Execute samtools idxstats and calculate mapping rate
    task = (taxa_id, bam_path)
    """
    taxa_id, bam_path = task
    
    if not os.path.exists(bam_path):
        return [taxa_id, bam_path, "FileMissing", "0", "0", "0.00"]

    try:
        proc = subprocess.run(
            ["samtools", "idxstats", bam_path],
            capture_output=True, text=True, check=True
        )
        
        total_mapped = 0
        total_unmapped = 0
        
        for line in proc.stdout.splitlines():
            parts = line.split('\t')
            if len(parts) >= 4:
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

def calculate_mapping_rate_idxstats(input_file, output_file, threads=32):
    print(f"[Idxstats] Reading input: {input_file}")
    tasks = parse_bam_map_file(input_file)
    if not tasks:
        return

    print(f"Total BAMs to process: {len(tasks)}")
    print(f"Running with {threads} parallel processes...")

    out_columns = ["TaxaID", "BamPath", "Mapped_Reads", "Unmapped_Reads", "Mapping_Rate_Pct"]

    with open(output_file, 'w') as f_out:
        f_out.write("\t".join(out_columns) + "\n")
        f_out.flush() 

        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            task_futures = {executor.submit(run_idxstats_task, t): t for t in tasks}
            
            for future in tqdm(concurrent.futures.as_completed(task_futures), total=len(tasks)):
                try:
                    res = future.result()
                    f_out.write("\t".join(res) + "\n")
                    f_out.flush()
                except Exception as exc:
                    print(f"Task generated an exception: {exc}")

    print(f"\nSuccess! Results saved to: {output_file}")

# ---------------------------------------------------------
# Flagstat Logic
# ---------------------------------------------------------

def run_flagstat_task(task):
    """
    Execute samtools flagstat and parse output
    task = (taxa_id, bam_path)
    """
    taxa_id, bam_path = task
    if not os.path.exists(bam_path):
        return [taxa_id, bam_path, "FileMissing", "0.00", "0", "0", "0"]

    try:
        proc = subprocess.run(
            ["samtools", "flagstat", bam_path],
            capture_output=True, text=True, check=True
        )
        stdout = proc.stdout

        m_pct = re.search(r"mapped \((\d+\.\d+)%", stdout)
        pp_pct = re.search(r"properly paired \((\d+\.\d+)%", stdout)
        pm_cnt = re.search(r"(\d+) \+ \d+ primary mapped", stdout)
        s_cnt = re.search(r"(\d+) \+ \d+ singletons", stdout)
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

def calculate_mapping_rate_flagstat(input_file, output_file, threads=32):
    print(f"[Flagstat] Reading input: {input_file}")
    tasks = parse_bam_map_file(input_file)
    if not tasks:
        return

    print(f"Total BAMs to process: {len(tasks)}")
    print(f"Running with {threads} parallel processes...")

    out_columns = [
        "TaxaID", "Mapped_Pct", "Properly_Paired_Pct", 
        "Primary_Mapped_Count", "Singletons_Count", "Diff_Chr_Count"
    ]

    with open(output_file, 'w') as f_out:
        f_out.write("\t".join(out_columns) + "\n")
        f_out.flush()

        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            task_futures = {executor.submit(run_flagstat_task, t): t for t in tasks}
            
            for future in tqdm(concurrent.futures.as_completed(task_futures), total=len(tasks)):
                try:
                    res = future.result()
                    f_out.write("\t".join(res) + "\n")
                    f_out.flush()
                except Exception as exc:
                    print(f"Task generated an exception: {exc}")

    print(f"\nSuccess! Results saved to: {output_file}")

# ---------------------------------------------------------
# Correlation Analysis Logic
# ---------------------------------------------------------

def reg_missing_idx_mapping(
    mapping_file, 
    smiss_file, 
    group_file=None, 
    tsv_file=None,
    output_prefix="missing_vs_mapping"
):
    print("Processing Missing Rate vs Mapping Rate Analysis...")

    if tsv_file and os.path.exists(tsv_file):
        print(f"Loading precomputed data from: {tsv_file}...")
        df_merged = load_df_from_tsv(tsv_file)
        if df_merged is None: return
    else:
        # 1. Read Mapping Rate Data
        df_idxstats = load_df_from_idxstats(mapping_file)
        df_map = df_idxstats[['TaxaID', 'Mapping_Rate_Pct']].rename(columns={'TaxaID': 'Sample'})
        # 2. Read .smiss Data using smiss_ana
        print(f"Reading .smiss file: {smiss_file}...")
        df_missing = load_df_from_plink2(smiss_file)
        if df_missing is None:
            return
        # 3. Merge Datasets
        print("Merging datasets...")
        df_merged = pd.merge(df_map, df_missing, on='Sample', how='inner')
        print(f"Matched samples: {len(df_merged)}")
        if len(df_merged) == 0:
            print("Error: No intersecting samples found between mapping file and smiss file.")
            return
        # 4. Integrate Group Info
        df_merged = anno_group(df_merged, group_file) if group_file else df_merged['Group'] = 'Unknown'
        # 5. Save merged data
        save_file = f"{output_prefix}.data.tsv"
        print(f"Saving merged data to {save_file}...")
        save_df_to_tsv(df_merged, save_file)
        
    # Print basic stats
    if 'Mapping_Rate_Pct' in df_merged.columns and 'Missing_Rate' in df_merged.columns:
         print(f"Correlation (Pearson): {df_merged['Mapping_Rate_Pct'].corr(df_merged['Missing_Rate']):.4f}")
    
    # Set style
    sns.set_theme(style="ticks")
    
    plot_joint_regression(
        df=df_merged, 
        x_col='Mapping_Rate_Pct', 
        y_col='Missing_Rate',
        group_col='Group',
        x_label='Mapping Rate (%)', 
        y_label='Missing Rate (F_MISS)', 
        filename=f"{output_prefix}_reg_miss_vs_map.png",
        title='Missing Rate vs Mapping Rate'
    )
                         
    print("Analysis Complete.")

# ---------------------------------------------------------
# Main CLI
# ---------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Alignment Analysis Tools: Idxstats, Flagstat, Correlation.")
    subparsers = parser.add_subparsers(dest='command', help='Task to perform')

    # 1. Idxstats
    p_idx = subparsers.add_parser('idxstats', help='Calculate mapping rate via samtools idxstats')
    p_idx.add_argument("-i", "--input", required=True, help="Input taxaBamMap.txt")
    p_idx.add_argument("-o", "--output", default="idxstats_mapping_rate.txt", help="Output TSV file")
    p_idx.add_argument("-t", "--threads", type=int, default=32, help="Parallel threads")

    # 2. Flagstat
    p_flag = subparsers.add_parser('flagstat', help='Calculate stats via samtools flagstat')
    p_flag.add_argument("-i", "--input", required=True, help="Input taxaBamMap.txt")
    p_flag.add_argument("-o", "--output", default="flagstat_summary.txt", help="Output TSV file")
    p_flag.add_argument("-t", "--threads", type=int, default=32, help="Parallel threads")

    # 3. Correlation
    p_corr = subparsers.add_parser('correlation', help='Correlate Mapping Rate vs Missing Rate')
    p_corr.add_argument("-m", "--mapping", help="Mapping rate file (from idxstats)")
    p_corr.add_argument("-s", "--smiss", help="PLINK .smiss file")
    p_corr.add_argument("-d", "--data", dest="tsv_file", help="Precomputed data TSV file")
    p_corr.add_argument("-g", "--group", help="Sample group file (optional)")
    p_corr.add_argument("-o", "--output_prefix", default="missing_vs_mapping", help="Output file prefix")

    args = parser.parse_args()

    if args.command == 'idxstats':
        calculate_mapping_rate_idxstats(args.input, args.output, args.threads)
    elif args.command == 'flagstat':
        calculate_mapping_rate_flagstat(args.input, args.output, args.threads)
    elif args.command == 'correlation':
        reg_missing_idx_mapping(args.mapping, args.smiss, args.group, args.tsv_file, args.output_prefix)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
