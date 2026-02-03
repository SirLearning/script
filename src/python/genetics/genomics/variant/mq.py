import os
import sys
import argparse
import pandas as pd
import numpy as np

# Ensure src/python is in path
current_dir = os.path.dirname(os.path.abspath(__file__))
src_python_dir = os.path.abspath(os.path.join(current_dir, "../../../.."))
if src_python_dir not in sys.path:
    sys.path.insert(0, src_python_dir)

from infra.file_utils import load_df_from_space_sep_no_header
from genetics.genomics.sample.smiss import load_smiss

def load_mq(filepath):
    """
    Loads Mapping Quality file.
    Assumes no header: Chrom, Position, MQ
    """
    print(f"[Info] Loading MQ: {filepath}")
    col_names = ['Chrom', 'Position', 'MQ']
    df = load_df_from_space_sep_no_header(filepath, col_names)
    if df is None: return None
    
    # Coerce MQ to numeric
    df['MQ'] = pd.to_numeric(df['MQ'], errors='coerce')
    df = df.dropna(subset=['MQ'])
    return df

# ---------------------------------------------------------
# BAM Selection Logic
# ---------------------------------------------------------

def parse_taxa_bam_map_with_coverage(filepath):
    """
    Parses the TaxaMap file and extracts coverage info.
    Returns:
        depth_data: DataFrame with columns ['Taxa', 'Coverage']
        bam_map: dict {Taxa: [bam_paths]}
    """
    depth_data = []
    bam_map = {}
    
    print(f"Parsing {filepath} line by line...")
    try:
        with open(filepath, 'r') as f:
            # Skip header logic: assuming first line is header if it starts with "Taxa"
            header = f.readline()
            if not header.strip().startswith("Taxa"):
                f.seek(0)
                
            for line in f:
                line = line.strip()
                if not line: continue
                
                parts = line.split('\t')
                # Need at least Taxa, Coverage, and preferably BAMs
                if len(parts) < 2:
                    continue
                    
                taxa = parts[0]
                try:
                    coverage = float(parts[1])
                except ValueError:
                    # Skipping malformed line (maybe header inside file)
                    continue
                
                # BAMs are from index 2 to end
                bams = [b.strip() for b in parts[2:] if b.strip()]
                
                depth_data.append({'Taxa': taxa, 'Coverage': coverage})
                bam_map[taxa] = bams
                
        return pd.DataFrame(depth_data), bam_map
    except FileNotFoundError:
        print(f"Error: Input file {filepath} not found.")
        return pd.DataFrame(), {}

def select_high_quality_bams(
    depth_file,
    missing_file,
    output_file,
    min_depth=5.5,
    max_missing=0.1,
    sample_size=500,
    seed=42
):
    print("Processing BAM Selection...")
    
    # 1. Read Depth and build BAM map
    df_depth, bam_map = parse_taxa_bam_map_with_coverage(depth_file)
    if df_depth.empty:
        print("Error: No depth info loaded.")
        return

    print(f"Loaded depth info for {len(df_depth)} samples.")
        
    # 2. Read Missing Rate
    print(f"Reading Missing Rate: {missing_file}")
    try:
        df_miss = load_smiss(missing_file)
        if df_miss is None: return
    except Exception as e:
        print(f"Error reading missing file: {e}")
        return

    # 3. Merge
    # df_depth has ['Taxa', 'Coverage']
    # df_miss has ['Sample', 'Missing_Rate']
    merged = pd.merge(df_depth, df_miss, left_on='Taxa', right_on='Sample')
    print(f"Merged samples: {len(merged)}")
    
    # 4. Filter
    mask = (merged['Coverage'] > min_depth) & (merged['Missing_Rate'] < max_missing)
    filtered = merged[mask].copy()
    
    print(f"Samples passing criteria (Depth > {min_depth}, Miss < {max_missing}): {len(filtered)}")
    
    if len(filtered) == 0:
        print("No samples found matching criteria.")
        return

    # 5. Sampling
    n_samples = min(len(filtered), sample_size)
    sampled = filtered.sample(n=n_samples, random_state=seed)
    print(f"Selected {len(sampled)} samples.")
    
    # 6. Output BAM paths
    count = 0
    try:
        # Create output directory if it doesn't exist
        out_dir = os.path.dirname(os.path.abspath(output_file))
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)

        with open(output_file, 'w') as f:
            for taxa in sampled['Taxa']:
                if taxa in bam_map:
                    for bam_path in bam_map[taxa]:
                        f.write(bam_path + '\n')
                        count += 1
                else:
                    print(f"Warning: Taxa {taxa} found in filtered list but missing from bam_map.")
        
        print(f"Successfully wrote {count} bam paths to {output_file}")
        print(f"Output file: {os.path.abspath(output_file)}")
    except IOError as e:
        print(f"Error writing output file: {e}")

# ---------------------------------------------------------
# Main CLI
# ---------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Mapping Quality (MQ) Analysis Tools")
    subparsers = parser.add_subparsers(dest='command', help='Task to perform')

    # Select BAMs Command
    p_sel = subparsers.add_parser('select_bams', help='Select high quality BAMs based on depth and missing rate')
    p_sel.add_argument("-d", "--depth", required=True, help="Input taxaBamMap.txt (for depth info)")
    p_sel.add_argument("-s", "--smiss", required=True, help="PLINK .smiss file")
    p_sel.add_argument("-o", "--output", default="selected_bams_paths.txt", help="Output text file with BAM paths")
    p_sel.add_argument("--min_depth", type=float, default=5.5, help="Minimum depth threshold (default: 5.5)")
    p_sel.add_argument("--max_missing", type=float, default=0.1, help="Maximum missing rate threshold (default: 0.1)")
    p_sel.add_argument("--sample_size", type=int, default=500, help="Number of samples to select (default: 500)")
    p_sel.add_argument("--seed", type=int, default=42, help="Random seed for sampling (default: 42)")

    args = parser.parse_args()

    if args.command == 'select_bams':
        select_high_quality_bams(
            args.depth, 
            args.smiss, 
            args.output, 
            args.min_depth, 
            args.max_missing, 
            args.sample_size,
            args.seed
        )
    else:
        # Default behavior if no command or just imported
        parser.print_help()

if __name__ == "__main__":
    main()
