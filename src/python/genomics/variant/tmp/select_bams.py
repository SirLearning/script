import pandas as pd
import numpy as np
import sys
import os

# Integers for reproducibility
SEED = 42
SAMPLE_SIZE = 500
OUTPUT_FILE = "selected_bams_paths.txt"

# Thresholds
MIN_DEPTH = 5.5
MAX_MISSING = 0.1

def parse_taxa_bam_map(filepath):
    """
    Manually parse the TaxaMap file to handle potentially variable number of BAM columns.
    Returns:
        depth_data: list of dicts [{'Taxa': t, 'Coverage': c}]
        bam_map: dict {Taxa: [bam_paths]}
    """
    depth_data = []
    bam_map = {}
    
    print(f"Parsing {filepath} line by line...")
    with open(filepath, 'r') as f:
        # Skip header logic: assuming first line is header if it starts with 'Taxa'
        header = f.readline()
        if not header.strip().startswith("Taxa"):
            # If no header, seek 0 (but example shows header)
            f.seek(0)
            
        for line in f:
            line = line.strip()
            if not line: continue
            
            parts = line.split('\t')
            if len(parts) < 3:
                # Need at least Taxa, Coverage, and one BAM
                # Although sometimes maybe 0 bams?
                continue
                
            taxa = parts[0]
            try:
                coverage = float(parts[1])
            except ValueError:
                # Header inside file or malformed line
                continue
                
            # BAMs are from index 2 to end
            bams = [b.strip() for b in parts[2:] if b.strip()]
            
            depth_data.append({'Taxa': taxa, 'Coverage-Of-All-Bams': coverage})
            bam_map[taxa] = bams
            
    return pd.DataFrame(depth_data), bam_map

def main():
    print("Processing BAM Selection...")
    
    depth_file = "/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap/all.A.taxaBamMap.txt"
    missing_file = "/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.missing.smiss"
    
    # 1. Read Depth and build BAM map
    try:
        df_depth, bam_map = parse_taxa_bam_map(depth_file)
        print(f"Loaded depth info for {len(df_depth)} samples.")
    except Exception as e:
        print(f"Error parsing depth file: {e}")
        sys.exit(1)
        
    # 2. Read Missing Rate
    print(f"Reading Missing Rate: {missing_file}")
    try:
        df_miss = pd.read_csv(missing_file, sep='\s+')
    except Exception as e:
        print(f"Error reading missing file: {e}")
        sys.exit(1)

    # 3. Merge
    merged = pd.merge(df_depth, df_miss, left_on='Taxa', right_on='#IID')
    print(f"Merged samples: {len(merged)}")
    
    # 4. Filter
    mask = (merged['Coverage-Of-All-Bams'] > MIN_DEPTH) & (merged['F_MISS'] < MAX_MISSING)
    filtered = merged[mask].copy()
    
    print(f"Samples passing criteria (Depth > {MIN_DEPTH}, Miss < {MAX_MISSING}): {len(filtered)}")
    
    if len(filtered) == 0:
        print("No samples found matching criteria.")
        sys.exit(0)

    # 5. Sampling
    n_samples = min(len(filtered), SAMPLE_SIZE)
    sampled = filtered.sample(n=n_samples, random_state=SEED)
    print(f"Selected {len(sampled)} samples.")
    
    # 6. Output BAM paths
    count = 0
    with open(OUTPUT_FILE, 'w') as f:
        for taxa in sampled['Taxa']:
            if taxa in bam_map:
                for bam_path in bam_map[taxa]:
                    f.write(bam_path + '\n')
                    count += 1
            else:
                print(f"Warning: Taxa {taxa} found in filtered list but missing from bam_map.")
                            
    print(f"Successfully wrote {count} bam paths to {OUTPUT_FILE}")
    print(f"Output file: {os.path.abspath(OUTPUT_FILE)}")

if __name__ == "__main__":
    main()
