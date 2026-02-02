import pandas as pd
import numpy as np
import os
import sys

def verify_file_integrity(
    target_file=None, 
    reference_file=DEFAULT_REF_MD5, 
    output_pass_file=None
):
    """
    Verifies file integrity against a reference (CSV or MD5/TXT).
    
    Args:
        target_file (str): Path to file containing calculated MD5s (format: 'hash filename').
        reference_file (str): Path to reference file. 
                              If .csv, parsing assumes specific columns.
                              If .md5/.txt, assumes 'hash filename' structure.
        output_pass_file (str): Optional path to write passed sample pairs.
    """
    print(f"[Info] Verifying {target_file} against {reference_file}")
    
    # 1. Load Reference Map
    ref_map = {}   # filename -> md5
    info_map = {}  # filename -> extra_info (e.g. sample title)

    try:
        if reference_file.endswith('.csv'):
            ref_df = pd.read_csv(reference_file)
            # Map based on columns 'Read filename 1/2' and 'Read file1/2 MD5'
            for idx, row in ref_df.iterrows():
                # Process File 1
                f1 = str(row.get('Read filename 1', '')).split()[0]
                if f1:
                    ref_map[f1] = str(row.get('Read file1 MD5', '')).strip()
                    info_map[f1] = row.get('Run title', 'Unknown')
                
                # Process File 2
                f2 = str(row.get('Read filename 2', '')).split()[0]
                if f2:
                    ref_map[f2] = str(row.get('Read file2 MD5', '')).strip()
                    info_map[f2] = row.get('Run title', 'Unknown')
        else:
            # Standard MD5 file
            ref_df = pd.read_csv(reference_file, sep=r'\s+', header=None, names=['md5', 'file'])
            ref_df['short'] = ref_df['file'].apply(os.path.basename)
            ref_map = dict(zip(ref_df['short'], ref_df['md5']))
            
    except Exception as e:
        print(f"[Error] Failed to load reference: {e}")
        return

    # 2. Load Target File
    try:
        target_df = pd.read_csv(target_file, sep=r'\s+', header=None, names=['md5', 'file'])
        target_df['filename'] = target_df['file'].apply(os.path.basename)
    except Exception as e:
        print(f"[Error] Failed to load target file: {e}")
        return

    # 3. Validation Logic
    passed_indices = set()
    
    for idx, row in target_df.iterrows():
        fname = row['filename']
        computed_md5 = str(row['md5']).strip()
        
        expected = ref_map.get(fname)
        if expected:
            if computed_md5 == expected:
                passed_indices.add(idx)
            else:
                print(f"[Fail] {fname} | Expected: {expected} | Got: {computed_md5}")
        else:
            print(f"[Warning] {fname} not found in reference.")

    # 4. Reporting & Output
    pass_rate = (len(passed_indices) / len(target_df)) * 100 if len(target_df) > 0 else 0
    print(f"Pass Rate: {pass_rate:.2f}% ({len(passed_indices)}/{len(target_df)})")

    if output_pass_file:
        try:
            with open(output_pass_file, 'w') as f:
                # Write pairs only if both f1 (i) and f2 (i+1) passed
                # Assumes strict ordering in target file
                written_count = 0
                for i in range(0, len(target_df) - 1, 2):
                    if i in passed_indices and (i+1) in passed_indices:
                        f1_name = target_df.iloc[i]['filename']
                        sample_id = f1_name.split('_')[0]
                        run_title = info_map.get(f1_name, "Unknown")
                        f.write(f"{sample_id}\t{run_title}\n")
                        written_count += 1
            print(f"[Info] Written {written_count} passed pairs to {output_pass_file}")
        except Exception as e:
            print(f"[Error] Failed to write output file: {e}")


if __name__ == '__main__':
    verify_file_integrity()

