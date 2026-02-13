from infra.utils.io import load_df_from_tsv, load_df_from_tsv_no_header
import pandas as pd
import numpy as np
import os

def load_df_from_plink2(input_file):
    """
    Loads a PLINK2 generated file (e.g., .scount) into a DataFrame.
    """
    print(f"[Info] Loading PLINK2 file: {input_file}")
    df = load_df_from_tsv(input_file)
    if df is None or df.empty:
        return None
    # Rename columns if needed, assuming #IID is standard
    if '#IID' in df.columns:
        df = df.rename(columns={'#IID': 'Sample'})
    return df

def load_df_from_king(input_file):
    """
    Loads KING output (.kin0, .king).
    """
    print(f"[Info] Loading KING file: {input_file}")
    df = load_df_from_tsv(input_file)
    if df is None: return None
    # Rename columns for consistency
    rename_map = {'#IID1': 'Sample1', 'IID1': 'Sample1', 'ID1': 'Sample1', 
                  'IID2': 'Sample2', 'ID2': 'Sample2', 'ID2': 'Sample2'}
    df = df.rename(columns=rename_map)
    return df

def load_ids_file(input_file):
    """
    Loads ID file (often .mibs.id or similar).
    Assuming space separated, usually 2 columns (FID IID).
    Returns list of IIDs.
    """
    print(f"[Info] Loading ID file: {input_file}")
    df = load_df_from_tsv(input_file)
    if df is not None:
        return df['IID'].astype(str).tolist()
    return []

def load_matrix_file(input_file):
    """
    Loads a matrix file (space separated, no header).
    """
    print(f"[Info] Loading Matrix file: {input_file}")
    df = load_df_from_tsv_no_header(input_file, col_names=None)
    return df

def load_group_file(input_file):
    """
    Loads group file. format: Sample Group
    Returns dictionary mapping Sample -> Group
    """
    if not os.path.exists(input_file):
        print(f"[Warning] Group file not found: {input_file}")
        return {}
    
    try:
        df = pd.read_csv(input_file, sep=r'\s+', header=None, names=['Sample', 'Group'])
        # Handle duplicates if any by taking first
        df = df.drop_duplicates(subset=['Sample'])
        return dict(zip(df['Sample'], df['Group']))
    except Exception as e:
        print(f"[Error] Failed to load group file: {e}")
        return {}