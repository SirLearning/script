import pandas as pd
from infra.utils.io import load_df_from_space_sep_no_header

def load_mq_data(filepath):
    """
    Loads Mapping Quality file.
    Assumes no header: Chrom, Position, MQ
    """
    print(f"[Info] Loading MQ: {filepath}")
    col_names = ['Chrom', 'Position', 'MQ']
    # Explicitly using load_df_from_space_sep_no_header which IO should provide or generic
    # If not present in io.py, we trust it will be added or we use generic. 
    # Based on inheritance, let's assume io.py logic.
    df = load_df_from_space_sep_no_header(filepath, col_names=col_names)
    if df is None: return None
    
    # Coerce MQ to numeric
    df['MQ'] = pd.to_numeric(df['MQ'], errors='coerce')
    df = df.dropna(subset=['MQ'])
    return df
