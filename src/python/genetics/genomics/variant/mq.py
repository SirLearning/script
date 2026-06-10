import pandas as pd
from infra.utils.io import load_df_from_space_sep_no_header

def load_mq_data(filepath, keep_missing=False):
    """
    Load site MQ (3 columns: Chrom, Position, MQ) from padded reference grid.
    """
    print(f"[Info] Loading MQ reference grid: {filepath}")
    col_names = ['Chrom', 'Position', 'MQ']
    df = load_df_from_space_sep_no_header(filepath, col_names=col_names)
    if df is None:
        return None
    df['MQ'] = pd.to_numeric(df['MQ'], errors='coerce')
    if not keep_missing:
        df = df.dropna(subset=['MQ'])
    return df


def load_mq_calls(filepath, keep_missing=False):
    """
    Load per-site mpileup MQ table (5 columns: Chrom, Position, REF, ALT, MQ).
    MQ is float mean MAPQ from bcftools mpileup INFO/I16 (one row per mpileup site).
    """
    print(f"[Info] Loading MQ calls: {filepath}")
    col_names = ['Chrom', 'Position', 'REF', 'ALT', 'MQ']
    df = load_df_from_space_sep_no_header(filepath, col_names=col_names)
    if df is None:
        return None
    df['MQ'] = pd.to_numeric(df['MQ'], errors='coerce')
    if not keep_missing:
        df = df.dropna(subset=['MQ'])
    return df
