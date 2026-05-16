import pandas as pd
from infra.utils.io import load_df_generic, save_df_to_tsv


def run_hapmap_build(input_file, output_prefix, window_size=100000):
    """Build a simple block-based HAPMAP table from marker genotypes."""
    df = load_df_generic(input_file)
    if df is None or df.empty:
        raise ValueError(f'No genotype records loaded from: {input_file}')

    for col in ('CHR', 'POS'):
        if col not in df.columns:
            raise ValueError(f'Missing required column: {col}')

    df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
    df = df.dropna(subset=['POS']).copy()
    df['BlockStart'] = (df['POS'] // int(window_size)) * int(window_size)
    df['BlockID'] = df['CHR'].astype(str) + ':' + df['BlockStart'].astype(int).astype(str)

    meta_cols = {'CHR', 'POS', 'ID', 'REF', 'ALT', 'BlockStart', 'BlockID'}
    sample_cols = [c for c in df.columns if c not in meta_cols]

    def consensus(series):
        s = series.dropna().astype(str)
        return s.mode().iloc[0] if not s.empty else 'NA'

    grouped = df.groupby(['BlockID', 'CHR', 'BlockStart'], as_index=False)
    out = grouped[sample_cols].agg(consensus)
    marker_counts = grouped.size().rename(columns={'size': 'MarkerCount'})
    out = out.merge(marker_counts, on=['BlockID', 'CHR', 'BlockStart'], how='left')

    save_df_to_tsv(out, f'{output_prefix}.hapmap.tsv')

