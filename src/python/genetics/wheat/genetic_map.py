import numpy as np
import pandas as pd
from infra.utils.io import load_df_generic, save_df_to_tsv
from infra.utils.graph import plot_multi_line_series


def _haldane_cm(r):
    r = min(max(float(r), 1e-6), 0.499999)
    return -50.0 * np.log(1.0 - 2.0 * r)


def run_genetic_map(input_file, output_prefix):
    """Build a simple genetic map from marker genotypes (WatSeqAnalysis ASMap-inspired Python flow)."""
    df = load_df_generic(input_file)
    if df is None or df.empty:
        raise ValueError(f'No marker table loaded from: {input_file}')

    for c in ('Marker', 'CHR', 'POS'):
        if c not in df.columns:
            raise ValueError(f'Missing required column: {c}')

    sample_cols = [c for c in df.columns if c not in {'Marker', 'CHR', 'POS'}]
    if not sample_cols:
        raise ValueError('No genotype sample columns found for map construction')

    df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
    df = df.dropna(subset=['POS']).sort_values(['CHR', 'POS']).copy()

    map_rows = []
    for chr_name, g in df.groupby('CHR'):
        g = g.reset_index(drop=True)
        cum_cm = 0.0
        map_rows.append({'CHR': chr_name, 'Marker': g.loc[0, 'Marker'], 'POS': g.loc[0, 'POS'], 'DistCM': 0.0, 'CumCM': cum_cm})

        for i in range(1, len(g)):
            prev = pd.to_numeric(g.loc[i - 1, sample_cols], errors='coerce')
            cur = pd.to_numeric(g.loc[i, sample_cols], errors='coerce')
            valid = (~prev.isna()) & (~cur.isna())
            if valid.sum() == 0:
                r = 0.0
            else:
                r = float((prev[valid] != cur[valid]).mean())
            d = _haldane_cm(r)
            cum_cm += d
            map_rows.append({'CHR': chr_name, 'Marker': g.loc[i, 'Marker'], 'POS': g.loc[i, 'POS'], 'DistCM': d, 'CumCM': cum_cm})

    out = pd.DataFrame(map_rows)
    save_df_to_tsv(out, f'{output_prefix}.map.tsv')

    chr_len = out.groupby('CHR', as_index=False)['CumCM'].max().rename(columns={'CumCM': 'MapLengthCM'})
    save_df_to_tsv(chr_len, f'{output_prefix}.chr_length.tsv')

    for chr_name, g in out.groupby('CHR'):
        g2 = g.sort_values('POS')
        plot_multi_line_series(
            data=g2,
            x_col='POS',
            y_specs=[{'y_col': 'CumCM', 'label': f'{chr_name} cumulative cM'}],
            title=f'Genetic Map ({chr_name})',
            filename=f'{output_prefix}.{chr_name}.png',
            x_label='Physical Position',
            y_label='Genetic Distance (cM)',
        )

