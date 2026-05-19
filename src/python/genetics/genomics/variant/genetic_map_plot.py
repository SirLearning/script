"""Plot precomputed genetic map tables (external ASMap / R/qtl output)."""

from infra.utils.graph import plot_multi_line_series
from infra.utils.io import load_df_generic


def plot_genetic_map(map_tsv_path, output_prefix):
    out = load_df_generic(map_tsv_path)
    if out is None or out.empty:
        raise ValueError(f'Empty genetic map table: {map_tsv_path}')
    for chr_name, g in out.groupby('CHR'):
        g2 = g.sort_values('POS')
        plot_multi_line_series(
            data=g2, x_col='POS',
            y_specs=[{'y_col': 'CumCM', 'label': f'{chr_name} cumulative cM'}],
            title=f'Genetic Map ({chr_name})', filename=f'{output_prefix}.{chr_name}.png',
            x_label='Physical Position', y_label='Genetic Distance (cM)',
        )
