"""Plot CNV calls from awk depth pipeline output."""

import pandas as pd

from infra.utils.graph import plot_distribution_with_stats
from infra.utils.io import load_df_generic


def plot_cnv_results(cnv_tsv_path, output_prefix, del_z=-2.0, dup_z=2.0):
    out = load_df_generic(cnv_tsv_path)
    zdf = out[['ZScore']].copy() if out is not None and not out.empty and 'ZScore' in out.columns else pd.DataFrame({'ZScore': [0.0]})
    plot_distribution_with_stats(
        data=zdf, col='ZScore', title='CNV Z-score Distribution',
        filename=f'{output_prefix}.zscore.png', x_label='Depth Z-score', y_label='CNV Count',
        thresholds=[
            {'value': float(del_z), 'label': f'DEL <= {del_z}', 'color': 'red', 'linestyle': '--'},
            {'value': float(dup_z), 'label': f'DUP >= {dup_z}', 'color': 'blue', 'linestyle': '--'},
        ],
    )
