import argparse
import numpy as np
import pandas as pd
from infra.utils.io import load_df_generic, save_df_to_tsv
from infra.utils.graph import plot_distribution_with_stats


def run_cnv_calling(input_file, output_prefix, del_z=-2.0, dup_z=2.0):
    """Call CNVs from normalized depth matrix (row=region, col=sample depth)."""
    df = load_df_generic(input_file)
    if df is None or df.empty:
        raise ValueError(f'No depth records loaded from: {input_file}')

    required = {'CHR', 'START', 'END'}
    if not required.issubset(df.columns):
        raise ValueError('Input requires CHR, START, END columns')

    sample_cols = [c for c in df.columns if c not in {'CHR', 'START', 'END'}]
    if not sample_cols:
        raise ValueError('No sample depth columns found')

    long_rows = []
    for sample in sample_cols:
        vals = pd.to_numeric(df[sample], errors='coerce')
        mu = vals.mean()
        sd = vals.std()
        z = (vals - mu) / sd if sd and np.isfinite(sd) and sd > 0 else pd.Series([0.0] * len(vals))
        for i, zv in enumerate(z):
            if pd.isna(zv):
                continue
            state = 'NORMAL'
            if zv <= float(del_z):
                state = 'DELETION'
            elif zv >= float(dup_z):
                state = 'DUPLICATION'
            if state != 'NORMAL':
                long_rows.append({
                    'Sample': sample,
                    'CHR': df.iloc[i]['CHR'],
                    'START': df.iloc[i]['START'],
                    'END': df.iloc[i]['END'],
                    'ZScore': float(zv),
                    'CNVType': state,
                })

    out = pd.DataFrame(long_rows)
    save_df_to_tsv(out, f'{output_prefix}.cnv.tsv')

    zdf = out[['ZScore']].copy() if not out.empty else pd.DataFrame({'ZScore': [0.0]})
    plot_distribution_with_stats(
        data=zdf,
        col='ZScore',
        title='CNV Z-score Distribution',
        filename=f'{output_prefix}.zscore.png',
        x_label='Depth Z-score',
        y_label='CNV Count',
        thresholds=[
            {'value': float(del_z), 'label': f'DEL <= {del_z}', 'color': 'red', 'linestyle': '--'},
            {'value': float(dup_z), 'label': f'DUP >= {dup_z}', 'color': 'blue', 'linestyle': '--'},
        ],
    )


def main():
    p = argparse.ArgumentParser(description='WWWG2B-style CNV calling from normalized depth')
    p.add_argument('--input', required=True)
    p.add_argument('--output-prefix', required=True)
    p.add_argument('--del-z', type=float, default=-2.0)
    p.add_argument('--dup-z', type=float, default=2.0)
    args = p.parse_args()
    run_cnv_calling(args.input, args.output_prefix, del_z=args.del_z, dup_z=args.dup_z)


if __name__ == '__main__':
    main()
