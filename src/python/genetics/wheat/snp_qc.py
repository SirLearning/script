import argparse
import pandas as pd
from infra.utils.io import load_df_generic, save_df_to_tsv
from infra.utils.graph import plot_distribution_with_stats


def run_snp_qc(input_file, output_prefix, maf=0.05, max_missing=0.1, min_qual=30.0):
    df = load_df_generic(input_file)
    if df is None or df.empty:
        raise ValueError(f'No SNP records loaded from: {input_file}')

    for col in ('MAF', 'MISSING_RATE', 'QUAL'):
        if col not in df.columns:
            raise ValueError(f'Missing required column: {col}')

    df['MAF'] = pd.to_numeric(df['MAF'], errors='coerce')
    df['MISSING_RATE'] = pd.to_numeric(df['MISSING_RATE'], errors='coerce')
    df['QUAL'] = pd.to_numeric(df['QUAL'], errors='coerce')

    qc = df[
        (df['MAF'] >= float(maf))
        & (df['MISSING_RATE'] <= float(max_missing))
        & (df['QUAL'] >= float(min_qual))
    ].copy()

    summary = pd.DataFrame([
        {'Metric': 'InputVariants', 'Value': len(df)},
        {'Metric': 'RetainedVariants', 'Value': len(qc)},
        {'Metric': 'RetentionRate', 'Value': (len(qc) / len(df)) if len(df) else 0.0},
        {'Metric': 'MAFThreshold', 'Value': float(maf)},
        {'Metric': 'MaxMissingRate', 'Value': float(max_missing)},
        {'Metric': 'MinQUAL', 'Value': float(min_qual)},
    ])

    save_df_to_tsv(qc, f'{output_prefix}.filtered.tsv')
    save_df_to_tsv(summary, f'{output_prefix}.summary.tsv')

    plot_distribution_with_stats(
        data=qc if not qc.empty else df,
        col='MAF',
        title='SNP MAF Distribution After QC',
        filename=f'{output_prefix}.maf.png',
        x_label='Minor Allele Frequency',
        y_label='Variant Count',
        thresholds=[{'value': float(maf), 'label': f'MAF >= {maf}', 'color': 'red', 'linestyle': '--'}],
    )


def main():
    p = argparse.ArgumentParser(description='WWWG2B-style SNP calling QC post-filter')
    p.add_argument('--input', required=True, help='SNP table with MAF/MISSING_RATE/QUAL columns')
    p.add_argument('--output-prefix', required=True)
    p.add_argument('--maf', type=float, default=0.05)
    p.add_argument('--max-missing', type=float, default=0.1)
    p.add_argument('--min-qual', type=float, default=30.0)
    args = p.parse_args()
    run_snp_qc(args.input, args.output_prefix, maf=args.maf, max_missing=args.max_missing, min_qual=args.min_qual)


if __name__ == '__main__':
    main()
