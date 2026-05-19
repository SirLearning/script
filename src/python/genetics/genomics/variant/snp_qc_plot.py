"""Summarize and plot SNP QC from PLINK2 --freq / --missing outputs."""

import pandas as pd

from infra.utils.graph import plot_distribution_with_stats
from infra.utils.io import load_df_generic, save_df_to_tsv


def _afreq_to_maf(afreq_df):
    df = afreq_df.copy()
    if 'ALT_FREQS' in df.columns:
        df['MAF'] = pd.to_numeric(df['ALT_FREQS'], errors='coerce')
        df['MAF'] = df['MAF'].apply(lambda v: min(v, 1.0 - v) if pd.notna(v) else v)
    elif 'MAF' in df.columns:
        df['MAF'] = pd.to_numeric(df['MAF'], errors='coerce')
    else:
        raise ValueError('afreq file missing ALT_FREQS or MAF column')
    if 'ID' not in df.columns and len(df.columns) > 2:
        df = df.rename(columns={df.columns[2]: 'ID'})
    return df


def _vmiss_to_rate(vmiss_df):
    df = vmiss_df.copy()
    miss_col = next((c for c in ('F_MISS', 'MISSING_RATE', 'MISSING_CT') if c in df.columns), None)
    if miss_col is None:
        raise ValueError('vmiss file missing F_MISS / MISSING_RATE column')
    df['MISSING_RATE'] = pd.to_numeric(df[miss_col], errors='coerce')
    if 'ID' not in df.columns and len(df.columns) > 2:
        df = df.rename(columns={df.columns[2]: 'ID'})
    return df


def summarize_and_plot_snp_qc(afreq_path, vmiss_path, output_prefix, maf=0.05, max_missing=0.1):
    afreq = _afreq_to_maf(load_df_generic(afreq_path))
    vmiss = _vmiss_to_rate(load_df_generic(vmiss_path))
    merged = afreq[['ID', 'MAF']].merge(vmiss[['ID', 'MISSING_RATE']], on='ID', how='inner')
    qc = merged[(merged['MAF'] >= float(maf)) & (merged['MISSING_RATE'] <= float(max_missing))].copy()
    summary = pd.DataFrame([
        {'Metric': 'InputVariants', 'Value': len(merged)},
        {'Metric': 'RetainedVariants', 'Value': len(qc)},
        {'Metric': 'RetentionRate', 'Value': (len(qc) / len(merged)) if len(merged) else 0.0},
        {'Metric': 'MAFThreshold', 'Value': float(maf)},
        {'Metric': 'MaxMissingRate', 'Value': float(max_missing)},
    ])
    save_df_to_tsv(qc, f'{output_prefix}.filtered.tsv')
    save_df_to_tsv(summary, f'{output_prefix}.summary.tsv')
    plot_distribution_with_stats(
        data=qc if not qc.empty else merged, col='MAF',
        title='SNP MAF Distribution After PLINK2 QC',
        filename=f'{output_prefix}.maf.png',
        x_label='Minor Allele Frequency', y_label='Variant Count',
        thresholds=[{'value': float(maf), 'label': f'MAF >= {maf}', 'color': 'red', 'linestyle': '--'}],
    )
