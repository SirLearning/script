"""Deprecated: genotype matrix export for sklearn-based wheat analytics.

Integrated wheat modes now use PLINK2 / GCTA in Nextflow for computation and Python for plotting only.
Kept for ad-hoc utilities; do not wire into integrated_wheat.nf.
"""

import os
import random

import pandas as pd

from genetics.wheat._deprecation import warn_deprecated_shim

warn_deprecated_shim("plink2_matrix", "genetics.genomics.plink (avoid sklearn matrix GWAS)")

from infra.utils.io import load_df_generic, save_df_to_tsv

_META_COLS = {'FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'}


def _read_pvar_variant_ids(pvar_path):
    """Collect variant IDs from a PLINK2 .pvar (skip ## and #CHROM header lines)."""
    variant_ids = []
    with open(pvar_path, encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('##') or line.startswith('#CHROM'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                variant_ids.append(parts[2])
    return variant_ids


def write_variant_extract_list(pvar_path, max_variants, seed):
    """Write variant IDs (pvar col 3) for plink2 --extract; subsample when over max_variants."""
    variant_ids = _read_pvar_variant_ids(pvar_path)
    if not variant_ids:
        raise ValueError(f'No variant IDs found in pvar: {pvar_path}')
    if max_variants > 0 and len(variant_ids) > max_variants:
        rng = random.Random(int(seed))
        variant_ids = rng.sample(variant_ids, int(max_variants))

    out_path = f'{os.path.splitext(pvar_path)[0]}.extract.id'
    with open(out_path, 'w', encoding='utf-8') as handle:
        for vid in variant_ids:
            handle.write(f'{vid}\n')
    return out_path, len(variant_ids)


def plink2_export_raw_to_geno_matrix(raw_path, output_path, sample_col='IID'):
    """Turn PLINK2 .raw (--export A) into TSV with Sample column and numeric variant columns."""
    df = load_df_generic(raw_path)
    if df is None or df.empty:
        raise ValueError(f'Empty PLINK2 export: {raw_path}')

    if sample_col not in df.columns:
        raise ValueError(f'Missing sample column {sample_col!r} in {raw_path}')

    df = df.rename(columns={sample_col: 'Sample'})
    drop_cols = [c for c in df.columns if c in _META_COLS]
    df = df.drop(columns=drop_cols, errors='ignore')

    feature_cols = [c for c in df.columns if c != 'Sample']
    df[feature_cols] = df[feature_cols].apply(pd.to_numeric, errors='coerce').fillna(0.0)
    save_df_to_tsv(df, output_path)
    return output_path
