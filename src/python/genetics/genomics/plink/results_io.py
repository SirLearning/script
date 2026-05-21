"""Load PLINK / PLINK2 / GCTA tabular outputs (format parsing only; no association math)."""

import pandas as pd
import numpy as np

from infra.utils.io import load_df_generic, save_df_to_tsv


def load_plink2_eigenvec(path):
    df = load_df_generic(path)
    if df is None or df.empty:
        raise ValueError(f'Empty eigenvec: {path}')
    if '#FID' in df.columns:
        df = df.rename(columns={'#FID': 'FID'})
    if '#IID' in df.columns and 'Sample' not in df.columns:
        df = df.rename(columns={'#IID': 'Sample'})
    elif 'IID' in df.columns and 'Sample' not in df.columns:
        df = df.rename(columns={'IID': 'Sample'})
    return df


def load_plink2_eigenval(path):
    vals = []
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals.append(float(line.split()[0]))
    if not vals:
        raise ValueError(f'No eigenvalues in {path}')
    total = float(sum(vals))
    ratios = [v / total if total > 0 else 0.0 for v in vals]
    return pd.DataFrame({
        'PC': [f'PC{i + 1}' for i in range(len(vals))],
        'Eigenvalue': vals,
        'ExplainedVarianceRatio': ratios,
    })


def load_plink2_glm(path):
    df = load_df_generic(path)
    if df is None or df.empty:
        raise ValueError(f'Empty GLM result: {path}')
    rename = {'#CHROM': 'CHR', 'POS': 'BP', 'ID': 'SNP', 'P': 'PVALUE', 'OBS_CT': 'N'}
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    if 'PVALUE' not in df.columns and 'P' in df.columns:
        df = df.rename(columns={'P': 'PVALUE'})
    df['PVALUE'] = pd.to_numeric(df['PVALUE'], errors='coerce')
    df = df.dropna(subset=['PVALUE']).sort_values('PVALUE').reset_index(drop=True)
    df['Index'] = np.arange(1, len(df) + 1)
    df['LOG10P'] = -np.log10(df['PVALUE'].clip(lower=1e-300))
    return df


def load_gcta_gwas(path):
    df = load_df_generic(path)
    if df is None or df.empty:
        raise ValueError(f'Empty GCTA GWAS result: {path}')
    if 'P' in df.columns and 'PVALUE' not in df.columns:
        df = df.rename(columns={'P': 'PVALUE'})
    df['PVALUE'] = pd.to_numeric(df['PVALUE'], errors='coerce')
    df = df.dropna(subset=['PVALUE']).sort_values('PVALUE').reset_index(drop=True)
    df['Index'] = np.arange(1, len(df) + 1)
    df['LOG10P'] = -np.log10(df['PVALUE'].clip(lower=1e-300))
    return df


def summarize_tagsnp_from_prune(prune_in_path, output_prefix, max_tags=1000):
    tags = []
    with open(prune_in_path, encoding='utf-8') as handle:
        for line in handle:
            vid = line.strip()
            if vid:
                tags.append(vid)
            if len(tags) >= int(max_tags):
                break
    out = pd.DataFrame({'TagSNP': tags})
    save_df_to_tsv(out, f'{output_prefix}.tagsnp.tsv')
    return out


def prepare_plink_phenotype_table(phenotype_file, trait, output_path):
    pheno = load_df_generic(phenotype_file)
    if pheno is None or pheno.empty:
        raise ValueError(f'Empty phenotype: {phenotype_file}')
    if 'Sample' not in pheno.columns:
        raise ValueError('Phenotype table requires Sample column')
    if trait not in pheno.columns:
        raise ValueError(f'Missing trait column: {trait}')
    out = pheno[['Sample', trait]].copy()
    out.columns = ['IID', trait]
    out.insert(0, '#FID', out['IID'].astype(str))
    save_df_to_tsv(out, output_path)
    return output_path
