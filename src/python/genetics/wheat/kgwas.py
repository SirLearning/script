import numpy as np
import pandas as pd
from infra.utils.io import load_df_generic, save_df_to_tsv
from infra.utils.graph import plot_gwas_qq


MIN_P = 1e-300


def _kmer_p(x, y):
    try:
        from scipy.stats import linregress
        return float(linregress(x, y).pvalue)
    except Exception:
        corr = np.corrcoef(x, y)[0, 1]
        if not np.isfinite(corr):
            return 1.0
        return max(MIN_P, 1.0 - min(1.0, abs(corr)))


def run_kgwas(kmer_matrix_file, phenotype_file, output_prefix, trait='Trait'):
    kmer = load_df_generic(kmer_matrix_file)
    pheno = load_df_generic(phenotype_file)
    if kmer is None or kmer.empty:
        raise ValueError('Empty k-mer matrix')
    if pheno is None or pheno.empty:
        raise ValueError('Empty phenotype table')

    if 'Sample' not in kmer.columns or 'Sample' not in pheno.columns:
        raise ValueError('Both files require Sample column')
    if trait not in pheno.columns:
        raise ValueError(f'Missing trait column: {trait}')

    merged = kmer.merge(pheno[['Sample', trait]], on='Sample', how='inner')
    y = pd.to_numeric(merged[trait], errors='coerce')

    rows = []
    for col in [c for c in kmer.columns if c != 'Sample']:
        x = pd.to_numeric(merged[col], errors='coerce')
        tmp = pd.DataFrame({'x': x, 'y': y}).dropna()
        if tmp.shape[0] < 3:
            continue
        p = _kmer_p(tmp['x'].values, tmp['y'].values)
        rows.append({'KMER': col, 'PVALUE': max(p, MIN_P)})

    if not rows:
        raise ValueError('No valid k-mer GWAS rows generated')

    out = pd.DataFrame(rows).sort_values('PVALUE').reset_index(drop=True)
    save_df_to_tsv(out, f'{output_prefix}.kgwas.tsv')

    exp = -np.log10((np.arange(1, len(out) + 1) - 0.5) / (len(out) + 0.5))
    obs = -np.log10(np.sort(out['PVALUE'].values))
    plot_gwas_qq(exp, obs, title='kGWAS QQ Plot', filename=f'{output_prefix}.qq.png')

