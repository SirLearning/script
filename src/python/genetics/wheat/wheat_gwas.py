import numpy as np
import pandas as pd
from infra.utils.io import load_df_generic, save_df_to_tsv
from infra.utils.graph import plot_gwas_qq, plot_scatter_with_thresholds


MIN_P = 1e-300


def _linreg_pvalue(x, y):
    try:
        from scipy.stats import linregress
        res = linregress(x, y)
        return float(res.slope), float(res.pvalue)
    except Exception:
        if len(x) < 3:
            return float('nan'), 1.0
        corr = np.corrcoef(x, y)[0, 1]
        if not np.isfinite(corr):
            return float('nan'), 1.0
        p_approx = max(MIN_P, 1.0 - min(1.0, abs(corr)))
        return float(corr), float(p_approx)


def run_gwas(genotype_file, phenotype_file, output_prefix, trait='Trait'):
    geno = load_df_generic(genotype_file)
    pheno = load_df_generic(phenotype_file)
    if geno is None or geno.empty:
        raise ValueError('Empty genotype table')
    if pheno is None or pheno.empty:
        raise ValueError('Empty phenotype table')
    if 'Sample' not in geno.columns or 'Sample' not in pheno.columns:
        raise ValueError('Both genotype and phenotype tables require Sample column')
    if trait not in pheno.columns:
        raise ValueError(f'Phenotype file missing trait column: {trait}')

    merged = geno.merge(pheno[['Sample', trait]], on='Sample', how='inner')
    y = pd.to_numeric(merged[trait], errors='coerce')
    marker_cols = [c for c in geno.columns if c != 'Sample']

    rows = []
    for mk in marker_cols:
        x = pd.to_numeric(merged[mk], errors='coerce')
        local = pd.DataFrame({'x': x, 'y': y}).dropna()
        if local.shape[0] < 3:
            continue
        slope, pval = _linreg_pvalue(local['x'].values, local['y'].values)
        rows.append({'SNP': mk, 'BETA': slope, 'PVALUE': max(float(pval), MIN_P)})

    if not rows:
        raise ValueError('No valid marker association results generated')

    out = pd.DataFrame(rows).sort_values('PVALUE').reset_index(drop=True)
    out['Index'] = np.arange(1, len(out) + 1)
    out['LOG10P'] = -np.log10(out['PVALUE'].clip(lower=MIN_P))

    save_df_to_tsv(out, f'{output_prefix}.gwas.tsv')

    exp = -np.log10((np.arange(1, len(out) + 1) - 0.5) / (len(out) + 0.5))
    obs = -np.log10(np.sort(out['PVALUE'].values))
    plot_gwas_qq(exp, obs, title='GWAS QQ Plot', filename=f'{output_prefix}.qq.png')

    plot_scatter_with_thresholds(
        data=out,
        x_col='Index',
        y_col='LOG10P',
        title='GWAS Manhattan-like Scatter',
        filename=f'{output_prefix}.manhattan.png',
        xlabel='Marker Index',
        ylabel='-log10(P)',
    )
