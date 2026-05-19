"""Plot GWAS association results from PLINK2 --glm or GCTA."""

import numpy as np

from genetics.genomics.plink.results_io import load_gcta_gwas, load_plink2_glm
from infra.utils.graph import plot_gwas_qq, plot_scatter_with_thresholds
from infra.utils.io import save_df_to_tsv


def plot_gwas_results(association_path, output_prefix, source='plink2'):
    if source == 'gcta':
        out = load_gcta_gwas(association_path)
    else:
        out = load_plink2_glm(association_path)
    save_df_to_tsv(out, f'{output_prefix}.gwas.tsv')
    exp = -np.log10((np.arange(1, len(out) + 1) - 0.5) / (len(out) + 0.5))
    obs = -np.log10(np.sort(out['PVALUE'].values))
    plot_gwas_qq(exp, obs, title='GWAS QQ Plot', filename=f'{output_prefix}.qq.png')
    plot_scatter_with_thresholds(
        data=out, x_col='Index', y_col='LOG10P',
        title='GWAS Manhattan-like Scatter', filename=f'{output_prefix}.manhattan.png',
        xlabel='Variant Index', ylabel='-log10(P)',
    )


def plot_kgwas_results(glm_path, output_prefix):
    plot_gwas_results(glm_path, output_prefix, source='plink2')
