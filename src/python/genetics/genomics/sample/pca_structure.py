"""Plot population structure from PLINK2 PCA outputs."""

from genetics.genomics.plink.results_io import load_plink2_eigenval, load_plink2_eigenvec
from infra.utils.graph import plot_bar_chart, plot_scatter_with_thresholds
from infra.utils.io import save_df_to_tsv


def plot_population_structure(eigenvec_path, eigenval_path, output_prefix):
    pcs = load_plink2_eigenvec(eigenvec_path)
    eigen = load_plink2_eigenval(eigenval_path)
    save_df_to_tsv(pcs, f'{output_prefix}.pca.tsv')
    save_df_to_tsv(eigen, f'{output_prefix}.eigen.tsv')
    if {'PC1', 'PC2'}.issubset(pcs.columns):
        plot_scatter_with_thresholds(
            data=pcs, x_col='PC1', y_col='PC2',
            title='Population Structure PCA (PLINK2)',
            filename=f'{output_prefix}.pca.png', xlabel='PC1', ylabel='PC2',
        )
    plot_bar_chart(
        names=eigen['PC'].tolist(), values=eigen['ExplainedVarianceRatio'].tolist(),
        title='PCA Explained Variance Ratio (PLINK2)', ylabel='Variance Ratio',
        filename=f'{output_prefix}.variance.png', ylim=(0.0, 1.0),
    )
