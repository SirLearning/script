import pandas as pd
from infra.stats import compute_pca_coordinates, compute_tsne_coordinates
from infra.utils.io import load_df_generic, save_df_to_tsv
from infra.utils.graph import plot_scatter_with_thresholds, plot_bar_chart


def run_population_structure(input_file, output_prefix, n_pcs=10, tsne_perplexity=30.0):
    df = load_df_generic(input_file)
    if df is None or df.empty:
        raise ValueError(f'No genotype matrix loaded from: {input_file}')
    if 'Sample' not in df.columns:
        raise ValueError('Input requires a Sample column')

    pcs, eigen = compute_pca_coordinates(df, n_components=n_pcs, id_col='Sample')
    tsne = compute_tsne_coordinates(
        pcs[['Sample'] + [c for c in pcs.columns if c.startswith('PC')]],
        n_components=2,
        perplexity=tsne_perplexity,
        id_col='Sample'
    )

    save_df_to_tsv(pcs, f'{output_prefix}.pca.tsv')
    save_df_to_tsv(eigen, f'{output_prefix}.eigen.tsv')
    save_df_to_tsv(tsne, f'{output_prefix}.tsne.tsv')

    if {'PC1', 'PC2'}.issubset(pcs.columns):
        plot_scatter_with_thresholds(
            data=pcs,
            x_col='PC1',
            y_col='PC2',
            title='Population Structure PCA',
            filename=f'{output_prefix}.pca.png',
            xlabel='PC1',
            ylabel='PC2',
        )

    if {'tSNE1', 'tSNE2'}.issubset(tsne.columns):
        plot_scatter_with_thresholds(
            data=tsne,
            x_col='tSNE1',
            y_col='tSNE2',
            title='Population Structure t-SNE',
            filename=f'{output_prefix}.tsne.png',
            xlabel='tSNE1',
            ylabel='tSNE2',
        )

    plot_bar_chart(
        names=eigen['PC'].tolist(),
        values=eigen['ExplainedVarianceRatio'].tolist(),
        title='PCA Explained Variance Ratio',
        ylabel='Variance Ratio',
        filename=f'{output_prefix}.variance.png',
        ylim=(0.0, 1.0),
    )
