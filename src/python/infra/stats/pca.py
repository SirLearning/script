import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


def compute_pca_coordinates(df, n_components=10, id_col='Sample'):
    """Compute PCA coordinates from a sample-by-feature table."""
    if id_col not in df.columns:
        raise ValueError(f"Missing id column: {id_col}")

    feature_df = df.drop(columns=[id_col]).apply(pd.to_numeric, errors='coerce').fillna(0.0)
    comp = max(1, min(int(n_components), feature_df.shape[0], feature_df.shape[1]))

    pca = PCA(n_components=comp)
    pcs = pca.fit_transform(feature_df.values)

    pc_cols = [f'PC{i}' for i in range(1, comp + 1)]
    pcs_df = pd.DataFrame(pcs, columns=pc_cols)
    pcs_df.insert(0, id_col, df[id_col].astype(str).values)

    eigen_df = pd.DataFrame({
        'PC': pc_cols,
        'ExplainedVarianceRatio': pca.explained_variance_ratio_,
        'ExplainedVariance': pca.explained_variance_,
    })
    return pcs_df, eigen_df


def compute_tsne_coordinates(df, n_components=2, random_state=42, perplexity=30.0, id_col='Sample'):
    """Compute t-SNE coordinates from PCA table (or any sample-by-feature table)."""
    if id_col not in df.columns:
        raise ValueError(f"Missing id column: {id_col}")

    feature_df = df.drop(columns=[id_col]).apply(pd.to_numeric, errors='coerce').fillna(0.0)
    n_samples = feature_df.shape[0]
    if n_samples < 3:
        raise ValueError('Need at least 3 samples for t-SNE')

    safe_perplexity = min(float(perplexity), max(2.0, (n_samples - 1) / 3.0))

    tsne = TSNE(
        n_components=int(n_components),
        random_state=int(random_state),
        perplexity=safe_perplexity,
        init='pca',
        learning_rate='auto'
    )
    embedding = tsne.fit_transform(feature_df.values)

    tsne_cols = [f'tSNE{i}' for i in range(1, int(n_components) + 1)]
    tsne_df = pd.DataFrame(embedding, columns=tsne_cols)
    tsne_df.insert(0, id_col, df[id_col].astype(str).values)
    return tsne_df
