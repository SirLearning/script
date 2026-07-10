"""Plot population structure from PLINK2 PCA outputs (PCA + t-SNE on PCs, optional sample groups)."""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
from sklearn.manifold import TSNE

from genetics.germplasm.sample.anno import anno_group
from genetics.genomics.plink.results_io import load_plink2_eigenval, load_plink2_eigenvec
from infra.utils.graph import plot_bar_chart, plot_scatter_with_thresholds
from infra.utils.io import save_df_to_tsv


def _ordered_pc_columns(df: pd.DataFrame) -> list[str]:
    """Return PC column names sorted as PC1, PC2, …"""

    def pc_key(name: str) -> int:
        if not isinstance(name, str) or not name.startswith("PC"):
            return 10**9
        try:
            return int(name[2:])
        except ValueError:
            return 10**9

    pcs = [c for c in df.columns if isinstance(c, str) and c.startswith("PC")]
    return sorted(pcs, key=pc_key)


def _attach_sample_groups(df: pd.DataFrame, group_file: str | None) -> pd.DataFrame:
    """Merge ``Group`` labels via ``anno_group`` (same path as ref_ibs / stats)."""

    if not group_file or not os.path.exists(group_file):
        if group_file:
            print(f"[Warning] Group file not found: {group_file}")
        return df
    annotated = anno_group(df, group_file=group_file, save_tsv=False)
    if annotated is None:
        return df
    if "Group" in annotated.columns:
        print(f"[Info] Sample groups: {annotated['Group'].nunique()} categories")
    return annotated


def _tsne_perplexity(n_samples: int) -> float:
    """Pick a valid perplexity for barnes_hut (must be well below n_samples)."""

    if n_samples < 3:
        return 1.0
    upper = max(2.0, (n_samples - 1) / 3.0 - 1e-6)
    return float(min(30.0, max(2.0, upper * 0.9)))


def _fit_tsne_2d(
    X: np.ndarray,
    *,
    max_iter: int,
    random_state: int,
    n_jobs: int,
) -> np.ndarray:
    """Run 2D t-SNE on row matrix X."""

    perp = _tsne_perplexity(X.shape[0])
    kw: dict = {
        "n_components": 2,
        "perplexity": perp,
        "max_iter": max_iter,
        "random_state": random_state,
        "init": "pca",
    }
    if n_jobs and n_jobs > 0:
        try:
            return TSNE(**kw, n_jobs=n_jobs).fit_transform(X)
        except TypeError:
            pass
    try:
        return TSNE(**kw, learning_rate="auto").fit_transform(X)
    except TypeError:
        return TSNE(**kw).fit_transform(X)


def plot_population_structure(
    eigenvec_path: str,
    eigenval_path: str,
    output_prefix: str,
    *,
    group_file: str | None = None,
    tsne_n_input_pcs: int = 10,
    tsne_max_iter: int = 1000,
    tsne_random_state: int = 42,
    tsne_n_jobs: int = 1,
    plot_tsne: bool = True,
    tsne_max_samples: int = 3000,
    pca_title: str = "Population Structure PCA (PLINK2)",
) -> None:
    """
    Load PLINK eigenvec/eigenval, attach optional sample ``Group`` labels,
    save tables, PCA scatter + variance bar, and optional t-SNE on leading PCs.

    Parameters
    ----------
    plot_tsne
        When False, skip t-SNE (recommended for large cohorts in PLINK1 partial stats).
    tsne_max_samples
        Skip t-SNE when sample count exceeds this threshold (avoids sklearn segfaults).
    pca_title
        Title for the PC1 vs PC2 scatter plot.
    """
    pcs = load_plink2_eigenvec(eigenvec_path)
    eigen = load_plink2_eigenval(eigenval_path)
    pcs = _attach_sample_groups(pcs, group_file)
    group_col = "Group" if "Group" in pcs.columns else None

    save_df_to_tsv(pcs, f"{output_prefix}.pca.tsv")
    save_df_to_tsv(eigen, f"{output_prefix}.eigen.tsv")
    if {"PC1", "PC2"}.issubset(pcs.columns):
        plot_scatter_with_thresholds(
            data=pcs,
            x_col="PC1",
            y_col="PC2",
            title=pca_title,
            filename=f"{output_prefix}.pca.png",
            xlabel="PC1",
            ylabel="PC2",
            group_col=group_col,
        )
    plot_bar_chart(
        names=eigen["PC"].tolist(),
        values=eigen["ExplainedVarianceRatio"].tolist(),
        title="PCA Explained Variance Ratio (PLINK2)",
        ylabel="Variance Ratio",
        filename=f"{output_prefix}.pca.variance.png",
        ylim=(0.0, 1.0),
    )

    pc_cols = _ordered_pc_columns(pcs)
    if not plot_tsne:
        print("[Info] plot_tsne=False; skipping t-SNE.")
        return
    if len(pc_cols) < 2:
        print("[Warning] Fewer than two PC columns; skipping t-SNE.")
        return
    k = max(2, min(int(tsne_n_input_pcs), len(pc_cols)))
    X = pcs[pc_cols[:k]].to_numpy(dtype=float)
    n = X.shape[0]
    if n > int(tsne_max_samples):
        print(f"[Info] Sample count {n} > tsne_max_samples={tsne_max_samples}; skipping t-SNE.")
        return
    if n < 3:
        print("[Warning] Fewer than three samples; skipping t-SNE.")
        return
    if not np.isfinite(X).all():
        print("[Warning] Non-finite PC values; skipping t-SNE.")
        return

    Z = _fit_tsne_2d(
        X,
        max_iter=int(tsne_max_iter),
        random_state=int(tsne_random_state),
        n_jobs=int(tsne_n_jobs),
    )
    id_col = "Sample" if "Sample" in pcs.columns else None
    if id_col is None:
        id_col = "IID" if "IID" in pcs.columns else None
    if id_col is None:
        tsne_df = pd.DataFrame({"idx": np.arange(n, dtype=int), "TSNE1": Z[:, 0], "TSNE2": Z[:, 1]})
    else:
        tsne_df = pd.DataFrame(
            {id_col: pcs[id_col].astype(str).values, "TSNE1": Z[:, 0], "TSNE2": Z[:, 1]},
        )
    if group_col:
        tsne_df[group_col] = pcs[group_col].astype(str).values
    save_df_to_tsv(tsne_df, f"{output_prefix}.tsne.tsv")
    plot_scatter_with_thresholds(
        data=tsne_df,
        x_col="TSNE1",
        y_col="TSNE2",
        title="t-SNE on PLINK2 PCs (sklearn)",
        filename=f"{output_prefix}.tsne.png",
        xlabel="t-SNE 1",
        ylabel="t-SNE 2",
        group_col=group_col,
    )
