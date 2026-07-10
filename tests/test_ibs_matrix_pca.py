"""Tests for IBS-matrix MDS and PLINK eigenvec loading."""

from __future__ import annotations

import numpy as np
import pandas as pd

from genetics.genomics.plink.results_io import load_plink_eigenvec
from genetics.genomics.sample.ibs import _classical_mds, plot_ibs_matrix_mds


def test_classical_mds_recovers_plane_structure(tmp_path):
    rng = np.random.default_rng(0)
    pts = rng.normal(size=(20, 2))
    dist = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(axis=2))
    coords, explained = _classical_mds(dist, n_components=2)
    assert coords.shape == (20, 2)
    assert explained.shape == (2,)
    assert explained.sum() <= 1.0 + 1e-9
    assert explained[0] >= explained[1]


def test_load_plink1_eigenvec_headerless(tmp_path):
    path = tmp_path / "tiny.eigenvec"
    path.write_text(
        "0 s1 0.1 0.2\n0 s2 -0.1 0.3\n",
        encoding="utf-8",
    )
    df = load_plink_eigenvec(str(path))
    assert list(df.columns[:4]) == ["FID", "Sample", "PC1", "PC2"]
    assert df["Sample"].tolist() == ["s1", "s2"]


def test_plot_ibs_matrix_mds_smoke(tmp_path):
    ids = [f"s{i}" for i in range(5)]
    sim = np.full((5, 5), 0.9)
    np.fill_diagonal(sim, 1.0)
    mat_path = tmp_path / "tiny.mibs"
    id_path = tmp_path / "tiny.mibs.id"
    pd.DataFrame(sim).to_csv(mat_path, sep="\t", header=False, index=False)
    pd.DataFrame({"FID": [0] * 5, "IID": ids}).to_csv(
        id_path, sep="\t", header=False, index=False
    )

    out_prefix = tmp_path / "out"
    plot_ibs_matrix_mds(str(mat_path), str(id_path), str(out_prefix), n_components=4)

    assert (tmp_path / "out.ibs_mds.tsv").exists()
    assert (tmp_path / "out.ibs_mds.png").exists()
    assert (tmp_path / "out.ibs_mds.variance.png").exists()
    assert (tmp_path / "out.ibs_mds_tsne.png").exists()
    assert (tmp_path / "out.ibs_mds_umap.png").exists()
    mds = pd.read_csv(tmp_path / "out.ibs_mds.tsv", sep="\t")
    assert {"MDS1", "MDS2"}.issubset(mds.columns)
