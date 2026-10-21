"""Unit tests for Wu 2025 LD block size parsing and plotting."""

from pathlib import Path

import matplotlib.image as mpimg
import numpy as np

from genetics.germplasm.meta.ld_block_size_plot import (
    load_wu_ld_block_table,
    run_ld_block_size_plot,
    summarize_ld_block_sizes,
)
from infra.utils.graph import local_xy_density

LD_BLOCK_CSV = """Table S26. Detailed information of LD blocks referred to Figure S2.,,,,,
Block number,Chromosome,Start (V2.1),End (V2.1),Number of Gene,Block size (Mb)
Block_1,Chr1A,1142753,8866181,315,7.723428
Block_2,Chr1A,8875965,12077625,180,3.20166
Block_3,Chr1A,12082077,13557312,77,1.475235
Block_4,Chr1A,13567659,15261628,99,0.5
"""


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_load_wu_ld_block_table(tmp_path):
    csv_path = _write(tmp_path / "blocks.csv", LD_BLOCK_CSV)
    df = load_wu_ld_block_table(csv_path)
    assert len(df) == 4
    assert list(df.columns) == [
        "block_number",
        "chromosome",
        "start_bp",
        "end_bp",
        "n_genes",
        "block_size_mb",
    ]
    assert df.loc[0, "block_size_mb"] == 7.723428


def test_summarize_ld_block_sizes(tmp_path):
    csv_path = _write(tmp_path / "blocks.csv", LD_BLOCK_CSV)
    block_df = load_wu_ld_block_table(csv_path)
    summary = summarize_ld_block_sizes(block_df)
    assert summary.loc[0, "n_blocks"] == 4
    assert summary.loc[0, "min_mb"] == 0.5
    assert summary.loc[0, "max_mb"] == 7.723428


def test_run_ld_block_size_plot_writes_outputs(tmp_path):
    csv_path = _write(tmp_path / "blocks.csv", LD_BLOCK_CSV)
    out_dir = tmp_path / "out"
    outputs = run_ld_block_size_plot(csv_path, out_dir)

    for key in (
        "blocks_tsv",
        "summary_tsv",
        "distribution_plot",
        "distribution_logy_plot",
        "size_vs_genes_plot",
    ):
        assert outputs[key].exists()

    img = mpimg.imread(str(outputs["size_vs_genes_plot"]))
    height, width = img.shape[:2]
    assert abs((width / height) - 2.0) < 0.08


def test_local_xy_density_cluster_exceeds_outlier():
    x = np.concatenate([np.zeros(40), [8.0]])
    y = np.concatenate([np.zeros(40), [8.0]])
    z = local_xy_density(x, y)
    assert z[0] > z[-1]
