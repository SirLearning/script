"""Local relatedness ranking from SINGER ARG texts."""

from pathlib import Path

import numpy as np

from genetics.genomics.variant.singer_relatedness import (
    plot_singer_arg_relatedness,
    rank_nearest_neighbors,
    sample_tmrca_from_arg,
)


def _write_two_sample_arg(tmp_path: Path) -> Path:
    """Four haplotypes: within-sample coalescence at 1 and 2, between at 10."""
    prefix = tmp_path / "singer_A.chr1A.0_100"
    (tmp_path / "singer_A.chr1A.0_100_nodes_0.txt").write_text(
        "0\n0\n0\n0\n1.0\n2.0\n10.0\n",
        encoding="utf-8",
    )
    (tmp_path / "singer_A.chr1A.0_100_branches_0.txt").write_text(
        "0 100 4 0\n"
        "0 100 4 1\n"
        "0 100 5 2\n"
        "0 100 5 3\n"
        "0 100 6 4\n"
        "0 100 6 5\n",
        encoding="utf-8",
    )
    return tmp_path


def test_sample_tmrca_mean_of_four_haplotype_pairs(tmp_path: Path):
    singer_dir = _write_two_sample_arg(tmp_path)
    matrix = sample_tmrca_from_arg(
        singer_dir,
        "A.chr1A.0_100",
        sample_ids=["S0", "S1"],
        n_positions=3,
        arg_sample=0,
    )
    assert matrix.shape == (2, 2)
    assert matrix[0, 0] == 0.0
    assert matrix[1, 1] == 0.0
    assert matrix[0, 1] == 10.0
    assert matrix[1, 0] == 10.0


def test_rank_nearest_neighbors_orders_by_smaller_tmrca():
    names = ["A", "B", "C"]
    matrix = np.array(
        [
            [0.0, 3.0, 10.0],
            [3.0, 0.0, 8.0],
            [10.0, 8.0, 0.0],
        ]
    )
    neighbors = rank_nearest_neighbors(matrix, names, top_k=1)
    by_id = {row["sample"]: row["neighbor"] for row in neighbors}
    assert by_id["A"] == "B"
    assert by_id["B"] == "A"
    assert by_id["C"] == "B"
    assert neighbors[0]["tmrca"] == 3.0


def test_plot_singer_arg_relatedness_writes_rank_tables(tmp_path: Path):
    singer_dir = _write_two_sample_arg(tmp_path)
    windows = tmp_path / "singer_windows.tsv"
    windows.write_text(
        "subgenome\tplink_chr\tref_chrom\tstart_bp\tend_bp\tn_snps\n"
        "A\t1\tchr1A\t0\t100\t10\n",
        encoding="utf-8",
    )
    keep = tmp_path / "WAP_China.keep"
    keep.write_text("0\tS0\n0\tS1\n", encoding="utf-8")
    prefix = tmp_path / "out" / "all.variant.singer_arg"
    written = plot_singer_arg_relatedness(
        windows,
        singer_dir,
        keep,
        prefix,
        n_positions=3,
        top_k=1,
        top_pairs=1,
    )
    assert written["n_samples"] == 2
    assert written["n_windows"] == 1
    assert (tmp_path / "out" / "all.variant.singer_arg.neighbor.info.tsv").is_file()
    assert (tmp_path / "out" / "all.variant.singer_arg.pair_rank.info.tsv").is_file()
    pngs = list((tmp_path / "out").glob("all.variant.singer_arg*.png"))
    assert pngs
