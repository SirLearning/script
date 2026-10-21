"""Tests for mapq_bin_popdep MQ30 population depth aggregation."""

import pandas as pd
import numpy as np
import pytest

from genetics.germplasm.sample.mapq_bin_popdep import (
    POPDEP_BIN_COLUMNS,
    add_sample_mq30_relative_depth,
    add_sample_relative_depth,
    build_mapq_bin_popdep,
    build_mapq_bin_popdep_from_dir,
    estimate_mq30_bin_depth,
    flag_bin_cnv_from_depth_rel,
    sample_bin_depth,
)
from genetics.germplasm.sample.mapq_bin_popdep_plot import (
    CNV_FRAC_BIN_LABELS,
    assign_cnv_frac_bin,
    plot_mapq_bin_popdep_qc_scatters,
)


def _sample_row(
    logical_chr: str,
    bin_idx: int,
    depth: float,
    mq30_frac: float,
    total_reads: int = 1000,
) -> dict:
    mq30_reads = int(total_reads * mq30_frac)
    mq30_each = mq30_reads // 3
    return {
        "logical_chr": logical_chr,
        "bin_idx": bin_idx,
        "bin_start": bin_idx * 1_000_000 + 1,
        "bin_end": (bin_idx + 1) * 1_000_000,
        "mq_0_10": 0,
        "mq_10_20": 0,
        "mq_20_30": total_reads - mq30_reads,
        "mq_30_40": mq30_each,
        "mq_40_50": mq30_each,
        "mq_50_60": mq30_reads - 2 * mq30_each,
        "total_reads": total_reads,
        "bin_mean_depth": depth,
    }


def test_estimate_mq30_bin_depth_fraction():
    df = pd.DataFrame([_sample_row("chr1A", 0, depth=10.0, mq30_frac=0.5)])
    est = estimate_mq30_bin_depth(df)
    assert est.iloc[0] == 5.0


def test_sample_bin_depth_all_vs_mq30():
    df = pd.DataFrame([_sample_row("chr1A", 0, depth=10.0, mq30_frac=0.5)])
    assert sample_bin_depth(df, depth_mode="all").iloc[0] == 10.0
    assert sample_bin_depth(df, depth_mode="mq30").iloc[0] == 5.0


def test_build_mapq_bin_popdep_aggregates_two_samples():
    rows = []
    for sample, depth in [("S1", 10.0), ("S2", 14.0)]:
        row = _sample_row("chr1A", 0, depth=depth, mq30_frac=1.0)
        row["sample"] = sample
        row = add_sample_mq30_relative_depth(pd.DataFrame([row]))
        rows.append(row)
    long = pd.concat(rows, ignore_index=True)
    grid = build_mapq_bin_popdep(long, min_samples=2)
    assert list(grid.columns) == list(POPDEP_BIN_COLUMNS)
    assert grid.iloc[0]["depth_mean"] == 12.0
    assert grid.iloc[0]["depth_sd"] == pytest.approx(2.828427, rel=1e-5)
    assert grid.iloc[0]["n_samples"] == 2
    assert abs(grid.iloc[0]["rel_mean"] - 1.0) < 1e-6


def test_flag_bin_cnv_from_depth_rel():
    rel = pd.Series([0.4, 0.9, 1.5, np.nan])
    flags = flag_bin_cnv_from_depth_rel(rel)
    assert flags.tolist() == [True, False, True, False]


def test_streaming_matches_batch_aggregation(tmp_path):
    from breachr.constants import MAPQ_DEPTH_BIN_COLUMNS

    rows = []
    for sample, depth in [("S1", 10.0), ("S2", 14.0), ("S3", 12.0)]:
        row = _sample_row("chr1A", 0, depth=depth, mq30_frac=1.0)
        row["sample"] = sample
        row = add_sample_mq30_relative_depth(pd.DataFrame([row]))
        rows.append(row)
    long = pd.concat(rows, ignore_index=True)
    batch = build_mapq_bin_popdep(long, min_samples=2)

    for sample, depth in [("S1", 10.0), ("S2", 14.0), ("S3", 12.0)]:
        base = _sample_row("chr1A", 0, depth=depth, mq30_frac=1.0)
        path = tmp_path / f"{sample}.mapq_depth_bin.tsv"
        pd.DataFrame([base])[list(MAPQ_DEPTH_BIN_COLUMNS)].to_csv(
            path, sep="\t", index=False
        )

    stream_grid, _ = build_mapq_bin_popdep_from_dir(
        tmp_path, max_workers=1, min_samples=2, keep_long=False, engine="streaming"
    )
    pd.testing.assert_frame_equal(
        batch.reset_index(drop=True),
        stream_grid.reset_index(drop=True),
        check_exact=False,
        rtol=1e-9,
    )


def test_polars_matches_pandas_on_sample_dir(tmp_path):
    from breachr.constants import MAPQ_DEPTH_BIN_COLUMNS

    depths = [("S1", 10.0), ("S2", 14.0), ("S3", 12.0)]
    for sample, depth in depths:
        base = _sample_row("chr1A", 0, depth=depth, mq30_frac=1.0)
        path = tmp_path / f"{sample}.mapq_depth_bin.tsv"
        pd.DataFrame([base])[list(MAPQ_DEPTH_BIN_COLUMNS)].to_csv(
            path, sep="\t", index=False
        )

    polars_grid, _ = build_mapq_bin_popdep_from_dir(
        tmp_path, max_workers=1, min_samples=2, keep_long=False, engine="polars"
    )
    pandas_grid, _ = build_mapq_bin_popdep_from_dir(
        tmp_path, max_workers=1, min_samples=2, keep_long=False, engine="pandas"
    )
    pd.testing.assert_frame_equal(
        polars_grid.reset_index(drop=True),
        pandas_grid.reset_index(drop=True),
        check_exact=False,
        rtol=1e-9,
    )


def test_assign_cnv_frac_bin_intervals():
    frac = pd.Series([0.0, 0.19, 0.2, 0.5, 0.99, 1.0])
    bins = assign_cnv_frac_bin(frac)
    assert bins.astype(str).tolist() == [
        "0.0-0.2",
        "0.0-0.2",
        "0.2-0.4",
        "0.4-0.6",
        "0.8-1.0",
        "0.8-1.0",
    ]


def test_plot_mapq_bin_popdep_qc_scatters_writes_pngs(tmp_path):
    grid = pd.DataFrame(
        [
            {
                "logical_chr": "chr1A",
                "bin_idx": 0,
                "bin_start": 1,
                "bin_end": 1_000_000,
                "n_samples": 100,
                "depth_mean": 10.0,
                "depth_sd": 1.0,
                "rel_mean": 1.0,
                "rel_sd": 0.1,
                "n_cnv": 5,
                "cnv_frac": 0.05,
            },
            {
                "logical_chr": "chr1A",
                "bin_idx": 1,
                "bin_start": 1_000_001,
                "bin_end": 2_000_000,
                "n_samples": 100,
                "depth_mean": 20.0,
                "depth_sd": 2.0,
                "rel_mean": 1.2,
                "rel_sd": 0.2,
                "n_cnv": 0,
                "cnv_frac": 0.0,
            },
        ]
    )
    prefix = tmp_path / "mapq30_1mb"
    paths = plot_mapq_bin_popdep_qc_scatters(grid, prefix)
    assert len(paths) == 8
    assert all(p.exists() and p.stat().st_size > 0 for p in paths)
    assert any("cv_vs_mean" in p.name for p in paths)
