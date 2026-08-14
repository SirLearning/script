"""Tests for breachr MAPQ-depth segment caller."""

from __future__ import annotations

import pandas as pd
import pytest

from breachr import CallConfig, call_mapq_depth_bins
from breachr.constants import LOGICAL_AUTOSOMES, MQ_BIN_LABELS, SEGMENT_OUTPUT_COLUMNS
from breachr.io import validate_mapq_depth_bin_df


def _empty_bin_table() -> pd.DataFrame:
    cols = [
        "logical_chr",
        "bin_idx",
        "bin_start",
        "bin_end",
        *MQ_BIN_LABELS,
        "total_reads",
        "bin_mean_depth",
    ]
    return pd.DataFrame(columns=cols)


def _make_bin(
    logical_chr: str,
    bin_idx: int,
    depth: float,
    *,
    total_reads: int = 10_000,
    mq_low: int = 500,
    mq_mid: int = 4000,
    mq_high: int = 4500,
) -> dict[str, object]:
    bin_size = 1_000_000
    start = bin_idx * bin_size + 1
    end = (bin_idx + 1) * bin_size
    mq = {
        "mq_0_10": mq_low // 2,
        "mq_10_20": mq_low - mq_low // 2,
        "mq_20_30": mq_mid // 2,
        "mq_30_40": mq_mid - mq_mid // 2,
        "mq_40_50": mq_high // 2,
        "mq_50_60": mq_high - mq_high // 2,
    }
    return {
        "logical_chr": logical_chr,
        "bin_idx": bin_idx,
        "bin_start": start,
        "bin_end": end,
        "total_reads": total_reads,
        "bin_mean_depth": depth,
        **mq,
    }


def _baseline_genome_df(base_depth: float = 10.0, *, bins_per_chr: int = 10) -> pd.DataFrame:
    """Disomic-like profile on all 21 logical chromosomes."""
    rows = []
    for logical_chr in LOGICAL_AUTOSOMES:
        for bin_idx in range(bins_per_chr):
            rows.append(_make_bin(logical_chr, bin_idx, base_depth))
    return pd.DataFrame(rows)


def test_validate_rejects_empty():
    with pytest.raises(ValueError, match="empty"):
        validate_mapq_depth_bin_df(_empty_bin_table())


def test_deletion_requires_low_depth_and_elevated_low_mq():
    df = _baseline_genome_df(base_depth=10.0)
    for bin_idx in (1, 2):
        mask = (df["logical_chr"] == "chr3B") & (df["bin_idx"] == bin_idx)
        df.loc[mask, "bin_mean_depth"] = 2.0
        df.loc[mask, "mq_0_10"] = 4000
        df.loc[mask, "mq_10_20"] = 3500
        df.loc[mask, "mq_40_50"] = 500
        df.loc[mask, "mq_50_60"] = 500

    calls = call_mapq_depth_bins(df, sample="test_del", config=CallConfig(min_segment_bins=2))
    del_calls = calls[calls["call_type"] == "deletion"]
    assert not del_calls.empty
    assert del_calls.iloc[0]["logical_chr"] == "chr3B"
    assert "misplaced" in del_calls.iloc[0]["evidence"]


def test_inter_homoeologue_disomic_depth_with_high_low_mq():
    df = _baseline_genome_df(base_depth=10.0)
    mask = (df["logical_chr"] == "chr2A") & (df["bin_idx"] == 1)
    df.loc[mask, "bin_mean_depth"] = 10.5
    df.loc[mask, "mq_0_10"] = 3500
    df.loc[mask, "mq_10_20"] = 3000
    df.loc[mask, "mq_40_50"] = 1000
    df.loc[mask, "mq_50_60"] = 1000

    config = CallConfig(min_segment_bins=1)
    calls = call_mapq_depth_bins(df, sample="test_ih", config=config)
    ih = calls[calls["call_type"] == "inter_homoeologue"]
    assert not ih.empty
    assert ih.iloc[0]["logical_chr"] == "chr2A"
    assert ih.iloc[0]["median_mq_low_frac"] > 0.3


def test_duplication_high_depth_with_normal_mq():
    df = _baseline_genome_df(base_depth=10.0)
    for bin_idx in (0, 1):
        mask = (df["logical_chr"] == "chr5D") & (df["bin_idx"] == bin_idx)
        df.loc[mask, "bin_mean_depth"] = 16.0
        df.loc[mask, "mq_0_10"] = 200
        df.loc[mask, "mq_10_20"] = 300

    calls = call_mapq_depth_bins(df, sample="test_dup", config=CallConfig(min_segment_bins=2))
    dup = calls[calls["call_type"] == "duplication"]
    assert not dup.empty
    assert dup.iloc[0]["logical_chr"] == "chr5D"
    assert dup.iloc[0]["median_depth_rel"] > 1.3


def test_aneuploidy_chr_loss_whole_chromosome():
    df = _baseline_genome_df(base_depth=10.0)
    chr_mask = df["logical_chr"] == "chr4A"
    df.loc[chr_mask, "bin_mean_depth"] = 4.5
    df.loc[chr_mask, "mq_0_10"] = 1500
    df.loc[chr_mask, "mq_10_20"] = 1200

    calls = call_mapq_depth_bins(df, sample="test_mono", config=CallConfig())
    loss = calls[calls["call_type"] == "aneuploidy_chr_loss"]
    assert len(loss) == 1
    assert loss.iloc[0]["logical_chr"] == "chr4A"
    assert loss.iloc[0]["n_bins"] == 10


def test_aneuploidy_chr_gain_whole_chromosome():
    df = _baseline_genome_df(base_depth=10.0)
    chr_mask = df["logical_chr"] == "chr6D"
    df.loc[chr_mask, "bin_mean_depth"] = 16.0
    df.loc[chr_mask, "mq_0_10"] = 200
    df.loc[chr_mask, "mq_10_20"] = 200

    calls = call_mapq_depth_bins(df, sample="test_tri", config=CallConfig())
    gain = calls[calls["call_type"] == "aneuploidy_chr_gain"]
    assert len(gain) == 1
    assert gain.iloc[0]["logical_chr"] == "chr6D"


def test_output_columns():
    df = _baseline_genome_df()
    calls = call_mapq_depth_bins(df, sample="test_cols")
    assert list(calls.columns) == list(SEGMENT_OUTPUT_COLUMNS)
