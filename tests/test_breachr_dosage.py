"""Tests for improved breachr dosage and subgenome calling."""

from __future__ import annotations

import pandas as pd
import pytest

from breachr import CallConfig, call_mapq_depth_bins
from breachr.caller import call_mapq_depth_bin_tsv
from breachr.constants import LOGICAL_AUTOSOMES, MQ_BIN_LABELS
from breachr.dosage import chr_dosage_profile, is_partial_cnv_pattern, is_uniform_monosomic


def _make_bin(logical_chr, bin_idx, depth, *, mq_low=500, mq_mid=4000, mq_high=4500, reads=10000):
    start = bin_idx * 1_000_000 + 1
    end = (bin_idx + 1) * 1_000_000
    return {
        "logical_chr": logical_chr,
        "bin_idx": bin_idx,
        "bin_start": start,
        "bin_end": end,
        "mq_0_10": mq_low // 2,
        "mq_10_20": mq_low - mq_low // 2,
        "mq_20_30": mq_mid // 2,
        "mq_30_40": mq_mid - mq_mid // 2,
        "mq_40_50": mq_high // 2,
        "mq_50_60": mq_high - mq_high // 2,
        "total_reads": reads,
        "bin_mean_depth": depth,
    }


def _baseline_df(base=10.0, bins_per_chr=10):
    rows = []
    for lc in LOGICAL_AUTOSOMES:
        for i in range(bins_per_chr):
            rows.append(_make_bin(lc, i, base))
    return pd.DataFrame(rows)


def test_uniform_monosomic_not_fragmented():
    df = _baseline_df(base=10.0, bins_per_chr=50)
    mask = df["logical_chr"] == "chr6B"
    df.loc[mask, "bin_mean_depth"] = 5.2
    calls = call_mapq_depth_bins(df, sample="test_mono_uniform", config=CallConfig(min_segment_bins=2))
    sub = calls[calls["logical_chr"] == "chr6B"]
    assert len(sub) == 1
    assert sub.iloc[0]["call_type"] == "aneuploidy_chr_loss"


def test_partial_cnv_not_whole_chr_monosomy():
    df = _baseline_df(base=10.0, bins_per_chr=40)
    mask = (df["logical_chr"] == "chr2B") & (df["bin_idx"].between(10, 28))
    df.loc[mask, "bin_mean_depth"] = 3.0
    calls = call_mapq_depth_bins(df, sample="test_partial", config=CallConfig(min_segment_bins=2))
    sub = calls[calls["logical_chr"] == "chr2B"]
    assert not (sub["call_type"] == "aneuploidy_chr_loss").any()
    assert (sub["call_type"] == "deletion").any()


@pytest.mark.skipif(
    not __import__("pathlib").Path(
        "/data1/dazheng_tusr1/vmap4.VCF.v1/benchmark/stats/aneuploidy_check/info/mapq_check/ABD_0249.mapq_depth_bin.tsv",
    ).exists(),
    reason="vmap4 benchmark MAPQ TSV not available",
)
def test_abd_0249_chr6b_monosomic_integration():
    path = "/data1/dazheng_tusr1/vmap4.VCF.v1/benchmark/stats/aneuploidy_check/info/mapq_check/ABD_0249.mapq_depth_bin.tsv"
    calls = call_mapq_depth_bin_tsv(path, sample="ABD_0249")
    sub = calls[calls["logical_chr"] == "chr6B"]
    assert len(sub) == 1
    assert sub.iloc[0]["call_type"] == "aneuploidy_chr_loss"
