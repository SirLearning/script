"""Tests for mosdepth vs breachr reconciliation."""

from __future__ import annotations

import pandas as pd

from breachr.constants import SEGMENT_OUTPUT_COLUMNS
from breachr.reconcile import (
    ReconcileConfig,
    parse_abnormal_chr_list,
    reconcile_flagged_chr,
    reconcile_flagged_samples,
)


def test_parse_abnormal_chr_list_multi():
    text = (
        "chr1A:tetrasomic+(2.32),chr1B:null(0.09),chr2A:trisomic(1.49),"
        "chr4B:monosomic(0.59)"
    )
    entries = parse_abnormal_chr_list(text)
    assert len(entries) == 4
    assert entries[0]["logical_chr"] == "chr1A"
    assert entries[0]["mosdepth_label"] == "tetrasomic+"
    assert entries[3]["mosdepth_label"] == "monosomic"


def test_consistent_whole_chr_aneuploidy():
    chr_segments = pd.DataFrame(
        [
            {
                "sample": "S1",
                "logical_chr": "chr6A",
                "bin_start": 1,
                "bin_end": 10_000_000,
                "length_bp": 10_000_000,
                "n_bins": 10,
                "median_depth_rel": 1.55,
                "median_mq_low_frac": 0.2,
                "call_type": "aneuploidy_chr_gain",
                "confidence": 0.8,
                "evidence": "gain",
            },
        ],
    )
    out = reconcile_flagged_chr(
        sample="S1",
        logical_chr="chr6A",
        mosdepth_label="trisomic",
        mosdepth_rel_depth=1.59,
        chr_segments=chr_segments,
    )
    assert out["reconciliation"] == "consistent_whole_chr_aneuploidy"


def test_likely_cnv_not_aneuploidy_duplication():
    chr_segments = pd.DataFrame(
        [
            {
                "sample": "S2",
                "logical_chr": "chr6A",
                "bin_start": 1,
                "bin_end": 4_000_000,
                "length_bp": 4_000_000,
                "n_bins": 4,
                "median_depth_rel": 1.8,
                "median_mq_low_frac": 0.15,
                "call_type": "duplication",
                "confidence": 0.7,
                "evidence": "high_depth",
            },
        ],
    )
    out = reconcile_flagged_chr(
        sample="S2",
        logical_chr="chr6A",
        mosdepth_label="trisomic",
        mosdepth_rel_depth=1.59,
        chr_segments=chr_segments,
        config=ReconcileConfig(min_cnv_length_bp=2_000_000, min_cnv_bins=2),
    )
    assert out["reconciliation"] == "likely_cnv_not_aneuploidy"


def test_missing_mapq_sample():
    flagged = pd.DataFrame(
        [
            {
                "sample": "S3",
                "aneuploid_flag": True,
                "abnormal_chr_list": "chr6A:trisomic(1.59)",
            },
        ],
    )
    out = reconcile_flagged_samples(
        flagged_df=flagged,
        segments_df=pd.DataFrame(columns=list(SEGMENT_OUTPUT_COLUMNS)),
        samples_with_mapq=set(),
    )
    assert len(out) == 1
    assert out.iloc[0]["reconciliation"] == "missing_mapq"
