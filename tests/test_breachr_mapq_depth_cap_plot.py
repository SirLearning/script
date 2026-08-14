"""Tests for breachr annotated MAPQ depth cap plots."""

from __future__ import annotations

import pandas as pd

from breachr.mapq_depth_cap_plot import build_chr_annotations


def test_build_chr_annotations_includes_reconcile_label():
    flagged = pd.DataFrame(
        [
            {
                "sample": "S1",
                "abnormal_chr_list": "chr6A:trisomic(1.59)",
                "baseline_depth": 2.0,
            },
        ],
    )
    reconcile = pd.DataFrame(
        [
            {
                "sample": "S1",
                "logical_chr": "chr6A",
                "reconciliation": "likely_cnv_not_aneuploidy",
                "breachr_call_type": "duplication",
                "breachr_length_bp": 4_000_000,
            },
        ],
    )
    segments = pd.DataFrame(
        [
            {
                "sample": "S1",
                "logical_chr": "chr6A",
                "call_type": "duplication",
                "length_bp": 4_000_000,
                "confidence": 0.8,
                "median_depth_rel": 1.8,
                "bin_start": 1,
                "bin_end": 4_000_000,
            },
        ],
    )
    ann = build_chr_annotations("S1", flagged, reconcile, segments)
    assert "chr6A" in ann
    assert "trisomic" in ann["chr6A"]["mosdepth_title"]
    assert ann["chr6A"]["reconciliation"] == "CNV not aneuploidy"
    assert "duplication" in ann["chr6A"]["breachr_call"]
