"""Tests for breachr reconciliation plotting."""

from __future__ import annotations

import pandas as pd

from breachr.reconcile_plot import plot_breachr_cnv_length_distribution, summarize_reconcile_by_sample


def test_summarize_reconcile_by_sample():
    reconcile = pd.DataFrame(
        [
            {"sample": "S1", "reconciliation": "likely_cnv_not_aneuploidy"},
            {"sample": "S1", "reconciliation": "consistent_whole_chr_aneuploidy"},
            {"sample": "S2", "reconciliation": "likely_cnv_not_aneuploidy"},
        ],
    )
    out = summarize_reconcile_by_sample(reconcile)
    assert out.loc[out["sample"] == "S1", "sample_class"].iloc[0] == "mixed_cnv_and_aneuploidy"
    assert out.loc[out["sample"] == "S2", "sample_class"].iloc[0] == "cnv_misclassification"


def test_plot_breachr_cnv_length_distribution(tmp_path):
    segments = pd.DataFrame(
        [
            {"sample": "S1", "call_type": "deletion", "length_bp": 5_000_000},
            {"sample": "S1", "call_type": "duplication", "length_bp": 12_000_000},
            {"sample": "S2", "call_type": "deletion", "length_bp": 2_000_000},
        ],
    )
    seg_path = tmp_path / "breachr_segments_all.tsv"
    segments.to_csv(seg_path, sep="\t", index=False)
    out = plot_breachr_cnv_length_distribution(seg_path, tmp_path / "plots")
    assert len(out) == 2
    assert (tmp_path / "plots" / "cnv_length_distribution.png").exists()
    assert (tmp_path / "plots" / "cnv_length_distribution_linear.png").exists()
