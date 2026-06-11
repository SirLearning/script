"""Unit tests for MAC computation from PLINK2 gcount columns."""

import pandas as pd
import pytest

from infra.utils.errors import DataLoadError
from genetics.genomics.variant.mac import (
    _compute_mac_miss_derived_metrics,
    _compute_mac_table,
    _mac_bin50_label,
    _mac_bin_maf_summary,
    _mac_bin_missing_summary,
    _mac_site_count_table,
    _plot_mac_maf_reg_by_mac_range,
    _plot_mac_missing_reg_by_mac_range,
    _stratified_sample_mac_bin50,
    _subset_by_mac_max,
    _summarize_mac_zero_sites,
)


def _sample_gcount_row(**overrides):
    base = {
        "HET_REF_ALT_CTS": 2,
        "TWO_ALT_GENO_CTS": 0,
        "HOM_REF_CT": 8,
        "HAP_REF_CT": 0,
        "HAP_ALT_CTS": 0,
    }
    base.update(overrides)
    return base


def test_compute_mac_table_het_site():
    df = pd.DataFrame([_sample_gcount_row()])
    out = _compute_mac_table(df)
    assert out.loc[0, "Alt_Count"] == 2
    assert out.loc[0, "Ref_Count"] == 18
    assert out.loc[0, "MAC"] == 2
    assert out.loc[0, "Het_Fraction"] == pytest.approx(1.0)


def test_compute_mac_table_mac1_singleton():
    df = pd.DataFrame([_sample_gcount_row(HET_REF_ALT_CTS=1, HOM_REF_CT=0, TWO_ALT_GENO_CTS=0)])
    out = _compute_mac_table(df)
    assert out.loc[0, "MAC"] == 1
    assert out.loc[0, "Alt_Count"] == 1
    assert out.loc[0, "Ref_Count"] == 1


def test_compute_mac_table_missing_columns_raises():
    df = pd.DataFrame([{"HOM_REF_CT": 1}])
    with pytest.raises(DataLoadError):
        _compute_mac_table(df)


def test_mac_site_count_table():
    df = pd.DataFrame(
        [
            _sample_gcount_row(),
            _sample_gcount_row(HET_REF_ALT_CTS=1, HOM_REF_CT=9),
        ]
    )
    mac_df = _compute_mac_table(df)
    counts = _mac_site_count_table(mac_df)
    assert set(counts["MAC"]) == {1, 2}
    assert counts.loc[counts["MAC"] == 2, "n_sites"].iloc[0] == 1
    assert counts.loc[counts["MAC"] == 1, "n_sites"].iloc[0] == 1


def test_compute_mac_miss_derived_metrics():
    df = pd.DataFrame({"MAC": [5], "F_MISS": [0.99], "N": [7675]})
    out = _compute_mac_miss_derived_metrics(df)
    assert out.loc[0, "AN"] == pytest.approx(2 * 7675 * 0.01)
    assert out.loc[0, "R_miss_boundary"] == pytest.approx(0.99 / (1 - 5 / 7675))


def test_mac_bin50_label():
    assert _mac_bin50_label(0) == "0-50"
    assert _mac_bin50_label(50) == "0-50"
    assert _mac_bin50_label(51) == "51-100"
    assert _mac_bin50_label(100) == "51-100"
    assert _mac_bin50_label(101) == "101-150"


def test_subset_by_mac_max():
    df = pd.DataFrame({"MAC": [0, 50, 51, 200], "F_MISS": [0.1, 0.2, 0.3, 0.4]})
    sub = _subset_by_mac_max(df, 100)
    assert len(sub) == 3
    assert sub["MAC"].max() == 51


def test_stratified_sample_mac_bin50_caps_per_bin():
    rows = []
    for mac in list(range(0, 120)) + [200]:
        rows.extend([{"MAC": mac, "F_MISS": 0.1} for _ in range(150)])
    df = pd.DataFrame(rows)
    out = _stratified_sample_mac_bin50(df, mac_max=100, per_bin=100, random_seed=1)
    assert len(out) == 200  # 0-50 and 51-100 bins => 100 each within mac_max=100
    assert out["MAC"].max() <= 100


def test_plot_mac_maf_reg_by_mac_range_filters(tmp_path):
    df = pd.DataFrame({"MAC": [0, 5, 50, 200, 800], "MAF": [0.0, 0.05, 0.1, 0.2, 0.4]})
    prefix = str(tmp_path / "out.variant.mac_maf")
    stats = _plot_mac_maf_reg_by_mac_range(df, prefix, mac_max=100, regression_max_points=100, random_seed=1)
    assert stats["N_MAC_0_100"] == 3
    assert (tmp_path / "out.variant.mac_maf.reg.mac0_100.png").exists()


def test_mac_bin_maf_summary():
    df = pd.DataFrame({"MAC": [0, 1, 1, 10], "MAF": [0.0, 0.01, 0.02, 0.1]})
    summary = _mac_bin_maf_summary(df)
    assert summary.loc[summary["mac_bucket"] == "1", "n_sites"].iloc[0] == 2


def test_plot_mac_missing_reg_by_mac_range_filters(tmp_path):
    df = pd.DataFrame({"MAC": [0, 5, 50, 200, 800], "F_MISS": [0.1, 0.2, 0.3, 0.4, 0.5]})
    prefix = str(tmp_path / "out.variant.mac_miss")
    stats = _plot_mac_missing_reg_by_mac_range(df, prefix, mac_max=100, regression_max_points=100, random_seed=1)
    assert stats["N_MAC_0_100"] == 3
    assert (tmp_path / "out.variant.mac_miss.reg.mac0_100.png").exists()


def test_mac_bin_missing_summary():
    df = pd.DataFrame(
        {
            "MAC": [0, 1, 1, 10, 100],
            "F_MISS": [0.2, 0.3, 0.5, 0.1, 0.05],
        }
    )
    summary = _mac_bin_missing_summary(df)
    assert list(summary["mac_bucket"]) == ["0", "1", "10-99", ">=100"]
    assert summary.loc[summary["mac_bucket"] == "1", "n_sites"].iloc[0] == 2
    assert summary.loc[summary["mac_bucket"] == "1", "mean_F_MISS"].iloc[0] == pytest.approx(0.4)


def test_summarize_mac_zero_sites_none():
    df = pd.DataFrame([_sample_gcount_row()])
    mac_df = _compute_mac_table(df)
    summary = _summarize_mac_zero_sites(mac_df)
    assert summary.loc[0, "category"] == "none"
    assert summary.loc[0, "n_sites"] == 0
