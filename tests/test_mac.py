"""Unit tests for MAC computation from PLINK2 gcount columns."""

import pandas as pd
import pytest

from infra.utils.errors import DataLoadError
from genetics.genomics.variant.mac import _compute_mac_table, _mac_site_count_table, _summarize_mac_zero_sites


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


def test_summarize_mac_zero_sites_none():
    df = pd.DataFrame([_sample_gcount_row()])
    mac_df = _compute_mac_table(df)
    summary = _summarize_mac_zero_sites(mac_df)
    assert summary.loc[0, "category"] == "none"
    assert summary.loc[0, "n_sites"] == 0
