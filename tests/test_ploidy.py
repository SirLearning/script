"""Corrected genome-formula ploidy from sample_group plus breachr flags."""

from pathlib import Path

import pandas as pd
import pytest

from genetics.germplasm.sample.ploidy import (
    assign_corrected_ploidy,
    formula_from_subgenomes,
    load_breachr_segments,
    nominal_ploidy_from_group,
    subgenome_flags_from_breachr,
    subgenome_ploidies,
)


def breachr_rows(sample: str, flags: dict[str, str]) -> list[dict]:
    """One whole-chromosome row per flagged subgenome (7 homoeologues)."""
    rows: list[dict] = []
    for homoe in range(1, 8):
        for sg, call_type in flags.items():
            rows.append(
                {
                    "sample": sample,
                    "logical_chr": f"chr{homoe}{sg}",
                    "bin_start": 1,
                    "bin_end": 1_000_000,
                    "length_bp": 1_000_000,
                    "n_bins": 1,
                    "median_depth_rel": 0.1,
                    "median_mq_low_frac": 0.4,
                    "call_type": call_type,
                    "confidence": 0.85,
                    "evidence": "test",
                }
            )
    return rows


def _seg_table(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def test_nominal_ploidy_follows_sample_group_records():
    assert nominal_ploidy_from_group("A") == "AA"
    assert nominal_ploidy_from_group("S") == "SS"
    assert nominal_ploidy_from_group("D") == "DD"
    assert nominal_ploidy_from_group("AB") == "AABB"
    for group in ("ABD", "WAP", "Watkins", "Nature", "HZNU"):
        assert nominal_ploidy_from_group(group) == "AABBDD"


def test_formula_from_present_subgenomes():
    assert formula_from_subgenomes(frozenset("A")) == "AA"
    assert formula_from_subgenomes(frozenset("B")) == "SS"
    assert formula_from_subgenomes(frozenset("D")) == "DD"
    assert formula_from_subgenomes(frozenset("AB")) == "AABB"
    assert formula_from_subgenomes(frozenset("ABD")) == "AABBDD"
    assert formula_from_subgenomes(frozenset("AD")) is None


def test_expected_absent_d_in_aabb_is_not_a_correction():
    segs = _seg_table(breachr_rows("AB_0001", {"D": "subgenome_absent"}))
    groups = pd.DataFrame({"Sample": ["AB_0001"], "Group": ["AB"]})
    out = assign_corrected_ploidy(groups, segs).set_index("Sample")
    assert out.loc["AB_0001", "ploidy"] == "AABB"
    assert out.loc["AB_0001", "nominal_ploidy"] == "AABB"
    assert out.loc["AB_0001", "reason"] == "agrees"


def test_hexaploid_d_reference_mismatch_is_reclassified_aabb():
    segs = _seg_table(breachr_rows("P1499A", {"D": "subgenome_reference_mismatch"}))
    groups = pd.DataFrame({"Sample": ["P1499A"], "Group": ["WAP"]})
    out = assign_corrected_ploidy(groups, segs).set_index("Sample")
    assert out.loc["P1499A", "nominal_ploidy"] == "AABBDD"
    assert out.loc["P1499A", "ploidy"] == "AABB"
    assert out.loc["P1499A", "reason"] == "breachr_ploidy_mismatch"
    assert out.loc["P1499A", "present_subgenomes"] == "AB"
    assert out.loc["P1499A", "breachr_drop"] == "D"


def test_ab_b_reference_mismatch_is_reclassified_aa():
    segs = _seg_table(
        breachr_rows("AB_0043", {"B": "subgenome_reference_mismatch", "D": "subgenome_absent"})
    )
    groups = pd.DataFrame({"Sample": ["AB_0043"], "Group": ["AB"]})
    out = assign_corrected_ploidy(groups, segs).set_index("Sample")
    assert out.loc["AB_0043", "ploidy"] == "AA"
    assert out.loc["AB_0043", "reason"] == "breachr_ploidy_mismatch"
    assert out.loc["AB_0043", "present_subgenomes"] == "A"


def test_a_diploid_expected_bd_absent_stays_aa():
    segs = _seg_table(
        breachr_rows("A_0001", {"B": "subgenome_absent", "D": "subgenome_absent"})
    )
    groups = pd.DataFrame({"Sample": ["A_0001"], "Group": ["A"]})
    out = assign_corrected_ploidy(groups, segs).set_index("Sample")
    assert out.loc["A_0001", "ploidy"] == "AA"
    assert out.loc["A_0001", "reason"] == "agrees"


def test_s_expected_ad_absent_stays_ss():
    segs = _seg_table(
        breachr_rows("S_0001", {"A": "subgenome_absent", "D": "subgenome_absent"})
    )
    groups = pd.DataFrame({"Sample": ["S_0001"], "Group": ["S"]})
    out = assign_corrected_ploidy(groups, segs).set_index("Sample")
    assert out.loc["S_0001", "ploidy"] == "SS"
    assert out.loc["S_0001", "reason"] == "agrees"


def test_focal_deletion_rows_do_not_change_ploidy():
    segs = pd.DataFrame(
        [
            {
                "sample": "ABD_0001",
                "logical_chr": "chr6B",
                "call_type": "deletion",
            }
        ]
    )
    groups = pd.DataFrame({"Sample": ["ABD_0001"], "Group": ["ABD"]})
    out = assign_corrected_ploidy(groups, segs).set_index("Sample")
    assert out.loc["ABD_0001", "ploidy"] == "AABBDD"
    assert out.loc["ABD_0001", "reason"] == "agrees"


def test_all_subgenomes_absent_keeps_nominal():
    segs = _seg_table(
        breachr_rows(
            "SAMEA1",
            {
                "A": "subgenome_absent",
                "B": "subgenome_absent",
                "D": "subgenome_absent",
            },
        )
    )
    groups = pd.DataFrame({"Sample": ["SAMEA1"], "Group": ["Nature"]})
    out = assign_corrected_ploidy(groups, segs).set_index("Sample")
    assert out.loc["SAMEA1", "ploidy"] == "AABBDD"
    assert out.loc["SAMEA1", "reason"] == "breachr_unmapped"


def test_join_key_strips_deduped_suffix():
    segs = _seg_table(breachr_rows("JZ174_deduped", {"D": "subgenome_reference_mismatch"}))
    groups = pd.DataFrame({"Sample": ["JZ174"], "Group": ["HZNU"]})
    out = assign_corrected_ploidy(groups, segs).set_index("Sample")
    assert out.loc["JZ174", "ploidy"] == "AABB"
    assert out.loc["JZ174", "sample_breachr"] == "JZ174_deduped"


def test_subgenome_ploidies_match_wheat_composition():
    assert subgenome_ploidies("A") == ("AA", "AABB", "AABBDD")
    assert subgenome_ploidies("B") == ("SS", "AABB", "AABBDD")
    assert subgenome_ploidies("D") == ("DD", "AABBDD")


def test_load_breachr_segments_and_flags(tmp_path: Path):
    path = tmp_path / "breachr_segments_all.tsv"
    _seg_table(breachr_rows("P1499A", {"D": "subgenome_reference_mismatch"})).to_csv(
        path, sep="\t", index=False
    )
    loaded = load_breachr_segments(path)
    flags = subgenome_flags_from_breachr(loaded).set_index("Sample")
    assert flags.loc["P1499A", "D"] == "subgenome_reference_mismatch"
    assert pd.isna(flags.loc["P1499A", "A"])
