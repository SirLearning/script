"""Analysis variant sets: coding mask x cohort frequency."""

from pathlib import Path

import pandas as pd
import pytest

from genetics.genomics.variant.coding_mask import (
    MASK_MIS_CONSENSUS,
    MASK_PLOF_HIGH,
    MASK_PROTEIN_ALTERING,
)
from genetics.genomics.variant.variant_sets import (
    SET_COLUMN_PREFIX,
    applicable_set_specs,
    build_variant_sets,
    join_gene_expression,
    join_mask_and_af,
    join_site_anno,
    normalize_variant_id,
    run_variant_sets,
    select_set_specs,
)


def _mask_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"Uploaded_variation": "1-10", "coding_mask": MASK_PLOF_HIGH, "Gene": "TraesCS1A02G000100"},
            {"Uploaded_variation": "1-20", "coding_mask": MASK_PLOF_HIGH, "Gene": "TraesCS1B02G000100"},
            {"Uploaded_variation": "1-30", "coding_mask": MASK_MIS_CONSENSUS, "Gene": "TraesCS1D02G000100"},
            {"Uploaded_variation": "1-40", "coding_mask": MASK_PROTEIN_ALTERING, "Gene": "TraesCS2A02G000100"},
        ]
    )


def _af_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"ID": "1-10", "MAF_wgs": 0.002, "MAF_wap": 0.003, "rare_wgs": True, "rare_wap": True, "rare_both": True},
            {"ID": "1-20", "MAF_wgs": 0.30, "MAF_wap": 0.28, "rare_wgs": False, "rare_wap": False, "rare_both": False},
            {"ID": "1-30", "MAF_wgs": 0.004, "MAF_wap": 0.20, "rare_wgs": True, "rare_wap": False, "rare_both": False},
            {"ID": "1-99", "MAF_wgs": 0.001, "MAF_wap": 0.001, "rare_wgs": True, "rare_wap": True, "rare_both": True},
        ]
    )


def test_normalize_variant_id_accepts_vep_column():
    out = normalize_variant_id(_mask_frame(), source="coding_mask")
    assert "ID" in out.columns
    assert list(out["ID"]) == ["1-10", "1-20", "1-30", "1-40"]


def test_join_keeps_only_shared_variants():
    merged = join_mask_and_af(_mask_frame(), _af_frame())
    assert set(merged["ID"]) == {"1-10", "1-20", "1-30"}
    assert "Gene" in merged.columns


def test_sets_intersect_mask_and_frequency():
    merged = join_mask_and_af(_mask_frame(), _af_frame())
    sets_df = build_variant_sets(merged).set_index("ID")
    assert bool(sets_df.loc["1-10", f"{SET_COLUMN_PREFIX}plof_rare_wgs"])
    assert not bool(sets_df.loc["1-20", f"{SET_COLUMN_PREFIX}plof_rare_wgs"])
    # 1-20 is a pLoF that is common in WAP (MAF 0.28 >= 0.05).
    assert bool(sets_df.loc["1-20", f"{SET_COLUMN_PREFIX}plof_common_wap"])
    # 1-30 is rare in WGS-core but common in WAP, so it fails the rare_both set.
    assert bool(sets_df.loc["1-30", f"{SET_COLUMN_PREFIX}mis_consensus_rare_wgs"])
    assert not bool(sets_df.loc["1-30", f"{SET_COLUMN_PREFIX}protein_altering_rare_both"])


def test_select_set_specs_rejects_unknown_name():
    assert [spec.name for spec in select_set_specs("plof_rare_wgs")] == ["plof_rare_wgs"]
    with pytest.raises(ValueError):
        select_set_specs("no_such_set")


def _site_anno_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "ID": "1-10",
                "degeneracy": "0fold",
                "element": "CDS_HC",
                "coding": True,
                "conserved": True,
                "regulatory": False,
                "neutral_4fold": False,
            },
            {
                "ID": "1-20",
                "degeneracy": "4fold",
                "element": "CDS_HC",
                "coding": True,
                "conserved": False,
                "regulatory": False,
                "neutral_4fold": True,
            },
            {
                "ID": "1-30",
                "degeneracy": "noncoding",
                "element": "intergenic",
                "coding": False,
                "conserved": True,
                "regulatory": True,
                "neutral_4fold": False,
            },
        ]
    )


def test_site_sets_are_skipped_without_a_site_anno_table():
    merged = join_mask_and_af(_mask_frame(), _af_frame())
    names = [spec.name for spec in applicable_set_specs(select_set_specs(None), merged)]
    assert "plof_rare_wgs" in names
    assert "neutral_4fold" not in names
    assert "plof_expressed_rare_wap" not in names


def test_site_sets_activate_once_site_classes_are_joined():
    merged = join_site_anno(join_mask_and_af(_mask_frame(), _af_frame()), _site_anno_frame())
    specs = applicable_set_specs(select_set_specs(None), merged)
    assert "neutral_4fold" in [spec.name for spec in specs]

    sets_df = build_variant_sets(merged, specs).set_index("ID")
    # 1-20 is 4-fold and unconstrained, and neutral_4fold ignores the coding mask.
    assert bool(sets_df.loc["1-20", f"{SET_COLUMN_PREFIX}neutral_4fold"])
    assert not bool(sets_df.loc["1-10", f"{SET_COLUMN_PREFIX}neutral_4fold"])
    # 1-10 is a conserved pLoF, but common in neither cohort check: it is rare in both.
    assert bool(sets_df.loc["1-10", f"{SET_COLUMN_PREFIX}plof_conserved_rare_wap"])
    assert not bool(sets_df.loc["1-20", f"{SET_COLUMN_PREFIX}plof_conserved_rare_wap"])


def test_variants_without_site_annotation_fall_out_of_site_sets():
    """A left join leaves NA for sites the featHighRes BED does not cover."""

    merged = join_site_anno(
        join_mask_and_af(_mask_frame(), _af_frame()),
        _site_anno_frame().iloc[:1],
    )
    specs = applicable_set_specs(select_set_specs("neutral_4fold"), merged)
    sets_df = build_variant_sets(merged, specs).set_index("ID")
    assert not sets_df[f"{SET_COLUMN_PREFIX}neutral_4fold"].any()


def test_gene_expression_joins_on_the_vep_gene_column():
    merged = join_mask_and_af(_mask_frame(), _af_frame())
    expression = pd.DataFrame(
        [
            {"gene_id": "TraesCS1A02G000100", "expressed": True, "max_tpm": 12.0},
            {"gene_id": "TraesCS1B02G000100", "expressed": False, "max_tpm": 0.0},
        ]
    )
    joined = join_gene_expression(merged, expression)
    specs = applicable_set_specs(select_set_specs("plof_expressed_rare_wap"), joined)
    sets_df = build_variant_sets(joined, specs).set_index("ID")
    assert bool(sets_df.loc["1-10", f"{SET_COLUMN_PREFIX}plof_expressed_rare_wap"])
    # 1-30's gene is absent from the atlas, so it is treated as not expressed.
    assert not bool(sets_df.loc["1-30", f"{SET_COLUMN_PREFIX}plof_expressed_rare_wap"])


def test_run_variant_sets_writes_table_and_extracts(tmp_path: Path):
    mask_path = tmp_path / "chr032.coding_mask.tsv"
    af_path = tmp_path / "chr032.dual_af.tsv"
    _mask_frame().to_csv(mask_path, sep="\t", index=False)
    _af_frame().to_csv(af_path, sep="\t", index=False)

    out = run_variant_sets(mask_path, af_path, tmp_path / "chr032", set_names="plof_rare_wgs")
    assert (tmp_path / "chr032.variant_sets.tsv").is_file()
    assert (tmp_path / "chr032.variant_sets.th.tsv").is_file()
    extract = tmp_path / "chr032.plof_rare_wgs.extract"
    assert extract.read_text(encoding="utf-8").split() == ["1-10"]
    assert int(out[f"{SET_COLUMN_PREFIX}plof_rare_wgs"].sum()) == 1
