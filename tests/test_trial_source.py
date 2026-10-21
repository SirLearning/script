"""Field-trial plot tables into WAP accession and check-cultivar sources."""

from pathlib import Path

import pandas as pd
import pytest

from genetics.phenotype.trial_source import (
    build_genotype_id_map,
    check_plot_slots,
    load_trial_plot_csv,
    normalize_genotype_id,
    parse_trial_sample_name,
    pivot_plot_means,
    run_trial_source_extract,
    split_check_rows,
)

CSV_HEADER = "sampleName,trialYear,trait,PH_1,PH_2,PH_3"
CSV_ROWS = [
    # P0001 and P0030 are real accessions; P0030 is genotyped as P0030A.
    "Z18-P0001,2018,PH_M_cm,50,52,",
    "Z19-P0001,2019,PH_M_cm,60,62,64",
    "Z18-P0030,2018,PH_M_cm,70,72,",
    # P0025 is a Jimai22 check slot: labelled in 2019, unlabelled in 2018.
    "Z18-P0025,2018,PH_M_cm,80,82,",
    "Z19-P0025-Jimai22,2019,PH_M_cm,76.5,79,80.5",
    # P0999 has phenotype but no genotype.
    "Z18-P0999,2018,PH_M_cm,90,,",
]


def _write_inputs(tmp_path: Path) -> tuple[Path, Path]:
    csv_path = tmp_path / "ph.csv"
    csv_path.write_text("\n".join([CSV_HEADER, *CSV_ROWS]) + "\n", encoding="utf-8")
    keep_path = tmp_path / "assoc.keep"
    keep_path.write_text("0\tP0001\n0\tP0030A\n", encoding="utf-8")
    return csv_path, keep_path


def test_parse_sample_name_splits_year_plot_and_check():
    assert parse_trial_sample_name("Z18-P0001") == ("18", "P0001", "")
    assert parse_trial_sample_name("Z19-P0025-Jimai22") == ("19", "P0025", "Jimai22")
    with pytest.raises(ValueError):
        parse_trial_sample_name("not-a-plot")


def test_genotype_id_normalisation_and_map():
    assert normalize_genotype_id("P0030A") == "P0030"
    assert normalize_genotype_id("P0001") == "P0001"
    assert build_genotype_id_map(["P0001", "P0030A"]) == {"P0001": "P0001", "P0030": "P0030A"}


def test_ambiguous_genotype_ids_rejected():
    with pytest.raises(ValueError):
        build_genotype_id_map(["P0030", "P0030A"])


def test_melt_drops_blank_replicates(tmp_path: Path):
    csv_path, _ = _write_inputs(tmp_path)
    tidy = load_trial_plot_csv(csv_path)
    # 2 + 3 + 2 + 2 + 3 + 1 measured plots, blanks dropped.
    assert len(tidy) == 13
    assert set(tidy["trial_year"].astype(int)) == {2018, 2019}


def test_unlabelled_rows_in_a_check_slot_are_checks(tmp_path: Path):
    csv_path, _ = _write_inputs(tmp_path)
    tidy = load_trial_plot_csv(csv_path)
    assert check_plot_slots(tidy) == {"P0025"}
    accessions, checks = split_check_rows(tidy)
    # Both the labelled 2019 row and the unlabelled 2018 row land in checks.
    assert set(checks["trial_year"].astype(int)) == {2018, 2019}
    assert set(checks["check_name"]) == {"Jimai22"}
    assert "P0025" not in set(accessions["plot_id"])


def test_pivot_gives_per_year_and_overall_means():
    matched = pd.DataFrame(
        [
            {"IID": "P0001", "trait": "PH", "trial_year": 2018, "value": 50.0},
            {"IID": "P0001", "trait": "PH", "trial_year": 2018, "value": 52.0},
            {"IID": "P0001", "trait": "PH", "trial_year": 2019, "value": 62.0},
        ]
    )
    wide = pivot_plot_means(matched).set_index("IID")
    assert float(wide.loc["P0001", "PH_2018"]) == 51.0
    assert float(wide.loc["P0001", "PH_2019"]) == 62.0
    assert float(wide.loc["P0001", "PH_mean"]) == pytest.approx((50 + 52 + 62) / 3)


def test_run_writes_accession_check_and_wide_sources(tmp_path: Path):
    csv_path, keep_path = _write_inputs(tmp_path)
    stats = run_trial_source_extract(
        csv_path,
        keep_path,
        tmp_path / "wap_ph",
        check_prefix=tmp_path / "jimai22_ph",
    )
    assert stats["n_check_slots"] == 1
    assert stats["n_genotyped_accessions"] == 2
    # P0999 has values but no genotype.
    assert stats["n_ungenotyped_plots"] == 1

    plots = pd.read_csv(tmp_path / "wap_ph.plot.tsv", sep="\t")
    assert set(plots["IID"]) == {"P0001", "P0030A"}
    checks = pd.read_csv(tmp_path / "jimai22_ph.plot.tsv", sep="\t")
    assert set(checks["check_name"]) == {"Jimai22"}

    wide = pd.read_csv(tmp_path / "wap_ph.pheno.tsv", sep="\t").set_index("IID")
    assert float(wide.loc["P0030A", "PH_M_cm_2018"]) == 71.0
    assert pd.isna(wide.loc["P0030A", "PH_M_cm_2019"])
