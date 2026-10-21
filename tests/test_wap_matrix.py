"""WAP phenotype matrix aligned to assoc.keep."""

from pathlib import Path

import pandas as pd
import pytest

from genetics.phenotype.wap_matrix import (
    FID_COLUMN,
    IID_COLUMN,
    align_phenotype_to_keep,
    load_keep_samples,
    run_wap_phenotype_matrix,
    select_traits,
    summarize_phenotype_matrix,
)


def _write_keep(path: Path, samples: list[str]) -> Path:
    path.write_text("".join(f"0\t{s}\n" for s in samples), encoding="utf-8")
    return path


def _pheno_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"IID": "P002", "HeadingDate": 120.0, "PlantHeight": 88.0},
            {"IID": "P001", "HeadingDate": 118.0, "PlantHeight": 91.0},
            {"IID": "P999", "HeadingDate": 130.0, "PlantHeight": 70.0},
        ]
    )


def test_keep_order_drives_matrix_order_and_membership(tmp_path: Path):
    keep = load_keep_samples(_write_keep(tmp_path / "assoc.keep", ["P001", "P002", "P003"]))
    matrix = align_phenotype_to_keep(_pheno_frame(), keep, ["HeadingDate"])
    assert list(matrix[IID_COLUMN]) == ["P001", "P002", "P003"]
    assert list(matrix[FID_COLUMN]) == ["0", "0", "0"]
    # P003 has no phenotype row; P999 is not in the keep list.
    assert pd.isna(matrix.loc[2, "HeadingDate"])
    assert "P999" not in set(matrix[IID_COLUMN])


def test_duplicate_samples_averaged_only_when_requested(tmp_path: Path):
    keep = load_keep_samples(_write_keep(tmp_path / "assoc.keep", ["P001"]))
    repeats = pd.DataFrame(
        [{"IID": "P001", "Yield": 4.0}, {"IID": "P001", "Yield": 6.0}]
    )
    first = align_phenotype_to_keep(repeats, keep, ["Yield"])
    assert float(first.loc[0, "Yield"]) == 4.0
    averaged = align_phenotype_to_keep(repeats, keep, ["Yield"], aggregate_duplicates=True)
    assert float(averaged.loc[0, "Yield"]) == 5.0


def test_select_traits_rejects_unknown_column():
    pheno = _pheno_frame()
    assert select_traits(pheno) == ["HeadingDate", "PlantHeight"]
    with pytest.raises(ValueError):
        select_traits(pheno, "NoSuchTrait")


def test_summary_counts_observed_values(tmp_path: Path):
    keep = load_keep_samples(_write_keep(tmp_path / "assoc.keep", ["P001", "P002", "P003"]))
    matrix = align_phenotype_to_keep(_pheno_frame(), keep, ["HeadingDate"])
    summary = summarize_phenotype_matrix(matrix, ["HeadingDate"]).set_index("trait")
    assert int(summary.loc["HeadingDate", "n_observed"]) == 2
    assert float(summary.loc["HeadingDate", "missing_rate"]) == pytest.approx(1 / 3)


def test_run_writes_pheno_and_qc(tmp_path: Path):
    keep_path = _write_keep(tmp_path / "assoc.keep", ["P001", "P002"])
    pheno_path = tmp_path / "wap.pheno.raw.tsv"
    _pheno_frame().to_csv(pheno_path, sep="\t", index=False)

    run_wap_phenotype_matrix(pheno_path, keep_path, tmp_path / "wap", traits="HeadingDate")
    assert (tmp_path / "wap.pheno.tsv").is_file()
    assert (tmp_path / "wap.pheno.qc.tsv").is_file()
    assert (tmp_path / "wap.pheno.th.tsv").is_file()
    written = pd.read_csv(tmp_path / "wap.pheno.tsv", sep="\t")
    assert list(written.columns) == [FID_COLUMN, IID_COLUMN, "HeadingDate"]
