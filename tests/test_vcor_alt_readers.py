"""PyArrow / Polars can read PLINK2 tab .vcor including the #CHROM_A header."""

from pathlib import Path

import pyarrow.csv as pacsv
import polars as pl

VCOR_HEADER = "#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2\n"
VCOR_ROW = "1\t1000\t1-1000\t1\t26000\t1-26000\t0.90\n"


def _tiny_vcor(tmp_path: Path) -> Path:
    path = tmp_path / "tiny.vcor"
    path.write_text(VCOR_HEADER + VCOR_ROW, encoding="utf-8")
    return path


def test_pyarrow_reads_vcor_hash_header_and_r2(tmp_path: Path):
    table = pacsv.read_csv(
        _tiny_vcor(tmp_path),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
        read_options=pacsv.ReadOptions(use_threads=False),
    )
    names = table.column_names
    assert names[0] == "#CHROM_A"
    assert names[-1] == "UNPHASED_R2"
    assert float(table.column("UNPHASED_R2")[0].as_py()) == 0.90


def test_polars_reads_vcor_hash_header_and_r2(tmp_path: Path):
    df = pl.read_csv(_tiny_vcor(tmp_path), separator="\t", infer_schema_length=10)
    assert df.columns[0] == "#CHROM_A"
    assert df.columns[-1] == "UNPHASED_R2"
    assert float(df["UNPHASED_R2"][0]) == 0.90
