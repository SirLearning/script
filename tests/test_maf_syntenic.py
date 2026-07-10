"""Tests for MAF syntenic site extraction."""

from __future__ import annotations

from pathlib import Path

from genetics.genomics.alignment.maf_syntenic import extract_syntenic_sites_from_maf


def _write_maf(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")


def test_extract_k2_filters_cross_chr_and_gaps(tmp_path: Path):
    maf = tmp_path / "1A.maf"
    _write_maf(
        maf,
        """##maf version=1 scoring=roast.v3.3
a score=1000.0
s traesA.chr1 100 3 + 1000 ACG
s trdicA.chr1 200 3 + 2000 ACG
s trdicB.chr6 300 3 + 3000 ACG

a score=500.0
s traesA.chr1 103 2 + 1000 GT
s aetau.chr1 400 1 + 4000 G
s aetau.chr1 401 1 + 4000 T

a score=-10.0
s traesA.chr1 200 1 + 1000 A
s trura.chr1 500 1 + 5000 A
""",
    )
    species_file = tmp_path / "species.txt"
    species_file.write_text("traesA\ntrdicA\naetau\n", encoding="utf-8")
    prefix = tmp_path / "1A.syntenic"
    summary = extract_syntenic_sites_from_maf(
        str(maf),
        str(prefix),
        "1A",
        species_list_path=str(species_file),
        min_species=2,
        min_score=0.0,
        filter_cross_chr=True,
    )

    assert summary["blocks_total"] == 3
    assert summary["blocks_processed"] == 2
    assert summary["blocks_skipped_score"] == 1
    assert summary["rows_dropped_cross_chr"] == 1
    assert summary["columns_emitted"] == 4
    assert summary["species_in_output"] == 3

    import gzip

    with gzip.open(tmp_path / "1A.syntenic.matrix.tsv.gz", "rt", encoding="utf-8") as handle:
        matrix = handle.read()
    lines = matrix.strip().splitlines()
    assert lines[0] == "ref_chr\tref_pos\tblock_score\tn_species\ttraesA\ttrdicA\taetau"
    assert lines[1] == "chr1\t100\t1000.000000\t2\tA\tA\t-"
    assert any(line.startswith("chr1\t103\t") for line in lines[1:])
    assert all("trdicB" not in line for line in lines[1:])


def test_chr_label_to_ref_mapping():
    from genetics.genomics.alignment.maf_syntenic import chr_label_to_ref

    assert chr_label_to_ref("1A") == ("traesA", "chr1")
    assert chr_label_to_ref("7D") == ("traesD", "chr7")
