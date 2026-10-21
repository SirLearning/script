"""RefSeq v1.0 asset inventory reports present FASTA and extra annotation paths."""

from pathlib import Path

from genetics.genomics.ref.v1_assets import REF_V1_RELATIVE_ASSETS, inventory_ref_v1_assets


def test_inventory_fasta_and_extra_annotation(tmp_path: Path):
    fasta = tmp_path / REF_V1_RELATIVE_ASSETS["fasta_ABD"]
    fasta.parent.mkdir(parents=True)
    fasta.write_text(">chr1A\nACGT\n", encoding="utf-8")
    gff = tmp_path / "wheat_v1.1_Lulab.gff3"
    gff.write_text("##gff-version 3\n", encoding="utf-8")
    table = inventory_ref_v1_assets(
        tmp_path,
        extra_assets={"gff_lulab": gff, "homoeolog_triads": tmp_path / "missing_triads.tsv"},
    ).set_index("asset")
    assert bool(table.loc["fasta_ABD", "present"])
    assert bool(table.loc["fasta_ABD", "required"])
    assert bool(table.loc["gff_lulab", "present"])
    assert not bool(table.loc["gff_lulab", "required"])
    assert not bool(table.loc["homoeolog_triads", "present"])
    assert "gff_hc" not in table.index
