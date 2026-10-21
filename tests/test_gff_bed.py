"""Merged GFF BED used to extract exon sites from main_raw."""

from pathlib import Path

import pytest

from genetics.genomics.gene.gff_bed import (
    bin_gff_feature_counts,
    parse_gff_feature_types,
    write_merged_gff_bed,
)


def test_parse_gff_feature_types_defaults_to_exon():
    assert parse_gff_feature_types(None) == ("exon",)
    assert parse_gff_feature_types("  ") == ("exon",)
    assert parse_gff_feature_types("CDS, exon") == ("CDS", "exon")


def test_write_merged_gff_bed_merges_overlapping_exons(tmp_path: Path):
    gff = tmp_path / "tiny.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "1\t.\tgene\t10\t100\t.\t+\t.\tID=g1\n"
        "1\t.\texon\t10\t20\t.\t+\t.\tParent=t1\n"
        "1\t.\texon\t18\t30\t.\t+\t.\tParent=t2\n"
        "1\t.\tCDS\t12\t20\t.\t+\t0\tParent=t1\n"
        "2\t.\texon\t5\t8\t.\t-\t.\tParent=t3\n",
        encoding="utf-8",
    )
    bed = tmp_path / "exon.bed"
    n = write_merged_gff_bed(gff, bed)
    assert n == 2
    assert bed.read_text(encoding="utf-8") == "1\t9\t30\n2\t4\t8\n"


def test_bin_gff_feature_counts_uses_gene_midpoints_and_strips_chr_prefix(tmp_path: Path):
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1A\t.\tgene\t100\t200\t.\t+\t.\tID=g1\n"
        "chr1A\t.\tgene\t10000100\t10000200\t.\t+\t.\tID=g2\n"
        "chr1A\t.\texon\t150\t180\t.\t+\t.\tParent=g1\n"
        "chr1B\t.\tgene\t50\t80\t.\t-\t.\tID=g3\n",
        encoding="utf-8",
    )
    bins = bin_gff_feature_counts(gff, bin_bp=10_000_000, feature_type="gene")
    a0 = bins[(bins["chrom"] == "1A") & (bins["bin_start"] == 0)].iloc[0]
    a1 = bins[(bins["chrom"] == "1A") & (bins["bin_start"] == 10_000_000)].iloc[0]
    b0 = bins[(bins["chrom"] == "1B") & (bins["bin_start"] == 0)].iloc[0]
    assert int(a0["value"]) == 1
    assert int(a1["value"]) == 1
    assert int(b0["value"]) == 1


def test_gene_lengths_from_gff_uses_longest_transcript_cds(tmp_path: Path):
    from genetics.genomics.gene.gff_bed import gene_lengths_from_gff

    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "1\t.\tgene\t1\t100\t.\t+\t.\tID=TraesCS1A02G000100\n"
        "1\t.\tmRNA\t1\t100\t.\t+\t.\tID=TraesCS1A02G000100.1;Parent=TraesCS1A02G000100\n"
        "1\t.\tCDS\t10\t19\t.\t+\t0\tParent=TraesCS1A02G000100.1\n"
        "1\t.\tCDS\t30\t49\t.\t+\t0\tParent=TraesCS1A02G000100.1\n"
        "1\t.\tmRNA\t1\t40\t.\t+\t.\tID=TraesCS1A02G000100.2;Parent=TraesCS1A02G000100\n"
        "1\t.\tCDS\t10\t19\t.\t+\t0\tParent=TraesCS1A02G000100.2\n"
        "1\t.\tgene\t1\t20\t.\t+\t.\tID=TraesCS1B02G000200\n"
        "1\t.\tmRNA\t1\t20\t.\t+\t.\tID=TraesCS1B02G000200.1;Parent=TraesCS1B02G000200\n"
        "1\t.\tCDS\t1\t5\t.\t+\t0\tParent=TraesCS1B02G000200.1\n"
        "1\t.\texon\t1\t5\t.\t+\t.\tParent=TraesCS1B02G000200.1\n"
        "1\t.\texon\t10\t20\t.\t+\t.\tParent=TraesCS1B02G000200.1\n",
        encoding="utf-8",
    )
    genes = gene_lengths_from_gff(gff).set_index("gene_id")
    assert int(genes.loc["TraesCS1A02G000100", "gene_span_bp"]) == 100
    assert int(genes.loc["TraesCS1A02G000100", "cds_bp"]) == 30
    assert int(genes.loc["TraesCS1B02G000200", "gene_span_bp"]) == 20
    assert int(genes.loc["TraesCS1B02G000200", "cds_bp"]) == 5
    assert int(genes.loc["TraesCS1B02G000200", "n_exons"]) == 2
    assert genes.loc["TraesCS1A02G000100", "subgenome"] == "A"
    assert genes.loc["TraesCS1B02G000200", "subgenome"] == "B"


def test_write_merged_gff_bed_rejects_empty_feature_set(tmp_path: Path):
    gff = tmp_path / "empty.gff3"
    gff.write_text("1\t.\tgene\t1\t10\t.\t+\t.\tID=g1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="No exon"):
        write_merged_gff_bed(gff, tmp_path / "out.bed")
