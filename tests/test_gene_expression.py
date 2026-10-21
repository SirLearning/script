"""Gene expression atlases normalized onto RefSeq v1.1 gene IDs."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from genetics.genomics.gene.expression import (
    GENE_ID_COLUMN,
    LEGACY_ID_COLUMN,
    load_gene_id_map,
    load_gene_loci,
    load_sample_groups,
    load_tpm_matrix,
    normalize_gene_ids,
    run_gene_expression,
    summarize_gene_expression,
    tissue_specificity_tau,
)

GFF_LINES = [
    "##gff-version 3",
    "1\tIWGSC\tgene\t40098\t70338\t33\t-\t.\tID=TraesCS1A02G000100;previous_id=TraesCS1A01G000100;primconf=HC",
    "1\tIWGSC\tmRNA\t40098\t70338\t.\t-\t.\tID=TraesCS1A02G000100.1;Parent=TraesCS1A02G000100",
    "1\tIWGSC\tgene\t70239\t72000\t20\t+\t.\tID=TraesCS1A02G000200;previous_id=TraesCS1A01G000200;primconf=HC",
    # New in v1.1: no previous_id.
    "1\tIWGSC\tgene\t80000\t81000\t20\t+\t.\tID=TraesCS1A02G000300;primconf=HC",
]


def _write_gff(tmp_path: Path) -> Path:
    path = tmp_path / "genes.gff3"
    path.write_text("\n".join(GFF_LINES) + "\n", encoding="utf-8")
    return path


def test_gene_id_map_reads_previous_id(tmp_path: Path):
    id_map = load_gene_id_map(_write_gff(tmp_path))
    assert list(id_map[GENE_ID_COLUMN]) == [
        "TraesCS1A02G000100",
        "TraesCS1A02G000200",
        "TraesCS1A02G000300",
    ]
    assert id_map.loc[0, LEGACY_ID_COLUMN] == "TraesCS1A01G000100"
    assert id_map.loc[2, LEGACY_ID_COLUMN] == ""


def test_gene_loci_keeps_named_chrom_and_midpoint(tmp_path: Path):
    gff = tmp_path / "hc.gff3"
    gff.write_text(
        "chr1A\tIWGSC\tgene\t100\t200\t.\t+\t.\tID=TraesCS1A02G000100\n",
        encoding="utf-8",
    )
    loci = load_gene_loci(gff)
    assert list(loci["chrom"]) == ["chr1A"]
    assert int(loci.loc[0, "pos"]) == 150


def test_gene_loci_linearizes_lulab_segment_ids(tmp_path: Path):
    """Segment 14 is the second arm of chr3A; the midpoint must be offset."""
    gff = tmp_path / "lulab.gff3"
    gff.write_text(
        "14\tIWGSC\tgene\t100\t200\t.\t+\t.\tID=TraesCS3A02G999900\n",
        encoding="utf-8",
    )
    loci = load_gene_loci(gff)
    assert list(loci["chrom"]) == ["chr3A"]
    assert int(loci.loc[0, "pos"]) > 200


def test_tpm_matrix_drops_the_length_column(tmp_path: Path):
    path = tmp_path / "tpm.tsv"
    pd.DataFrame(
        {
            "Geneid": ["TraesCS1A02G000100", "TraesCS1A02G000200"],
            "Length": [5268, 927],
            "S1": [1.0, 0.0],
            "S2": [3.0, 2.0],
        }
    ).to_csv(path, sep="\t", index=False)
    matrix = load_tpm_matrix(path)
    assert list(matrix.columns) == [GENE_ID_COLUMN, "S1", "S2"]


def test_v10_ids_map_through_previous_id_and_lc_is_not_stripped(tmp_path: Path):
    """HC and LC share numbering in v1.0, so stripping the suffix would merge genes."""

    id_map = load_gene_id_map(_write_gff(tmp_path))
    matrix = pd.DataFrame(
        {
            GENE_ID_COLUMN: [
                "TraesCS1A01G000100",
                "TraesCS1A01G000100LC",
                "TraesCS1A01G000200",
            ],
            "S1": [1.0, 9.0, 2.0],
        }
    )
    mapped, report = normalize_gene_ids(matrix, id_map)
    assert list(mapped[GENE_ID_COLUMN]) == ["TraesCS1A02G000100", "TraesCS1A02G000200"]
    # The LC row is dropped rather than collapsed onto the HC gene, so its TPM never leaks.
    assert list(mapped["S1"]) == [1.0, 2.0]
    assert report["n_mapped_from_v10"] == 2
    assert report["n_low_confidence_dropped"] == 1


def test_v11_ids_pass_through_unchanged(tmp_path: Path):
    id_map = load_gene_id_map(_write_gff(tmp_path))
    matrix = pd.DataFrame({GENE_ID_COLUMN: ["TraesCS1A02G000300"], "S1": [4.0]})
    mapped, report = normalize_gene_ids(matrix, id_map)
    assert list(mapped[GENE_ID_COLUMN]) == ["TraesCS1A02G000300"]
    assert report["n_direct_v11"] == 1
    assert report["n_mapped_from_v10"] == 0


def test_tau_is_zero_for_uniform_and_one_for_single_group():
    uniform = np.array([[5.0, 5.0, 5.0]], dtype="float32")
    specific = np.array([[10.0, 0.0, 0.0]], dtype="float32")
    assert tissue_specificity_tau(uniform)[0] == pytest.approx(0.0)
    assert tissue_specificity_tau(specific)[0] == pytest.approx(1.0)


def test_summary_uses_group_means_when_metadata_is_given():
    matrix = pd.DataFrame(
        {
            GENE_ID_COLUMN: ["g1", "g2"],
            "S1": [10.0, 0.0],
            "S2": [8.0, 0.0],
            "S3": [0.0, 6.0],
        }
    )
    groups = pd.Series(["leaf", "leaf", "root"], index=["S1", "S2", "S3"])
    summary = summarize_gene_expression(matrix, sample_groups=groups).set_index(GENE_ID_COLUMN)
    assert summary.loc["g1", "tpm_leaf"] == pytest.approx(9.0)
    assert summary.loc["g1", "tpm_root"] == pytest.approx(0.0)
    assert summary.loc["g1", "breadth"] == 1
    assert summary.loc["g1", "tau"] == pytest.approx(1.0)
    assert bool(summary.loc["g2", "expressed"])


def test_load_sample_groups_rejects_a_missing_column(tmp_path: Path):
    path = tmp_path / "meta.txt"
    pd.DataFrame({"run_accession": ["S1"], "Tissue": ["leaf"]}).to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError):
        load_sample_groups(path, sample_col="run_accession", group_col="High level tissue")


def test_run_gene_expression_writes_summary(tmp_path: Path):
    gff = _write_gff(tmp_path)
    tpm = tmp_path / "Development_tpm.tsv"
    pd.DataFrame(
        {
            "gene": ["TraesCS1A01G000100", "TraesCS1A01G000200", "TraesCS1A01G000900LC"],
            "Sample_1": [10.0, 0.0, 5.0],
            "Sample_2": [8.0, 0.2, 5.0],
        }
    ).to_csv(tpm, sep="\t", index=False)
    metadata = tmp_path / "default_metadata.txt"
    pd.DataFrame(
        {
            "run_accession": ["Sample_1", "Sample_2"],
            "High level tissue": ["leaf", "root"],
        }
    ).to_csv(metadata, sep="\t", index=False)

    summary = run_gene_expression(tpm, gff, tmp_path / "dev", metadata_path=metadata)
    assert (tmp_path / "dev.gene_expr.tsv").is_file()
    assert (tmp_path / "dev.gene_expr.th.tsv").is_file()
    assert set(summary[GENE_ID_COLUMN]) == {"TraesCS1A02G000100", "TraesCS1A02G000200"}
    thresholds = pd.read_csv(tmp_path / "dev.gene_expr.th.tsv", sep="\t")
    assert int(thresholds.loc[0, "n_low_confidence_dropped"]) == 1
    assert int(thresholds.loc[0, "n_expressed_genes"]) == 1
