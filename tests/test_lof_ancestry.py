"""Polarize LoF alleles against Triticeae lineage outgroups."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from genetics.genomics.variant.lof_ancestry import (
    CLASS_ANCESTRAL,
    CLASS_DERIVED,
    CLASS_UNKNOWN,
    classify_lof_ancestry,
    infer_ancestral_allele,
    lookup_ancestral_alleles,
    outgroups_for_subgenome,
    run_lof_ancestry,
    segment_to_matrix_key,
)


def test_classify_lof_ancestry_three_classes():
    assert classify_lof_ancestry("T", "C") == CLASS_DERIVED
    assert classify_lof_ancestry("T", "T") == CLASS_ANCESTRAL
    assert classify_lof_ancestry("t", "T") == CLASS_ANCESTRAL
    assert classify_lof_ancestry("T", "") == CLASS_UNKNOWN
    assert classify_lof_ancestry("T", "-") == CLASS_UNKNOWN


def test_infer_ancestral_allele_consensus_conflict_and_gaps():
    outgroups = outgroups_for_subgenome("A")
    assert "traesA" not in outgroups
    consensus = infer_ancestral_allele(
        {"traesA": "C", "trdurA": "t", "trboe": "T", "trmon": "-"},
        outgroups,
    )
    assert consensus == "T"
    conflict = infer_ancestral_allele(
        {"trdurA": "T", "trboe": "C"},
        outgroups,
    )
    assert conflict == ""
    missing = infer_ancestral_allele({"traesA": "C", "trdurA": "-"}, outgroups)
    assert missing == ""


def test_lookup_and_run_on_tiny_matrix(tmp_path: Path):
    matrix = tmp_path / "1A.syntenic.matrix.tsv.gz"
    header = "ref_chr\tref_pos\tblock_score\tn_species\ttraesA\ttrboe\ttrdurA\n"
    rows = [
        header,
        "chr1\t40104\t100\t3\tC\tT\tT\n",
        "chr1\t40113\t100\t2\tG\t-\t-\n",
        "chr1\t40116\t100\t3\tA\tA\tC\n",
    ]
    with gzip.open(matrix, "wt", encoding="utf-8") as handle:
        handle.writelines(rows)

    found = lookup_ancestral_alleles(
        matrix,
        {40104, 40113, 40116},
        outgroups_for_subgenome("A"),
    )
    assert found[40104] == "T"
    assert found[40113] == ""
    assert found[40116] == ""

    chr_label, ref_pos = segment_to_matrix_key(1, 40105)
    assert chr_label == "1A"
    assert ref_pos == 40104

    vep = tmp_path / "A.vep.tsv"
    vep.write_text(
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n"
        "1-40105\t1:40105\tT\tTraesCS1A02G000100\tTraesCS1A02G000100.1\tTranscript\tstop_gained\n"
        "1-40114\t1:40114\tA\tTraesCS1A02G000100\tTraesCS1A02G000100.1\tTranscript\tstop_gained\n"
        "1-40117\t1:40117\tG\tTraesCS1A02G000100\tTraesCS1A02G000100.1\tTranscript\tsynonymous_variant\n"
        "1-40105\t1:40105\tT\tTraesCS1A02G000100\tTraesCS1A02G000100.1\tTranscript\tmissense_variant\n",
        encoding="utf-8",
    )
    prefix = tmp_path / "lof_anc"
    sites = run_lof_ancestry([vep], tmp_path, prefix)
    by_id = {row["id"]: row for row in sites.to_dict("records")}
    assert by_id["1-40105"]["ancestry_class"] == CLASS_ANCESTRAL
    assert by_id["1-40114"]["ancestry_class"] == CLASS_UNKNOWN
    assert "1-40117" not in by_id
    assert (prefix.with_name(prefix.name + ".lof_ancestry.th.tsv")).exists()


def test_lof_count_scales_with_cds_length(tmp_path: Path):
    from genetics.genomics.gene.gff_bed import gene_lengths_from_gff
    from genetics.genomics.variant.lof_ancestry import run_lof_vs_gene_length

    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "1\t.\tgene\t1\t300\t.\t+\t.\tID=TraesCS1A02G000100\n"
        "1\t.\tmRNA\t1\t300\t.\t+\t.\tID=t1;Parent=TraesCS1A02G000100\n"
        "1\t.\tCDS\t1\t100\t.\t+\t0\tParent=t1\n"
        "1\t.\tgene\t1\t50\t.\t+\t.\tID=TraesCS1A02G000200\n"
        "1\t.\tmRNA\t1\t50\t.\t+\t.\tID=t2;Parent=TraesCS1A02G000200\n"
        "1\t.\tCDS\t1\t10\t.\t+\t0\tParent=t2\n",
        encoding="utf-8",
    )
    lof = tmp_path / "lof.tsv"
    lof.write_text(
        "id\tgene\n"
        "1-1\tTraesCS1A02G000100\n"
        "1-2\tTraesCS1A02G000100\n"
        "1-3\tTraesCS1A02G000100\n"
        "1-4\tTraesCS1A02G000200\n",
        encoding="utf-8",
    )
    prefix = tmp_path / "vs_len"
    table = run_lof_vs_gene_length(lof, gff, prefix)
    by_gene = table.set_index("gene_id")
    assert int(by_gene.loc["TraesCS1A02G000100", "n_lof"]) == 3
    assert int(by_gene.loc["TraesCS1A02G000200", "n_lof"]) == 1
    assert int(by_gene.loc["TraesCS1A02G000100", "cds_bp"]) == 100
    assert Path(f"{prefix}.lof_vs_cds.png").exists()
    assert Path(f"{prefix}.lof_vs_gene_length.th.tsv").exists()


def test_rare_wap_count_vs_cds_length(tmp_path: Path):
    from genetics.genomics.variant.lof_ancestry import (
        collect_flagged_sites,
        count_flagged_per_gene,
        run_lof_vs_gene_length,
    )

    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "1\t.\tgene\t1\t300\t.\t+\t.\tID=TraesCS1A02G000100\n"
        "1\t.\tmRNA\t1\t300\t.\t+\t.\tID=t1;Parent=TraesCS1A02G000100\n"
        "1\t.\tCDS\t1\t100\t.\t+\t0\tParent=t1\n"
        "1\t.\tgene\t1\t50\t.\t+\t.\tID=TraesCS1A02G000200\n"
        "1\t.\tmRNA\t1\t50\t.\t+\t.\tID=t2;Parent=TraesCS1A02G000200\n"
        "1\t.\tCDS\t1\t10\t.\t+\t0\tParent=t2\n",
        encoding="utf-8",
    )
    sets = tmp_path / "A.variant_sets.tsv"
    sets.write_text(
        "ID\tGene\trare_wap\n"
        "1-1\tTraesCS1A02G000100\tTrue\n"
        "1-2\tTraesCS1A02G000100\tTrue\n"
        "1-3\tTraesCS1A02G000100\tFalse\n"
        "1-4\tTraesCS1A02G000200\tTrue\n",
        encoding="utf-8",
    )
    sites = collect_flagged_sites([sets], flag_col="rare_wap")
    assert len(sites) == 3
    counted = count_flagged_per_gene([sets], flag_col="rare_wap", count_col="n_rare")
    assert int(counted.set_index("gene").loc["TraesCS1A02G000100", "n_rare"]) == 2
    assert int(counted.set_index("gene").loc["TraesCS1A02G000200", "n_rare"]) == 1
    prefix = tmp_path / "rare"
    table = run_lof_vs_gene_length(
        sites,
        gff,
        prefix,
        count_col="n_rare",
        y_noun="Rare (WAP MAF<1%)",
        file_stem="rare",
    )
    by_gene = table.set_index("gene_id")
    assert int(by_gene.loc["TraesCS1A02G000100", "n_rare"]) == 2
    assert int(by_gene.loc["TraesCS1A02G000200", "n_rare"]) == 1
    assert Path(f"{prefix}.rare_vs_cds.png").exists()
    assert Path(f"{prefix}.rare_per_kb_vs_cds.png").exists()
    pre = tmp_path / "rare_pre"
    pre_table = run_lof_vs_gene_length(
        pd.DataFrame({"gene": ["TraesCS1A02G000100", "TraesCS1A02G000200"], "n_rare": [2, 1]}),
        gff,
        pre,
        count_col="n_rare",
        y_noun="Rare (WAP MAF<1%)",
        file_stem="rare",
    )
    assert int(pre_table.set_index("gene_id").loc["TraesCS1A02G000100", "n_rare"]) == 2
