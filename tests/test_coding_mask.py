"""Project coding masks from VEP-like consequence rows."""

from pathlib import Path

import pandas as pd

from genetics.genomics.variant.coding_mask import (
    MASK_MIS_CONSENSUS,
    MASK_MIS_SEVERE,
    MASK_PLOF_HIGH,
    MASK_PROTEIN_ALTERING,
    classify_vep_table,
    load_vep_tab_table,
    most_severe_consequence,
    run_coding_mask,
    subgenome_from_gene_id,
)

VEP_TAB_FILE = "\n".join(
    [
        "## ENSEMBL VARIANT EFFECT PREDICTOR v116.1",
        "## Output produced at 2026-09-03 00:00:00",
        "## Column descriptions:",
        "#Uploaded_variation\tLocation\tAllele\tGene\tConsequence",
        "1-275\t1:275\tT\tTraesCS1A02G000100\tmissense_variant",
        "0-41\tchrUn:41\tT\tTraesCSU02G000100\tdownstream_gene_variant",
    ]
)


def test_vep_tab_reader_skips_metadata_and_strips_header_hash(tmp_path: Path):
    path = tmp_path / "x.vep.tsv"
    path.write_text(VEP_TAB_FILE + "\n", encoding="utf-8")
    frame = load_vep_tab_table(path)
    assert list(frame.columns) == [
        "Uploaded_variation", "Location", "Allele", "Gene", "Consequence",
    ]
    assert len(frame) == 2
    assert frame["Uploaded_variation"].tolist() == ["1-275", "0-41"]


def test_most_severe_picks_frameshift_before_missense():
    assert most_severe_consequence("missense_variant,frameshift_variant") == "frameshift_variant"


def test_subgenome_from_iwgsc_gene_id():
    assert subgenome_from_gene_id("TraesCS1A02G000100") == "A"
    assert subgenome_from_gene_id("TraesCS3B02G123400") == "B"
    assert subgenome_from_gene_id("not_a_wheat_gene") == ""


def test_mask_assignment_and_terminal_demotion(tmp_path: Path):
    raw = pd.DataFrame(
        [
            {
                "Uploaded_variation": "1-100",
                "Gene": "TraesCS1A02G000100",
                "Consequence": "stop_gained",
                "Protein_position": "40/400",
                "SIFT": "-",
                "PolyPhen": "-",
            },
            {
                "Uploaded_variation": "1-200",
                "Gene": "TraesCS1B02G000100",
                "Consequence": "stop_gained",
                "Protein_position": "392/400",
                "SIFT": "-",
                "PolyPhen": "-",
            },
            {
                "Uploaded_variation": "1-300",
                "Gene": "TraesCS1D02G000100",
                "Consequence": "missense_variant",
                "Protein_position": "10/200",
                "SIFT": "deleterious(0.01)",
                "PolyPhen": "probably_damaging(0.98)",
            },
            {
                "Uploaded_variation": "1-400",
                "Gene": "TraesCS2A02G000100",
                "Consequence": "missense_variant",
                "Protein_position": "11/200",
                "SIFT": "deleterious(0.02)",
                "PolyPhen": "benign(0.01)",
            },
        ]
    )
    classified = classify_vep_table(raw)
    masks = classified.set_index("Uploaded_variation")["coding_mask"]
    assert masks["1-100"] == MASK_PLOF_HIGH
    assert masks["1-200"] == MASK_PROTEIN_ALTERING
    assert masks["1-300"] == MASK_MIS_CONSENSUS
    assert masks["1-400"] == MASK_MIS_SEVERE
    assert classified.set_index("Uploaded_variation").loc["1-100", "subgenome"] == "A"

    vep_path = tmp_path / "vep.tsv"
    raw.to_csv(vep_path, sep="\t", index=False)
    out = run_coding_mask(vep_path, tmp_path / "chr_pilot")
    assert (tmp_path / "chr_pilot.coding_mask.tsv").is_file()
    assert int((out["coding_mask"] == MASK_PLOF_HIGH).sum()) == 1
