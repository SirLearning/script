"""Site annotation from the featHighRes conservation / epigenome BED."""

from pathlib import Path

import pandas as pd

from genetics.genomics.variant.site_anno import (
    COL_CODING,
    COL_CONSERVED,
    COL_DEGENERACY,
    COL_ELEMENT,
    COL_NEUTRAL,
    COL_REGULATORY,
    EPI_FEATURE_COLUMNS,
    SITE_ANNO_COLUMNS,
    add_site_classes,
    classify_degeneracy,
    load_site_anno_chunks,
    run_site_anno,
    summarize_site_anno,
)

DEFAULTS = {name: 0 for name in EPI_FEATURE_COLUMNS}
DEFAULTS.update(
    {
        "heptamer_USscore": 1.0,
        "degeneration": 0,
        "triticeae_phyloP": 0.0,
        "poaceae_phyloP": 0.0,
        "monocots_phyloP": 0.0,
        "triticeae_phastCons": 0.0,
        "poaceae_phastCons": 0.0,
        "monocots_phastCons": 0.0,
        "GERP29": 0.0,
        "PPH2": 0.0,
        "LIST2": 0.0,
        "SIFT": 1.0,
    }
)


def _site(variant_id: str, pos: int, **overrides) -> dict:
    row = {"ID": variant_id, "chrom": "chr1A", "pos": pos}
    row.update(DEFAULTS)
    row.update(overrides)
    return row


def _frame(rows) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=list(SITE_ANNO_COLUMNS))


def _write_chunk(path: Path, rows) -> Path:
    _frame(rows).to_csv(path, sep="\t", header=False, index=False)
    return path


def test_degeneracy_codes_follow_the_verified_mapping():
    """2 is non-degenerate and 1 is 4-fold, which is the inverse of the intuitive read."""

    table = _frame(
        [
            _site("1-10", 10, CDS_HC=1, degeneration=2),
            _site("1-20", 20, CDS_HC=1, degeneration=1),
            _site("1-30", 30, CDS_HC=1, degeneration=0),
            _site("1-40", 40, CDS_HC=0, degeneration=1),
        ]
    )
    assert list(classify_degeneracy(table)) == ["0fold", "4fold", "2fold", "noncoding"]


def test_element_takes_the_most_specific_flag():
    table = _frame(
        [
            _site("1-10", 10, CDS_HC=1, TE=1),
            _site("1-20", 20, ncRNA=1, CGI=1),
            _site("1-30", 30, TE=1),
            _site("1-40", 40),
        ]
    )
    assert list(add_site_classes(table)[COL_ELEMENT]) == ["CDS_HC", "ncRNA", "TE", "intergenic"]


def test_conserved_needs_both_alignment_scores_or_gerp_alone():
    table = _frame(
        [
            _site("1-10", 10, monocots_phyloP=3.0, monocots_phastCons=0.95),
            # High phyloP but weak phastCons and low GERP is not enough.
            _site("1-20", 20, monocots_phyloP=3.0, monocots_phastCons=0.10),
            _site("1-30", 30, GERP29=4.0),
            _site("1-40", 40),
        ]
    )
    assert list(add_site_classes(table)[COL_CONSERVED]) == [True, False, True, False]


def test_regulatory_is_the_union_of_the_chromatin_tracks():
    table = _frame(
        [
            _site("1-10", 10, ATAC=1),
            _site("1-20", 20, epiRE=1, DHS=1),
            _site("1-30", 30),
        ]
    )
    assert list(add_site_classes(table)[COL_REGULATORY]) == [True, True, False]


def test_neutral_set_is_fourfold_and_unconstrained():
    table = _frame(
        [
            _site("1-10", 10, CDS_HC=1, degeneration=1),
            # 4-fold but conserved, so not neutral.
            _site("1-20", 20, CDS_HC=1, degeneration=1, GERP29=5.0),
            # 0-fold.
            _site("1-30", 30, CDS_HC=1, degeneration=2),
        ]
    )
    classed = add_site_classes(table)
    assert list(classed[COL_NEUTRAL]) == [True, False, False]
    assert list(classed[COL_CODING]) == [True, True, True]


def test_summary_counts_every_class():
    classed = add_site_classes(
        _frame(
            [
                _site("1-10", 10, CDS_HC=1, degeneration=1),
                _site("1-20", 20, CDS_HC=1, degeneration=2, GERP29=3.0),
                _site("1-30", 30, ATAC=1),
            ]
        )
    )
    stats = summarize_site_anno(classed)
    assert stats["n_sites"] == 3
    assert stats["n_4fold"] == 1
    assert stats["n_0fold"] == 1
    assert stats["n_noncoding"] == 1
    assert stats["n_conserved"] == 1
    assert stats["n_regulatory"] == 1
    assert stats["n_neutral_4fold"] == 1


def test_chunks_are_headerless_and_concatenate(tmp_path: Path):
    _write_chunk(tmp_path / "chr1A.tsv", [_site("1-10", 10), _site("1-20", 20)])
    _write_chunk(tmp_path / "chr2A.tsv", [_site("7-10", 10)])
    merged = load_site_anno_chunks([tmp_path / "chr1A.tsv", tmp_path / "chr2A.tsv"])
    assert list(merged["ID"]) == ["1-10", "1-20", "7-10"]
    assert list(merged.columns) == list(SITE_ANNO_COLUMNS)


def test_run_site_anno_writes_table_and_thresholds(tmp_path: Path):
    chunk = _write_chunk(
        tmp_path / "chunk.tsv",
        [
            _site("1-10", 10, CDS_HC=1, degeneration=1),
            _site("1-20", 20, CDS_HC=1, degeneration=2, monocots_phyloP=4.0, monocots_phastCons=0.9),
            _site("1-30", 30, CDS_HC=1, degeneration=0),
            _site("1-40", 40, ATAC=1),
        ],
    )
    out = run_site_anno([chunk], tmp_path / "A")
    assert (tmp_path / "A.site_anno.tsv").is_file()
    assert (tmp_path / "A.site_anno.th.tsv").is_file()
    assert int(out[COL_NEUTRAL].sum()) == 1
    written = pd.read_csv(tmp_path / "A.site_anno.tsv", sep="\t")
    assert COL_DEGENERACY in written.columns
    assert "monocots_phyloP" in written.columns
