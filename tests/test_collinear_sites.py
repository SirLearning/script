"""Tests for collinear syntenic site overlap and sampling."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from genetics.genomics.alignment.collinear_sites import (
    _count_overlaps,
    _reservoir_sample,
    build_common_variant_sites,
    collinear_layout,
    extract_common_variant_sites_chr,
    plot_full_lineage_subgenome_fraction,
    prepare_collinear_test_dataset,
    proportional_sample_candidates,
    CollinearCandidate,
)


def _write_pvar(path: Path, plink_chr: int, rows: list[tuple[int, str, float]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tINFO\n")
        for pos, vid, maf in rows:
            handle.write(f"{plink_chr}\t{pos}\t{vid}\tA\tT\tMAF={maf}\n")


def _write_full_lineage(path: Path, chr_label: str, positions: list[int]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("chr_label\tref_pos_1based\n")
        for pos in positions:
            handle.write(f"{chr_label}\t{pos}\n")


def test_plot_full_lineage_subgenome_fraction(tmp_path: Path) -> None:
    summary = tmp_path / "summary.tsv"
    summary.write_text(
        "scope\tn_total\tn_full_lineage\tfraction_full_lineage\n"
        "A\t100\t33\t0.33\n"
        "B\t100\t24\t0.24\n"
        "D\t100\t16\t0.16\n"
        "ALL\t300\t73\t0.243333\n",
        encoding="utf-8",
    )
    out_png = tmp_path / "plot.png"
    df = plot_full_lineage_subgenome_fraction(summary, out_png)
    assert len(df) == 3
    assert out_png.exists()


def test_extract_common_variant_sites_chr(tmp_path: Path) -> None:
    main_raw = tmp_path / "main_raw"
    main_raw.mkdir()
    _write_pvar(main_raw / "chr001.plink2.pvar", 1, [(10, "1-10", 0.10), (20, "1-20", 0.01)])
    out = tmp_path / "common_variant_sites" / "chr001.tsv.gz"
    n = extract_common_variant_sites_chr("chr001", main_raw, out, maf_threshold=0.05)
    assert n == 1
    with gzip.open(out, "rt", encoding="utf-8") as handle:
        lines = handle.readlines()
    assert lines[0].startswith("chr_label")
    assert "1A\t10\t1-10" in lines[1]


def test_join_and_reservoir_sample(tmp_path: Path) -> None:
    main_raw = tmp_path / "main_raw"
    main_raw.mkdir()
    layout = collinear_layout(tmp_path / "collinear")
    _write_full_lineage(layout["full_lineage_sites"] / "1A.tsv.gz", "1A", [10, 20, 30])
    _write_pvar(main_raw / "chr001.plink2.pvar", 1, [(10, "1-10", 0.10), (30, "1-30", 0.20)])
    build_common_variant_sites(tmp_path / "main_raw", layout["common_variant_sites"], workers=1)

    totals, per_chr = _count_overlaps(layout["full_lineage_sites"], layout["common_variant_sites"])
    assert totals["n_full_lineage"] == 3
    assert totals["n_overlap_maf"] == 2

    selected = _reservoir_sample(
        layout["full_lineage_sites"],
        layout["common_variant_sites"],
        target_n=1,
        seed=1,
        totals=totals,
        per_chr_rows=per_chr,
    )
    assert len(selected) == 1


def test_proportional_sample_candidates() -> None:
    candidates = [
        CollinearCandidate("1A", "A", 1, "chr001", "1-1", 0.1),
        CollinearCandidate("1B", "B", 2, "chr003", "3-2", 0.2),
    ]
    selected = proportional_sample_candidates(candidates, target_n=1, seed=1)
    assert len(selected) == 1


def test_prepare_collinear_test_dataset_smoke(tmp_path: Path) -> None:
    layout = collinear_layout(tmp_path / "collinear")
    matrix_dir = tmp_path / "matrices"
    main_raw = tmp_path / "main_raw"
    matrix_dir.mkdir()
    main_raw.mkdir()

    _write_pvar(main_raw / "chr001.plink2.pvar", 1, [(5, "1-5", 0.10)])
    for label in ["1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D",
                  "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"]:
        _write_full_lineage(layout["full_lineage_sites"] / f"{label}.tsv.gz", label, [5] if label == "1A" else [])

    build_common_variant_sites(main_raw, layout["common_variant_sites"], workers=1)

    summary = tmp_path / "full_lineage_summary.tsv"
    summary.write_text(
        "scope\tn_total\tn_full_lineage\tfraction_full_lineage\n"
        "A\t1\t1\t1.0\nB\t1\t0\t0.0\nD\t1\t0\t0.0\nALL\t3\t1\t0.333333\n",
        encoding="utf-8",
    )

    # Skip matrix mawk step — site caches already populated.
    result = prepare_collinear_test_dataset(
        syntenic_matrix_dir=matrix_dir,
        main_raw_dir=main_raw,
        full_lineage_summary_tsv=summary,
        output_dir=layout["root"],
        target_n=1,
    )
    assert result["n_selected"] == 1
    assert (layout["full_lineage_sites"] / "1A.tsv.gz").exists()
    assert (layout["common_variant_sites"] / "chr001.tsv.gz").exists()
    assert (layout["extract"] / "chr001.collinear.txt").exists()
