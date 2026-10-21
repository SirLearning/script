"""Cohort roles: WAP is the phenotyped association panel inside WGS-core."""

from pathlib import Path

import pytest

from genetics.germplasm.sample.cohort_manifest import (
    KEEP_ASSOC,
    KEEP_WGS_CORE,
    ROLE_ASSOC,
    ROLE_WGS_CORE,
    ROLE_WGS_OTHER,
    build_sample_manifest,
    run_sample_manifest,
    summarize_cohort_qc,
)


def _write_group(path: Path, rows: list[tuple[str, str]]) -> Path:
    path.write_text("".join(f"{sample} {group}\n" for sample, group in rows), encoding="utf-8")
    return path


def test_wap_is_assoc_and_inside_wgs_core(tmp_path: Path):
    group = _write_group(
        tmp_path / "sample_group.txt",
        [("P0001", "WAP"), ("ABD_0001", "ABD"), ("WAT_1", "Watkins"), ("D_1", "D")],
    )
    manifest = build_sample_manifest(group)
    by_sample = manifest.set_index("Sample")
    assert by_sample.loc["P0001", "cohort_role"] == ROLE_ASSOC
    assert bool(by_sample.loc["P0001", "in_assoc"])
    assert bool(by_sample.loc["P0001", "in_wgs_core"])
    assert by_sample.loc["ABD_0001", "cohort_role"] == ROLE_WGS_CORE
    assert not bool(by_sample.loc["ABD_0001", "in_assoc"])
    assert by_sample.loc["WAT_1", "cohort_role"] == ROLE_WGS_OTHER
    assert by_sample.loc["D_1", "cohort_role"] == ROLE_WGS_OTHER


def test_keep_files_and_qc_counts(tmp_path: Path):
    group = _write_group(
        tmp_path / "sample_group.txt",
        [("P0001", "WAP"), ("P0002", "WAP"), ("ABD_0001", "ABD"), ("WAT_1", "Watkins")],
    )
    prefix = tmp_path / "cohort"
    manifest = run_sample_manifest(group, prefix)
    assoc_keep = (tmp_path / KEEP_ASSOC).read_text(encoding="utf-8").splitlines()
    wgs_keep = (tmp_path / KEEP_WGS_CORE).read_text(encoding="utf-8").splitlines()
    assert [line.split("\t")[1] for line in assoc_keep] == ["P0001", "P0002"]
    assert {line.split("\t")[1] for line in wgs_keep} == {"P0001", "P0002", "ABD_0001"}
    summary = summarize_cohort_qc(manifest)
    totals = summary.set_index("key").loc["total"]
    assert int(totals["n"]) == 4
    assert int(totals["n_assoc"]) == 2
    assert int(totals["n_wgs_core"]) == 3
    assert (tmp_path / "cohort.manifest.tsv").is_file()
    assert (tmp_path / "cohort.qc.th.tsv").is_file()


def test_taxa_map_coverage_join(tmp_path: Path):
    group = _write_group(tmp_path / "sample_group.txt", [("P0001", "WAP"), ("ABD_0001", "ABD")])
    taxa_dir = tmp_path / "maps"
    taxa_dir.mkdir()
    (taxa_dir / "WAP.taxaBamMap.txt").write_text("P0001\t3.10\t/tmp/p.bam\n", encoding="utf-8")
    manifest = build_sample_manifest(group, taxa_bam_map_dir=taxa_dir)
    by_sample = manifest.set_index("Sample")
    assert bool(by_sample.loc["P0001", "in_taxa_map"])
    assert float(by_sample.loc["P0001", "Coverage"]) == pytest.approx(3.10)
    assert not bool(by_sample.loc["ABD_0001", "in_taxa_map"])
