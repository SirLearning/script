"""WGS-core vs WAP allele-frequency join."""

from pathlib import Path

from genetics.genomics.variant.dual_af import merge_cohort_af, run_dual_af


def _write_afreq(path: Path, rows: list[tuple[str, float]]) -> Path:
    lines = ["ID\tALT_FREQS"]
    lines.extend(f"{vid}\t{af}" for vid, af in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_dual_af_marks_shared_rare_and_discordant(tmp_path: Path):
    wgs = _write_afreq(
        tmp_path / "wgs.afreq",
        [("1-10", 0.002), ("1-20", 0.20), ("1-30", 0.008)],
    )
    wap = _write_afreq(
        tmp_path / "wap.afreq",
        [("1-10", 0.003), ("1-20", 0.22), ("1-30", 0.04), ("1-99", 0.5)],
    )
    merged = merge_cohort_af(wgs, wap, rare_max=0.01)
    by_id = merged.set_index("ID")
    assert set(merged["ID"]) == {"1-10", "1-20", "1-30"}
    assert bool(by_id.loc["1-10", "rare_both"])
    assert not bool(by_id.loc["1-20", "rare_wgs"])
    assert bool(by_id.loc["1-30", "rare_discordant"])


def test_run_dual_af_writes_tables(tmp_path: Path):
    wgs = _write_afreq(tmp_path / "wgs.afreq", [("2-5", 0.1)])
    wap = _write_afreq(tmp_path / "wap.afreq", [("2-5", 0.12)])
    run_dual_af(wgs, wap, tmp_path / "af", rare_max=0.01)
    assert (tmp_path / "af.dual_af.tsv").is_file()
    assert (tmp_path / "af.dual_af.th.tsv").is_file()
