"""WAP China keep list and SINGER 2 Mb window selection."""

import random
from pathlib import Path

from genetics.germplasm.sample.wap_china import (
    china_wap_taxa,
    resolve_r1_id_to_taxa,
    write_plink2_keep,
)
from genetics.genomics.variant.singer_arg import (
    count_pvar_window_bins,
    list_windows_missing_arg,
    plot_singer_arg_preview,
    random_singer_seeds,
    select_spread_windows,
    singer_window_id,
    summarize_singer_windows,
    window_overlaps_centromere,
)


def test_resolve_r1_id_exact_and_a_suffix():
    available = {"P0001", "P0030A", "P2094A"}
    assert resolve_r1_id_to_taxa("P0001", available) == "P0001"
    assert resolve_r1_id_to_taxa("P0030", available) == "P0030A"
    assert resolve_r1_id_to_taxa("P9999", available) is None


def test_china_wap_taxa_intersects_psam(tmp_path: Path):
    tsv = tmp_path / "germplasm.tsv"
    tsv.write_text(
        "r1_id\tcountry\n"
        "P0001\tChina\n"
        "P0030\tChina\n"
        "P0007\tAustria\n"
        "P0099\tChina\n",
        encoding="utf-8",
    )
    taxa = china_wap_taxa(tsv, {"P0001", "P0030A", "P0007"}, country="China")
    assert taxa == ["P0001", "P0030A"]


def test_write_plink2_keep_fid_iid(tmp_path: Path):
    path = tmp_path / "china.keep"
    write_plink2_keep(["P0001", "P0030A"], path)
    assert path.read_text(encoding="utf-8") == "0\tP0001\n0\tP0030A\n"


def test_centromere_overlap_chr1a():
    assert window_overlaps_centromere("chr1A", 210_000_000, 212_000_000)
    assert not window_overlaps_centromere("chr1A", 2_000_000, 4_000_000)


def test_spread_windows_round_robin(tmp_path: Path):
    pvar = tmp_path / "A.pvar"
    lines = ["#CHROM\tPOS\tID\tREF\tALT\n"]
    for pos in range(2_000_100, 2_000_100 + 80 * 100, 100):
        lines.append(f"1\t{pos}\ts\tA\tT\n")
    for pos in range(10_000_100, 10_000_100 + 60 * 100, 100):
        lines.append(f"1\t{pos}\ts\tA\tT\n")
    pvar.write_text("".join(lines), encoding="utf-8")
    bins = count_pvar_window_bins(pvar, subgenome="A", window_bp=2_000_000)
    assert any(row["start_bp"] == 2_000_000 and row["n_snps"] == 80 for row in bins)

    mixed = [
        {"subgenome": "A", "plink_chr": "1", "ref_chrom": "chr1A", "start_bp": 2_000_000, "end_bp": 4_000_000, "n_snps": 80},
        {"subgenome": "A", "plink_chr": "1", "ref_chrom": "chr1A", "start_bp": 10_000_000, "end_bp": 12_000_000, "n_snps": 60},
        {"subgenome": "B", "plink_chr": "3", "ref_chrom": "chr1B", "start_bp": 2_000_000, "end_bp": 4_000_000, "n_snps": 200},
        {"subgenome": "B", "plink_chr": "3", "ref_chrom": "chr1B", "start_bp": 8_000_000, "end_bp": 10_000_000, "n_snps": 90},
        {"subgenome": "D", "plink_chr": "5", "ref_chrom": "chr1D", "start_bp": 2_000_000, "end_bp": 4_000_000, "n_snps": 55},
    ]
    picked = select_spread_windows(mixed, n_windows=3, min_snps=50)
    chroms = [row["ref_chrom"] for row in picked]
    assert chroms == ["chr1A", "chr1B", "chr1D"]
    assert picked[1]["n_snps"] == 200


def test_list_windows_missing_arg_skips_complete(tmp_path: Path):
    windows = tmp_path / "singer_windows.tsv"
    windows.write_text(
        "subgenome\tplink_chr\tref_chrom\tstart_bp\tend_bp\tn_snps\n"
        "A\t1\tchr1A\t2000000\t4000000\t100\n"
        "B\t3\tchr1B\t2000000\t4000000\t200\n",
        encoding="utf-8",
    )
    singer_dir = tmp_path / "singer"
    singer_dir.mkdir()
    ok_id = singer_window_id(
        {
            "subgenome": "A",
            "ref_chrom": "chr1A",
            "start_bp": 2_000_000,
            "end_bp": 4_000_000,
        }
    )
    (singer_dir / f"singer_{ok_id}_nodes_0.txt").write_text("ok\n", encoding="utf-8")
    missing = list_windows_missing_arg(windows, singer_dir)
    assert [row["window_id"] for row in missing] == ["B.chr1B.2000000_4000000"]


def test_random_singer_seeds_are_ours_not_default_42():
    seeds = random_singer_seeds(5, rng=random.Random(7))
    assert len(seeds) == 5
    assert len(set(seeds)) == 5
    assert 42 not in seeds


def _write_preview_fixtures(tmp_path: Path) -> tuple[Path, Path]:
    windows = tmp_path / "singer_windows.tsv"
    windows.write_text(
        "subgenome\tplink_chr\tref_chrom\tstart_bp\tend_bp\tn_snps\n"
        "A\t1\tchr1A\t2000000\t4000000\t100\n"
        "B\t3\tchr1B\t2000000\t4000000\t200\n",
        encoding="utf-8",
    )
    singer_dir = tmp_path / "singer"
    singer_dir.mkdir()
    ok_id = "A.chr1A.2000000_4000000"
    (singer_dir / f"singer_{ok_id}_nodes_0.txt").write_text("0\n0\n1.5\n", encoding="utf-8")
    (singer_dir / f"singer_{ok_id}_recombs_0.txt").write_text(
        "10 1 2 0.1\n20 3 4 0.2\n",
        encoding="utf-8",
    )
    (singer_dir / f"singer_{ok_id}_muts_0.txt").write_text("1 2 3 1\n", encoding="utf-8")
    (singer_dir / f"singer_{ok_id}.run.log").write_text(
        "estimated_Ne=12.5\nvalid mutations: 3\n",
        encoding="utf-8",
    )
    return windows, singer_dir


def test_summarize_singer_windows_marks_ok_and_fail(tmp_path: Path):
    windows, singer_dir = _write_preview_fixtures(tmp_path)
    rows = summarize_singer_windows(windows, singer_dir)
    assert [row["window_id"] for row in rows] == [
        "A.chr1A.2000000_4000000",
        "B.chr1B.2000000_4000000",
    ]
    ok, fail = rows
    assert ok["arg_ok"] is True
    assert ok["estimated_ne"] == 12.5
    assert ok["n_nodes"] == 3
    assert ok["n_leaves"] == 2
    assert ok["t_max"] == 1.5
    assert ok["n_recomb"] == 2
    assert ok["n_muts"] == 1
    assert ok["n_arg_samples"] == 1
    assert fail["arg_ok"] is False
    assert fail["estimated_ne"] is None


def test_plot_singer_arg_preview_writes_info_and_plots(tmp_path: Path):
    windows, singer_dir = _write_preview_fixtures(tmp_path)
    prefix = tmp_path / "out" / "all.variant.singer_arg"
    written = plot_singer_arg_preview(windows, singer_dir, prefix)
    assert (tmp_path / "out" / "all.variant.singer_arg.window.info.tsv").is_file()
    assert (tmp_path / "out" / "all.variant.singer_arg.window.th.tsv").is_file()
    pngs = list((tmp_path / "out").glob("all.variant.singer_arg*.png"))
    assert pngs
    assert written["n_ok"] == 1
    assert written["n_fail"] == 1
