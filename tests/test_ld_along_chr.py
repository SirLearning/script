"""Evenly sampled along-chromosome LD: distance where mean r2 falls to 0.2."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from genetics.genomics.variant.ld_along_chr import (
    aggregate_decay_by_spatial_bin,
    ana_ld_r2_decay_along_chr,
    decay_plot_ylim,
    evenly_sample_variants_along_chromosome,
    ld_along_chr_output_stem,
    load_pvar_ref_variants,
    mean_r2_decay_distance_kb,
    plink_chr_ids_for_ref_chrom,
    plot_decay_along_subgenome,
    sample_even_focals_from_pvar,
    sample_focals_per_spatial_bin,
    subgenome_ref_chroms,
    summarize_decay_distance_by_chromosome,
    write_all_ref_merge_update_maps,
    write_ref_merge_update_maps,
)
from infra.wheat.ref_v1 import plink_segment_to_ref_position


VCOR_HEADER = "#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2\n"
PVAR_HEADER = "#CHROM\tPOS\tID\tREF\tALT\n"


def test_evenly_sample_picks_nearest_to_even_targets():
    variants = pd.DataFrame(
        {
            "ID": [f"s-{p}" for p in range(0, 101, 10)],
            "ref_pos": list(range(0, 101, 10)),
        }
    )
    sampled = evenly_sample_variants_along_chromosome(variants, n=5)
    assert list(sampled["ref_pos"]) == [0, 20, 50, 70, 100]


def test_evenly_sample_does_not_reuse_the_same_variant():
    variants = pd.DataFrame({"ID": ["a", "b", "c"], "ref_pos": [0, 50, 100]})
    sampled = evenly_sample_variants_along_chromosome(variants, n=5)
    assert len(sampled) == 3
    assert sampled["ID"].is_unique


def test_mean_r2_decay_interpolates_crossing_of_threshold():
    # Bin mids 25 / 75 / 125 kb; means 0.80 / 0.50 / 0.10. Threshold 0.2
    # crosses between 75 (0.50) and 125 (0.10) at 112.5 kb.
    distance_kb = []
    r2 = []
    for start, mean in ((0, 0.80), (50, 0.50), (100, 0.10)):
        for offset in (5, 15, 25, 35, 45):
            distance_kb.append(start + offset)
            r2.append(mean)
    decay_kb, reached, n_pairs = mean_r2_decay_distance_kb(
        distance_kb,
        r2,
        threshold=0.2,
        bin_size_kb=50.0,
        min_pairs_per_bin=5,
        max_distance_kb=150.0,
    )
    assert reached
    assert n_pairs == 15
    assert decay_kb == pytest.approx(112.5)


def test_summarize_decay_distance_by_chromosome_uses_reached_focals():
    focals = pd.DataFrame(
        {
            "ref_chrom": ["1A", "1A", "1A", "1B"],
            "reached": [True, True, False, True],
            "decay_kb": [100.0, 300.0, 20000.0, 50.0],
        }
    )
    out = summarize_decay_distance_by_chromosome(focals, r2_threshold=0.5)
    by_chr = out.set_index("chromosome")
    assert by_chr.loc["chr1A", "n_focals"] == 3
    assert by_chr.loc["chr1A", "n_reached"] == 2
    assert by_chr.loc["chr1A", "median_decay_kb"] == pytest.approx(200.0)
    assert by_chr.loc["chr1A", "mean_decay_kb"] == pytest.approx(200.0)
    assert by_chr.loc["chr1B", "median_decay_kb"] == pytest.approx(50.0)
    assert by_chr.loc["chr1A", "r2_threshold"] == pytest.approx(0.5)


def test_ld_along_chr_output_stem_keeps_default_r2_0p2_name():
    assert ld_along_chr_output_stem(0.2) == "ld_r2_decay_along"
    assert ld_along_chr_output_stem(0.5) == "ld_r2_decay_along.r2_0p5"


def test_mean_r2_decay_returns_nan_when_never_crossing():
    distance_kb = np.repeat(np.arange(25, 500, 50), 6)
    r2 = np.full(distance_kb.shape, 0.6)
    decay_kb, reached, n_pairs = mean_r2_decay_distance_kb(
        distance_kb,
        r2,
        threshold=0.2,
        bin_size_kb=50.0,
        min_pairs_per_bin=5,
        max_distance_kb=500.0,
    )
    assert not reached
    assert n_pairs == len(distance_kb)
    assert np.isnan(decay_kb)


def test_decay_plot_ylim_follows_data_not_window_cap():
    lo, hi = decay_plot_ylim([50.0, 192.2686, 120.1])
    assert lo == 0.0
    assert hi == pytest.approx(192.2686 * 1.15)
    lo_empty, hi_empty = decay_plot_ylim([np.nan, np.nan])
    assert lo_empty == 0.0
    assert hi_empty == 50.0


def test_plink_chr_ids_for_chr3_subgenomes():
    assert plink_chr_ids_for_ref_chrom("3A") == ["13", "14"]
    assert plink_chr_ids_for_ref_chrom("chr3B") == ["15", "16"]
    assert plink_chr_ids_for_ref_chrom("3D") == ["17", "18"]


def test_plink_segment_to_ref_position_joins_chr3_arms():
    assert plink_segment_to_ref_position(13, 100) == ("3A", 100)
    assert plink_segment_to_ref_position(14, 100) == ("3A", 454103970 + 100)
    assert plink_segment_to_ref_position(15, 50) == ("3B", 50)
    assert plink_segment_to_ref_position(16, 50) == ("3B", 448155269 + 50)
    assert plink_segment_to_ref_position(17, 10) == ("3D", 10)
    assert plink_segment_to_ref_position(18, 10) == ("3D", 476235359 + 10)


def test_write_ref_merge_update_maps_collapses_both_arms(tmp_path: Path):
    pvar = tmp_path / "A_test.plink2.pvar"
    pvar.write_text(
        PVAR_HEADER
        + "13\t100\t13-100\tA\tT\n"
        + "14\t50\t14-50\tA\tT\n"
        + "15\t9\t15-9\tA\tT\n",
        encoding="utf-8",
    )
    chr_map = tmp_path / "update.chr"
    pos_map = tmp_path / "update.map"
    sub = write_ref_merge_update_maps(str(pvar), "3A", str(chr_map), str(pos_map))
    assert set(sub["ID"]) == {"13-100", "14-50"}
    chr_rows = {
        line.split("\t")[0]: line.split("\t")[1]
        for line in chr_map.read_text(encoding="utf-8").splitlines()
    }
    pos_rows = {
        line.split("\t")[0]: int(line.split("\t")[1])
        for line in pos_map.read_text(encoding="utf-8").splitlines()
    }
    assert chr_rows == {"13-100": "chr3A", "14-50": "chr3A"}
    assert pos_rows == {"13-100": 100, "14-50": 454103970 + 50}


def test_sample_even_focals_from_pvar_filters_chromosome(tmp_path: Path):
    pvar = tmp_path / "A_test.plink2.pvar"
    lines = [PVAR_HEADER]
    for pos in range(1_000_000, 10_000_001, 1_000_000):
        lines.append(f"13\t{pos}\t13-{pos}\tA\tT\n")
    for pos in range(1_000_000, 4_000_001, 1_000_000):
        lines.append(f"14\t{pos}\t14-{pos}\tA\tT\n")
    for pos in range(1_000_000, 4_000_001, 1_000_000):
        lines.append(f"15\t{pos}\t15-{pos}\tA\tT\n")
    pvar.write_text("".join(lines), encoding="utf-8")

    focals = sample_even_focals_from_pvar(
        str(pvar),
        chromosomes=["3A"],
        n_per_chr=3,
        output_prefix=str(tmp_path / "A.3A"),
    )
    assert set(focals["ref_chrom"]) == {"3A"}
    assert len(focals) == 3
    assert {13, 14} <= set(focals["plink_chr"])
    ids = Path(tmp_path / "A.3A.focals.id").read_text(encoding="utf-8").splitlines()
    assert ids == list(focals["ID"])
    assert Path(tmp_path / "A.3A.plink_chr_ids.txt").read_text(encoding="utf-8") == "13,14"


def test_sample_even_focals_writes_merged_label_for_recoded_pvar(tmp_path: Path):
    """A ref-chr recoded pvar names the whole chromosome, so --chr gets one token."""
    pvar = tmp_path / "A_test.plink2.pvar"
    lines = [PVAR_HEADER]
    for pos in (100, 454103970 + 50):
        lines.append(f"chr3A\t{pos}\t13-{pos}\tA\tT\n")
    pvar.write_text("".join(lines), encoding="utf-8")
    sample_even_focals_from_pvar(
        str(pvar),
        chromosomes=["3A"],
        n_per_chr=2,
        output_prefix=str(tmp_path / "A.3A"),
    )
    text = Path(tmp_path / "A.3A.plink_chr_ids.txt").read_text(encoding="utf-8")
    assert text == "chr3A"


def test_recoded_pvar_keeps_linear_pos_without_re_offsetting(tmp_path: Path):
    """POS in a recoded pvar is already whole-chromosome, so no arm offset is added."""
    pvar = tmp_path / "A_test.plink2.pvar"
    pvar.write_text(
        PVAR_HEADER + f"chr3A\t{454103970 + 50}\t14-50\tA\tT\n",
        encoding="utf-8",
    )
    variants = load_pvar_ref_variants(str(pvar))
    assert list(variants["ref_chrom"]) == ["3A"]
    assert list(variants["ref_pos"]) == [454103970 + 50]


def test_legacy_first_arm_coded_merged_pvar_keeps_linear_pos(tmp_path: Path):
    """Ref-chr pfiles made before the naming change still read correctly.

    They carry the first arm's id with whole-chromosome POS. A first arm starts at
    ref_start 0, so the segment mapping is the identity and POS survives unchanged.
    """
    pvar = tmp_path / "A_test.plink2.pvar"
    pvar.write_text(
        PVAR_HEADER + f"13\t{454103970 + 50}\t14-50\tA\tT\n",
        encoding="utf-8",
    )
    variants = load_pvar_ref_variants(str(pvar))
    assert list(variants["ref_chrom"]) == ["3A"]
    assert list(variants["ref_pos"]) == [454103970 + 50]


def test_write_all_ref_merge_update_maps_names_every_chromosome(tmp_path: Path):
    """Every arm-numbered row is renamed; only second-arm rows also move POS."""
    pvar = tmp_path / "A_test.plink2.pvar"
    pvar.write_text(
        PVAR_HEADER
        + "13\t100\t13-100\tA\tT\n"
        + "14\t50\t14-50\tA\tT\n"
        + "1\t20\t1-20\tA\tT\n"
        + "2\t10\t2-10\tA\tT\n"
        + "0\t5\t0-5\tA\tT\n"
        + "43\t7\t43-7\tA\tT\n",
        encoding="utf-8",
    )
    chr_map = tmp_path / "update.chr"
    pos_map = tmp_path / "update.map"
    updated = write_all_ref_merge_update_maps(str(pvar), str(chr_map), str(pos_map))
    assert set(updated["ID"]) == {"13-100", "14-50", "1-20", "2-10", "0-5", "43-7"}
    chr_rows = {
        line.split("\t")[0]: line.split("\t")[1]
        for line in chr_map.read_text(encoding="utf-8").splitlines()
    }
    pos_rows = {
        line.split("\t")[0]: int(line.split("\t")[1])
        for line in pos_map.read_text(encoding="utf-8").splitlines()
    }
    assert chr_rows == {
        "13-100": "chr3A",
        "14-50": "chr3A",
        "1-20": "chr1A",
        "2-10": "chr1A",
        "0-5": "chrUn",
        "43-7": "Mit",
    }
    # First arms and organelles keep POS; second arms gain the first-arm length.
    assert pos_rows["13-100"] == 100
    assert pos_rows["1-20"] == 20
    assert pos_rows["0-5"] == 5
    assert pos_rows["43-7"] == 7
    assert pos_rows["14-50"] == 454103970 + 50
    assert pos_rows["2-10"] == 471304005 + 10


def test_write_all_ref_merge_update_maps_is_a_no_op_on_recoded_pvar(tmp_path: Path):
    """Re-running the recode on an already-named pvar must not emit any update rows."""
    pvar = tmp_path / "A_test.plink2.pvar"
    pvar.write_text(
        PVAR_HEADER
        + "chr3A\t100\t13-100\tA\tT\n"
        + f"chr3A\t{454103970 + 50}\t14-50\tA\tT\n",
        encoding="utf-8",
    )
    chr_map = tmp_path / "update.chr"
    pos_map = tmp_path / "update.map"
    updated = write_all_ref_merge_update_maps(str(pvar), str(chr_map), str(pos_map))
    assert updated.empty
    assert chr_map.read_text(encoding="utf-8") == ""


def test_ana_ld_r2_decay_along_chr_longer_near_high_ld_focal(tmp_path: Path):
    """A focal whose partners stay at high r2 farther away gets a longer decay distance."""
    vcor = tmp_path / "A.3A.focal.ld.vcor"
    rows = [VCOR_HEADER]
    short_focal = ("13", 10_000_000, "13-10000000")
    long_focal = ("13", 317_000_000, "13-317000000")
    for focal, hold_kb in ((short_focal, 100), (long_focal, 2000)):
        chrom, pos, vid = focal
        for dist_bp in range(5_000, 3_000_001, 5_000):
            r2 = 0.85 if dist_bp <= hold_kb * 1000 else 0.05
            partner_pos = pos + dist_bp
            rows.append(
                f"{chrom}\t{pos}\t{vid}\t{chrom}\t{partner_pos}\t"
                f"{chrom}-{partner_pos}\t{r2:.4f}\n"
            )
    vcor.write_text("".join(rows), encoding="utf-8")

    focals = pd.DataFrame(
        {
            "ID": ["13-10000000", "13-317000000"],
            "plink_chr": [13, 13],
            "pos": [10_000_000, 317_000_000],
            "ref_chrom": ["3A", "3A"],
            "ref_pos": [10_000_000, 317_000_000],
        }
    )
    focals_path = tmp_path / "A.3A.focals.info.tsv"
    focals.to_csv(focals_path, sep="\t", index=False)

    prefix = str(tmp_path / "A.variant.ld_r2_decay_along.3A")
    table = ana_ld_r2_decay_along_chr(
        input_file=str(vcor),
        focals_file=str(focals_path),
        output_prefix=prefix,
        r2_threshold=0.2,
        max_distance_kb=5000,
        bin_size_kb=50,
        min_pairs_per_bin=5,
    )
    by_id = table.set_index("ID")
    assert by_id.loc["13-317000000", "reached"]
    assert by_id.loc["13-10000000", "reached"]
    assert (
        by_id.loc["13-317000000", "decay_kb"]
        > by_id.loc["13-10000000", "decay_kb"] * 5
    )
    assert Path(f"{prefix}.info.tsv").is_file()
    assert Path(f"{prefix}.png").is_file()
    assert Path(f"{prefix}.th.tsv").is_file()


def test_sample_focals_per_spatial_bin_caps_at_n_per_bin():
    rows = []
    for bin_i in range(3):
        start = bin_i * 5_000_000 + 1
        for k in range(120):
            pos = start + k * 40_000
            rows.append({"ID": f"s-{bin_i}-{k}", "ref_pos": pos, "ref_chrom": "3A"})
    variants = pd.DataFrame(rows)
    sampled = sample_focals_per_spatial_bin(
        variants, n_per_bin=100, bin_mb=5.0, chrom_length_bp=15_000_000
    )
    assert sampled["bin_idx"].nunique() == 3
    assert (sampled.groupby("bin_idx").size() == 100).all()
    assert sampled["ID"].is_unique


def test_sample_focals_per_spatial_bin_keeps_all_when_bin_is_small():
    variants = pd.DataFrame(
        {
            "ID": [f"s-{i}" for i in range(40)],
            "ref_pos": [1 + i * 1000 for i in range(40)],
            "ref_chrom": "3A",
        }
    )
    sampled = sample_focals_per_spatial_bin(
        variants, n_per_bin=100, bin_mb=5.0, chrom_length_bp=5_000_000
    )
    assert len(sampled) == 40
    assert sampled["n_in_bin"].iloc[0] == 40


def test_aggregate_decay_by_spatial_bin_means_focal_distances():
    focals = pd.DataFrame(
        {
            "ID": ["a", "b", "c", "d"],
            "ref_chrom": ["3A", "3A", "3A", "3A"],
            "bin_idx": [0, 0, 0, 1],
            "bin_start_bp": [1, 1, 1, 5_000_001],
            "bin_end_bp": [5_000_000, 5_000_000, 5_000_000, 10_000_000],
            "bin_mid_mb": [2.5, 2.5, 2.5, 7.5],
            "n_in_bin": [200, 200, 200, 200],
            "decay_kb": [10.0, 20.0, 30.0, 100.0],
            "reached": [True, True, True, True],
        }
    )
    binned = aggregate_decay_by_spatial_bin(focals)
    by_bin = binned.set_index("bin_idx")
    assert by_bin.loc[0, "mean_decay_kb"] == pytest.approx(20.0)
    assert by_bin.loc[0, "n_focals"] == 3
    assert by_bin.loc[1, "mean_decay_kb"] == pytest.approx(100.0)
    assert by_bin.loc[0, "plot_decay_kb"] == pytest.approx(20.0)
    assert bool(by_bin.loc[0, "reached"])


def test_aggregate_unreached_bin_plots_at_window_cap():
    focals = pd.DataFrame(
        {
            "ID": ["a", "b", "c"],
            "ref_chrom": ["3A", "3A", "3A"],
            "bin_idx": [0, 0, 1],
            "bin_start_bp": [1, 1, 5_000_001],
            "bin_end_bp": [5_000_000, 5_000_000, 10_000_000],
            "bin_mid_mb": [2.5, 2.5, 7.5],
            "n_in_bin": [10, 10, 8],
            "decay_kb": [50.0, 80.0, np.nan],
            "reached": [True, True, False],
        }
    )
    binned = aggregate_decay_by_spatial_bin(focals, max_distance_kb=20000.0)
    by_bin = binned.set_index("bin_idx")
    assert by_bin.loc[0, "plot_decay_kb"] == pytest.approx(65.0)
    assert bool(by_bin.loc[0, "reached"])
    assert not bool(by_bin.loc[1, "reached"])
    assert np.isnan(by_bin.loc[1, "mean_decay_kb"])
    assert by_bin.loc[1, "plot_decay_kb"] == pytest.approx(20000.0)


def test_subgenome_ref_chroms_are_homoeologs_1_to_7():
    assert subgenome_ref_chroms("A") == ["1A", "2A", "3A", "4A", "5A", "6A", "7A"]
    assert subgenome_ref_chroms("B") == ["1B", "2B", "3B", "4B", "5B", "6B", "7B"]
    assert subgenome_ref_chroms("D") == ["1D", "2D", "3D", "4D", "5D", "6D", "7D"]


def test_plot_decay_along_subgenome_stacks_chr1_to_chr7(tmp_path: Path):
    rows = []
    for i in range(1, 8):
        chrom = f"{i}A"
        for mid, decay in ((10.0, 80.0 * i), (80.0, 120.0 * i)):
            rows.append(
                {
                    "ref_chrom": chrom,
                    "bin_mid_mb": mid,
                    "ref_pos_mb": mid,
                    "plot_decay_kb": decay,
                    "decay_kb": decay,
                    "mean_decay_kb": decay,
                    "reached": True,
                    "bin_start_bp": int(mid * 1e6 - 2.5e6),
                    "bin_end_bp": int(mid * 1e6 + 2.5e6),
                }
            )
    prefix = str(tmp_path / "A.variant.ld_r2_decay_along")
    fig = plot_decay_along_subgenome(
        pd.DataFrame(rows), "A", prefix, max_distance_kb=20000.0
    )
    assert Path(f"{prefix}.png").is_file()
    assert Path(f"{prefix}.png").stat().st_size > 1000
    axes = getattr(fig, "axes", [])
    assert len(axes) == 7
    titles = [ax.get_title() for ax in axes]
    assert titles == [f"chr{i}A" for i in range(1, 8)]

