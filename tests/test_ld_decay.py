"""LD decay summaries: per-chromosome mean r2 curves plus overall mean."""

from pathlib import Path

import pandas as pd
import pytest

from genetics.genomics.variant.ld_decay import (
    ana_ld_decay,
    ana_ld_decay_mean_r2_threshold,
    distance_mean_r2_to_threshold,
    homoeolog_subgenome_ld_table,
    ld_decay_chrom_label,
    plot_ld_homoeolog_subgenome_compare,
    summarize_ld_decay_chromosomes,
)


VCOR_HEADER = "#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2\n"


def test_ld_decay_chrom_label_maps_vmap4_ids_to_ref_names():
    assert ld_decay_chrom_label(1) == "chr1A"
    assert ld_decay_chrom_label(2) == "chr1A"
    assert ld_decay_chrom_label(7) == "chr2A"
    assert ld_decay_chrom_label("13") == "chr3A"
    assert ld_decay_chrom_label(0) == "chrUn"
    assert ld_decay_chrom_label(43) == "Mit"


def test_ana_ld_decay_writes_per_chromosome_curves_and_overall_mean(tmp_path: Path):
    vcor = tmp_path / "A.info.ld.vcor"
    vcor.write_text(
        VCOR_HEADER
        + "1\t1000\ts1a\t1\t26000\ts1b\t0.80\n"
        + "1\t1000\ts1a\t1\t51000\ts1c\t0.70\n"
        + "7\t1000\ts2a\t7\t26000\ts2b\t0.20\n"
        + "7\t1000\ts2a\t7\t51000\ts2c\t0.10\n",
        encoding="utf-8",
    )
    prefix = str(tmp_path / "A.variant.ld_decay")
    ana_ld_decay(
        input_file=str(vcor),
        output_prefix=prefix,
        max_distance_kb=5000,
        bin_size_kb=50,
        scatter_max_points=10,
        chunksize=10,
        median_sample_per_bin=10,
        workers=1,
    )

    overall = pd.read_csv(f"{prefix}.info.tsv", sep="\t")
    by_chr = pd.read_csv(f"{prefix}.by_chr.info.tsv", sep="\t")
    assert Path(f"{prefix}.png").is_file()

    first = overall.loc[overall["bin_start_kb"] == 0].iloc[0]
    assert first["n_pairs"] == 2
    assert first["mean_r2"] == pytest.approx(0.5)

    chr1 = by_chr[(by_chr["chromosome"] == "chr1A") & (by_chr["bin_start_kb"] == 0)].iloc[0]
    chr2 = by_chr[(by_chr["chromosome"] == "chr2A") & (by_chr["bin_start_kb"] == 0)].iloc[0]
    assert chr1["mean_r2"] == pytest.approx(0.80)
    assert chr2["mean_r2"] == pytest.approx(0.20)
    assert list(by_chr["chromosome"].unique()) == ["chr1A", "chr2A"]


def _by_chr_rows(chrom: str, mean_r2: float, n_pairs: int = 10) -> str:
    return (
        f"{chrom}\t0.0\t50.0\t25.0\t{mean_r2}\t{n_pairs}\n"
        f"{chrom}\t50.0\t100.0\t75.0\t{mean_r2}\t{n_pairs}\n"
    )


def test_homoeolog_compare_ranks_by_weighted_mean_r2(tmp_path: Path):
    header = "chromosome\tbin_start_kb\tbin_end_kb\tbin_mid_kb\tmean_r2\tn_pairs\n"
    a = tmp_path / "A.by_chr.tsv"
    b = tmp_path / "B.by_chr.tsv"
    d = tmp_path / "D.by_chr.tsv"
    a.write_text(header + _by_chr_rows("chr1A", 0.80) + _by_chr_rows("chr2A", 0.20), encoding="utf-8")
    b.write_text(header + _by_chr_rows("chr1B", 0.70) + _by_chr_rows("chr2B", 0.30), encoding="utf-8")
    d.write_text(header + _by_chr_rows("chr1D", 0.10) + _by_chr_rows("chr2D", 0.40), encoding="utf-8")

    summary_a = summarize_ld_decay_chromosomes(pd.read_csv(a, sep="\t"))
    assert summary_a.loc[summary_a["chromosome"] == "chr1A", "mean_r2"].iloc[0] == pytest.approx(0.80)

    table = homoeolog_subgenome_ld_table(
        {
            "A": summary_a,
            "B": summarize_ld_decay_chromosomes(pd.read_csv(b, sep="\t")),
            "D": summarize_ld_decay_chromosomes(pd.read_csv(d, sep="\t")),
        }
    )
    row1 = table.loc[table["homoeolog"] == 1].iloc[0]
    row2 = table.loc[table["homoeolog"] == 2].iloc[0]
    assert row1["rank_A"] == 1
    assert row1["rank_B"] == 1
    assert row1["rank_D"] == 2
    assert row2["rank_D"] == 1

    prefix = str(tmp_path / "ABD.variant.ld_decay.homoeolog_compare")
    plot_ld_homoeolog_subgenome_compare(str(a), str(b), str(d), prefix)
    assert Path(f"{prefix}.png").is_file()
    assert Path(f"{prefix}.rank.png").is_file()
    published = pd.read_csv(f"{prefix}.info.tsv", sep="\t")
    assert list(published["homoeolog"]) == [1, 2]


def test_distance_mean_r2_to_threshold_interpolates_crossing():
    by_chr = pd.DataFrame(
        {
            "chromosome": ["chr1A", "chr1A", "chr1A", "chr2A"],
            "bin_start_kb": [0.0, 1.0, 2.0, 0.0],
            "bin_end_kb": [1.0, 2.0, 3.0, 1.0],
            "bin_mid_kb": [0.5, 1.5, 2.5, 0.5],
            "mean_r2": [0.80, 0.60, 0.40, 0.30],
            "n_pairs": [10, 10, 10, 8],
        }
    )
    out = distance_mean_r2_to_threshold(by_chr, threshold=0.5)
    by = out.set_index("chromosome")
    assert by.loc["chr1A", "decay_kb"] == pytest.approx(2.0)
    assert bool(by.loc["chr1A", "reached"])
    assert by.loc["chr1A", "first_bin_mean_r2"] == pytest.approx(0.80)
    assert by.loc["chr2A", "decay_kb"] == pytest.approx(0.5)
    assert bool(by.loc["chr2A", "reached"])
    assert by.loc["chr2A", "first_bin_mean_r2"] == pytest.approx(0.30)


def test_ana_ld_decay_mean_r2_threshold_writes_chr_summary(tmp_path: Path):
    vcor = tmp_path / "A.info.ld.vcor"
    vcor.write_text(
        VCOR_HEADER
        + "1\t1000\ts1a\t1\t1400\ts1b\t0.80\n"
        + "1\t1000\ts1a\t1\t2500\ts1c\t0.60\n"
        + "1\t1000\ts1a\t1\t3500\ts1d\t0.40\n"
        + "7\t1000\ts2a\t7\t1400\ts2b\t0.30\n",
        encoding="utf-8",
    )
    prefix = str(tmp_path / "A.variant.ld_decay.r2_0p5.1kb")
    summary = ana_ld_decay_mean_r2_threshold(
        input_file=str(vcor),
        output_prefix=prefix,
        r2_threshold=0.5,
        max_distance_kb=10,
        bin_size_kb=1,
        scatter_max_points=10,
        chunksize=10,
        workers=1,
    )
    by = summary.set_index("chromosome")
    assert by.loc["chr1A", "decay_kb"] == pytest.approx(2.0)
    assert by.loc["chr2A", "decay_kb"] == pytest.approx(0.5)
    assert Path(f"{prefix}.chr_summary.info.tsv").is_file()
    assert Path(f"{prefix}.by_chr.info.tsv").is_file()
    assert Path(f"{prefix}.png").is_file()

