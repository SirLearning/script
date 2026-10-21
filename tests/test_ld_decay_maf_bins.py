"""Within-bin MAF LD decay: pairs only when both sites share one MAF interval."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from genetics.genomics.variant.ld_decay import (
    MAF_LD_BIN_LABELS,
    ana_ld_decay_maf_bins,
    combine_ld_maf_bin_decay_tables,
    maf_ld_bin_index,
    plot_ld_maf_bin_decay_from_table,
    within_maf_bin_masks,
)


VCOR_HEADER = "#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2\n"
AFREQ_HEADER = "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n"


def test_maf_ld_bin_index_edges():
    maf = np.array([0.0, 0.005, 0.01, 0.011, 0.05, 0.051, 0.10, 0.15, 0.20, 0.21, 0.50, 0.51])
    idx = maf_ld_bin_index(maf)
    assert list(idx) == [-1, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, -1]


def test_within_maf_bin_masks_keep_same_interval_only():
    maf_a = np.array([0.005, 0.005, 0.08, 0.25, 0.25])
    maf_b = np.array([0.009, 0.04, 0.09, 0.40, 0.08])
    masks = within_maf_bin_masks(maf_a, maf_b)
    assert list(masks["(0, 0.01]"]) == [True, False, False, False, False]
    assert list(masks["(0.01, 0.05]"]) == [False, False, False, False, False]
    assert list(masks["(0.05, 0.1]"]) == [False, False, True, False, False]
    assert list(masks["(0.2, 0.5]"]) == [False, False, False, True, False]


def test_ana_ld_decay_maf_bins_uses_within_bin_pairs_only(tmp_path: Path):
    vcor = tmp_path / "A.info.ld.vcor"
    afreq = tmp_path / "A.info.afreq"
    vcor.write_text(
        VCOR_HEADER
        + "1\t1000\t1-1000\t1\t26000\t1-26000\t0.90\n"
        + "1\t1000\t1-1000\t1\t51000\t1-51000\t0.30\n"
        + "1\t51000\t1-51000\t1\t76000\t1-76000\t0.50\n"
        + "1\t80000\t1-80000\t1\t130000\t1-130000\t0.40\n",
        encoding="utf-8",
    )
    afreq.write_text(
        AFREQ_HEADER
        + "1\t1-1000\tA\tT\t0.005\t100\n"
        + "1\t1-26000\tA\tT\t0.008\t100\n"
        + "1\t1-51000\tA\tT\t0.20\t100\n"
        + "1\t1-76000\tA\tT\t0.30\t100\n"
        + "1\t1-80000\tA\tT\t0.08\t100\n"
        + "1\t1-130000\tA\tT\t0.09\t100\n",
        encoding="utf-8",
    )
    prefix = str(tmp_path / "A.variant.ld_decay.maf_bin")
    ana_ld_decay_maf_bins(
        input_file=str(vcor),
        afreq_file=str(afreq),
        output_prefix=prefix,
        max_distance_kb=5000,
        bin_size_kb=50,
        chunksize=10,
    )
    info = pd.read_csv(f"{prefix}.info.tsv", sep="\t")
    th = pd.read_csv(f"{prefix}.th.tsv", sep="\t").iloc[0].to_dict()
    by_class = info.set_index(["pair_class", "bin_mid_kb"])
    assert float(by_class.loc[("(0, 0.01]", 25.0), "mean_r2"]) == 0.90
    assert int(by_class.loc[("(0, 0.01]", 25.0), "n_pairs"]) == 1
    assert float(by_class.loc[("(0.05, 0.1]", 75.0), "mean_r2"]) == 0.40
    assert "(0.1, 0.2]" not in set(info["pair_class"])
    assert int(th["Pairs_(0, 0.01]"]) == 1
    assert int(th["Pairs_Dropped_Cross_Bin"]) == 2
    assert Path(f"{prefix}.png").is_file()


def test_plot_ld_maf_bin_decay_from_table_writes_png(tmp_path: Path):
    info = tmp_path / "A.variant.ld_decay.maf_bin.info.tsv"
    rows = ["pair_class\tbin_start_kb\tbin_end_kb\tbin_mid_kb\tmean_r2\tmedian_r2\tn_pairs"]
    for label in MAF_LD_BIN_LABELS:
        rows.append(f"{label}\t0.0\t50.0\t25.0\t0.40\t0.30\t10")
        rows.append(f"{label}\t50.0\t100.0\t75.0\t0.20\t0.10\t8")
    info.write_text("\n".join(rows) + "\n", encoding="utf-8")
    prefix = str(tmp_path / "A.variant.ld_decay.maf_bin")
    plot_ld_maf_bin_decay_from_table(str(info), prefix, max_distance_kb=5000)
    assert Path(f"{prefix}.png").is_file()


def test_plot_ld_maf_bin_decay_from_table_uses_r_squared_superscript(
    tmp_path: Path, monkeypatch
):
    captured: dict = {}

    def _capture_plot(**kwargs):
        captured.update(kwargs)
        Path(kwargs["filename"]).touch()

    monkeypatch.setattr(
        "genetics.genomics.variant.ld_decay.plot_multi_line_series",
        _capture_plot,
    )
    info = tmp_path / "A.variant.ld_decay.maf_bin.info.tsv"
    info.write_text(
        "pair_class\tbin_start_kb\tbin_end_kb\tbin_mid_kb\tmean_r2\tmedian_r2\tn_pairs\n"
        "(0, 0.01]\t0.0\t50.0\t25.0\t0.40\t0.30\t10\n",
        encoding="utf-8",
    )
    plot_ld_maf_bin_decay_from_table(
        str(info),
        str(tmp_path / "A.variant.ld_decay.maf_bin"),
        max_distance_kb=5000,
    )
    assert captured["y_label"] == r"LD ($r^2$)"
    labels = [spec["label"] for spec in captured["y_specs"]]
    assert labels
    assert all(r"$r^2$" in label for label in labels)


def test_combine_ld_maf_bin_decay_tables_weights_by_pair_count(tmp_path: Path):
    cols = "pair_class\tbin_start_kb\tbin_end_kb\tbin_mid_kb\tmean_r2\tmedian_r2\tn_pairs"
    a = tmp_path / "A.variant.ld_decay.maf_bin.info.tsv"
    b = tmp_path / "B.variant.ld_decay.maf_bin.info.tsv"
    d = tmp_path / "D.variant.ld_decay.maf_bin.info.tsv"
    a.write_text(
        cols
        + "\n(0, 0.01]\t0.0\t50.0\t25.0\t0.20\t0.10\t3\n"
        + "(0.2, 0.5]\t0.0\t50.0\t25.0\t0.80\t0.70\t1\n",
        encoding="utf-8",
    )
    b.write_text(
        cols + "\n(0, 0.01]\t0.0\t50.0\t25.0\t0.80\t0.70\t1\n",
        encoding="utf-8",
    )
    d.write_text(
        cols + "\n(0, 0.01]\t50.0\t100.0\t75.0\t0.10\t0.10\t4\n",
        encoding="utf-8",
    )
    prefix = str(tmp_path / "all.variant.ld_decay.maf_bin")
    grouped = combine_ld_maf_bin_decay_tables(
        [a, b, d],
        prefix,
        max_distance_kb=5000,
    )
    by_class = grouped.set_index(["pair_class", "bin_mid_kb"])
    assert float(by_class.loc[("(0, 0.01]", 25.0), "mean_r2"]) == pytest.approx(0.35)
    assert int(by_class.loc[("(0, 0.01]", 25.0), "n_pairs"]) == 4
    assert float(by_class.loc[("(0, 0.01]", 75.0), "mean_r2"]) == pytest.approx(0.10)
    assert float(by_class.loc[("(0.2, 0.5]", 25.0), "mean_r2"]) == pytest.approx(0.80)
    th = pd.read_csv(f"{prefix}.th.tsv", sep="\t").iloc[0].to_dict()
    assert int(th["N_Subgenome_Tables"]) == 3
    assert int(th["Pairs_(0, 0.01]"]) == 8
    assert Path(f"{prefix}.png").is_file()
