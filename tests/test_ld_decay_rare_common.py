"""Rare-rare vs rare-common LD pair classification and decay plots."""

from pathlib import Path

import pandas as pd
import pytest

from infra.utils.graph import SEQUENTIAL_WHITE_RED, log10_offset, plot_grouped_density_lines
from genetics.genomics.variant.ld_decay import (
    ALL_BIN_LABEL,
    DPRIME_CLASS_NEG1,
    DPRIME_CLASS_NEGATIVE,
    DPRIME_CLASS_POS1,
    DPRIME_CLASS_POSITIVE,
    DPRIME_CLASS_ZERO,
    PAIR_CLASS_LINE_ALPHA,
    PAIR_CLASS_DECAY_COLORS,
    PAIR_COMMON_COMMON,
    PAIR_RARE_COMMON,
    PAIR_RARE_RARE,
    _iter_plink2_vcor_chunks,
    _iter_plink2_vcor_chunks_polars,
    _reservoir_add,
    ana_ld_decay_rare_common,
    annotate_pairs_pairwise_maf,
    bin_signed_dprime_pplus_vs_q,
    bin_signed_dprime_windows,
    classify_signed_dprime,
    pair_class_masks,
    pair_type_masks,
    pair_haplotype_background_roles,
    pairwise_complete_allele_stats,
    pairwise_complete_ld_stats,
    plot_ld_maf_r2_scatter_from_table,
    plot_ld_dprime_decay_from_table,
    plot_ld_rare_common_decay_from_table,
    plot_signed_dprime_pplus_vs_q_overall,
    abs_dprime_from_maf_r2,
    pair_distance_kb_from_ids,
    r2_max_contour_branches,
    r2_max_from_maf,
    r2_max_scatter_line_overlays,
)


VCOR_HEADER = "#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2\n"
AFREQ_HEADER = "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n"


def test_pair_class_masks_keep_rare_pairs_and_drop_common_common():
    import numpy as np

    maf_a = np.array([0.01, 0.01, 0.20, 0.049, np.nan])
    maf_b = np.array([0.02, 0.20, 0.30, 0.05, 0.01])
    rare_rare, rare_common = pair_class_masks(maf_a, maf_b, maf_threshold=0.05)
    assert list(rare_rare) == [True, False, False, False, False]
    assert list(rare_common) == [False, True, False, True, False]
    types = pair_type_masks(maf_a, maf_b, maf_threshold=0.05)
    assert list(types[PAIR_RARE_RARE]) == [True, False, False, False, False]
    assert list(types[PAIR_RARE_COMMON]) == [False, True, False, True, False]
    assert list(types[PAIR_COMMON_COMMON]) == [False, False, True, False, False]


def test_ana_ld_decay_rare_common_splits_three_pair_classes(tmp_path: Path):
    vcor = tmp_path / "A.info.ld.vcor"
    afreq = tmp_path / "A.info.afreq"
    vcor.write_text(
        VCOR_HEADER
        + "1\t1000\t1-1000\t1\t26000\t1-26000\t0.90\n"
        + "1\t1000\t1-1000\t1\t2000000\t1-2000000\t0.10\n"
        + "1\t1000\t1-1000\t1\t51000\t1-51000\t0.30\n"
        + "1\t26000\t1-26000\t1\t51000\t1-51000\t0.25\n"
        + "1\t51000\t1-51000\t1\t76000\t1-76000\t0.99\n",
        encoding="utf-8",
    )
    afreq.write_text(
        AFREQ_HEADER
        + "1\t1-1000\tA\tT\t0.01\t100\n"
        + "1\t1-26000\tA\tT\t0.02\t100\n"
        + "1\t1-2000000\tA\tT\t0.01\t100\n"
        + "1\t1-51000\tA\tT\t0.20\t100\n"
        + "1\t1-76000\tA\tT\t0.30\t100\n",
        encoding="utf-8",
    )
    prefix = str(tmp_path / "A.variant.ld_decay.rare_common")
    ana_ld_decay_rare_common(
        input_file=str(vcor),
        afreq_file=str(afreq),
        output_prefix=prefix,
        maf_threshold=0.05,
        max_distance_kb=5000,
        bin_size_kb=50,
        scatter_max_points=100,
        chunksize=10,
        workers=1,
        random_seed=42,
    )

    info = pd.read_csv(f"{prefix}.info.tsv", sep="\t")
    th = pd.read_csv(f"{prefix}.th.tsv", sep="\t")
    th_map = th.iloc[0].to_dict()

    assert set(info["pair_class"]) == {PAIR_RARE_RARE, PAIR_RARE_COMMON, PAIR_COMMON_COMMON}
    rr = info[info["pair_class"] == PAIR_RARE_RARE].set_index("bin_mid_kb")
    rc = info[info["pair_class"] == PAIR_RARE_COMMON].set_index("bin_mid_kb")
    cc = info[info["pair_class"] == PAIR_COMMON_COMMON].set_index("bin_mid_kb")
    assert rr.loc[25.0, "mean_r2"] == 0.90
    assert rr.loc[1975.0, "mean_r2"] == 0.10
    assert int(rc.loc[25.0, "n_pairs"]) == 1
    assert int(rc.loc[75.0, "n_pairs"]) == 1
    assert float(cc.loc[25.0, "mean_r2"]) == 0.99
    assert int(th_map["Pairs_Rare_Rare"]) == 2
    assert int(th_map["Pairs_Rare_Common"]) == 2
    assert int(th_map["Pairs_Common_Common"]) == 1
    assert int(th_map["Maf_Scatter_Pairs"]) == 5
    assert int(th_map["Maf_Scatter_Pairs_R2_gt_0.2"]) == 4
    assert int(th_map["Pairs_R2_gt_R2max_GlobalMAF"]) >= 1

    scatter_tbl = pd.read_csv(f"{prefix}.maf_scatter.info.tsv", sep="\t")
    assert {"ID_A", "ID_B", "MAF_A", "MAF_B", "R2", "R2_MAX", "OBS_CT_A", "OBS_CT_B", "DISTANCE_KB"} <= set(
        scatter_tbl.columns
    )
    cc = scatter_tbl[(scatter_tbl["ID_A"] == "1-51000") & (scatter_tbl["ID_B"] == "1-76000")]
    assert len(cc) == 1
    assert float(cc["R2"].iloc[0]) == pytest.approx(0.99, rel=0, abs=1e-5)
    assert float(cc["R2_MAX"].iloc[0]) == pytest.approx(0.20 * 0.70 / (0.30 * 0.80), abs=1e-3)
    assert float(cc["R2"].iloc[0]) > float(cc["R2_MAX"].iloc[0])
    assert "DISTANCE_KB" in scatter_tbl.columns
    dprime_info = pd.read_csv(f"{prefix}.dprime.info.tsv", sep="\t")
    assert {"pair_class", "bin_mid_kb", "mean_dprime", "n_pairs"} <= set(dprime_info.columns)

    for suffix in (
        "",
        ".scatter",
        ".r2_dist",
        ".r2_dist.log",
        ".maf_scatter",
        ".maf_scatter_0.2",
        ".maf_scatter_rare",
        ".maf_scatter_0.2_rare",
        ".dprime",
        ".dprime.scatter",
        ".dprime.dist",
        ".dprime.dist.log",
        ".dprime.signed_prop",
        ".dprime.signed_boundary_stack",
        ".dprime.signed_interior_dist",
    ):
        assert Path(f"{prefix}{suffix}.png").is_file()


def test_plot_ld_rare_common_decay_from_table_writes_png(tmp_path: Path):
    info = tmp_path / "A.variant.ld_decay.rare_common.info.tsv"
    info.write_text(
        "pair_class\tbin_start_kb\tbin_end_kb\tbin_mid_kb\tmean_r2\tmedian_r2\tn_pairs\n"
        "rare-rare\t0.0\t50.0\t25.0\t0.40\t0.30\t10\n"
        "rare-common\t0.0\t50.0\t25.0\t0.40\t0.30\t10\n"
        "common-common\t0.0\t50.0\t25.0\t0.20\t0.10\t10\n"
        "rare-rare\t50.0\t100.0\t75.0\t0.20\t0.10\t8\n"
        "rare-common\t50.0\t100.0\t75.0\t0.21\t0.10\t8\n"
        "common-common\t50.0\t100.0\t75.0\t0.05\t0.00\t8\n",
        encoding="utf-8",
    )
    prefix = str(tmp_path / "A.variant.ld_decay.rare_common")
    plot_ld_rare_common_decay_from_table(str(info), prefix, max_distance_kb=5000)
    assert Path(f"{prefix}.png").is_file()
    assert PAIR_CLASS_LINE_ALPHA[PAIR_RARE_RARE] == 1.0
    assert PAIR_CLASS_LINE_ALPHA[PAIR_RARE_COMMON] < 1.0
    assert PAIR_CLASS_DECAY_COLORS[PAIR_RARE_COMMON] == "#ff4d4d"


def test_abs_dprime_from_maf_r2_complete_coupling_example():
    import numpy as np

    p = np.array([0.01])
    q = np.array([0.20])
    r2 = p * (1.0 - q) / (q * (1.0 - p))
    got = abs_dprime_from_maf_r2(p, q, r2)
    assert float(got[0]) == pytest.approx(1.0, abs=1e-6)
    rmax = r2_max_from_maf(p, q)
    assert float(r2[0]) == pytest.approx(float(rmax[0]), abs=1e-6)
    equal = abs_dprime_from_maf_r2(np.array([0.2]), np.array([0.2]), np.array([0.25]))
    assert float(equal[0]) == pytest.approx(0.5)
    clipped = abs_dprime_from_maf_r2(
        np.array([0.20]), np.array([0.30]), np.array([0.99])
    )
    assert float(clipped[0]) == pytest.approx(1.0)


def test_r2_only_dprime_wrong_for_complete_repulsion():
    import numpy as np

    p = np.array([0.01])
    q = np.array([0.01])
    r2 = p * q / ((1.0 - p) * (1.0 - q))
    got = abs_dprime_from_maf_r2(p, q, r2)
    assert float(got[0]) == pytest.approx(0.01 / 0.99, abs=1e-4)
    p2 = np.array([0.01])
    q2 = np.array([0.20])
    r2b = p2 * q2 / ((1.0 - p2) * (1.0 - q2))
    gotb = abs_dprime_from_maf_r2(p2, q2, r2b)
    assert float(gotb[0]) == pytest.approx(0.20 / 0.80, abs=1e-4)


def test_pairwise_complete_ld_stats_sign_aware_denominator():
    import numpy as np

    n = 10000

    def inbred(alt_idx):
        x = np.zeros(n, dtype=np.float64)
        x[list(alt_idx)] = 2.0
        return x

    maf_a, maf_b, n_ok, d, dp = pairwise_complete_ld_stats(
        inbred(range(0, 100)), inbred(range(100, 200))
    )
    assert n_ok == n
    assert maf_a == pytest.approx(0.01)
    assert maf_b == pytest.approx(0.01)
    assert d < 0
    assert dp == pytest.approx(-1.0, abs=1e-6)

    maf_a, maf_b, n_ok, d, dp = pairwise_complete_ld_stats(
        inbred(range(0, 100)), inbred(range(0, 2000))
    )
    assert maf_a == pytest.approx(0.01)
    assert maf_b == pytest.approx(0.20)
    assert d > 0
    assert dp == pytest.approx(1.0, abs=1e-6)

    maf_a, maf_b, n_ok, d, dp = pairwise_complete_ld_stats(
        inbred(range(0, 100)), inbred(range(100, 2100))
    )
    assert maf_a == pytest.approx(0.01)
    assert maf_b == pytest.approx(0.20)
    assert d < 0
    assert dp == pytest.approx(-1.0, abs=1e-6)

    rng = np.random.default_rng(0)
    da = (rng.random(n) < 0.01).astype(np.float64) * 2.0
    db = (rng.random(n) < 0.20).astype(np.float64) * 2.0
    _maf_a, _maf_b, _n_ok, d, _dp = pairwise_complete_ld_stats(da, db)
    assert abs(d) < 0.003


def test_log10_offset_maps_zero_to_zero():
    import numpy as np

    got = log10_offset(np.array([0.0, 9.0, 99.0]), offset=1.0)
    assert got[0] == 0.0
    assert got[1] == pytest.approx(1.0)
    assert got[2] == pytest.approx(2.0)
    assert np.isfinite(got).all()


def test_plot_grouped_density_lines_log_keeps_empty_bins(tmp_path: Path):
    import numpy as np

    df = pd.DataFrame(
        {
            "x": np.concatenate([np.zeros(80), np.ones(20)]),
            "hue": ["a"] * 100,
        }
    )
    out = tmp_path / "dens_log.png"
    plot_grouped_density_lines(
        data=df,
        x_col="x",
        hue_col="hue",
        hue_order=["a"],
        palette_map={"a": "#4c72b0"},
        title="density log10+1",
        filename=str(out),
        xlim=(0, 1),
        n_bins=40,
        density=True,
        log_scale=True,
    )
    assert out.is_file()


def test_plot_ld_dprime_writes_density_lines(tmp_path: Path):
    import numpy as np

    prefix = str(tmp_path / "x.variant.ld_decay.rare_common")
    rng = np.random.default_rng(0)
    n_rr, n_rc = 80, 40
    dprime = np.concatenate([rng.beta(8.0, 1.2, n_rr), rng.beta(2.0, 2.0, n_rc)])
    dprime = np.clip(dprime, 0.0, 1.0)
    df = pd.DataFrame(
        {
            "ID_A": ["1-1000"] * (n_rr + n_rc),
            "ID_B": [f"1-{2000 + 1000 * i}" for i in range(n_rr + n_rc)],
            "MAF_A": [0.01] * (n_rr + n_rc),
            "MAF_B": [0.02] * n_rr + [0.20] * n_rc,
            "R2": [0.10] * (n_rr + n_rc),
            "D_PRIME": dprime,
            "D_PRIME_SIGNED": np.concatenate(
                [np.full(40, 1.0), np.full(40, -1.0), np.full(20, 0.4), np.full(20, -0.4)]
            ),
        }
    )
    plot_ld_dprime_decay_from_table(df, prefix, max_distance_kb=5000, bin_size_kb=50)
    assert Path(f"{prefix}.dprime.dist.png").is_file()
    assert Path(f"{prefix}.dprime.dist.log.png").is_file()
    assert Path(f"{prefix}.dprime.signed_prop.png").is_file()
    assert Path(f"{prefix}.dprime.signed_boundary_stack.png").is_file()
    assert Path(f"{prefix}.dprime.signed_interior_dist.png").is_file()
    assert Path(f"{prefix}.dprime.signed.info.tsv").is_file()
    assert Path(f"{prefix}.dprime.signed_pplus_vs_q.png").is_file()
    assert Path(f"{prefix}.dprime.signed_pplus_vs_q.info.tsv").is_file()
    assert Path(f"{prefix}.dprime.signed_pplus_vs_q.rare_common.mac_dist.png").is_file()
    assert Path(f"{prefix}.dprime.signed_pplus_vs_q.rare_rare.mac_dist.png").is_file()


def test_classify_signed_dprime_five_classes():
    import numpy as np

    s = np.array([1.0, -1.0, 0.4, -0.3, 0.0, 1.0 - 1e-7, np.nan])
    got = classify_signed_dprime(s)
    assert list(got[:5]) == [
        DPRIME_CLASS_POS1,
        DPRIME_CLASS_NEG1,
        DPRIME_CLASS_POSITIVE,
        DPRIME_CLASS_NEGATIVE,
        DPRIME_CLASS_ZERO,
    ]
    assert got[5] == DPRIME_CLASS_POS1
    assert got[6] == ""


def test_bin_signed_dprime_windows_proportions():
    df = pd.DataFrame(
        {
            "MAF_A": [0.01] * 8,
            "MAF_B": [0.02, 0.02, 0.02, 0.20, 0.20, 0.20, 0.02, 0.02],
            "D_PRIME_SIGNED": [1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 0.4, -1.0],
            "DISTANCE_KB": [25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 75.0, 75.0],
        }
    )
    got = bin_signed_dprime_windows(
        df, maf_threshold=0.05, max_distance_kb=100, bin_size_kb=50
    )
    rr_near = got[(got["pair_class"] == PAIR_RARE_RARE) & (got["bin_mid_kb"] == 25.0)].iloc[0]
    rc_near = got[(got["pair_class"] == PAIR_RARE_COMMON) & (got["bin_mid_kb"] == 25.0)].iloc[0]
    rr_far = got[(got["pair_class"] == PAIR_RARE_RARE) & (got["bin_mid_kb"] == 75.0)].iloc[0]
    assert int(rr_near["n_all"]) == 3
    assert int(rr_near["n_pos1"]) == 2
    assert int(rr_near["n_neg1"]) == 1
    assert float(rr_near["p_plus1"]) == pytest.approx(2 / 3)
    assert float(rr_near["p_minus1"]) == pytest.approx(1 / 3)
    assert float(rr_near["p_plus1_given_boundary"]) == pytest.approx(2 / 3)
    assert float(rr_near["asymmetry"]) == pytest.approx((2 - 1) / 3)
    assert int(rc_near["n_pos1"]) == 1
    assert int(rc_near["n_neg1"]) == 2
    assert float(rc_near["p_plus1_given_boundary"]) == pytest.approx(1 / 3)
    assert int(rr_far["n_other"]) == 1
    assert int(rr_far["n_neg1"]) == 1
    assert float(rr_far["p_other"]) == pytest.approx(0.5)


def test_pair_distance_kb_from_ids_same_chrom_only():
    import numpy as np

    got = pair_distance_kb_from_ids(["19-196669359"], ["19-196719359"])
    assert float(got[0]) == pytest.approx(50.0)
    cross = pair_distance_kb_from_ids(["19-1"], ["20-1"])
    assert np.isnan(cross[0])


def test_plot_ld_maf_r2_scatter_writes_ultrarare_zooms(tmp_path: Path):
    prefix = str(tmp_path / "x.variant.ld_decay.rare_common")
    df = pd.DataFrame(
        {
            "MAF_A": [0.005, 0.004, 0.02],
            "MAF_B": [0.006, 0.003, 0.03],
            "R2": [0.90, 0.05, 0.50],
        }
    )
    plot_ld_maf_r2_scatter_from_table(df, prefix, pairwise_maf=False)
    assert Path(f"{prefix}.maf_scatter_ultrarare.png").is_file()
    assert Path(f"{prefix}.maf_scatter_0.2_ultrarare.png").is_file()


def test_sequential_white_red_cmap_is_colorless_to_red():
    rgba0 = SEQUENTIAL_WHITE_RED(0.0)
    rgba1 = SEQUENTIAL_WHITE_RED(1.0)
    assert rgba0[0] > 0.95 and rgba0[1] > 0.95 and rgba0[2] > 0.95
    assert rgba1[0] > rgba1[2]
    assert rgba1[1] < 0.35


def test_r2_max_from_maf_matches_allele_frequency_bound():
    import numpy as np

    p = np.array([0.01, 0.01, 0.01, 0.17, 0.20])
    q = np.array([0.01, 0.05, 0.20, 0.47, 0.30])
    got = r2_max_from_maf(p, q)
    want = p * (1.0 - q) / (q * (1.0 - p))
    assert got == pytest.approx(want)
    assert float(got[2]) == pytest.approx(0.040404, abs=1e-5)
    swapped = r2_max_from_maf(q, p)
    assert swapped == pytest.approx(got)
    assert float(r2_max_from_maf(np.array([0.0]), np.array([0.0]))[0]) == 1.0


def test_r2_max_contour_branches_match_formula_and_mirror_across_diagonal():
    import numpy as np

    r2 = 0.2
    p, q = r2_max_contour_branches(r2, n=50)
    assert np.all(p <= q)
    assert float(q[-1]) == pytest.approx(0.5, abs=1e-6)
    assert float(p[-1]) == pytest.approx(r2 / (1.0 + r2), abs=1e-6)
    got = r2_max_from_maf(p, q)
    assert got == pytest.approx(np.full_like(got, r2), abs=1e-5)
    overlays = r2_max_scatter_line_overlays((0.8, 0.6, 0.2))
    assert len(overlays) == 6
    labels = [spec.get("label") for spec in overlays]
    assert labels.count(None) == 3
    first = overlays[0]
    mirror = overlays[1]
    assert first["x"] == pytest.approx(mirror["y"])
    assert first["y"] == pytest.approx(mirror["x"])


def test_pairwise_complete_allele_stats_use_intersection_not_global_called():
    import numpy as np

    # Site A called in 6/8 samples (MAF ~0.5); site B missing exactly where A is 0.
    dosage_a = np.array([0.0, 0.0, 0.0, 2.0, 2.0, 2.0, np.nan, np.nan])
    dosage_b = np.array([np.nan, np.nan, np.nan, 2.0, 2.0, 2.0, 0.0, 0.0])
    maf_a_global = min(np.nanmean(dosage_a) / 2.0, 1.0 - np.nanmean(dosage_a) / 2.0)
    maf_b_global = min(np.nanmean(dosage_b) / 2.0, 1.0 - np.nanmean(dosage_b) / 2.0)
    assert maf_a_global == pytest.approx(0.5)
    assert maf_b_global == pytest.approx(0.4)
    maf_a, maf_b, n = pairwise_complete_allele_stats(dosage_a, dosage_b)
    assert n == 3
    assert maf_a == pytest.approx(0.0)
    assert maf_b == pytest.approx(0.0)


def test_annotate_pairs_pairwise_maf_adds_intersection_columns():
    import numpy as np

    pairs = pd.DataFrame(
        {
            "ID_A": ["s1", "s1"],
            "ID_B": ["s2", "s3"],
            "MAF_A": [0.50, 0.50],
            "MAF_B": [0.40, 0.10],
            "R2": [0.95, 0.10],
            "R2_MAX": [0.67, 0.11],
        }
    )
    dosage_by_id = {
        "s1": np.array([0.0, 0.0, 0.0, 2.0, 2.0, 2.0, np.nan, np.nan]),
        "s2": np.array([np.nan, np.nan, np.nan, 2.0, 2.0, 2.0, 0.0, 0.0]),
        "s3": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0]),
    }
    out = annotate_pairs_pairwise_maf(pairs, dosage_by_id)
    assert float(out.loc[0, "MAF_A_GLOBAL"]) == pytest.approx(0.50)
    assert float(out.loc[0, "MAF_A"]) == pytest.approx(0.0)
    assert float(out.loc[0, "MAF_B"]) == pytest.approx(0.0)
    assert int(out.loc[0, "N_PAIRWISE"]) == 3
    assert float(out.loc[0, "R2_MAX"]) == pytest.approx(1.0)
    assert int(out.loc[1, "N_PAIRWISE"]) == 6


def test_iter_plink2_vcor_chunks_reads_tab_header_not_hash_comment(tmp_path: Path):
    vcor = tmp_path / "tiny.vcor"
    vcor.write_text(
        VCOR_HEADER + "1\t1000\t1-1000\t1\t26000\t1-26000\t0.90\n",
        encoding="utf-8",
    )
    chunks = list(_iter_plink2_vcor_chunks(str(vcor), chunksize=10))
    assert len(chunks) == 1
    df = chunks[0]
    assert list(df.columns)[:3] == ["#CHROM_A", "POS_A", "ID_A"]
    assert float(df["UNPHASED_R2"].iloc[0]) == pytest.approx(0.90)


def test_iter_plink2_vcor_chunks_polars_byte_windows_cover_all_rows(tmp_path: Path):
    rows = [
        f"1\t{1000 + i}\t1-{1000 + i}\t1\t{2000 + i}\t1-{2000 + i}\t0.{i}\n"
        for i in range(5)
    ]
    vcor = tmp_path / "chunk.vcor"
    vcor.write_text(VCOR_HEADER + "".join(rows), encoding="utf-8")
    chunks = list(_iter_plink2_vcor_chunks_polars(str(vcor), chunksize=2))
    assert len(chunks) >= 2
    df = pd.concat(chunks, ignore_index=True)
    assert len(df) == 5
    assert list(df.columns)[:3] == ["#CHROM_A", "POS_A", "ID_A"]
    assert list(df["POS_A"].astype(int)) == [1000, 1001, 1002, 1003, 1004]
    assert float(df["UNPHASED_R2"].iloc[-1]) == pytest.approx(0.4)
    subset = list(
        _iter_plink2_vcor_chunks_polars(
            str(vcor),
            chunksize=2,
            usecols=["POS_A", "POS_B", "UNPHASED_R2"],
        )
    )
    sub = pd.concat(subset, ignore_index=True)
    assert list(sub.columns) == ["POS_A", "POS_B", "UNPHASED_R2"]
    assert len(sub) == 5


def test_reservoir_add_caps_without_concatenating_all_pairs():
    import numpy as np

    rng = np.random.default_rng(0)
    reservoir = np.empty((0, 3), dtype=np.float32)
    n_seen = 0
    for _ in range(5):
        incoming = rng.random((10000, 3)).astype(np.float32)
        reservoir, n_seen = _reservoir_add(reservoir, incoming, n_seen, 200, rng)
    assert n_seen == 50000
    assert len(reservoir) == 200


def test_pair_haplotype_background_roles_rare_common_uses_common_maf():
    import numpy as np

    maf_a = np.array([0.01, 0.20, 0.01])
    maf_b = np.array([0.25, 0.02, 0.03])
    n_pair = np.array([100.0, 100.0, 50.0])
    roles = pair_haplotype_background_roles(
        maf_a, maf_b, maf_threshold=0.05, n_pairwise=n_pair
    )
    assert list(roles["pair_class"]) == [
        PAIR_RARE_COMMON,
        PAIR_RARE_COMMON,
        PAIR_RARE_RARE,
    ]
    assert float(roles.loc[0, "background_maf"]) == pytest.approx(0.25)
    assert float(roles.loc[0, "rare_maf"]) == pytest.approx(0.01)
    assert float(roles.loc[0, "rare_mac"]) == pytest.approx(0.01 * 2 * 100)
    assert float(roles.loc[1, "background_maf"]) == pytest.approx(0.20)
    assert float(roles.loc[2, "background_maf"]) == pytest.approx(0.03)
    assert float(roles.loc[2, "rare_maf"]) == pytest.approx(0.01)
    assert float(roles.loc[2, "rare_mac"]) == pytest.approx(0.01 * 2 * 50)


def test_bin_pplus_vs_q_matches_common_maf_on_boundary():
    rows = []
    for _ in range(1):
        rows.append({"MAF_A": 0.01, "MAF_B": 0.10, "D_PRIME_SIGNED": 1.0, "DISTANCE_KB": 10.0, "N_PAIRWISE": 200})
    for _ in range(9):
        rows.append({"MAF_A": 0.01, "MAF_B": 0.10, "D_PRIME_SIGNED": -1.0, "DISTANCE_KB": 10.0, "N_PAIRWISE": 200})
    for _ in range(4):
        rows.append({"MAF_A": 0.01, "MAF_B": 0.40, "D_PRIME_SIGNED": 1.0, "DISTANCE_KB": 10.0, "N_PAIRWISE": 200})
    for _ in range(6):
        rows.append({"MAF_A": 0.01, "MAF_B": 0.40, "D_PRIME_SIGNED": -1.0, "DISTANCE_KB": 10.0, "N_PAIRWISE": 200})
    got = bin_signed_dprime_pplus_vs_q(pd.DataFrame(rows), maf_threshold=0.05)
    overall = got[
        (got["pair_class"] == PAIR_RARE_COMMON)
        & (got["mac_bin"] == ALL_BIN_LABEL)
        & (got["dist_bin"] == ALL_BIN_LABEL)
    ]
    by_q = overall.set_index(overall["q_mean"].round(2))
    assert float(by_q.loc[0.10, "p_plus"]) == pytest.approx(0.10)
    assert float(by_q.loc[0.40, "p_plus"]) == pytest.approx(0.40)
    assert int(by_q.loc[0.10, "n_bound"]) == 10
    assert int(by_q.loc[0.40, "n_bound"]) == 10


def test_bin_pplus_vs_q_drops_pairs_beyond_max_distance():
    rows = [
        {"MAF_A": 0.01, "MAF_B": 0.20, "D_PRIME_SIGNED": 1.0, "DISTANCE_KB": 500.0, "N_PAIRWISE": 200},
        {"MAF_A": 0.01, "MAF_B": 0.20, "D_PRIME_SIGNED": -1.0, "DISTANCE_KB": 500.0, "N_PAIRWISE": 200},
        {"MAF_A": 0.01, "MAF_B": 0.20, "D_PRIME_SIGNED": 1.0, "DISTANCE_KB": 2000.0, "N_PAIRWISE": 200},
        {"MAF_A": 0.01, "MAF_B": 0.20, "D_PRIME_SIGNED": -1.0, "DISTANCE_KB": 2000.0, "N_PAIRWISE": 200},
    ]
    got = bin_signed_dprime_pplus_vs_q(
        pd.DataFrame(rows), maf_threshold=0.05, max_distance_kb=1000
    )
    overall = got[
        (got["pair_class"] == PAIR_RARE_COMMON)
        & (got["mac_bin"] == ALL_BIN_LABEL)
        & (got["dist_bin"] == ALL_BIN_LABEL)
    ]
    assert int(overall["n_bound"].sum()) == 2


def test_plot_signed_dprime_pplus_vs_q_overall_1mb(tmp_path: Path):
    prefix = str(tmp_path / "x.variant.ld_decay.rare_common")
    df = pd.DataFrame(
        {
            "MAF_A": [0.01, 0.01, 0.01, 0.01],
            "MAF_B": [0.20, 0.20, 0.20, 0.20],
            "D_PRIME_SIGNED": [1.0, -1.0, 1.0, -1.0],
            "DISTANCE_KB": [100.0, 200.0, 2500.0, 3000.0],
            "N_PAIRWISE": [200, 200, 200, 200],
        }
    )
    plot_signed_dprime_pplus_vs_q_overall(
        df, prefix, maf_threshold=0.05, max_distance_kb=1000, filename_suffix="1mb"
    )
    assert Path(f"{prefix}.dprime.signed_pplus_vs_q.1mb.png").is_file()
    assert Path(f"{prefix}.dprime.signed_pplus_vs_q.1mb.info.tsv").is_file()
    assert not Path(f"{prefix}.dprime.signed_pplus_vs_q.png").is_file()
    info = pd.read_csv(f"{prefix}.dprime.signed_pplus_vs_q.1mb.info.tsv", sep="\t")
    overall = info[
        (info["mac_bin"] == ALL_BIN_LABEL) & (info["dist_bin"] == ALL_BIN_LABEL)
    ]
    assert int(overall["n_bound"].sum()) == 2
