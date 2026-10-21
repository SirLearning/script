"""Trusted-cell selection on a RelativeDepth mean-vs-SD grid."""

from pathlib import Path

import gzip
import numpy as np
import pandas as pd
import pytest

from genetics.genomics.variant.popdep_density_filter import (
    DENSITY_DROP_LABEL,
    DENSITY_KEEP_LABEL,
    KEEP_FRAC_LOOSE,
    KEEP_FRAC_STRICT,
    KEEP_FRAC_WIDE,
    LINEAR_RELDEPTH_PLOT_RANGE,
    LINEAR_SITE_KEEP_FRACS,
    LINEAR_SITE_RANGE,
    LOG_N_CELLS,
    MEAN_RANGE_0_6,
    CV_RANGE_0_25,
    CV_RELIABLE_BINS,
    CV_RELIABLE_KEEP_FRACS,
    CV_RELIABLE_N_CELLS,
    accumulate_rel_mean_sd_histogram,
    accumulate_rel_mean_cv_histogram_grid,
    ana_test_popdep_cv_density,
    ana_test_popdep_density_keep_plots,
    ana_test_popdep_linear_density,
    ana_test_popdep_log_density,
    assign_bin_indices,
    concat_gzip_beds,
    histogram_log_log,
    histogram_rel_mean_sd,
    keep_frac_pct,
    label_density_filter,
    load_histogram_npz,
    plot_density_keep_drop_scatters,
    reliable_cv_hist_ranges,
    site_coverage_curve,
    sites_in_trusted_cells,
    square_bins_for_n_cells,
    trusted_mask_from_counts,
    trusted_mask_from_site_coverage,
    write_chr_cv_reliable_beds,
    write_chr_mean_cv_histogram,
)


def test_histogram_counts_points_inside_range():
    x = np.array([0.5, 0.5, 4.9, 9.0])
    y = np.array([0.5, 0.6, 4.9, 0.5])
    counts, n_valid, n_oob = histogram_rel_mean_sd(
        x, y, x_range=(0.0, 5.0), y_range=(0.0, 5.0), bins_x=5, bins_y=5
    )
    assert counts.shape == (5, 5)
    assert int(counts.sum()) == 3
    assert n_valid == 3
    assert n_oob == 1


def test_trusted_mask_keeps_top_fraction_of_nonempty_cells():
    counts = np.array(
        [
            [10, 9, 8, 0],
            [7, 6, 5, 0],
            [4, 3, 2, 0],
            [1, 0, 0, 0],
        ],
        dtype=np.int64,
    )
    mask_80 = trusted_mask_from_counts(counts, keep_frac=0.8)
    mask_60 = trusted_mask_from_counts(counts, keep_frac=0.6)
    assert int(mask_60.sum()) <= int(mask_80.sum())
    assert counts[mask_60].min() >= counts[mask_80].min()
    assert not bool(mask_80[counts == 0].any())


def test_site_coverage_mask_stops_at_target_fraction_of_sites():
    counts = np.array([[100, 50, 30, 0], [10, 5, 4, 1]], dtype=np.int64)
    mask_80, info_80 = trusted_mask_from_site_coverage(counts, site_frac=0.8)
    mask_60, info_60 = trusted_mask_from_site_coverage(counts, site_frac=0.6)
    mask_90, _info_90 = trusted_mask_from_site_coverage(counts, site_frac=0.9)
    total = int(counts.sum())
    assert int(counts[mask_80].sum()) >= int(np.ceil(0.8 * total))
    assert int(counts[mask_60].sum()) >= int(np.ceil(0.6 * total))
    assert int(counts[mask_90].sum()) >= int(np.ceil(0.9 * total))
    assert info_80["n_cells_kept"] == 3
    assert info_60["n_cells_kept"] == 2
    assert info_80["n_nonempty"] == 7
    assert info_80["cell_top_frac"] == pytest.approx(3 / 7)
    assert not bool(mask_80[counts == 0].any())
    assert int(mask_60.sum()) < int(mask_80.sum())
    assert int(mask_80.sum()) <= int(mask_90.sum())


def test_site_coverage_curve_tracks_densest_cell_prefix():
    counts = np.array([[100, 50, 30, 0], [10, 5, 4, 1]], dtype=np.int64)
    curve = site_coverage_curve(counts)
    assert float(curve.iloc[0]["site_frac_kept"]) == 0.0
    assert float(curve.iloc[0]["cell_top_frac"]) == 0.0
    assert float(curve.iloc[-1]["site_frac_kept"]) == pytest.approx(1.0)
    assert float(curve.iloc[-1]["cell_top_frac"]) == pytest.approx(1.0)
    _, info_80 = trusted_mask_from_site_coverage(counts, site_frac=0.8)
    row = curve.loc[curve["n_cells_kept"] == info_80["n_cells_kept"]].iloc[0]
    assert int(row["n_sites_kept"]) == info_80["n_sites_kept"]
    assert float(row["cell_top_frac"]) == pytest.approx(info_80["cell_top_frac"])


def test_sites_in_trusted_cells_follow_the_mask():
    mask = np.array([[True, False], [False, False]])
    sites = pd.DataFrame(
        {
            "ID": ["dense", "sparse", "oob"],
            "RelativeDepth_Mean": [0.25, 0.75, 9.0],
            "RelativeDepth_SD": [0.25, 0.75, 0.25],
        }
    )
    kept = sites_in_trusted_cells(
        sites, mask, x_range=(0.0, 1.0), y_range=(0.0, 1.0)
    )
    assert list(kept["ID"]) == ["dense"]


def test_accumulate_histogram_from_tsv(tmp_path: Path):
    table = tmp_path / "chr000.popdep.txt"
    table.write_text(
        "Chrom\tPosition\tDepth_Mean\tDepth_SD\tRelativeDepth_Mean\tRelativeDepth_SD\n"
        "0\t1\t10\t1\t0.5\t0.4\n"
        "0\t2\t10\t1\t0.5\t0.4\n"
        "0\t3\t10\t1\t9.0\t0.4\n"
        "0\t4\t10\t1\tnan\t0.4\n"
    )
    counts, stats = accumulate_rel_mean_sd_histogram(
        table, x_range=(0.0, 5.0), y_range=(0.0, 5.0), bins_x=5, bins_y=5
    )
    assert int(stats["n_rows"]) == 4
    assert int(stats["n_valid"]) == 2
    assert int(stats["n_oob"]) == 1
    assert int(stats["n_nan"]) == 1
    assert int(counts.sum()) == 2
    i, j = assign_bin_indices(
        np.array([0.5]), np.array([0.4]), x_range=(0.0, 5.0), y_range=(0.0, 5.0), bins_x=5, bins_y=5
    )
    assert counts[i[0], j[0]] == 2
    assert KEEP_FRAC_LOOSE == 0.8
    assert KEEP_FRAC_STRICT == 0.6
    assert KEEP_FRAC_WIDE == 0.9


def _mean_sd_sites() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Depth_Mean": [10.0, 40.0, 8.0],
            "Depth_SD": [2.0, 20.0, 1.5],
            "RelativeDepth_Mean": [0.25, 0.75, 9.0],
            "RelativeDepth_SD": [0.25, 0.75, 0.25],
        }
    )


def test_label_density_filter_marks_oob_as_dropped():
    mask = np.array([[True, False], [False, False]])
    labeled = label_density_filter(
        _mean_sd_sites(), mask, x_range=(0.0, 1.0), y_range=(0.0, 1.0)
    )
    assert list(labeled["density_filter"]) == [
        DENSITY_KEEP_LABEL,
        DENSITY_DROP_LABEL,
        DENSITY_DROP_LABEL,
    ]


def test_ana_test_popdep_density_keep_plots_writes_mean_sd_pngs(tmp_path: Path):
    mask = np.array([[True, False], [False, False]])
    masks = tmp_path / "masks.npz"
    np.savez_compressed(masks, mask_60=mask, mask_80=mask)
    sites = tmp_path / "sites.tsv"
    _mean_sd_sites().to_csv(sites, sep="\t", index=False)
    prefix = tmp_path / "test_popdep.popdep.density"
    ana_test_popdep_density_keep_plots(
        sites,
        masks,
        prefix,
        x_range=(0.0, 1.0),
        y_range=(0.0, 1.0),
    )
    for pct in (60, 80):
        assert (tmp_path / f"test_popdep.popdep.density.keep{pct}.counts.tsv").exists()
        for pair in (
            "Depth_SD_vs_Depth_Mean",
            "RelativeDepth_SD_vs_RelativeDepth_Mean",
        ):
            assert (tmp_path / f"test_popdep.popdep.density.keep{pct}.scatter_{pair}.png").exists()
            assert (tmp_path / f"test_popdep.popdep.density.keep{pct}.scatter_{pair}.log.png").exists()


def test_square_bins_for_90000_cells_is_300_by_300():
    assert LOG_N_CELLS == 90000
    assert square_bins_for_n_cells(90000) == (300, 300)


def test_log_log_bins_split_low_relmean_that_share_a_linear_cell():
    x = np.array([0.01, 0.9])
    y = np.array([0.01, 0.9])
    linear, _, _ = histogram_rel_mean_sd(
        x, y, x_range=(0.0, 5.0), y_range=(0.0, 5.0), bins_x=5, bins_y=5
    )
    assert int((linear > 0).sum()) == 1
    log_counts, n_valid, n_oob = histogram_log_log(
        x, y, bins_x=8, bins_y=8
    )
    assert n_valid == 2
    assert n_oob == 0
    assert int((log_counts > 0).sum()) == 2


def test_ana_test_popdep_log_density_writes_heatmap_and_reldepth_scatters(tmp_path: Path):
    rng = np.random.default_rng(0)
    n = 80
    sites = pd.DataFrame(
        {
            "RelativeDepth_Mean": np.concatenate(
                [rng.uniform(0.02, 0.08, 20), rng.uniform(0.8, 1.4, n - 20)]
            ),
            "RelativeDepth_SD": np.concatenate(
                [rng.uniform(0.02, 0.08, 20), rng.uniform(0.4, 1.2, n - 20)]
            ),
            "Depth_Mean": rng.uniform(5.0, 20.0, n),
            "Depth_SD": rng.uniform(1.0, 8.0, n),
        }
    )
    table = tmp_path / "chr000.popdep.sites.tsv"
    sites.to_csv(table, sep="\t", index=False)
    prefix = tmp_path / "test_popdep.popdep.log_density"
    ana_test_popdep_log_density(table, prefix, n_cells=16)
    assert (tmp_path / "test_popdep.popdep.log_density.rel_mean_sd.log.density.png").exists()
    assert (tmp_path / "test_popdep.popdep.log_density.site_coverage_curve.tsv").exists()
    assert (tmp_path / "test_popdep.popdep.log_density.site_coverage_curve.png").exists()
    assert (tmp_path / "test_popdep.popdep.log_density.site_coverage_curve.log.png").exists()
    assert (tmp_path / "test_popdep.popdep.log_density.keep60.counts.tsv").exists()
    assert (
        tmp_path
        / "test_popdep.popdep.log_density.keep80.scatter_RelativeDepth_SD_vs_RelativeDepth_Mean.log.png"
    ).exists()
    assert (
        tmp_path
        / "test_popdep.popdep.log_density.keep80.scatter_RelativeDepth_SD_vs_RelativeDepth_Mean.png"
    ).exists()
    assert (
        tmp_path
        / "test_popdep.popdep.log_density.keep90.scatter_RelativeDepth_SD_vs_RelativeDepth_Mean.png"
    ).exists()
    assert (tmp_path / "test_popdep.popdep.log_density.keep90.counts.tsv").exists()
    assert LINEAR_RELDEPTH_PLOT_RANGE == (0.0, 10.0)


def test_linear_reldepth_keep_scatter_clips_to_0_10(tmp_path, monkeypatch):
    captured = []

    def fake_plot_joint_regression(df, *args, **kwargs):
        captured.append({"df": df, **kwargs})
        Path(str(kwargs["filename"])).touch()

    monkeypatch.setattr(
        "genetics.genomics.variant.popdep_density_filter.plot_joint_regression",
        fake_plot_joint_regression,
    )
    mask = np.array([[True, False], [False, False]])
    sites = pd.concat(
        [
            _mean_sd_sites(),
            pd.DataFrame(
                {
                    "Depth_Mean": [80.0],
                    "Depth_SD": [4.0],
                    "RelativeDepth_Mean": [20.0],
                    "RelativeDepth_SD": [1.0],
                }
            ),
        ],
        ignore_index=True,
    )
    plot_density_keep_drop_scatters(
        sites,
        mask,
        tmp_path / "out",
        90,
        x_range=(0.0, 1.0),
        y_range=(0.0, 1.0),
        scatter_pairs=(
            (
                "RelativeDepth_Mean",
                "RelativeDepth_SD",
                "Relative mean depth",
                "Relative depth SD",
            ),
        ),
        linear_x_lim=LINEAR_RELDEPTH_PLOT_RANGE,
        linear_y_lim=LINEAR_RELDEPTH_PLOT_RANGE,
    )
    linear = [row for row in captured if not row.get("log_x")]
    logged = [row for row in captured if row.get("log_x")]
    assert linear
    assert logged
    for row in linear:
        assert row["x_lim"] == (0.0, 10.0)
        assert row["y_lim"] == (0.0, 10.0)
        assert float(row["df"]["RelativeDepth_Mean"].max()) <= 10.0
        assert float(row["df"]["RelativeDepth_SD"].max()) <= 10.0
        assert not bool((row["df"]["RelativeDepth_Mean"] > 10.0).any())
    for row in logged:
        assert row.get("x_lim") is None
        assert row.get("y_lim") is None
        assert float(row["df"]["RelativeDepth_Mean"].max()) == pytest.approx(20.0)


def test_ana_test_popdep_linear_density_keeps_sites_inside_0p02_to_6(tmp_path):
    rng = np.random.default_rng(1)
    n_in = 70
    sites = pd.DataFrame(
        {
            "RelativeDepth_Mean": np.concatenate(
                [
                    np.full(5, 0.005),
                    rng.uniform(0.8, 1.4, n_in),
                    np.full(5, 7.0),
                ]
            ),
            "RelativeDepth_SD": np.concatenate(
                [
                    np.full(5, 1.0),
                    rng.uniform(0.4, 1.2, n_in),
                    np.full(5, 1.0),
                ]
            ),
            "Depth_Mean": rng.uniform(5.0, 20.0, n_in + 10),
            "Depth_SD": rng.uniform(1.0, 8.0, n_in + 10),
        }
    )
    table = tmp_path / "chr000.popdep.sites.tsv"
    sites.to_csv(table, sep="\t", index=False)
    prefix = tmp_path / "test_popdep.popdep.lin_density"
    masks = ana_test_popdep_linear_density(table, prefix, n_cells=16)
    assert LINEAR_SITE_RANGE == (0.02, 6.0)
    assert (tmp_path / "test_popdep.popdep.lin_density.rel_mean_sd.density.png").exists()
    assert (tmp_path / "test_popdep.popdep.lin_density.site_coverage_curve.tsv").exists()
    packed = np.load(tmp_path / "test_popdep.popdep.lin_density.masks.npz")
    assert "mask_90" in packed.files
    assert "mask_99" in packed.files
    pair = "RelativeDepth_SD_vs_RelativeDepth_Mean"
    assert LINEAR_SITE_KEEP_FRACS == (0.6, 0.8, 0.9, 0.95, 0.97, 0.98, 0.99)
    for frac in LINEAR_SITE_KEEP_FRACS:
        pct = keep_frac_pct(frac)
        assert (
            tmp_path / f"test_popdep.popdep.lin_density.keep{pct}.scatter_{pair}.png"
        ).exists()
        assert not (
            tmp_path / f"test_popdep.popdep.lin_density.keep{pct}.scatter_{pair}.log.png"
        ).exists()
    labeled = label_density_filter(
        sites,
        masks["mask_90"],
        x_range=tuple(packed["x_range"]),
        y_range=tuple(packed["y_range"]),
    )
    oob = labeled["RelativeDepth_Mean"].isin([0.005, 7.0])
    assert set(labeled.loc[oob, "density_filter"]) == {DENSITY_DROP_LABEL}
    assert DENSITY_KEEP_LABEL in set(labeled.loc[~oob, "density_filter"])
    th = pd.read_csv(tmp_path / "test_popdep.popdep.lin_density.th.tsv", sep="\t")
    assert float(th["n_oob"].iloc[0]) == 10
    assert not (
        tmp_path / "test_popdep.popdep.lin_density.rel_mean_sd.log.density.png"
    ).exists()


def test_ana_test_popdep_cv_density_keeps_mean_0_to_6(tmp_path):
    rng = np.random.default_rng(2)
    n_in = 80
    sites = pd.DataFrame(
        {
            "RelativeDepth_Mean": np.concatenate(
                [
                    np.zeros(4),
                    rng.uniform(0.3, 1.4, n_in),
                    np.full(4, 7.0),
                ]
            ),
            "RelativeDepth_SD": np.concatenate(
                [
                    np.full(4, 0.1),
                    rng.uniform(0.2, 0.6, n_in),
                    np.full(4, 1.0),
                ]
            ),
            "region_class": np.concatenate(
                [
                    np.full(4, "gap"),
                    rng.choice(["repeat", "normal", "gc"], n_in),
                    np.full(4, "cnv"),
                ]
            ),
        }
    )
    table = tmp_path / "chr000.popdep.sites.tsv"
    sites.to_csv(table, sep="\t", index=False)
    prefix = tmp_path / "test_popdep.popdep.cv_density"
    masks = ana_test_popdep_cv_density(table, prefix, n_cells=16)
    assert MEAN_RANGE_0_6 == (0.0, 6.0)
    assert (tmp_path / "test_popdep.popdep.cv_density.rel_mean_cv.density.png").exists()
    assert (
        tmp_path / "test_popdep.popdep.cv_density.reldepth.cv_vs_mean_by_class.png"
    ).exists()
    packed = np.load(tmp_path / "test_popdep.popdep.cv_density.masks.npz")
    assert "mask_70" in packed.files
    assert "mask_85" in packed.files
    assert "mask_99" in packed.files
    pair = "RelativeDepth_CV_vs_RelativeDepth_Mean"
    for pct in (60, 70, 80, 85, 90, 95, 97, 98, 99):
        assert (
            tmp_path / f"test_popdep.popdep.cv_density.keep{pct}.scatter_{pair}.png"
        ).exists()
        assert not (
            tmp_path / f"test_popdep.popdep.cv_density.keep{pct}.scatter_{pair}.log.png"
        ).exists()
    labeled = label_density_filter(
        sites.assign(
            RelativeDepth_CV=np.where(
                sites["RelativeDepth_Mean"] > 0,
                sites["RelativeDepth_SD"] / sites["RelativeDepth_Mean"],
                np.nan,
            )
        ),
        masks["mask_90"],
        x_col="RelativeDepth_Mean",
        y_col="RelativeDepth_CV",
        x_range=tuple(packed["x_range"]),
        y_range=tuple(packed["y_range"]),
    )
    oob = labeled["RelativeDepth_Mean"] >= 7.0
    assert set(labeled.loc[oob, "density_filter"]) == {DENSITY_DROP_LABEL}
    assert DENSITY_KEEP_LABEL in set(labeled["density_filter"])


def _tiger_grid_text() -> str:
    return (
        "Chrom\tPosition\tDepth_Mean\tDepth_SD\tRelativeDepth_Mean\tRelativeDepth_SD\n"
        "0\t1\t10\t1\t1.0\t0.5\n"
        "0\t2\t10\t1\t1.0\t0.5\n"
        "0\t3\t10\t1\t1.0\t0.5\n"
        "0\t4\t10\t1\t5.0\t20.0\n"
        "0\t5\t10\t1\t7.0\t1.0\n"
        "0\t6\t10\t1\t1.0\t30.0\n"
    )


def test_n500_cv_reliable_grid_is_500_by_500_mean_0_6_cv_0_25():
    assert MEAN_RANGE_0_6 == (0.0, 6.0)
    assert CV_RANGE_0_25 == (0.0, 25.0)
    assert CV_RELIABLE_BINS == 500
    assert CV_RELIABLE_N_CELLS == 250_000
    assert square_bins_for_n_cells(CV_RELIABLE_N_CELLS) == (500, 500)
    assert CV_RELIABLE_KEEP_FRACS == (0.6, 0.8)
    x_range, y_range = reliable_cv_hist_ranges()
    assert x_range[0] == 0.0
    assert x_range[1] > 6.0
    assert y_range[0] == 0.0
    assert y_range[1] > 25.0


def test_accumulate_rel_mean_cv_histogram_grid_drops_oob(tmp_path: Path):
    table = tmp_path / "chr000.popdep.txt"
    table.write_text(_tiger_grid_text())
    x_range, y_range = reliable_cv_hist_ranges()
    counts, stats = accumulate_rel_mean_cv_histogram_grid(
        table, x_range=x_range, y_range=y_range, bins_x=2, bins_y=2
    )
    assert counts.shape == (2, 2)
    assert stats["n_rows"] == 6
    assert stats["n_valid"] == 4
    assert stats["n_oob"] == 2
    assert int(counts.sum()) == 4
    assert int(counts[0, 0]) == 3


def test_write_chr_cv_reliable_beds_merge_keep60_inside_keep80(tmp_path: Path):
    table = tmp_path / "chr000.popdep.txt"
    table.write_text(_tiger_grid_text())
    x_range, y_range = reliable_cv_hist_ranges()
    npz = tmp_path / "chr000.rel_mean_cv.counts.npz"
    write_chr_mean_cv_histogram(
        table, npz, x_range=x_range, y_range=y_range, bins_x=2, bins_y=2
    )
    counts, _stats, hist_x, hist_y = load_histogram_npz(npz)
    mask_60, _info_60 = trusted_mask_from_site_coverage(counts, 0.6)
    mask_80, _info_80 = trusted_mask_from_site_coverage(counts, 0.8)
    masks_path = tmp_path / "masks.npz"
    np.savez_compressed(
        masks_path,
        mask_60=mask_60,
        mask_80=mask_80,
        x_range=np.asarray(hist_x, dtype=float),
        y_range=np.asarray(hist_y, dtype=float),
    )
    prefix = tmp_path / "chr000.popdep.cv_reliable"
    stats = write_chr_cv_reliable_beds(table, masks_path, prefix, batch_size=2)
    assert stats["n_kept60"] == 3
    assert stats["n_kept80"] == 4
    assert stats["n_intervals60"] == 1
    bed60 = tmp_path / "chr000.popdep.cv_reliable.keep60.bed.gz"
    bed80 = tmp_path / "chr000.popdep.cv_reliable.keep80.bed.gz"
    with gzip.open(bed60, "rt", encoding="utf-8") as handle:
        rows60 = [line.rstrip("\n").split("\t") for line in handle]
    with gzip.open(bed80, "rt", encoding="utf-8") as handle:
        rows80 = [line.rstrip("\n").split("\t") for line in handle]
    assert rows60 == [["0", "0", "3"]]
    assert rows80 == [["0", "0", "4"]]
    assert stats["n_intervals80"] == 1
    merged = tmp_path / "all.site.popdep.cv_reliable.keep80.bed.gz"
    n_lines = concat_gzip_beds([bed80], merged)
    assert n_lines == len(rows80)
    assert (tmp_path / "chr000.popdep.cv_reliable.lib.info.tsv").exists()


