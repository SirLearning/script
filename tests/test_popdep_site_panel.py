"""Tests for test_site random genome popdep panel."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from genetics.genomics.variant import popdep_site_panel as tsp


def test_sample_random_genome_sites_count_and_columns():
    df = tsp.sample_random_genome_sites(
        1000,
        seed=7,
        chr_file_ids=["chr043", "chr044"],
    )
    assert len(df) == 1000
    assert set(df.columns) >= {"Chr_id", "CHROM", "Position"}
    assert df["Chr_id"].isin(["chr043", "chr044"]).all()


def test_allocate_proportional_sums_to_target():
    counts = {"a": 100, "b": 300}
    alloc = tsp._allocate_proportional(counts, 40)
    assert sum(alloc.values()) == 40
    assert alloc["b"] > alloc["a"]


def test_scan_variant_positions(tmp_path: Path):
    pvar = tmp_path / "chr001.plink2.pvar"
    pvar.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "1\t100\tv100\tA\tG\n"
        "1\t200\tv200\tC\tT\n"
        "1\t300\tv300\tG\tA\n",
        encoding="utf-8",
    )
    found = tsp.scan_variant_positions(tmp_path, "chr001", {100, 300, 999})
    assert found == {100: "v100", 300: "v300"}

    # Early stop: positions above max(wanted) are never scanned.
    pvar.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\n"
        + "\n".join(f"1\t{pos}\tv{pos}\tA\tG" for pos in (10, 20, 500, 600))
        + "\n",
        encoding="utf-8",
    )
    assert tsp.scan_variant_positions(tmp_path, "chr001", {10, 20}) == {10: "v10", 20: "v20"}


def test_prepare_test_site_panel_smoke(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    popdep_dir = tmp_path / "popdep"
    main_raw = tmp_path / "main_raw"
    (popdep_dir / "variant").mkdir(parents=True)
    main_raw.mkdir()

    pvar = main_raw / "chr043.plink2.pvar"
    pvar.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "43\t10\tv10\tA\tG\n"
        "43\t20\tv20\tC\tT\n",
        encoding="utf-8",
    )

    def fake_sample(*args, **kwargs):
        return pd.DataFrame(
            {
                "Chr_id": ["chr043", "chr043", "chr044"],
                "CHROM": [43, 43, 44],
                "Position": [10, 20, 5],
            }
        )

    def fake_tabix(chrom, positions, popdep_dir_arg):
        return {
            pos: (1.0 + pos * 0.01, 0.1, 0.9 + pos * 0.001, 0.05)
            for pos in positions
        }

    monkeypatch.setattr(tsp, "sample_random_genome_sites", fake_sample)
    monkeypatch.setattr(tsp, "_tabix_lookup_chunked", lambda c, p, d: fake_tabix(c, p, d))
    monkeypatch.setattr(tsp, "list_plink_chr_file_ids", lambda **kwargs: ["chr043", "chr044"])

    out = tmp_path / "process" / "test_site"
    summary = tsp.prepare_test_site_panel(
        popdep_dir=popdep_dir,
        main_raw_dir=main_raw,
        output_dir=out,
        target_bp=3,
        seed=1,
        force=True,
    )
    panel = pd.read_csv(out / "sites" / "panel.sites.tsv", sep="\t")
    assert len(panel) == 3
    assert int(panel["Is_Variant"].sum()) == 2
    assert summary["n_variant_sites"] == 2
    assert (out / "sites" / "extract" / "chr043.variant.txt").exists()


def test_ana_test_site_popdep_plots(tmp_path: Path):
    proc = tmp_path / "test_site"
    sites = proc / "sites"
    sites.mkdir(parents=True)
    panel = pd.DataFrame(
        {
            "Chr_id": ["chr043"] * 4,
            "CHROM": [43] * 4,
            "Position": [1, 2, 3, 4],
            "Depth_Mean": [1.0, 2.0, 3.0, 4.0],
            "Depth_SD": [0.1, 0.2, 0.3, 0.4],
            "Depth_CV": [0.1, 0.1, 0.1, 0.1],
            "RelativeDepth_Mean": [0.9, 1.0, 1.1, 1.2],
            "RelativeDepth_SD": [0.05, 0.05, 0.05, 0.05],
            "RelativeDepth_CV": [0.05, 0.05, 0.05, 0.05],
            "Is_Variant": [False, True, False, True],
            "Variant_ID": [".", "v2", ".", "v4"],
        }
    )
    panel.to_csv(sites / "panel.sites.tsv", sep="\t", index=False)

    out_prefix = tmp_path / "stats" / "test_site.popdep"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tsp.ana_test_site_popdep_plots(proc, str(out_prefix), max_along_chr_sites=10)

    pngs = list((tmp_path / "stats").glob("*.png"))
    assert pngs, "expected popdep plot PNGs"
    summary = pd.read_csv(f"{out_prefix}.plot_summary.tsv", sep="\t")
    assert int(summary.iloc[0]["n_variant_sites"]) == 2


def test_allocate_test_popdep_sites_sums_to_target():
    alloc = tsp.allocate_test_popdep_sites(
        1000,
        chr_file_ids=["chr043", "chr044"],
    )
    assert int(alloc["n_sites"].sum()) == 1000
    assert set(alloc["Chr_id"]) == {"chr043", "chr044"}


def test_prepare_test_popdep_chr_panel(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    from genetics.genomics.variant.popdep_region_class import write_segment_bed

    anno = tmp_path / "annotation"
    write_segment_bed({43: [(0, 15)]}, anno / "gap.seg.bed")
    write_segment_bed({43: []}, anno / "repeat.seg.bed")
    write_segment_bed({43: []}, anno / "gc.seg.bed")
    write_segment_bed({43: []}, anno / "cnv.seg.bed")

    def fake_sample(chr_file_id, n, *, seed):
        return pd.DataFrame(
            {
                "Chr_id": [chr_file_id, chr_file_id],
                "CHROM": [43, 43],
                "Position": [10, 20],
            }
        )

    def fake_tabix(chrom, positions, popdep_dir_arg):
        return {pos: (2.0, 0.5, 1.0, 0.1) for pos in positions}

    monkeypatch.setattr(tsp, "sample_chr_sites", fake_sample)
    monkeypatch.setattr(tsp, "_tabix_lookup_chunked", fake_tabix)

    out = tmp_path / "chr043.popdep.sites.tsv"
    panel = tsp.prepare_test_popdep_chr_panel(
        "chr043",
        2,
        popdep_dir=tmp_path / "popdep",
        annotation_dir=anno,
        output_tsv=out,
        seed=1,
    )
    assert list(panel["region_class"]) == ["gap", "normal"]
    assert out.exists()


def test_ana_test_popdep_class_plots(tmp_path: Path):
    panel = pd.DataFrame(
        {
            "Chr_id": ["chr043"] * 5,
            "CHROM": [43] * 5,
            "Position": [1, 2, 3, 4, 5],
            "Depth_Mean": [1.0, 2.0, 3.0, 4.0, 5.0],
            "Depth_SD": [0.2, 0.3, 0.4, 0.5, 0.6],
            "Depth_CV": [0.2, 0.15, 0.13, 0.12, 0.12],
            "RelativeDepth_Mean": [0.8, 0.9, 1.0, 1.1, 1.2],
                "RelativeDepth_SD": [0.16, 0.25, 0.10, 0.40, 0.18],
            "region_class": ["gap", "repeat", "gc", "cnv", "normal"],
        }
    )
    prefix = tmp_path / "test_popdep.popdep"
    tsp.ana_test_popdep_class_plots(panel, str(prefix), scatter_cap=5)
    assert (tmp_path / "test_popdep.popdep.class_counts.tsv").exists()
    pngs = list(tmp_path.glob("test_popdep.popdep*.png"))
    assert pngs, "expected class-coloured popdep PNGs"
    assert (tmp_path / "test_popdep.popdep.reldepth.sd_vs_mean_by_class.png").exists()
    assert (tmp_path / "test_popdep.popdep.reldepth.cv_vs_mean_by_class.png").exists()
    assert (tmp_path / "test_popdep.popdep.scatter.sites.tsv").exists()


def test_ana_test_popdep_class_plots_joints_only_skips_mean_hist(tmp_path: Path):
    panel = pd.DataFrame(
        {
            "Depth_Mean": [1.0, 2.0, 3.0, 4.0, 5.0],
            "Depth_SD": [0.2, 0.3, 0.4, 0.5, 0.6],
            "RelativeDepth_Mean": [0.8, 0.9, 1.0, 1.1, 1.2],
            "RelativeDepth_SD": [0.16, 0.25, 0.10, 0.40, 0.18],
            "region_class": ["gap", "repeat", "gc", "cnv", "normal"],
        }
    )
    prefix = tmp_path / "joints_only.popdep"
    tsp.ana_test_popdep_class_plots(panel, str(prefix), scatter_cap=5, joints_only=True)
    assert (tmp_path / "joints_only.popdep.reldepth.sd_vs_mean_by_class.png").exists()
    assert not (tmp_path / "joints_only.popdep.depth.mean_dist_by_class.png").exists()


def test_hist_counts_by_group_clips_into_bins():
    edges = tsp.hist_bin_edges((0.0, 10.0), 2)
    counts = tsp.hist_counts_by_group(
        [1.0, 1.0, 9.0, 100.0],
        ["gap", "gap", "repeat", "repeat"],
        edges,
        group_order=("gap", "repeat"),
    )
    assert counts["gap"].tolist() == [2, 0]
    assert counts["repeat"].tolist() == [0, 2]


def test_ana_test_popdep_class_plots_from_chr_paths(tmp_path: Path):
    rows = (
        "Chr_id\tCHROM\tPosition\tDepth_Mean\tDepth_SD\tregion_class\n"
        "chr001\t1\t1\t1.0\t0.2\tgap\n"
        "chr001\t1\t2\t2.0\t0.3\trepeat\n"
        "chr001\t1\t3\t3.0\t0.4\tgc\n"
        "chr001\t1\t4\t4.0\t0.5\tcnv\n"
        "chr001\t1\t5\t5.0\t0.6\tnormal\n"
    )
    src = tmp_path / "chr001.popdep.sites.tsv"
    src.write_text(rows, encoding="utf-8")
    prefix = tmp_path / "from_paths"
    tsp.ana_test_popdep_class_plots([src], str(prefix), scatter_cap=5, hist_bins=4)
    summary = pd.read_csv(tmp_path / "from_paths.plot_summary.tsv", sep="\t")
    assert int(summary.loc[0, "n_sites"]) == 5
    assert (tmp_path / "from_paths.depth.mean_dist_by_class.png").exists()

def test_rescale_zero_padded_mean_sd_recovers_unpadded_panel():
    # Two taxa at depth 10, two padded zeros → n_old=4.
    mean_old, sd_old = 5.0, 5.0
    mean_new, sd_new = tsp.rescale_zero_padded_mean_sd(
        mean_old, sd_old, n_old=4, n_new=2
    )
    assert mean_new == pytest.approx(10.0)
    assert sd_new == pytest.approx(0.0, abs=1e-12)


def test_n500_crosschr_ntaxa_used_skips_chr032():
    assert tsp.n500_crosschr_ntaxa_used(32) == tsp.N500_PANEL_N
    assert tsp.n500_crosschr_ntaxa_used(31) == 7675
    assert tsp.n500_crosschr_ntaxa_used(3) == 7594
    assert tsp.n500_crosschr_ntaxa_used(5) == 7738
    assert tsp.n500_crosschr_ntaxa_used(0) == 8285


def test_rescale_test_popdep_sites_chr032_unchanged(tmp_path: Path):
    src = tmp_path / "chr032.popdep.sites.tsv"
    src.write_text(
        "Chr_id\tCHROM\tPosition\tDepth_Mean\tDepth_SD\tDepth_CV\t"
        "RelativeDepth_Mean\tRelativeDepth_SD\tRelativeDepth_CV\tregion_class\n"
        "chr032\t32\t10\t11.0\t5.0\t0.4545\t1.0\t0.4\t0.4\tnormal\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.tsv"
    df = tsp.rescale_test_popdep_sites(src, out)
    assert df["Depth_Mean"].iloc[0] == pytest.approx(11.0)
    assert df["Depth_SD"].iloc[0] == pytest.approx(5.0)
    assert out.exists()


def test_rescale_test_popdep_sites_a_subgenome_n_to_500(tmp_path: Path):
    src = tmp_path / "chr031.popdep.sites.tsv"
    n_old = 7675
    n_new = 500
    mean_old = 0.8
    sd_old = 3.5
    src.write_text(
        "Chr_id\tCHROM\tPosition\tDepth_Mean\tDepth_SD\tDepth_CV\t"
        "RelativeDepth_Mean\tRelativeDepth_SD\tRelativeDepth_CV\tregion_class\n"
        f"chr031\t31\t10\t{mean_old}\t{sd_old}\t4.375\t0.06\t0.26\t4.33\trepeat\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.tsv"
    df = tsp.rescale_test_popdep_sites(src, out)
    mean_new, sd_new = tsp.rescale_zero_padded_mean_sd(
        mean_old, sd_old, n_old=n_old, n_new=n_new
    )
    assert df["Depth_Mean"].iloc[0] == pytest.approx(mean_new)
    assert df["Depth_SD"].iloc[0] == pytest.approx(sd_new)
    assert df["Depth_CV"].iloc[0] == pytest.approx(sd_new / mean_new)
    rel_mean, rel_sd = tsp.rescale_zero_padded_mean_sd(
        0.06, 0.26, n_old=n_old, n_new=n_new
    )
    assert df["RelativeDepth_Mean"].iloc[0] == pytest.approx(rel_mean)
    assert df["RelativeDepth_SD"].iloc[0] == pytest.approx(rel_sd)
