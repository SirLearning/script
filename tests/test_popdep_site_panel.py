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
