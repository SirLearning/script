"""Wheat selfing-model excess-heterozygosity tail probabilities and retain sweeps."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from genetics.genomics.variant.mac_het_excess import (
    ana_mac_het_excess_sweep,
    max_allowed_het,
    median_inbreeding_from_het_table,
    p_excess_het_binomial,
    site_excess_het_table,
    sweep_excess_het_thresholds,
)


def test_all_het_mac5_matches_selfing_tail():
    p = p_excess_het_binomial(h_obs=5, mac=5, f=0.95)
    assert p == pytest.approx(0.05**5)


def test_heatmap_x_tick_step_widens_with_mac_max():
    from genetics.genomics.variant.mac_het_excess import _heatmap_x_tick_step

    assert _heatmap_x_tick_step(100) == 10
    assert _heatmap_x_tick_step(200) == 20
    assert _heatmap_x_tick_step(500) == 50
    assert _heatmap_x_tick_step(1000) == 100


def test_singleton_het_is_not_extreme():
    p = p_excess_het_binomial(h_obs=1, mac=1, f=0.95)
    assert p == pytest.approx(0.05)


def test_zero_het_always_retained():
    p = p_excess_het_binomial(h_obs=0, mac=10, f=0.99)
    assert p == pytest.approx(1.0)


def test_max_allowed_het_widens_with_mac_but_stays_far_below_mac():
    allowed = [max_allowed_het(mac, f=0.95, alpha=1e-6) for mac in (1, 5, 20, 50)]
    assert allowed[0] == 1
    assert allowed == sorted(allowed)
    assert allowed[-1] < 50


def test_sweep_retains_more_sites_as_alpha_gets_smaller():
    sites = site_excess_het_table(
        pd.DataFrame(
            [
                {"HET_REF_ALT_CTS": 1, "TWO_ALT_GENO_CTS": 0, "HOM_REF_CT": 100, "HAP_REF_CT": 0, "HAP_ALT_CTS": 0},
                {"HET_REF_ALT_CTS": 5, "TWO_ALT_GENO_CTS": 0, "HOM_REF_CT": 100, "HAP_REF_CT": 0, "HAP_ALT_CTS": 0},
            ]
        ),
        f_panel=0.95,
    )
    sweep = sweep_excess_het_thresholds(sites, f_panel=0.95, thresholds=(1e-2, 1e-6, 1e-12))
    retain = sweep.set_index("p_threshold")["n_retain"]
    assert retain.loc[1e-12] >= retain.loc[1e-6] >= retain.loc[1e-2]
    assert int(retain.loc[1e-2]) == 1
    assert int(retain.loc[1e-6]) == 1
    assert int(retain.loc[1e-12]) == 2


def test_median_inbreeding_uses_sample_f_not_mean():
    het = pd.DataFrame({"IID": ["a", "b", "c"], "F": [0.99, 0.98, 0.10]})
    assert median_inbreeding_from_het_table(het) == pytest.approx(0.98)


def test_ana_writes_sweep_not_a_filtered_pfile(tmp_path: Path):
    gcount = tmp_path / "A.info.gcount"
    gcount.write_text(
        "CHROM\tID\tREF\tALT\tHOM_REF_CT\tHET_REF_ALT_CTS\tTWO_ALT_GENO_CTS\tHAP_REF_CT\tHAP_ALT_CTS\tMISSING_CT\n"
        "1\t1-1\tA\tG\t100\t1\t0\t0\t0\t0\n"
        "1\t1-2\tA\tG\t100\t5\t0\t0\t0\t0\n"
    )
    het = tmp_path / "A.info.het"
    het.write_text("#FID\tIID\tO_HOM\tE_HOM\tO_HET\tE_HET\tOBS_CT\tF\n0\ts1\t10\t8\t1\t4\t11\t0.95\n")
    prefix = str(tmp_path / "A.variant.mac_het")
    ana_mac_het_excess_sweep(str(gcount), prefix, het_path=str(het), f_grid=(0.95,), thresholds=(1e-6, 1e-12))
    sweep = pd.read_csv(f"{prefix}.sweep.tsv", sep="\t")
    assert (tmp_path / "A.variant.mac_het.sweep.tsv").exists()
    assert not list(tmp_path.glob("*.pgen"))
    assert set(sweep["n_tested"]) == {2}
    data_f_rows = sweep[sweep["f_source"] == "het_median"]
    assert not data_f_rows.empty
    assert np.isclose(data_f_rows["f_panel"].iloc[0], 0.95)


def test_retain_drops_mac5_all_het_at_f095():
    sites = site_excess_het_table(
        pd.DataFrame(
            [
                {"HET_REF_ALT_CTS": 1, "TWO_ALT_GENO_CTS": 0, "HOM_REF_CT": 100, "HAP_REF_CT": 0, "HAP_ALT_CTS": 0},
                {"HET_REF_ALT_CTS": 5, "TWO_ALT_GENO_CTS": 0, "HOM_REF_CT": 100, "HAP_REF_CT": 0, "HAP_ALT_CTS": 0},
            ]
        ),
        f_panel=0.95,
    )
    from genetics.genomics.variant.mac_het_excess import drop_sites_by_p_excess, retain_sites_by_p_excess

    kept = retain_sites_by_p_excess(sites, alpha=1e-6)
    dropped = drop_sites_by_p_excess(sites, alpha=1e-6)
    assert len(kept) == 1
    assert int(kept["MAC"].iloc[0]) == 1
    assert len(dropped) == 1
    assert int(dropped["MAC"].iloc[0]) == 5


def test_retained_heatmap_uses_distinct_filename(tmp_path: Path):
    from genetics.genomics.variant.mac_het_excess import HEATMAP_MAC_MAXES, ana_mac_hetfrac_heatmap_retained

    rows = [
        "CHROM\tID\tREF\tALT\tHOM_REF_CT\tHET_REF_ALT_CTS\tTWO_ALT_GENO_CTS\tHAP_REF_CT\tHAP_ALT_CTS\tMISSING_CT"
    ]
    for i in range(12):
        rows.append(f"1\t1-{i}\tA\tG\t100\t1\t0\t0\t0\t0")
    for i in range(12):
        rows.append(f"1\t10-{i}\tA\tG\t100\t10\t0\t0\t0\t0")
    gcount = tmp_path / "A.info.gcount"
    gcount.write_text("\n".join(rows) + "\n")
    het = tmp_path / "A.info.het"
    het.write_text("#FID\tIID\tO_HOM\tE_HOM\tO_HET\tE_HET\tOBS_CT\tF\n0\ts1\t10\t8\t1\t4\t11\t0.95\n")
    prefix = str(tmp_path / "A.variant.mac")
    ana_mac_hetfrac_heatmap_retained(str(gcount), prefix, het_path=str(het), alpha=1e-6)
    for mac_max in HEATMAP_MAC_MAXES:
        for kind in ("hetfrac", "h"):
            for status in ("", ".retained", ".dropped"):
                stem = f"heatmap_mac_{kind}{status}.0_{mac_max}.png"
                assert (tmp_path / f"A.variant.mac.{stem}").exists(), stem
    th = pd.read_csv(f"{prefix}.het_retain.th.tsv", sep="\t")
    assert int(th.loc[0, "N_retain"]) == 12
    assert not list(tmp_path.glob("*.pgen"))
    from genetics.genomics.variant.mac_het_excess import HEATMAP_MAC_FILTER_OVERLAYS

    for mac_max in HEATMAP_MAC_FILTER_OVERLAYS:
        stem = f"A.variant.mac.heatmap_mac_h.filter.0_{mac_max}.png"
        assert (tmp_path / stem).exists(), stem
    assert not (tmp_path / "A.variant.mac.heatmap_mac_h.filter.0_1000.png").exists()


def test_excess_het_cutoff_heatmap_coords_sit_on_retained_h_top_edge():
    from genetics.genomics.variant.mac_het_excess import (
        excess_het_cutoff_heatmap_coords,
        max_allowed_het,
    )

    mac_max = 100
    f_panel = 0.95
    alpha = 1e-6
    xs, ys = excess_het_cutoff_heatmap_coords(mac_max, f_panel, alpha=alpha)
    assert xs[0] == 1.0
    assert xs[1] == 2.0
    h5 = max_allowed_het(5, f_panel, alpha)
    # MAC=5 occupies heatmap x [5, 6]; two vertices per MAC starting at MAC=1.
    i5 = (5 - 1) * 2
    assert xs[i5 : i5 + 2] == [5.0, 6.0]
    assert ys[i5] == pytest.approx(mac_max - h5)
    assert ys[i5] == ys[i5 + 1]


def test_popdep_het_filter_scatter_colors_pass_and_fail(tmp_path: Path):
    from genetics.genomics.variant.mac_het_excess import (
        HET_FILTER_FAIL,
        HET_FILTER_PASS,
        _join_popdep_het_filter,
        ana_popdep_scatter_by_het_filter,
        site_excess_het_table,
    )

    gcount = tmp_path / "A.info.gcount"
    gcount.write_text(
        "CHROM\tID\tREF\tALT\tHOM_REF_CT\tHET_REF_ALT_CTS\tTWO_ALT_GENO_CTS\tHAP_REF_CT\tHAP_ALT_CTS\tMISSING_CT\n"
        "1\t1-1\tA\tG\t100\t1\t0\t0\t0\t0\n"
        "1\t1-2\tA\tG\t100\t5\t0\t0\t0\t0\n"
        "1\t1-3\tA\tG\t100\t0\t0\t0\t0\t0\n"
    )
    het = tmp_path / "A.info.het"
    het.write_text("#FID\tIID\tO_HOM\tE_HOM\tO_HET\tE_HET\tOBS_CT\tF\n0\ts1\t10\t8\t1\t4\t11\t0.95\n")
    popdep = tmp_path / "A.popdep.info.tsv"
    popdep.write_text(
        "#CHROM\tID\tPOS\tDepth_Mean\tDepth_SD\tDepth_CV\tRelativeDepth_Mean\tRelativeDepth_SD\tRelativeDepth_CV\n"
        "1\t1-1\t1\t10.0\t2.0\t0.2\t1.0\t0.2\t0.2\n"
        "1\t1-2\t2\t20.0\t4.0\t0.2\t2.0\t0.4\t0.2\n"
        "1\t1-3\t3\t30.0\t6.0\t0.2\t3.0\t0.6\t0.2\n"
    )
    sites = site_excess_het_table(
        pd.read_csv(gcount, sep="\t"),
        f_panel=0.95,
    )
    merged = _join_popdep_het_filter(sites, pd.read_csv(popdep, sep="\t"), alpha=1e-6)
    assert set(merged["ID"]) == {"1-1", "1-2"}
    by_id = merged.set_index("ID")["het_filter"]
    assert int(by_id.loc["1-1"]) == HET_FILTER_PASS
    assert int(by_id.loc["1-2"]) == HET_FILTER_FAIL

    prefix = str(tmp_path / "A.variant.popdep.het_filter")
    ana_popdep_scatter_by_het_filter(
        str(gcount),
        str(popdep),
        prefix,
        het_path=str(het),
        alpha=1e-6,
    )
    assert (tmp_path / "A.variant.popdep.het_filter.scatter_Depth_SD_vs_Depth_Mean.png").exists()
    assert (tmp_path / "A.variant.popdep.het_filter.scatter_Depth_SD_vs_Depth_Mean.log.png").exists()
    assert (
        tmp_path / "A.variant.popdep.het_filter.scatter_RelativeDepth_SD_vs_RelativeDepth_Mean.png"
    ).exists()
    assert (
        tmp_path
        / "A.variant.popdep.het_filter.scatter_RelativeDepth_SD_vs_RelativeDepth_Mean.log.png"
    ).exists()
    assert (tmp_path / "A.variant.popdep.het_filter.scatter_Depth_CV_vs_Depth_Mean.png").exists()
    assert (tmp_path / "A.variant.popdep.het_filter.scatter_Depth_CV_vs_Depth_Mean.log.png").exists()
    assert (
        tmp_path / "A.variant.popdep.het_filter.scatter_RelativeDepth_CV_vs_RelativeDepth_Mean.png"
    ).exists()
    assert (
        tmp_path
        / "A.variant.popdep.het_filter.scatter_RelativeDepth_CV_vs_RelativeDepth_Mean.log.png"
    ).exists()
    assert not list(tmp_path.glob("*.pgen"))
