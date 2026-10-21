"""R_H = H/MAC band classification and robust inbreeding for wheat selfing QC."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from genetics.genomics.variant.mac_het_rh import (
    RH_BOTTOM,
    RH_DIAGONAL,
    RH_MAC_LT_MIN,
    RH_MIDDLE,
    allhet_fraction_by_mac,
    ana_allhet_sample_burden,
    ana_mac_het_rh_diagnostics,
    annotate_rh_class,
    inbreeding_from_het_fraction,
    sample_het_concentration,
)


def _gcount_frame():
    # MAC = H + 2A when ALT is minor. Bottom: H=0,A=5 → MAC=10.
    # Mixed: H=4,A=3 → MAC=10. Diagonal: H=10,A=0 → MAC=10.
    # Low MAC kept out of the bands: H=2,A=0 → MAC=2.
    return pd.DataFrame(
        [
            {
                "CHROM": 1,
                "ID": "1-1",
                "REF": "A",
                "ALT": "G",
                "HOM_REF_CT": 100,
                "HET_REF_ALT_CTS": 0,
                "TWO_ALT_GENO_CTS": 5,
                "HAP_REF_CT": 0,
                "HAP_ALT_CTS": 0,
                "MISSING_CT": 2,
            },
            {
                "CHROM": 1,
                "ID": "1-2",
                "REF": "A",
                "ALT": "G",
                "HOM_REF_CT": 100,
                "HET_REF_ALT_CTS": 4,
                "TWO_ALT_GENO_CTS": 3,
                "HAP_REF_CT": 0,
                "HAP_ALT_CTS": 0,
                "MISSING_CT": 1,
            },
            {
                "CHROM": 1,
                "ID": "1-3",
                "REF": "A",
                "ALT": "G",
                "HOM_REF_CT": 100,
                "HET_REF_ALT_CTS": 10,
                "TWO_ALT_GENO_CTS": 0,
                "HAP_REF_CT": 0,
                "HAP_ALT_CTS": 0,
                "MISSING_CT": 0,
            },
            {
                "CHROM": 1,
                "ID": "1-4",
                "REF": "A",
                "ALT": "G",
                "HOM_REF_CT": 100,
                "HET_REF_ALT_CTS": 2,
                "TWO_ALT_GENO_CTS": 0,
                "HAP_REF_CT": 0,
                "HAP_ALT_CTS": 0,
                "MISSING_CT": 0,
            },
        ]
    )


def test_rh_class_splits_bottom_middle_diagonal_and_low_mac():
    from genetics.genomics.variant.mac import _compute_mac_table

    sites = annotate_rh_class(_compute_mac_table(_gcount_frame()), mac_min=5)
    by_id = sites.set_index("ID")
    assert by_id.loc["1-1", "rh_class"] == RH_BOTTOM
    assert by_id.loc["1-1", "R_H"] == pytest.approx(0.0)
    assert by_id.loc["1-2", "rh_class"] == RH_MIDDLE
    assert 0.0 < float(by_id.loc["1-2", "R_H"]) < 1.0
    assert by_id.loc["1-3", "rh_class"] == RH_DIAGONAL
    assert by_id.loc["1-3", "R_H"] == pytest.approx(1.0)
    assert bool(by_id.loc["1-3", "all_het"])
    assert by_id.loc["1-4", "rh_class"] == RH_MAC_LT_MIN
    assert not bool(by_id.loc["1-4", "all_het"])


def test_excluding_diagonal_raises_site_moment_f():
    from genetics.genomics.variant.mac import _compute_mac_table

    sites = annotate_rh_class(_compute_mac_table(_gcount_frame()), mac_min=5)
    f_all = inbreeding_from_het_fraction(sites, mac_min=5, exclude_diagonal=False)
    f_ex = inbreeding_from_het_fraction(sites, mac_min=5, exclude_diagonal=True)
    f_bottom = inbreeding_from_het_fraction(sites, mac_min=5, rh_classes=(RH_BOTTOM,))
    assert f_ex["F_mom"] > f_all["F_mom"]
    assert f_bottom["F_mom"] > f_ex["F_mom"]
    assert f_bottom["F_mom"] == pytest.approx(1.0)


def test_allhet_fraction_at_mac10_is_one_third_in_fixture():
    from genetics.genomics.variant.mac import _compute_mac_table

    sites = annotate_rh_class(_compute_mac_table(_gcount_frame()), mac_min=5)
    by_mac = allhet_fraction_by_mac(sites, mac_max=20, f_panel=0.95, n_diploid=110.0)
    row = by_mac.set_index("MAC").loc[10]
    assert int(row["n_sites"]) == 3
    assert int(row["n_all_het"]) == 1
    assert float(row["frac_all_het"]) == pytest.approx(1.0 / 3.0)
    assert float(row["expected_frac_all_het"]) < 1e-10


def test_sample_het_concentration_flags_one_dominant_sample():
    het = pd.DataFrame(
        {
            "IID": ["a", "b", "c", "d"],
            "O_HET": [97, 1, 1, 1],
        }
    )
    stats = sample_het_concentration(het)
    assert stats["n_samples"] == 4
    assert stats["n_for_half_hets"] == 1
    assert stats["max_frac"] == pytest.approx(0.97)


def test_sample_het_concentration_derives_ohet_from_ohom():
    het = pd.DataFrame(
        {
            "IID": ["a", "b"],
            "O(HOM)": [10, 90],
            "OBS_CT": [100, 100],
            "F": [-1.0, 0.9],
        }
    )
    stats = sample_het_concentration(het)
    assert stats["sum_O_HET"] == pytest.approx(100.0)
    assert stats["max_frac"] == pytest.approx(0.9)


def test_ana_writes_rh_tables_plots_and_extract(tmp_path: Path):
    gcount = tmp_path / "A.info.gcount"
    rows = [
        "CHROM\tID\tREF\tALT\tHOM_REF_CT\tHET_REF_ALT_CTS\tTWO_ALT_GENO_CTS\tHAP_REF_CT\tHAP_ALT_CTS\tMISSING_CT"
    ]
    for i in range(12):
        rows.append(f"1\tb{i}\tA\tG\t100\t0\t5\t0\t0\t1")
    for i in range(12):
        rows.append(f"1\tm{i}\tA\tG\t100\t4\t3\t0\t0\t2")
    for i in range(12):
        rows.append(f"1\td{i}\tA\tG\t100\t10\t0\t0\t0\t0")
    gcount.write_text("\n".join(rows) + "\n")
    het = tmp_path / "A.info.het"
    het.write_text("#FID\tIID\tO_HOM\tE_HOM\tO_HET\tE_HET\tOBS_CT\tF\n0\ts1\t10\t8\t1\t4\t11\t0.86\n")
    popdep = tmp_path / "A.popdep.info.tsv"
    mq30 = tmp_path / "A.mapq30.popdep.info.tsv"
    mq = tmp_path / "A.mq.info.tsv"
    vmiss = tmp_path / "A.info.vmiss"
    pop_rows = [
        "#CHROM\tID\tPOS\tDepth_Mean\tDepth_SD\tDepth_CV\tRelativeDepth_Mean\tRelativeDepth_SD\tRelativeDepth_CV"
    ]
    mq30_rows = list(pop_rows)
    mq_rows = ["#CHROM\tID\tPOS\tMQ"]
    vmiss_rows = ["#CHROM\tID\tMISSING_CT\tOBS_CT\tF_MISS"]
    for prefix, depth, depth30, mq_val, miss in (
        ("b", 8.0, 7.5, 55.0, 0.01),
        ("m", 10.0, 9.0, 50.0, 0.02),
        ("d", 25.0, 12.0, 25.0, 0.08),
    ):
        for i in range(12):
            vid = f"{prefix}{i}"
            pop_rows.append(f"1\t{vid}\t{i}\t{depth}\t2.0\t0.2\t1.0\t0.2\t0.2")
            mq30_rows.append(f"1\t{vid}\t{i}\t{depth30}\t2.0\t0.2\t1.0\t0.2\t0.2")
            mq_rows.append(f"1\t{vid}\t{i}\t{mq_val}")
            vmiss_rows.append(f"1\t{vid}\t1\t110\t{miss}")
    popdep.write_text("\n".join(pop_rows) + "\n")
    mq30.write_text("\n".join(mq30_rows) + "\n")
    mq.write_text("\n".join(mq_rows) + "\n")
    vmiss.write_text("\n".join(vmiss_rows) + "\n")

    prefix = str(tmp_path / "A.variant.mac_het.rh")
    ana_mac_het_rh_diagnostics(
        str(gcount),
        prefix,
        het_path=str(het),
        popdep_file=str(popdep),
        popdep_mapq30_file=str(mq30),
        mq_file=str(mq),
        vmiss_file=str(vmiss),
    )
    th = pd.read_csv(f"{prefix}.th.tsv", sep="\t")
    assert float(th.loc[0, "F_het_median"]) == pytest.approx(0.86)
    assert float(th.loc[0, "F_mom_exclude_diagonal"]) > float(th.loc[0, "F_mom_all"])
    assert int(th.loc[0, "n_diagonal"]) == 12
    classes = pd.read_csv(f"{prefix}.rh_class.info.tsv", sep="\t")
    assert set(classes["rh_class"]) >= {RH_BOTTOM, RH_MIDDLE, RH_DIAGONAL}
    extract = Path(f"{prefix}.all_het.extract")
    assert extract.exists()
    ids = [line for line in extract.read_text().splitlines() if line]
    assert len(ids) == 12
    assert all(i.startswith("d") for i in ids)
    enrich = pd.read_csv(f"{prefix}.rh_enrich.info.tsv", sep="\t")
    by_class = enrich.set_index("rh_class")
    assert float(by_class.loc[RH_DIAGONAL, "median_Depth_Mean"]) > float(
        by_class.loc[RH_BOTTOM, "median_Depth_Mean"]
    )
    assert (tmp_path / "A.variant.mac_het.rh.rh_hist.png").exists()
    assert (tmp_path / "A.variant.mac_het.rh.frac_allhet_by_mac.png").exists()
    assert not list(tmp_path.glob("*.pgen"))


def test_ana_allhet_sample_burden_writes_concentration(tmp_path: Path):
    het = tmp_path / "A.all_het.info.het"
    het.write_text(
        "#FID\tIID\tO_HOM\tE_HOM\tO_HET\tE_HET\tOBS_CT\tF\n"
        "0\ta\t10\t8\t90\t4\t100\t-2.0\n"
        "0\tb\t10\t8\t5\t4\t100\t0.2\n"
        "0\tc\t10\t8\t5\t4\t100\t0.2\n"
    )
    prefix = str(tmp_path / "A.variant.mac_het.rh_sample")
    ana_allhet_sample_burden(str(het), prefix)
    th = pd.read_csv(f"{prefix}.th.tsv", sep="\t")
    assert int(th.loc[0, "n_for_half_hets"]) == 1
    assert (tmp_path / "A.variant.mac_het.rh_sample.ohet_hist.png").exists()
