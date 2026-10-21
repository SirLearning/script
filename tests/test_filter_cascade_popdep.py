"""PopDep-density step on the test_thin filter cascade."""

from pathlib import Path

import pandas as pd

from genetics.genomics.variant.filter_cascade_maf_plot import (
    estimate_full_library_from_thin,
    filter_cascade_test_thin,
    run_cv_reliable_three_step,
    run_filter_cascade_cv_reliable,
    variant_ids_in_reliable_beds,
    variant_ids_in_reliable_beds_from_ids,
)


def _write_gcount(path: Path) -> None:
    path.write_text(
        "CHROM\tID\tREF\tALT\tHOM_REF_CT\tHET_REF_ALT_CTS\tTWO_ALT_GENO_CTS\tHAP_REF_CT\tHAP_ALT_CTS\tMISSING_CT\n"
        "1\ts_keep\tA\tG\t90\t1\t10\t0\t0\t0\n"
        "1\ts_hetfail\tA\tG\t90\t20\t0\t0\t0\t0\n"
        "1\ts_popdep_fail\tA\tG\t90\t1\t10\t0\t0\t0\n"
        "1\ts_mac_fail\tA\tG\t97\t1\t1\t0\t0\t0\n"
    )


def test_popdep_then_het_then_mac5_records_each_step(tmp_path: Path):
    variant = tmp_path / "variant"
    variant.mkdir()
    _write_gcount(variant / "A.info.gcount")
    het = tmp_path / "A.info.het"
    het.write_text("IID\tF\na\t0.95\nb\t0.95\nc\t0.95\n")

    _, by_sg, totals = filter_cascade_test_thin(
        tmp_path,
        het_paths={"A": het},
        wap_keep_path=None,
        subgenomes=("A",),
        skip_wap=True,
        popdep_keep_ids={"s_keep", "s_hetfail", "s_mac_fail"},
        popdep_keep_ids_60={"s_keep"},
    )
    row = by_sg.iloc[0]
    assert int(row["n_raw"]) == 4
    assert int(row["n_popdep_60"]) == 1
    assert int(row["n_popdep_80"]) == 3
    assert int(row["n_het_pass"]) == 2
    assert int(row["n_mac_ge5"]) == 1
    assert int(row["n_wap_mac_ge2"]) == 1
    by_step = totals.set_index("step")["n_sites"].to_dict()
    assert by_step["raw"] == 4
    assert by_step["popdep_80"] == 3
    assert by_step["het_pass"] == 2
    assert by_step["mac_ge5"] == 1
    estimate = estimate_full_library_from_thin(totals, thin_rate=0.001)
    assert estimate.set_index("step").loc["wap_mac_ge2", "n_estimate_full"] == 1000.0


def test_reliable_bed_keeps_1based_position_inside_half_open_interval(tmp_path: Path):
    bed = tmp_path / "chr001.keep80.bed.gz"
    import gzip

    with gzip.open(bed, "wt", encoding="utf-8") as handle:
        handle.write("1\t9\t11\n")
    pvar = tmp_path / "A_test.plink2.pvar"
    pvar.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "1\t10\ts_in\tA\tG\n"
        "1\t11\ts_edge\tA\tG\n"
        "1\t12\ts_out\tA\tG\n"
        "2\t10\ts_other_chr\tA\tG\n"
    )
    assert variant_ids_in_reliable_beds([pvar], [bed]) == {"s_in", "s_edge"}


def test_cv_reliable_cascade_names_the_popdep_step(tmp_path: Path):
    import gzip

    variant = tmp_path / "variant"
    variant.mkdir()
    _write_gcount(variant / "A.info.gcount")
    het = tmp_path / "A.info.het"
    het.write_text("IID\tF\na\t0.95\nb\t0.95\nc\t0.95\n")
    pvar = tmp_path / "A_test.plink2.pvar"
    pvar.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "1\t10\ts_keep\tA\tG\n"
        "1\t10\ts_hetfail\tA\tG\n"
        "1\t10\ts_mac_fail\tA\tG\n"
        "1\t50\ts_popdep_fail\tA\tG\n"
    )
    bed = tmp_path / "chr001.keep60.bed.gz"
    with gzip.open(bed, "wt", encoding="utf-8") as handle:
        handle.write("1\t9\t11\n")
    prefix = tmp_path / "all.variant.cv_reliable.keep60"
    run_filter_cascade_cv_reliable(
        tmp_path,
        prefix,
        bed_paths=[bed],
        het_paths={"A": het},
        wap_keep_path=tmp_path / "unused.keep",
        keep_pct=60,
        subgenomes=("A",),
        thin_rate=0.001,
        skip_wap=True,
    )
    totals = pd.read_csv(tmp_path / "all.variant.cv_reliable.keep60.filter_steps_total.tsv", sep="\t")
    by_step = totals.set_index("step")["n_sites"].to_dict()
    assert "popdep_80" not in by_step
    assert by_step["raw"] == 4
    assert by_step["popdep_60"] == 3
    assert by_step["het_pass"] == 2
    assert by_step["mac_ge5"] == 1
    estimate = pd.read_csv(
        tmp_path / "all.variant.cv_reliable.keep60.filter_steps_genome_estimate.tsv", sep="\t"
    )
    assert float(estimate.set_index("step").loc["mac_ge5", "n_estimate_full"]) == 1000.0


def test_segment_id_bed_ignores_linear_refchr_coordinates(tmp_path: Path):
    import gzip

    bed = tmp_path / "chr002.keep80.bed.gz"
    with gzip.open(bed, "wt", encoding="utf-8") as handle:
        handle.write("2\t9\t11\n")
    pvar = tmp_path / "A_test.plink2.pvar"
    pvar.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "chr1A\t471305444\t2-10\tA\tG\n"
        "chr1A\t10\t1-10\tA\tG\n"
        "chr1A\t471305446\t2-12\tA\tG\n",
        encoding="utf-8",
    )
    assert variant_ids_in_reliable_beds_from_ids([pvar], [bed]) == {"2-10"}
    assert variant_ids_in_reliable_beds([pvar], [bed]) == set()


def test_three_step_filter_uses_segment_ids_and_skips_wap(tmp_path: Path):
    import gzip

    variant = tmp_path / "variant"
    variant.mkdir()
    (variant / "A.info.gcount").write_text(
        "CHROM\tID\tREF\tALT\tHOM_REF_CT\tHET_REF_ALT_CTS\tTWO_ALT_GENO_CTS\tHAP_REF_CT\tHAP_ALT_CTS\tMISSING_CT\n"
        "1\t1-10\tA\tG\t90\t1\t10\t0\t0\t0\n"
        "1\t1-11\tA\tG\t90\t20\t0\t0\t0\t0\n"
        "1\t1-12\tA\tG\t90\t1\t10\t0\t0\t0\n"
        "1\t1-40\tA\tG\t97\t1\t1\t0\t0\t0\n",
        encoding="utf-8",
    )
    pvar = tmp_path / "A_test.plink2.pvar"
    pvar.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "chr1A\t999\t1-10\tA\tG\n"
        "chr1A\t999\t1-11\tA\tG\n"
        "chr1A\t999\t1-12\tA\tG\n"
        "chr1A\t999\t1-40\tA\tG\n",
        encoding="utf-8",
    )
    bed = tmp_path / "chr001.keep80.bed.gz"
    with gzip.open(bed, "wt", encoding="utf-8") as handle:
        handle.write("1\t9\t11\n1\t39\t40\n")
    het = tmp_path / "A.info.het"
    het.write_text("IID\tF\na\t0.95\nb\t0.95\nc\t0.95\n", encoding="utf-8")

    final = run_cv_reliable_three_step(
        tmp_path,
        tmp_path / "all.variant.cv_reliable.keep80",
        bed_paths=[bed],
        keep_pct=80,
        het_paths={"A": het},
        subgenomes=("A",),
    )
    assert set(final["ID"].astype(str)) == {"1-10"}
    totals = pd.read_csv(tmp_path / "all.variant.cv_reliable.keep80.filter_steps_total.tsv", sep="\t")
    by_step = totals.set_index("step")["n_sites"].to_dict()
    assert list(by_step) == ["raw", "popdep_80", "het_pass", "mac_ge5"]
    assert by_step["raw"] == 4
    assert by_step["popdep_80"] == 3
    assert by_step["het_pass"] == 2
    assert by_step["mac_ge5"] == 1
    assert (tmp_path / "all.variant.cv_reliable.keep80.filter_steps.png").is_file()
    assert (tmp_path / "all.variant.cv_reliable.keep80.maf_rare_common.png").is_file()
    assert not (tmp_path / "all.variant.cv_reliable.keep80.filter_steps_genome_estimate.tsv").exists()

    final_mac1 = run_cv_reliable_three_step(
        tmp_path,
        tmp_path / "all.variant.cv_reliable.keep80.mac1",
        bed_paths=[bed],
        keep_pct=80,
        het_paths={"A": het},
        subgenomes=("A",),
        min_mac=1,
    )
    assert set(final_mac1["ID"].astype(str)) == {"1-10", "1-40"}
    mac1_steps = pd.read_csv(
        tmp_path / "all.variant.cv_reliable.keep80.mac1.filter_steps_total.tsv",
        sep="\t",
    )
    mac1_by_step = mac1_steps.set_index("step")["n_sites"].to_dict()
    assert list(mac1_by_step) == ["raw", "popdep_80", "het_pass", "mac_ge1"]
    assert mac1_by_step["mac_ge1"] == 2
