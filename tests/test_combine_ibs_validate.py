"""Combine-name duplicate validation helpers."""

import numpy as np
import pandas as pd
import pytest

from infra.utils.graph import plot_grouped_histogram

from genetics.genomics.sample.combine_ibs_validate import (
    _load_highlight_pairs,
    _random_group_stats,
    add_dedup_fields,
    assign_ibs_tiers,
    assign_random_classes,
    build_combine_pairs,
    load_final_merge_edges,
    load_taxa_coverage,
    map_to_genotype_name,
    parse_name_list,
    plot_combine_ibs_low_coverage,
    plot_merge_coverage_hist,
    summarize_by_ploidy,
)


def test_parse_name_list_splits_pipes_and_strips():
    assert parse_name_list("P0025|P0050") == ["P0025", "P0050"]
    assert parse_name_list("") == []
    assert parse_name_list(None) == []


def test_map_to_genotype_name_uses_trailing_a():
    genotype = {"P0030A", "ABD_0001"}
    assert map_to_genotype_name("P0030", genotype) == "P0030A"
    assert map_to_genotype_name("ABD_0001", genotype) == "ABD_0001"
    assert map_to_genotype_name("P0001", genotype) is None


def test_build_combine_pairs_excludes_controls(tmp_path):
    combine = tmp_path / "combine.tsv"
    rows = [
        {
            "germplasm_name": "A",
            "wap_pq_combine_names": "P0030|P0050",
            "vmap3_combine_names": "ABD_0001",
        },
        {
            "germplasm_name": "B",
            "wap_pq_combine_names": "P0025|P0050|P0075",
            "vmap3_combine_names": "",
        },
        {
            "germplasm_name": "C",
            "wap_pq_combine_names": "",
            "vmap3_combine_names": "ABD_0001|ABD_0002",
        },
        {
            "germplasm_name": "D",
            "wap_pq_combine_names": "P0100",
            "vmap3_combine_names": "",
        },
        {
            "germplasm_name": "E",
            "wap_pq_combine_names": "P0200|P0210",
            "vmap3_combine_names": "",
        },
    ]
    pd.DataFrame(rows).to_csv(combine, sep="\t", index=False)
    genotype = {
        "P0030A",
        "P0050A",
        "P0025A",
        "P0075A",
        "P0100A",
        "P0200A",
        "P0210A",
        "ABD_0001",
        "ABD_0002",
    }
    pairs, stats = build_combine_pairs(
        combine, genotype, excluded_names={"P0050", "P0200", "P0210"}
    )
    assert len(pairs) == 3
    assert not pairs["sample1"].str.contains("P0050").any()
    assert not pairs["sample2"].str.contains("P0050").any()
    assert set(pairs["pair_type"]) == {
        "cross",
        "within_wap",
        "within_vmap3",
    }
    assert stats["germplasm_rows"] == 5
    assert stats["single_name_rows_raw"] == 1
    assert stats["multi_name_rows_raw"] == 4
    assert stats["rows_with_zero_names_after_exclusion"] == 1
    assert stats["rows_with_ge2_genotype_names"] == 3
    assert stats["excluded_name_occurrences"] == 4


def test_build_combine_pairs_multi_cohort(tmp_path):
    combine = tmp_path / "combine_multi.tsv"
    rows = [
        {
            "germplasm_name": "A",
            "ipk_combine_names": "SAMEA1|SAMEA2",
            "hzau_combine_names": "",
            "watkins_combine_names": "",
        },
        {
            "germplasm_name": "B",
            "ipk_combine_names": "SAMEA3",
            "hzau_combine_names": "JZ1",
            "watkins_combine_names": "",
        },
        {
            "germplasm_name": "C",
            "ipk_combine_names": "",
            "hzau_combine_names": "",
            "watkins_combine_names": "CRR1|CRR2|CRR3",
        },
        {
            "germplasm_name": "D",
            "ipk_combine_names": "",
            "hzau_combine_names": "JZ2",
            "watkins_combine_names": "",
        },
    ]
    pd.DataFrame(rows).to_csv(combine, sep="\t", index=False)
    genotype = {
        "SAMEA1",
        "SAMEA2",
        "SAMEA3",
        "JZ1",
        "JZ2",
        "CRR1",
        "CRR2",
        "CRR3",
    }
    pairs, stats = build_combine_pairs(
        combine,
        genotype,
        name_columns=[
            "ipk_combine_names",
            "hzau_combine_names",
            "watkins_combine_names",
        ],
        cohort_labels=["ipk", "hzau", "watkins"],
    )
    assert len(pairs) == 5
    counts = pairs["pair_type"].value_counts().to_dict()
    assert counts == {"within_watkins": 3, "within_ipk": 1, "cross": 1}
    assert stats["germplasm_rows"] == 4
    assert stats["rows_with_ge2_genotype_names"] == 3


def test_random_classes_cover_non_aabbdd_groups(tmp_path):
    random_path = tmp_path / "random.tsv"
    bounds = {
        "AABBDD": (0.95, 0.99),
        "AABB": (0.90, 0.95),
        "DD": (0.80, 0.90),
        "AA": (0.80, 0.90),
        "SS": (0.85, 0.90),
    }
    rows = []
    for group, (low, high) in bounds.items():
        for value in np.linspace(low, high, 20):
            rows.append(
                {
                    "ploidy": group,
                    "plink_ibs": value,
                    "hom_ibs": min(value + 0.005, 1.0),
                }
            )
    pd.DataFrame(rows).to_csv(random_path, sep="\t", index=False)
    _, by_group = _random_group_stats(random_path)
    pairs = pd.DataFrame(
        [
            {
                "ploidy1": "AABBDD",
                "ploidy2": "AABBDD",
                "plink_ibs": 0.999,
                "hom_ibs": 0.999,
                "common_subgenomes": "ABD",
                "hom_by_sg": "A:0.999;B:0.999;D:0.999",
            },
            {
                "ploidy1": "DD",
                "ploidy2": "DD",
                "plink_ibs": 0.95,
                "hom_ibs": 0.98,
                "common_subgenomes": "D",
                "hom_by_sg": "D:0.98",
            },
            {
                "ploidy1": "AABB",
                "ploidy2": "AABB",
                "plink_ibs": 0.96,
                "hom_ibs": 0.99,
                "common_subgenomes": "AB",
                "hom_by_sg": "A:0.99;B:0.99",
            },
            {
                "ploidy1": "AABB",
                "ploidy2": "AABBDD",
                "plink_ibs": 0.97,
                "hom_ibs": 0.99,
                "common_subgenomes": "AB",
                "hom_by_sg": "A:0.99;B:0.99",
            },
        ]
    )
    out = assign_random_classes(pairs, by_group)
    assert out.loc[0, "random_group"] == "AABBDD"
    assert out.loc[0, "plink_random_class"] == "high"
    assert out.loc[1, "random_group"] == "DD"
    assert out.loc[1, "plink_random_class"] == "high"
    assert out.loc[2, "random_group"] == "AABB"
    assert out.loc[2, "hom_random_class"] == "high"
    assert out.loc[3, "random_group"] == "NA"

    summary = summarize_by_ploidy(out, by_group)
    indexed = summary.set_index("ploidy_group")
    assert {"AABBDD", "AABB", "DD"}.issubset(indexed.index)
    assert indexed.loc["DD", "subgenomes"] == "D"
    assert indexed.loc["AABB", "subgenomes"] == "AB"
    assert indexed.loc["AABBDD", "subgenomes"] == "ABD"
    assert indexed.loc["DD", "dup_hom_D_median"] == pytest.approx(0.98)


def test_plot_grouped_histogram_writes_png(tmp_path):
    frame = pd.DataFrame(
        {
            "value": [0.90, 0.92, 0.94, 0.96, 0.98, 0.99],
            "group": ["random"] * 3 + ["duplicate"] * 3,
        }
    )
    output = tmp_path / "distribution.hist.png"
    plot_grouped_histogram(
        data=frame,
        x_col="value",
        hue_col="group",
        hue_order=["random", "duplicate"],
        palette_map={"random": "#2a78d6", "duplicate": "#eb6834"},
        title="test",
        filename=str(output),
        xlim=(0.89, 1.0),
        bins=20,
    )
    assert output.exists()
    assert output.stat().st_size > 0


def _tier_frame():
    rng = np.random.default_rng(0)
    low = 0.965 + rng.normal(0, 0.001, 30)
    mid = 0.982 + rng.normal(0, 0.001, 30)
    high = 0.994 + rng.normal(0, 0.001, 30)
    plink = np.concatenate([low, mid, high])
    hom = np.concatenate([low + 0.008, mid + 0.008, high + 0.005])
    frame = pd.DataFrame(
        {
            "ploidy1": ["AABBDD"] * 90 + ["AABB", "DD"],
            "ploidy2": ["AABBDD"] * 90 + ["AABB", "DD"],
            "plink_ibs": np.concatenate([plink, [0.99, 0.99]]),
            "hom_ibs": np.concatenate([hom, [0.999, 0.999]]),
        }
    )
    return frame


def test_assign_ibs_tiers_separates_three_modes():
    frame = assign_ibs_tiers(_tier_frame())
    abd = frame.iloc[:90]
    assert set(abd["tier"]) == {"low", "mid", "high"}
    assert (abd.iloc[:30]["tier"] == "low").mean() > 0.9
    assert (abd.iloc[30:60]["tier"] == "mid").mean() > 0.9
    assert (abd.iloc[60:]["tier"] == "high").mean() > 0.9
    assert set(frame.iloc[90:]["tier"]) == {"NA"}


def test_assign_ibs_tiers_uses_hom_metric_when_requested():
    frame = _tier_frame()
    frame["plink_ibs"] = 0.5
    out = assign_ibs_tiers(
        frame,
        metric="hom_ibs",
        low_max=0.983,
        high_min=0.995,
    )
    abd = out.iloc[:90]
    assert (abd.iloc[:30]["tier"] == "low").mean() > 0.9
    assert (abd.iloc[30:60]["tier"] == "mid").mean() > 0.9
    assert (abd.iloc[60:]["tier"] == "high").mean() > 0.9


def test_plot_merge_coverage_hist_uses_multi_member_rows(tmp_path):
    merge_map = tmp_path / "taxaMergeMap.txt"
    merge_map.write_text(
        "merged_ab\t.\ta\tb\n"
        "c\t.\tc\n"
    )
    combine = tmp_path / "combine.tsv"
    pd.DataFrame(
        [
            {"members": "a|b", "n_members": "2", "sources": "WAP:2"},
            {"members": "c", "n_members": "1", "sources": "WAP:1"},
        ]
    ).to_csv(combine, sep="\t", index=False)
    taxa = tmp_path / "taxa.tsv"
    pd.DataFrame(
        {
            "Taxa": ["a", "b", "c"],
            "Coverage-Of-All-Bams": ["4", "2", "10"],
            "Bams": ["a.bam", "b.bam", "c.bam"],
        }
    ).to_csv(taxa, sep="\t", index=False)

    stats = plot_merge_coverage_hist(
        merge_map, taxa, tmp_path / "coverage_hist", bins=10, x_max=20
    )

    assert (tmp_path / "coverage_hist.png").exists()
    assert (tmp_path / "coverage_hist.hist.png").exists()
    assert stats["rows"] == 1
    assert stats["member_samples"] == 2
    assert stats["before_median"] == pytest.approx(3.0)
    assert stats["after_median"] == pytest.approx(6.0)

    all_stats = plot_merge_coverage_hist(
        merge_map,
        taxa,
        tmp_path / "coverage_hist_all",
        bins=10,
        x_max=20,
        combine_map_path=combine,
        sources_pattern="WAP",
        include_singletons=True,
    )
    assert all_stats["rows"] == 2
    assert all_stats["member_samples"] == 3
    assert all_stats["after_median"] == pytest.approx(8.0)


def test_load_final_merge_edges_excludes_watkins_and_blocked(tmp_path):
    decisions = tmp_path / "pair_decisions.tsv"
    pd.DataFrame(
        [
            {
                "sample1": "a",
                "sample2": "b",
                "hom_ibs": "0.9995",
                "sources": "VMap3+WAP",
                "merges": "True",
            },
            {
                "sample1": "c",
                "sample2": "d",
                "hom_ibs": "0.9000",
                "sources": "WAP",
                "merges": "False",
            },
            {
                "sample1": "e",
                "sample2": "f",
                "hom_ibs": "0.9999",
                "sources": "Watkins",
                "merges": "True",
            },
            {
                "sample1": "g",
                "sample2": "h",
                "hom_ibs": "0.9999",
                "sources": "IPK",
                "merges": "True",
            },
        ]
    ).to_csv(decisions, sep="\t", index=False)

    edges = load_final_merge_edges(decisions)

    assert len(edges) == 1
    assert edges.iloc[0]["sample1"] == "a"
    assert edges.iloc[0]["sources"] == "VMap3+WAP"


def test_add_dedup_fields_fills_missing_coverage_from_taxa_bam_map(tmp_path):
    dedup = tmp_path / "dedup.tsv"
    pd.DataFrame(
        [
            {
                "Sample": "s1",
                "Group": "G",
                "Coverage": "5.0000",
                "Dup_Group_ID": "2",
                "Dup_Status": "Unique",
                "Clean_CN": "",
                "Clean_Acc": "",
            }
        ]
    ).to_csv(dedup, sep="\t", index=False)
    coverage = tmp_path / "taxaBamMap.txt"
    pd.DataFrame(
        {
            "Taxa": ["s2"],
            "Coverage-Of-All-Bams": ["0.4000"],
            "Bams": ["s2.bam"],
        }
    ).to_csv(coverage, sep="\t", index=False)
    pairs = pd.DataFrame(
        [
            {"sample1": "s1", "sample2": "s2"},
            {"sample1": "s2", "sample2": "s3"},
        ]
    )

    out = add_dedup_fields(pairs, dedup, load_taxa_coverage(coverage))

    assert out.loc[0, "coverage1"] == "5.0000"
    assert out.loc[0, "coverage2"] == "0.4000"
    assert out.loc[0, "coverage_min"] == pytest.approx(0.4)
    # s3 is in neither the dedup table nor the map, so the pair stays unknown.
    assert pd.isna(out.loc[1, "coverage_min"])


def test_plot_combine_ibs_low_coverage_uses_lower_depth(tmp_path):
    pairs = tmp_path / "pairs.tsv"
    depths = [
        (0.10, 5.00, 0.9900),
        (0.12, 4.50, 0.9902),
        (0.20, 4.00, 0.9910),
        (0.50, 3.00, 0.9940),
        (1.00, 2.50, 0.9980),
        (2.50, 6.00, 0.9994),
        (8.00, 9.00, 0.9996),
    ]
    rows = [
        {
            "sample1": f"a{index}",
            "sample2": f"b{index}",
            "pair_type": "within_watkins",
            "coverage1": first,
            "coverage2": second,
            "coverage_min": min(first, second),
            "hom_ibs": hom,
            "union_n": 48 if index < 3 else 4,
        }
        for index, (first, second, hom) in enumerate(depths)
    ]
    rows.append(
        {
            "sample1": "ipk1",
            "sample2": "ipk2",
            "pair_type": "within_ipk",
            "coverage1": 10.0,
            "coverage2": 10.0,
            "coverage_min": 10.0,
            "hom_ibs": 0.9999,
        }
    )
    rows.append(
        {
            "sample1": "nocov1",
            "sample2": "nocov2",
            "pair_type": "within_watkins",
            "coverage1": "",
            "coverage2": "",
            "coverage_min": "",
            "hom_ibs": 0.9990,
        }
    )
    pd.DataFrame(rows).to_csv(pairs, sep="\t", index=False)
    miss_rows = []
    for index in range(len(depths)):
        miss_rows.append({"Sample": f"a{index}", "F_MISS": 0.05 + 0.01 * index})
        miss_rows.append({"Sample": f"b{index}", "F_MISS": 0.07 + 0.01 * index})
    missing_rates = []
    for subgenome in "ABD":
        path = tmp_path / f"{subgenome}.sample.miss.rate.info.tsv"
        pd.DataFrame(miss_rows).to_csv(path, sep="\t", index=False)
        missing_rates.append(path)

    summary = plot_combine_ibs_low_coverage(
        pairs,
        tmp_path / "watkins_lowcov",
        pair_type="within_watkins",
        bins=(0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0, 20.0),
        highlight_union_n=48,
        missing_rate_paths=missing_rates,
        y_min=0.99,
    )
    assert (tmp_path / "watkins_lowcov.png").exists()
    assert (tmp_path / "watkins_lowcov.png").stat().st_size > 0
    assert (tmp_path / "watkins_lowcov.linear.png").exists()
    assert (tmp_path / "watkins_lowcov.linear.png").stat().st_size > 0
    for stem in ("high_missing_rate", "missing_sum"):
        for suffix in ("", ".linear"):
            path = tmp_path / f"watkins_lowcov.{stem}{suffix}.png"
            assert path.exists()
            assert path.stat().st_size > 0
    # The within_ipk pair and the pair without coverage do not enter the plot.
    assert summary["n_pairs"].sum() == len(depths)
    low_bin = summary.iloc[0]
    assert low_bin["coverage_low"] == pytest.approx(0.10)
    assert low_bin["n_pairs"] == 2
    assert low_bin["median_low_coverage"] == pytest.approx(0.11)
    assert low_bin["median_hom_ibs"] == pytest.approx(0.9901)


def test_plot_combine_ibs_low_coverage_fills_from_taxa_bam_map(tmp_path):
    pairs = tmp_path / "pairs.tsv"
    pd.DataFrame(
        [
            {
                "sample1": "a1",
                "sample2": "b1",
                "pair_type": "within_watkins",
                "coverage1": "",
                "coverage2": "",
                "hom_ibs": 0.999,
            },
            {
                "sample1": "a2",
                "sample2": "b2",
                "pair_type": "within_watkins",
                "coverage1": "",
                "coverage2": "",
                "hom_ibs": 0.995,
            },
        ]
    ).to_csv(pairs, sep="\t", index=False)
    coverage_map = tmp_path / "taxaBamMap.txt"
    pd.DataFrame(
        {
            "Taxa": ["a1", "b1", "a2", "b2"],
            "Coverage-Of-All-Bams": ["4.0", "0.5", "3.0", "1.5"],
            "Bams": ["a1.bam", "b1.bam", "a2.bam", "b2.bam"],
        }
    ).to_csv(coverage_map, sep="\t", index=False)

    summary = plot_combine_ibs_low_coverage(
        pairs,
        tmp_path / "filled",
        coverage_map_path=coverage_map,
        bins=(0.05, 1.0, 10.0),
    )
    # Without the map both pairs would be dropped; the map supplies both ends.
    assert summary["n_pairs"].sum() == 2
    assert summary.iloc[0]["median_low_coverage"] == pytest.approx(0.5)


def test_load_highlight_pairs_thresholds_and_filters(tmp_path):
    source = tmp_path / "inter_pairs.tsv"
    pd.DataFrame(
        [
            {"sample1": "a1", "sample2": "b1", "pair_type": "within_watkins", "hom_ibs": 0.9995},
            {"sample1": "a2", "sample2": "b2", "pair_type": "within_watkins", "hom_ibs": 0.998},
            {"sample1": "a3", "sample2": "b3", "pair_type": "within_ipk", "hom_ibs": 0.9999},
        ]
    ).to_csv(source, sep="\t", index=False)

    selected = _load_highlight_pairs(source, "hom_ibs", 0.999, pair_type="within_watkins")
    # Strictly above the threshold, same pair_type only, order-insensitive.
    assert selected == {("a1", "b1")}
    assert _load_highlight_pairs(source, "hom_ibs", 0.999, pair_type="all") == {
        ("a1", "b1"),
        ("a3", "b3"),
    }


def test_plot_combine_ibs_low_coverage_highlights_threshold_pairs(tmp_path):
    pairs = tmp_path / "pairs.tsv"
    pd.DataFrame(
        [
            {"sample1": "a1", "sample2": "b1", "pair_type": "within_watkins",
             "coverage1": 4.0, "coverage2": 5.0, "union_n": 48, "hom_ibs": 0.9995},
            {"sample1": "a2", "sample2": "b2", "pair_type": "within_watkins",
             "coverage1": 0.5, "coverage2": 0.6, "union_n": 2, "hom_ibs": 0.970},
            {"sample1": "a3", "sample2": "b3", "pair_type": "within_watkins",
             "coverage1": 3.0, "coverage2": 3.5, "union_n": 2, "hom_ibs": 0.9997},
        ]
    ).to_csv(pairs, sep="\t", index=False)
    highlight = tmp_path / "inter_pairs.tsv"
    pd.DataFrame(
        [
            {"sample1": "a1", "sample2": "b1", "pair_type": "within_watkins", "hom_ibs": 0.9995},
            {"sample1": "a2", "sample2": "b2", "pair_type": "within_watkins", "hom_ibs": 0.998},
        ]
    ).to_csv(highlight, sep="\t", index=False)

    summary = plot_combine_ibs_low_coverage(
        pairs,
        tmp_path / "highlighted",
        highlight_union_n=48,
        highlight_pairs_path=highlight,
        highlight_min_ibs=0.999,
        highlight_label="IBS>0.999 on test_common_inter",
        bins=(0.05, 1.0, 10.0),
    )
    assert summary["n_pairs"].sum() == 3
    assert (tmp_path / "highlighted.png").stat().st_size > 0
    assert (tmp_path / "highlighted.linear.png").stat().st_size > 0
