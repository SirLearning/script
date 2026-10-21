"""Identity / duplicate metrics: rare-allele Jaccard, discordance, somalier-like."""

import numpy as np
import pandas as pd
import pytest

from genetics.genomics.sample.identity import (
    _plink1_path,
    ana_hom_ibs_one_vs_plink,
    ana_hom_snp_ibs_from_mibs,
    ana_ibs_sampled_pair_dist,
    cluster_undirected_pairs,
    compare_square_ibs,
    consensus_pairs,
    discord_rate_pairs,
    histogram_peak_centers,
    hom_snp_ibs_from_bed,
    hom_snp_ibs_from_gt,
    is_snp_alleles,
    sample_square_ibs_pairs,
    sample_square_ibs_pairs_by_label,
    ana_ibs_ploidy_pair_dist,
    ana_ibs_ploidy_group_pair_dist,
    sample_combined_subgenome_ibs_pairs,
    jaccard_from_carriers,
    nearest_neighbors_from_square,
    pairs_from_square_ibs,
    pairs_hom_ibs_one,
    somalier_relatedness_pairs,
    write_square_mibs,
)


def test_jaccard_identical_rare_carriers_is_one():
    carriers = {"s1": {"v1", "v2", "v3"}, "s2": {"v1", "v2", "v3"}, "s3": {"v9"}}
    out = jaccard_from_carriers(carriers, min_shared=1, min_jaccard=0.0)
    pair = out[(out["Sample1"] == "s1") & (out["Sample2"] == "s2")].iloc[0]
    assert pair["n_shared"] == 3
    assert pair["jaccard"] == pytest.approx(1.0)


def test_jaccard_omits_pairs_below_min_shared():
    carriers = {"s1": {"v1"}, "s2": {"v1", "v2", "v3"}}
    out = jaccard_from_carriers(carriers, min_shared=2, min_jaccard=0.0)
    assert out.empty


def test_discord_rate_identical_genotypes_is_zero():
    ids = ["a", "b", "c"]
    gt = np.array(
        [
            [0, 1, 2],
            [0, 1, 2],
            [2, 2, 0],
        ],
        dtype=float,
    )
    out = discord_rate_pairs(ids, gt, max_discord=1.0, min_nobs=1)
    ab = out[(out["Sample1"] == "a") & (out["Sample2"] == "b")].iloc[0]
    assert ab["n_obs"] == 3
    assert ab["n_diff"] == 0
    assert ab["discord_rate"] == pytest.approx(0.0)


def test_discord_rate_opposite_homozygotes():
    ids = ["a", "b"]
    gt = np.array([[0.0, 0.0, np.nan], [2.0, 2.0, 1.0]])
    out = discord_rate_pairs(ids, gt, max_discord=1.0, min_nobs=1)
    row = out.iloc[0]
    assert row["n_obs"] == 2
    assert row["n_diff"] == 2
    assert row["discord_rate"] == pytest.approx(1.0)


def test_somalier_relatedness_is_one_for_identical_hets():
    ids = ["a", "b"]
    gt = np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0]])
    out = somalier_relatedness_pairs(ids, gt, min_relatedness=0.0, min_nobs=1)
    row = out.iloc[0]
    assert row["shared_hets"] == 2
    assert row["ibs0"] == 0
    assert row["relatedness"] == pytest.approx(1.0)


def test_somalier_relatedness_drops_with_ibs0():
    ids = ["a", "b"]
    gt = np.array([[0.0, 1.0], [2.0, 1.0]])
    out = somalier_relatedness_pairs(ids, gt, min_relatedness=-2.0, min_nobs=1)
    assert out.iloc[0]["ibs0"] == 1
    assert out.iloc[0]["relatedness"] < 1.0


def test_cluster_undirected_pairs_connects_overlap():
    pairs = pd.DataFrame(
        {"Sample1": ["a", "b"], "Sample2": ["b", "c"], "score": [1.0, 1.0]}
    )
    clusters = cluster_undirected_pairs(pairs)
    assert len(clusters) == 1
    assert set(clusters[0]) == {"a", "b", "c"}


def test_consensus_pairs_requires_min_methods():
    king = pd.DataFrame({"Sample1": ["a"], "Sample2": ["b"]})
    ibd = pd.DataFrame({"Sample1": ["a"], "Sample2": ["b"]})
    jac = pd.DataFrame({"Sample1": ["a"], "Sample2": ["c"]})
    out = consensus_pairs({"king": king, "ibd": ibd, "jaccard": jac}, min_methods=2)
    assert list(out["n_methods"]) == [2]
    assert out.iloc[0]["Sample1"] == "a"
    assert out.iloc[0]["Sample2"] == "b"
    assert "king" in out.iloc[0]["methods"]
    assert "ibd" in out.iloc[0]["methods"]


def test_hom_snp_ibs_identical_homozygotes_is_one():
    ids = ["s1", "s2"]
    gt = np.array([[0.0, 2.0, 0.0], [0.0, 2.0, 0.0]])
    ibs, n_comp, n_match = hom_snp_ibs_from_gt(ids, gt)
    assert n_comp[0, 1] == 3
    assert n_match[0, 1] == 3
    assert ibs[0, 1] == pytest.approx(1.0)


def test_hom_snp_ibs_ignores_hets_and_missing():
    ids = ["s1", "s2"]
    # Site0 both hom match; site1 het vs hom (skip); site2 missing (skip);
    # site3 opposite homs (mismatch).
    gt = np.array([[0.0, 1.0, np.nan, 0.0], [0.0, 0.0, 2.0, 2.0]])
    ibs, n_comp, n_match = hom_snp_ibs_from_gt(ids, gt)
    assert n_comp[0, 1] == 2
    assert n_match[0, 1] == 1
    assert ibs[0, 1] == pytest.approx(0.5)


def test_hom_snp_ibs_opposite_homs_are_zero():
    ids = ["s1", "s2"]
    gt = np.array([[0.0, 0.0], [2.0, 2.0]])
    ibs, n_comp, _n_match = hom_snp_ibs_from_gt(ids, gt)
    assert n_comp[0, 1] == 2
    assert ibs[0, 1] == pytest.approx(0.0)


def test_is_snp_alleles_rejects_indels():
    assert is_snp_alleles("A", "T")
    assert not is_snp_alleles("A", "AT")
    assert not is_snp_alleles("-", "A")


def test_plink1_path_appends_after_dot_plink():
    assert str(_plink1_path("A_test.plink", ".fam")).endswith("A_test.plink.fam")
    assert str(_plink1_path("tiny", ".bed")).endswith("tiny.bed")


def test_pairs_from_square_ibs_keeps_high_pairs():
    ids = ["a", "b", "c"]
    ibs = np.array([[1.0, 0.99, 0.5], [0.99, 1.0, 0.8], [0.5, 0.8, 1.0]])
    n_comp = np.full((3, 3), 100)
    n_match = (ibs * n_comp).astype(int)
    out = pairs_from_square_ibs(ids, ibs, n_comp, n_match, min_ibs=0.9, min_nobs=10)
    assert len(out) == 1
    assert out.iloc[0]["Sample1"] == "a"
    assert out.iloc[0]["Sample2"] == "b"
    assert out.iloc[0]["hom_ibs"] == pytest.approx(0.99)


def test_compare_square_ibs_keeps_pairs_high_on_either():
    ids = ["a", "b", "c"]
    plink = np.array([[1.0, 0.96, 0.4], [0.96, 1.0, 0.99], [0.4, 0.99, 1.0]])
    hom = np.array([[1.0, 1.0, 0.4], [1.0, 1.0, 0.5], [0.4, 0.5, 1.0]])
    n_comp = np.full((3, 3), 50)
    out = compare_square_ibs(ids, plink, hom, n_comp, min_either=0.95)
    pairs = set(zip(out["Sample1"], out["Sample2"]))
    assert ("a", "b") in pairs
    assert ("b", "c") in pairs
    assert ("a", "c") not in pairs
    ab = out[(out["Sample1"] == "a") & (out["Sample2"] == "b")].iloc[0]
    assert ab["plink_ibs"] == pytest.approx(0.96)
    assert ab["hom_ibs"] == pytest.approx(1.0)


def test_sample_square_ibs_pairs_keeps_low_scores():
    ids = ["a", "b", "c"]
    plink = np.array([[1.0, 0.96, 0.40], [0.96, 1.0, 0.50], [0.40, 0.50, 1.0]])
    hom = np.array([[1.0, 0.97, 0.41], [0.97, 1.0, 0.52], [0.41, 0.52, 1.0]])
    out = sample_square_ibs_pairs(ids, plink, hom, n_pairs=10, seed=1)
    assert len(out) == 3
    ac = out[(out["Sample1"] == "a") & (out["Sample2"] == "c")].iloc[0]
    assert ac["plink_ibs"] == pytest.approx(0.40)
    capped = sample_square_ibs_pairs(ids, plink, hom, n_pairs=2, seed=1)
    assert len(capped) == 2


def test_nearest_neighbors_from_square_orders_by_score():
    ids = ["a", "b", "c"]
    mat = np.array([[1.0, 0.9, 0.7], [0.9, 1.0, 0.8], [0.7, 0.8, 1.0]])
    out = nearest_neighbors_from_square(ids, mat, k=1)
    a_nb = out[out["Sample"] == "a"].iloc[0]
    assert a_nb["Neighbor"] == "b"
    assert a_nb["rank"] == 1


def _write_tiny_plink_bed(tmp_path, ids, dosages, alleles):
    """Write a SNP-major PLINK1 bed/bim/fam for tests. dosages: (n_var, n_samples)."""
    prefix = tmp_path / "tiny"
    n_var, n_samples = dosages.shape
    n_bytes = (n_samples + 3) // 4
    body = np.zeros((n_var, n_bytes), dtype=np.uint8)
    for v in range(n_var):
        codes = np.full(n_samples, 1, dtype=np.uint8)
        row = dosages[v]
        codes[row == 0] = 0
        codes[row == 1] = 2
        codes[row == 2] = 3
        for i, code in enumerate(codes):
            body[v, i // 4] |= np.uint8(int(code) << (2 * (i % 4)))
    magic = np.array([0x6C, 0x1B, 0x01], dtype=np.uint8)
    prefix.with_suffix(".bed").write_bytes(magic.tobytes() + body.tobytes())
    bim_lines = []
    for i, (a1, a2) in enumerate(alleles):
        bim_lines.append(f"1\tsnp{i}\t0\t{i + 1}\t{a1}\t{a2}\n")
    prefix.with_suffix(".bim").write_text("".join(bim_lines), encoding="utf-8")
    fam_lines = [f"{sid}\t{sid}\t0\t0\t0\t-9\n" for sid in ids]
    prefix.with_suffix(".fam").write_text("".join(fam_lines), encoding="utf-8")
    return str(prefix)


def test_hom_snp_ibs_from_bed_skips_indels(tmp_path):
    ids = ["s1", "s2"]
    # SNP0 match; indel would mismatch if counted; SNP1 opposite homs.
    dosages = np.array(
        [
            [0.0, 0.0],
            [0.0, 2.0],
            [2.0, 0.0],
        ],
        dtype=float,
    )
    alleles = [("A", "T"), ("A", "AT"), ("C", "G")]
    prefix = _write_tiny_plink_bed(tmp_path, ids, dosages, alleles)
    got_ids, ibs, n_comp, n_match = hom_snp_ibs_from_bed(prefix, chunk_sites=1)
    assert got_ids == ids
    assert n_comp[0, 1] == 2
    assert n_match[0, 1] == 1
    assert ibs[0, 1] == pytest.approx(0.5)


def test_ana_hom_snp_ibs_from_mibs_writes_compare_tables(tmp_path):
    ids = ["a", "b", "c"]
    hom = np.array([[1.0, 1.0, 0.4], [1.0, 1.0, 0.5], [0.4, 0.5, 1.0]])
    plink = np.array([[1.0, 0.96, 0.4], [0.96, 1.0, 0.99], [0.4, 0.99, 1.0]])
    write_square_mibs(tmp_path / "hom.mibs", tmp_path / "hom.mibs.id", ids, hom)
    write_square_mibs(tmp_path / "plink.mibs", tmp_path / "plink.mibs.id", ids, plink)
    prefix = tmp_path / "out"
    pairs = ana_hom_snp_ibs_from_mibs(
        str(tmp_path / "hom.mibs"),
        str(tmp_path / "hom.mibs.id"),
        str(prefix),
        plink_mibs=str(tmp_path / "plink.mibs"),
        plink_mibs_id=str(tmp_path / "plink.mibs.id"),
        min_pair_ibs=0.99,
        min_compare_ibs=0.95,
        min_nobs=0,
    )
    assert not pairs.empty
    cmp = pd.read_csv(f"{prefix}.vs_plink.pairs.tsv", sep="\t")
    assert {"plink_ibs", "hom_ibs", "ibs_delta"}.issubset(cmp.columns)
    assert (tmp_path / "out.vs_plink.scatter.png").exists()
    assert (tmp_path / "out.dist.png").exists()


def test_pairs_hom_ibs_one_keeps_only_exact_one_and_plink_score():
    ids = ["a", "b", "c"]
    hom = np.array([[1.0, 1.0, 0.99], [1.0, 1.0, 1.0], [0.99, 1.0, 1.0]])
    plink = np.array([[1.0, 0.96, 0.90], [0.96, 1.0, 0.99], [0.90, 0.99, 1.0]])
    out = pairs_hom_ibs_one(ids, hom, plink)
    pairs = {(r.Sample1, r.Sample2): r for r in out.itertuples()}
    assert set(pairs) == {("a", "b"), ("b", "c")}
    assert pairs[("a", "b")].hom_ibs == pytest.approx(1.0)
    assert pairs[("a", "b")].plink_ibs == pytest.approx(0.96)
    assert pairs[("a", "b")].ibs_delta == pytest.approx(0.04)
    assert pairs[("b", "c")].plink_ibs == pytest.approx(0.99)
    assert ("a", "c") not in pairs


def test_ana_hom_ibs_one_vs_plink_writes_pair_table_and_plink_dist(tmp_path):
    ids = ["a", "b", "c"]
    hom = np.array([[1.0, 1.0, 0.4], [1.0, 1.0, 1.0], [0.4, 1.0, 1.0]])
    plink = np.array([[1.0, 0.96, 0.4], [0.96, 1.0, 0.88], [0.4, 0.88, 1.0]])
    write_square_mibs(tmp_path / "hom.mibs", tmp_path / "hom.mibs.id", ids, hom)
    write_square_mibs(tmp_path / "plink.mibs", tmp_path / "plink.mibs.id", ids, plink)
    prefix = tmp_path / "out"
    out = ana_hom_ibs_one_vs_plink(
        str(tmp_path / "hom.mibs"),
        str(tmp_path / "hom.mibs.id"),
        str(tmp_path / "plink.mibs"),
        str(tmp_path / "plink.mibs.id"),
        str(prefix),
    )
    assert len(out) == 2
    table = pd.read_csv(f"{prefix}.eq1.pairs.tsv", sep="\t")
    assert {"Sample1", "Sample2", "hom_ibs", "plink_ibs", "ibs_delta"}.issubset(table.columns)
    assert set(zip(table["Sample1"], table["Sample2"])) == {("a", "b"), ("b", "c")}
    assert (tmp_path / "out.eq1.plink.dist.png").exists()
    th = pd.read_csv(f"{prefix}.eq1.th.tsv", sep="\t")
    assert int(th.loc[0, "n_pairs"]) == 2
    assert float(th.loc[0, "min_plink_ibs"]) == pytest.approx(0.88)


def test_histogram_peak_centers_finds_two_separated_modes():
    rng = np.random.default_rng(0)
    scores = np.concatenate(
        [rng.normal(0.85, 0.008, 4000), rng.normal(0.99, 0.002, 4000)]
    )
    peaks = histogram_peak_centers(scores, n_bins=60, min_rel_height=0.15)
    assert len(peaks) == 2
    assert min(peaks) < 0.90
    assert max(peaks) > 0.96


def test_histogram_peak_centers_ignores_a_tiny_high_tail():
    rng = np.random.default_rng(0)
    scores = np.concatenate([rng.normal(0.95, 0.008, 20000), np.full(15, 0.999)])
    peaks = histogram_peak_centers(scores, n_bins=50, min_rel_height=0.05)
    assert len(peaks) == 1
    assert abs(peaks[0] - 0.95) < 0.04


def test_ana_ibs_sampled_pair_dist_plots_full_range(tmp_path):
    table = pd.DataFrame(
        {
            "Sample1": ["a", "a", "b"],
            "Sample2": ["b", "c", "c"],
            "plink_ibs": [0.70, 0.95, 0.99],
            "hom_ibs": [0.72, 0.96, 1.00],
        }
    )
    pairs = tmp_path / "pairs.tsv"
    table.to_csv(pairs, sep="\t", index=False)
    prefix = tmp_path / "out"
    out = ana_ibs_sampled_pair_dist(str(pairs), str(prefix))
    assert len(out) == 3
    assert (tmp_path / "out.plink.dist.png").is_file()
    assert (tmp_path / "out.hom.dist.png").is_file()
    assert (tmp_path / "out.plink.dist_log.png").is_file()
    assert (tmp_path / "out.hom.dist_log.png").is_file()
    th = pd.read_csv(f"{prefix}.th.tsv", sep="\t")
    assert float(th.loc[0, "min_plink_ibs"]) == pytest.approx(0.70)
    assert float(th.loc[0, "min_hom_ibs"]) == pytest.approx(0.72)
    assert int(th.loc[0, "n_pairs"]) == 3


def test_sample_square_ibs_pairs_by_label_keeps_within_class_pairs_only():
    ids = ["a1", "a2", "h1", "h2"]
    plink = np.array(
        [
            [1.00, 0.99, 0.80, 0.81],
            [0.99, 1.00, 0.82, 0.83],
            [0.80, 0.82, 1.00, 0.95],
            [0.81, 0.83, 0.95, 1.00],
        ]
    )
    hom = plink + 0.01
    np.fill_diagonal(hom, 1.0)
    labels = {"a1": "AA", "a2": "AA", "h1": "AABBDD", "h2": "AABBDD"}
    out = sample_square_ibs_pairs_by_label(
        ids,
        plink,
        hom,
        labels,
        allowed_labels=("AA", "AABBDD"),
        n_pairs_per_class=100,
        seed=1,
    )
    pairs = set(zip(out["Sample1"], out["Sample2"]))
    assert pairs == {("a1", "a2"), ("h1", "h2")}
    aa = out[out["ploidy"] == "AA"].iloc[0]
    assert aa["plink_ibs"] == pytest.approx(0.99)
    hex_row = out[out["ploidy"] == "AABBDD"].iloc[0]
    assert hex_row["plink_ibs"] == pytest.approx(0.95)
    assert hex_row["hom_ibs"] == pytest.approx(0.96)


def test_ana_ibs_ploidy_pair_dist_writes_faceted_plots(tmp_path):
    ids = ["a1", "a2", "h1", "h2"]
    plink = np.array(
        [
            [1.00, 0.99, 0.80, 0.81],
            [0.99, 1.00, 0.82, 0.83],
            [0.80, 0.82, 1.00, 0.95],
            [0.81, 0.83, 0.95, 1.00],
        ]
    )
    hom = plink + 0.01
    np.fill_diagonal(hom, 1.0)
    write_square_mibs(tmp_path / "plink.mibs", tmp_path / "plink.mibs.id", ids, plink)
    write_square_mibs(tmp_path / "hom.mibs", tmp_path / "hom.mibs.id", ids, hom)
    groups = tmp_path / "sample_group.txt"
    groups.write_text("a1 A\na2 A\nh1 ABD\nh2 ABD\n", encoding="utf-8")
    seg_rows = []
    for sample in ("a1", "a2"):
        for homoe in range(1, 8):
            for sg in ("B", "D"):
                seg_rows.append(
                    {
                        "sample": sample,
                        "logical_chr": f"chr{homoe}{sg}",
                        "call_type": "subgenome_absent",
                    }
                )
    seg_path = tmp_path / "breachr_segments_all.tsv"
    pd.DataFrame(seg_rows).to_csv(seg_path, sep="\t", index=False)
    prefix = tmp_path / "A.sample.identity.ploidy_ibs"
    out = ana_ibs_ploidy_pair_dist(
        str(tmp_path / "plink.mibs"),
        str(tmp_path / "plink.mibs.id"),
        str(tmp_path / "hom.mibs"),
        str(tmp_path / "hom.mibs.id"),
        str(groups),
        str(seg_path),
        str(prefix),
        subgenome="A",
    )
    assert set(out["ploidy"]) == {"AA", "AABBDD"}
    assert (tmp_path / "A.sample.identity.ploidy_ibs.plink.dist.png").is_file()
    assert (tmp_path / "A.sample.identity.ploidy_ibs.hom.dist.png").is_file()
    assert (tmp_path / "A.sample.identity.ploidy_ibs.plink.dist_log.png").is_file()
    assert (tmp_path / "A.sample.identity.ploidy_ibs.hom.dist_log.png").is_file()
    assert (tmp_path / "A.sample.identity.ploidy_ibs.plink.density.png").is_file()
    assert (tmp_path / "A.sample.identity.ploidy_ibs.hom.density.png").is_file()
    samples = pd.read_csv(f"{prefix}.samples.tsv", sep="\t")
    assert set(zip(samples["Sample"], samples["ploidy"])) == {
        ("a1", "AA"),
        ("a2", "AA"),
        ("h1", "AABBDD"),
        ("h2", "AABBDD"),
    }


def test_sample_combined_subgenome_ibs_pairs_means_required_subgenomes():
    a_ids = ["aa1", "aa2", "t1", "t2", "h1", "h2"]
    a_plink = np.full((6, 6), 0.50)
    np.fill_diagonal(a_plink, 1.0)
    a_plink[0, 1] = a_plink[1, 0] = 0.99
    a_plink[2, 3] = a_plink[3, 2] = 0.90
    a_plink[4, 5] = a_plink[5, 4] = 0.70
    a_hom = a_plink + 0.01
    np.fill_diagonal(a_hom, 1.0)

    b_ids = ["ss1", "ss2", "t1", "t2", "h1", "h2"]
    b_plink = np.full((6, 6), 0.50)
    np.fill_diagonal(b_plink, 1.0)
    b_plink[0, 1] = b_plink[1, 0] = 0.88
    b_plink[2, 3] = b_plink[3, 2] = 0.80
    b_plink[4, 5] = b_plink[5, 4] = 0.80
    b_hom = b_plink + 0.01
    np.fill_diagonal(b_hom, 1.0)

    d_ids = ["dd1", "dd2", "h1", "h2"]
    d_plink = np.full((4, 4), 0.50)
    np.fill_diagonal(d_plink, 1.0)
    d_plink[0, 1] = d_plink[1, 0] = 0.77
    d_plink[2, 3] = d_plink[3, 2] = 0.90
    d_hom = d_plink + 0.01
    np.fill_diagonal(d_hom, 1.0)

    labels = {
        "aa1": "AA",
        "aa2": "AA",
        "ss1": "SS",
        "ss2": "SS",
        "dd1": "DD",
        "dd2": "DD",
        "t1": "AABB",
        "t2": "AABB",
        "h1": "AABBDD",
        "h2": "AABBDD",
    }
    out = sample_combined_subgenome_ibs_pairs(
        {
            "A": (a_ids, a_plink, a_hom),
            "B": (b_ids, b_plink, b_hom),
            "D": (d_ids, d_plink, d_hom),
        },
        labels,
        n_pairs_per_class=100,
        seed=1,
    )
    by_class = {row["ploidy"]: row for row in out.to_dict(orient="records")}
    assert set(by_class) == {"AA", "SS", "DD", "AABB", "AABBDD"}
    assert by_class["AA"]["plink_ibs"] == pytest.approx(0.99)
    assert by_class["SS"]["plink_ibs"] == pytest.approx(0.88)
    assert by_class["DD"]["plink_ibs"] == pytest.approx(0.77)
    assert by_class["AABB"]["plink_ibs"] == pytest.approx(0.85)
    assert by_class["AABB"]["plink_ibs_A"] == pytest.approx(0.90)
    assert by_class["AABB"]["plink_ibs_B"] == pytest.approx(0.80)
    assert by_class["AABBDD"]["plink_ibs"] == pytest.approx(0.80)
    assert by_class["AABBDD"]["hom_ibs"] == pytest.approx(0.81)


def test_ana_ibs_ploidy_group_pair_dist_writes_combined_plots(tmp_path):
    def _write(sg, ids, plink):
        hom = plink + 0.01
        np.fill_diagonal(hom, 1.0)
        write_square_mibs(tmp_path / f"{sg}.plink.mibs", tmp_path / f"{sg}.plink.mibs.id", ids, plink)
        write_square_mibs(tmp_path / f"{sg}.hom.mibs", tmp_path / f"{sg}.hom.mibs.id", ids, hom)

    a_ids = ["aa1", "aa2", "t1", "t2", "h1", "h2"]
    a_plink = np.full((6, 6), 0.50)
    np.fill_diagonal(a_plink, 1.0)
    a_plink[0, 1] = a_plink[1, 0] = 0.99
    a_plink[2, 3] = a_plink[3, 2] = 0.90
    a_plink[4, 5] = a_plink[5, 4] = 0.70
    _write("A", a_ids, a_plink)

    b_ids = ["ss1", "ss2", "t1", "t2", "h1", "h2"]
    b_plink = np.full((6, 6), 0.50)
    np.fill_diagonal(b_plink, 1.0)
    b_plink[0, 1] = b_plink[1, 0] = 0.88
    b_plink[2, 3] = b_plink[3, 2] = 0.80
    b_plink[4, 5] = b_plink[5, 4] = 0.80
    _write("B", b_ids, b_plink)

    d_ids = ["dd1", "dd2", "h1", "h2"]
    d_plink = np.full((4, 4), 0.50)
    np.fill_diagonal(d_plink, 1.0)
    d_plink[0, 1] = d_plink[1, 0] = 0.77
    d_plink[2, 3] = d_plink[3, 2] = 0.90
    _write("D", d_ids, d_plink)

    groups = tmp_path / "sample_group.txt"
    groups.write_text(
        "aa1 A\naa2 A\nss1 S\nss2 S\ndd1 D\ndd2 D\nt1 AB\nt2 AB\nh1 ABD\nh2 ABD\n",
        encoding="utf-8",
    )
    seg_rows = []
    for sample in ("aa1", "aa2"):
        for homoe in range(1, 8):
            for sg in ("B", "D"):
                seg_rows.append(
                    {
                        "sample": sample,
                        "logical_chr": f"chr{homoe}{sg}",
                        "call_type": "subgenome_absent",
                    }
                )
    for sample in ("ss1", "ss2"):
        for homoe in range(1, 8):
            for sg in ("A", "D"):
                seg_rows.append(
                    {
                        "sample": sample,
                        "logical_chr": f"chr{homoe}{sg}",
                        "call_type": "subgenome_absent",
                    }
                )
    for sample in ("dd1", "dd2"):
        for homoe in range(1, 8):
            for sg in ("A", "B"):
                seg_rows.append(
                    {
                        "sample": sample,
                        "logical_chr": f"chr{homoe}{sg}",
                        "call_type": "subgenome_absent",
                    }
                )
    for sample in ("t1", "t2"):
        for homoe in range(1, 8):
            seg_rows.append(
                {
                    "sample": sample,
                    "logical_chr": f"chr{homoe}D",
                    "call_type": "subgenome_absent",
                }
            )
    seg_path = tmp_path / "breachr_segments_all.tsv"
    pd.DataFrame(seg_rows).to_csv(seg_path, sep="\t", index=False)
    prefix = tmp_path / "all.sample.identity.ploidy_group_ibs"
    out = ana_ibs_ploidy_group_pair_dist(
        str(tmp_path / "A.plink.mibs"),
        str(tmp_path / "A.plink.mibs.id"),
        str(tmp_path / "A.hom.mibs"),
        str(tmp_path / "A.hom.mibs.id"),
        str(tmp_path / "B.plink.mibs"),
        str(tmp_path / "B.plink.mibs.id"),
        str(tmp_path / "B.hom.mibs"),
        str(tmp_path / "B.hom.mibs.id"),
        str(tmp_path / "D.plink.mibs"),
        str(tmp_path / "D.plink.mibs.id"),
        str(tmp_path / "D.hom.mibs"),
        str(tmp_path / "D.hom.mibs.id"),
        str(groups),
        str(seg_path),
        str(prefix),
    )
    assert set(out["ploidy"]) == {"AA", "SS", "DD", "AABB", "AABBDD"}
    aabb = out[out["ploidy"] == "AABB"].iloc[0]
    assert aabb["plink_ibs"] == pytest.approx(0.85)
    hex_row = out[out["ploidy"] == "AABBDD"].iloc[0]
    assert hex_row["plink_ibs"] == pytest.approx(0.80)
    assert (tmp_path / "all.sample.identity.ploidy_group_ibs.plink.dist.png").is_file()
    assert (tmp_path / "all.sample.identity.ploidy_group_ibs.hom.dist.png").is_file()
    assert (tmp_path / "all.sample.identity.ploidy_group_ibs.plink.density.png").is_file()
    samples = pd.read_csv(f"{prefix}.samples.tsv", sep="\t")
    assert set(zip(samples["Sample"], samples["ploidy"])) >= {
        ("aa1", "AA"),
        ("ss1", "SS"),
        ("dd1", "DD"),
        ("t1", "AABB"),
        ("h1", "AABBDD"),
    }
