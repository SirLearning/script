"""Subgenome VEP atlas: unique sites, consequence groups, composition plots."""

from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.variant.coding_mask import VEP_META_PREFIX
from genetics.genomics.variant.vep_atlas import (
    CONSEQUENCE_GROUPS,
    add_rare_enrichment,
    attach_frequency,
    chromosome_region_intervals,
    collapse_vep_sites,
    consequence_group,
    frequency_class,
    impact_from_consequence,
    plot_common_rare_mirror,
    plot_vep_atlas_circos,
    plot_window_subgenome_violins,
    run_vep_atlas,
    subgenome_count_table,
)


def _write_vep(path: Path, rows: list[tuple[str, str, str, str]]):
    lines = [
        f"{VEP_META_PREFIX} ENSEMBL VARIANT EFFECT PREDICTOR v116.1\n",
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t"
        "Consequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\n",
    ]
    for var_id, location, gene, consequence in rows:
        lines.append(
            f"{var_id}\t{location}\tA\t{gene}\t{gene}.1\tTranscript\t"
            f"{consequence}\t-\t-\t-\t-\t-\n"
        )
    path.write_text("".join(lines), encoding="utf-8")


def test_collapse_vep_sites_keeps_one_row_with_the_most_severe_consequence(tmp_path: Path):
    vep = tmp_path / "A.vep.tsv"
    _write_vep(
        vep,
        [
            ("1-100", "1:100", "TraesCS1A02G000100", "intron_variant"),
            ("1-100", "1:100", "TraesCS1A02G000100", "stop_gained"),
            ("1-200", "1:200", "-", "intergenic_variant"),
        ],
    )
    sites = collapse_vep_sites(vep)
    assert list(sites["id"]) == ["1-100", "1-200"]
    assert list(sites["consequence"]) == ["stop_gained", "intergenic_variant"]
    assert list(sites["group"]) == ["predicted LoF", "intergenic"]
    assert list(sites["impact"]) == ["HIGH", "MODIFIER"]


def test_consequence_group_and_impact_cover_the_atlas_categories():
    assert consequence_group("splice_acceptor_variant") == "predicted LoF"
    assert consequence_group("splice_region_variant") == "splice-related"
    assert consequence_group("missense_variant") == "missense"
    assert consequence_group("3_prime_UTR_variant") == "UTR"
    assert consequence_group("upstream_gene_variant") == "upstream/downstream"
    assert all(name in CONSEQUENCE_GROUPS for name in (
        "predicted LoF", "missense", "synonymous", "UTR", "intergenic", "other",
    ))
    assert impact_from_consequence("frameshift_variant") == "HIGH"
    assert impact_from_consequence("missense_variant") == "MODERATE"
    assert impact_from_consequence("synonymous_variant") == "LOW"
    assert impact_from_consequence("intergenic_variant") == "MODIFIER"


def test_subgenome_count_table_normalises_by_reference_length_not_callable_bases():
    sites = {
        "A": pd.DataFrame({"id": ["1-1", "1-2"], "group": ["missense", "intergenic"]}),
        "B": pd.DataFrame({"id": ["2-1"], "group": ["synonymous"]}),
        "D": pd.DataFrame({"id": ["3-1", "3-2", "3-3"], "group": ["missense"] * 3}),
    }
    table = subgenome_count_table(sites)
    assert list(table.index) == ["A", "B", "D"]
    assert list(table["n_variants"]) == [2, 1, 3]
    assert (table["n_per_ref_gb"] > 0).all()
    assert table.loc["D", "n_per_ref_gb"] > table.loc["B", "n_per_ref_gb"]


def test_frequency_class_uses_mac_then_maf_bins():
    assert frequency_class(mac=1, maf=0.0003) == "singleton"
    assert frequency_class(mac=2, maf=0.0005) == "doubleton"
    assert frequency_class(mac=5, maf=0.002) == "MAC 3-10"
    assert frequency_class(mac=20, maf=0.0004) == "MAF <0.1%"
    assert frequency_class(mac=50, maf=0.03) == "1%-5%"
    assert frequency_class(mac=200, maf=0.08) == ">=5%"


def test_attach_frequency_joins_acount_and_afreq_on_variant_id(tmp_path: Path):
    sites = pd.DataFrame({"id": ["1-100", "1-200"]})
    afreq = tmp_path / "A.info.afreq"
    acount = tmp_path / "A.assoc.acount"
    afreq.write_text(
        "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n"
        "1\t1-100\tA\tT\t0.08\t2000\n"
        "1\t1-200\tA\tT\t0.002\t2000\n",
        encoding="utf-8",
    )
    acount.write_text(
        "#CHROM\tID\tREF\tALT\tALT_CTS\tOBS_CT\n"
        "1\t1-100\tA\tT\t160\t2000\n"
        "1\t1-200\tA\tT\t4\t2000\n",
        encoding="utf-8",
    )
    got = attach_frequency(sites, afreq, acount)
    assert list(got["freq_class"]) == [">=5%", "MAC 3-10"]
    assert got.loc[got["id"] == "1-100", "is_common"].iloc[0]
    assert got.loc[got["id"] == "1-200", "is_rare"].iloc[0]


def test_run_vep_atlas_writes_composition_and_heatmap(tmp_path: Path):
    genome = {}
    exon = {}
    afreq = {}
    acount = {}
    for sg, extra in (("A", "stop_gained"), ("B", "missense_variant"), ("D", "intergenic_variant")):
        genome[sg] = tmp_path / f"{sg}.genome.vep.tsv"
        exon[sg] = tmp_path / f"{sg}.exon.vep.tsv"
        _write_vep(
            genome[sg],
            [
                (f"{sg}-1", "1:100", "TraesCS1A02G000100", extra),
                (f"{sg}-2", "1:200", "-", "intergenic_variant"),
                (f"{sg}-3", "1:300", "TraesCS1A02G000200", "synonymous_variant"),
            ],
        )
        _write_vep(
            exon[sg],
            [
                (f"{sg}-e1", "1:110", "TraesCS1A02G000100", "missense_variant"),
                (f"{sg}-e2", "1:120", "TraesCS1A02G000100", "3_prime_UTR_variant"),
            ],
        )
        afreq[sg] = tmp_path / f"{sg}.info.afreq"
        acount[sg] = tmp_path / f"{sg}.assoc.acount"
        afreq[sg].write_text(
            "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n"
            f"1\t{sg}-1\tA\tT\t0.001\t2000\n"
            f"1\t{sg}-2\tA\tT\t0.10\t2000\n"
            f"1\t{sg}-3\tA\tT\t0.02\t2000\n",
            encoding="utf-8",
        )
        acount[sg].write_text(
            "#CHROM\tID\tREF\tALT\tALT_CTS\tOBS_CT\n"
            f"1\t{sg}-1\tA\tT\t2\t2000\n"
            f"1\t{sg}-2\tA\tT\t200\t2000\n"
            f"1\t{sg}-3\tA\tT\t40\t2000\n",
            encoding="utf-8",
        )

    prefix = tmp_path / "all.variant.vep_atlas"
    run_vep_atlas(
        {sg: str(path) for sg, path in genome.items()},
        {sg: str(path) for sg, path in exon.items()},
        str(prefix),
        genome_afreq={sg: str(path) for sg, path in afreq.items()},
        genome_acount={sg: str(path) for sg, path in acount.items()},
    )
    assert (tmp_path / "all.variant.vep_atlas.composition.png").is_file()
    assert (tmp_path / "all.variant.vep_atlas.freq_function.png").is_file()
    counts = pd.read_csv(tmp_path / "all.variant.vep_atlas.counts.tsv", sep="\t")
    assert set(counts["subgenome"]) == {"A", "B", "D"}
    assert int(counts["n_variants"].sum()) == 9


def test_run_vep_atlas_rare_lof_uses_exon_vep_not_genome(tmp_path: Path):
    genome = tmp_path / "A.genome.vep.tsv"
    exon = tmp_path / "A.exon.vep.tsv"
    _write_vep(genome, [("A-g1", "1:100", "TraesCS1A02G000100", "stop_gained")])
    _write_vep(
        exon,
        [
            ("A-e1", "1:110", "TraesCS1A02G000100", "stop_gained"),
            ("A-e2", "1:120", "TraesCS1A02G000100", "frameshift_variant"),
            ("A-e3", "1:130", "TraesCS1A02G000100", "missense_variant"),
        ],
    )
    genome_afreq = tmp_path / "A.genome.afreq"
    genome_acount = tmp_path / "A.genome.acount"
    exon_afreq = tmp_path / "A.exon.afreq"
    exon_acount = tmp_path / "A.exon.acount"
    genome_afreq.write_text(
        "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n1\tA-g1\tA\tT\t0.001\t2000\n",
        encoding="utf-8",
    )
    genome_acount.write_text(
        "#CHROM\tID\tREF\tALT\tALT_CTS\tOBS_CT\n1\tA-g1\tA\tT\t2\t2000\n",
        encoding="utf-8",
    )
    exon_afreq.write_text(
        "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n"
        "1\tA-e1\tA\tT\t0.001\t2000\n"
        "1\tA-e2\tA\tT\t0.002\t2000\n"
        "1\tA-e3\tA\tT\t0.10\t2000\n",
        encoding="utf-8",
    )
    exon_acount.write_text(
        "#CHROM\tID\tREF\tALT\tALT_CTS\tOBS_CT\n"
        "1\tA-e1\tA\tT\t2\t2000\n"
        "1\tA-e2\tA\tT\t4\t2000\n"
        "1\tA-e3\tA\tT\t200\t2000\n",
        encoding="utf-8",
    )
    empty = tmp_path / "empty.vep.tsv"
    _write_vep(empty, [])
    prefix = tmp_path / "atlas"
    run_vep_atlas(
        {"A": str(genome), "B": str(empty), "D": str(empty)},
        {"A": str(exon), "B": str(empty), "D": str(empty)},
        str(prefix),
        genome_afreq={"A": str(genome_afreq)},
        genome_acount={"A": str(genome_acount)},
        exon_afreq={"A": str(exon_afreq)},
        exon_acount={"A": str(exon_acount)},
    )
    counts = pd.read_csv(tmp_path / "atlas.counts.tsv", sep="\t").set_index("subgenome")
    assert int(counts.loc["A", "n_rare_functional"]) == 2


def test_window_mirror_and_violin_plots_write_pngs(tmp_path: Path):
    windows = pd.DataFrame(
        {
            "chrom": ["1A", "2A", "1B", "1D"],
            "bin_start": [0, 0, 0, 0],
            "bin_end": [10_000_000] * 4,
            "subgenome": ["A", "A", "B", "D"],
            "n_variants": [20, 10, 15, 8],
            "n_common": [12, 4, 6, 3],
            "n_rare": [5, 4, 7, 2],
            "n_rare_functional": [1, 0, 2, 0],
        }
    )
    windows = add_rare_enrichment(windows)
    assert "rare_enrichment" in windows.columns
    assert np.isfinite(windows["rare_enrichment"]).any()
    plot_common_rare_mirror(windows, str(tmp_path / "mirror.png"))
    plot_window_subgenome_violins(windows, str(tmp_path / "violin.png"))
    assert (tmp_path / "mirror.png").is_file()
    assert (tmp_path / "violin.png").is_file()


def test_chromosome_region_intervals_flank_iwgsc_centromere_with_pericentromere():
    from infra.wheat.ref_v1 import get_centromere_end, get_centromere_start, get_ref_chromosome_length

    regions = chromosome_region_intervals(flank_bp=50_000_000)
    one_a = regions.loc[regions["chrom"] == "1A"].sort_values("start_pos")
    labels = list(one_a["label"])
    assert labels[0] == "Arm"
    assert "Pericentromere" in labels
    assert "Centromere" in labels
    centro = one_a.loc[one_a["label"] == "Centromere"].iloc[0]
    assert int(centro["start_pos"]) == int(get_centromere_start("1A"))
    assert int(centro["end_pos"]) == int(get_centromere_end("1A"))
    peri = one_a.loc[one_a["label"] == "Pericentromere"]
    assert peri["start_pos"].min() == int(get_centromere_start("1A")) - 50_000_000
    assert peri["end_pos"].max() == int(get_centromere_end("1A")) + 50_000_000
    assert int(one_a["end_pos"].max()) == int(get_ref_chromosome_length("1A"))
    assert one_a["start_pos"].min() == 0


def test_plot_vep_atlas_circos_adds_ld_blocks_regions_and_gene_histogram(tmp_path: Path):
    windows = pd.DataFrame(
        {
            "chrom": ["1A", "1A", "1B"],
            "bin_start": [0, 10_000_000, 0],
            "bin_end": [10_000_000, 20_000_000, 10_000_000],
            "n_common": [12, 4, 6],
            "n_rare": [5, 4, 7],
            "n_rare_functional": [1, 0, 2],
        }
    )
    ld_blocks = pd.DataFrame(
        {
            "chrom": ["1A", "1B"],
            "start_pos": [1_000_000, 500_000],
            "end_pos": [3_000_000, 2_000_000],
        }
    )
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1A\t.\tgene\t100\t200\t.\t+\t.\tID=g1\n"
        "chr1A\t.\tgene\t15000000\t15000100\t.\t+\t.\tID=g2\n"
        "chr1B\t.\tgene\t50\t80\t.\t-\t.\tID=g3\n",
        encoding="utf-8",
    )
    out = tmp_path / "all.variant.vep_atlas.circos_windows.png"
    plot_vep_atlas_circos(
        windows,
        str(out),
        bin_bp=10_000_000,
        ld_blocks=ld_blocks,
        gene_gff=gff,
    )
    assert out.is_file()
    assert out.stat().st_size > 0


def test_run_vep_atlas_keep_ids_drop_sites_before_counts(tmp_path: Path):
    genome = tmp_path / "A.genome.vep.tsv"
    exon = tmp_path / "A.exon.vep.tsv"
    _write_vep(
        genome,
        [
            ("A-1", "1:100", "TraesCS1A02G000100", "stop_gained"),
            ("A-2", "1:200", "-", "intergenic_variant"),
        ],
    )
    _write_vep(
        exon,
        [
            ("A-e1", "1:110", "TraesCS1A02G000100", "stop_gained"),
            ("A-e2", "1:120", "TraesCS1A02G000100", "missense_variant"),
        ],
    )
    (tmp_path / "genome.afreq").write_text(
        "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n"
        "1\tA-1\tA\tT\t0.001\t2000\n"
        "1\tA-2\tA\tT\t0.10\t2000\n",
        encoding="utf-8",
    )
    (tmp_path / "exon.afreq").write_text(
        "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n"
        "1\tA-e1\tA\tT\t0.001\t2000\n"
        "1\tA-e2\tA\tT\t0.001\t2000\n",
        encoding="utf-8",
    )
    (tmp_path / "genome.acount").write_text(
        "#CHROM\tID\tREF\tALT\tALT_CTS\tOBS_CT\n"
        "1\tA-1\tA\tT\t2\t2000\n"
        "1\tA-2\tA\tT\t200\t2000\n",
        encoding="utf-8",
    )
    (tmp_path / "exon.acount").write_text(
        "#CHROM\tID\tREF\tALT\tALT_CTS\tOBS_CT\n"
        "1\tA-e1\tA\tT\t2\t2000\n"
        "1\tA-e2\tA\tT\t2\t2000\n",
        encoding="utf-8",
    )
    prefix = tmp_path / "all.variant.vep_atlas.cv_reliable.keep80"
    run_vep_atlas(
        {"A": str(genome)},
        {"A": str(exon)},
        str(prefix),
        genome_afreq={"A": str(tmp_path / "genome.afreq")},
        genome_acount={"A": str(tmp_path / "genome.acount")},
        exon_afreq={"A": str(tmp_path / "exon.afreq")},
        exon_acount={"A": str(tmp_path / "exon.acount")},
        genome_keep_ids={"A-1"},
        exon_keep_ids={"A-e1"},
    )
    counts = pd.read_csv(prefix.parent / "all.variant.vep_atlas.cv_reliable.keep80.counts.tsv", sep="\t")
    assert int(counts["n_variants"].sum()) == 1
    assert int(counts["n_rare_functional"].sum()) == 1


def test_rare_common_share_is_within_subgenome(tmp_path: Path):
    from genetics.genomics.variant.vep_atlas import (
        plot_subgenome_rare_common_share,
        rare_common_share_table,
    )

    counts = pd.DataFrame(
        {
            "subgenome": ["A", "B", "D"],
            "n_rare": [90, 80, 99],
            "n_common": [10, 20, 1],
        }
    )
    prop = rare_common_share_table(counts)
    assert list(prop.index) == ["A", "B", "D"]
    assert prop.loc["D", "Common"] == 0.01
    assert prop.loc["A", "Rare"] == 0.9
    out = tmp_path / "share.png"
    plot_subgenome_rare_common_share(counts, str(out))
    assert out.is_file()
    assert out.stat().st_size > 0
