"""Rare/common density and r2 / D-prime LD blocks for a three-ring Circos map."""

from pathlib import Path

import pandas as pd

from genetics.genomics.variant.circos_ld_map import (
    COMMON_TRACK_LABEL,
    LD_BLOCK_TRACK_LABEL,
    RARE_TRACK_LABEL,
    ana_circos_ld_map,
    bin_variant_density,
    ld_block_coverage_intervals,
    ld_blocks_from_pairs,
    load_sites_from_pvar_afreq,
    merge_adjacent_blocks_by_dprime_q1,
    plot_circos_ld_map_from_tables,
)


PVAR_HEADER = "#CHROM\tPOS\tID\tREF\tALT\n"
AFREQ_HEADER = "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n"
LD_HEADER = "CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2 DP\n"


def test_bin_variant_density_splits_rare_and_common_per_mb():
    sites = pd.DataFrame(
        {
            "chrom": ["1A", "1A", "1A", "1B"],
            "pos": [1_000, 1_500_000, 1_600_000, 100],
            "maf": [0.01, 0.20, 0.04, 0.30],
        }
    )
    bins = bin_variant_density(sites, bin_bp=1_000_000, maf_threshold=0.05)
    a0 = bins[(bins["chrom"] == "1A") & (bins["bin_start"] == 0)].iloc[0]
    a1 = bins[(bins["chrom"] == "1A") & (bins["bin_start"] == 1_000_000)].iloc[0]
    b0 = bins[(bins["chrom"] == "1B") & (bins["bin_start"] == 0)].iloc[0]
    assert int(a0["n_rare"]) == 1
    assert int(a0["n_common"]) == 0
    assert int(a1["n_rare"]) == 1
    assert int(a1["n_common"]) == 1
    assert int(b0["n_rare"]) == 0
    assert int(b0["n_common"]) == 1


def test_ld_blocks_merge_adjacent_when_dprime_q1_exceeds_0_98():
    snps = pd.DataFrame(
        {
            "chrom": ["1A", "1A", "1A", "1A"],
            "id": ["s1", "s2", "s3", "s4"],
            "pos": [100, 200, 300, 500_000],
        }
    )
    pairs = pd.DataFrame(
        {
            "id_a": ["s1", "s2", "s3"],
            "id_b": ["s2", "s3", "s4"],
            "r2": [0.90, 0.10, 0.05],
            "dprime": [0.995, 0.992, 0.20],
        }
    )
    initial = ld_blocks_from_pairs(snps, pairs, r2_threshold=0.8)
    assert list(zip(initial["start_id"], initial["end_id"])) == [("s1", "s2"), ("s3", "s3"), ("s4", "s4")]
    merged = merge_adjacent_blocks_by_dprime_q1(initial, snps, pairs, q1_threshold=0.98)
    assert list(zip(merged["start_id"], merged["end_id"])) == [("s1", "s3"), ("s4", "s4")]
    assert int(merged.iloc[0]["end_pos"]) == 300
    assert int(merged.iloc[0]["n_snps"]) == 3


def test_ana_circos_ld_map_writes_three_ring_png(tmp_path: Path):
    pvar = tmp_path / "A_test.plink2.pvar"
    afreq = tmp_path / "A.info.afreq"
    ld = tmp_path / "A.info.ld_block.ld"
    pvar.write_text(
        PVAR_HEADER
        + "chr1A\t1000\ts1\tA\tT\n"
        + "chr1A\t2000\ts2\tA\tT\n"
        + "chr1A\t3000\ts3\tA\tT\n"
        + "chr1B\t500\ts4\tA\tT\n"
        + "chr1B\t1500000\ts5\tA\tT\n",
        encoding="utf-8",
    )
    afreq.write_text(
        AFREQ_HEADER
        + "1A\ts1\tA\tT\t0.01\t100\n"
        + "1A\ts2\tA\tT\t0.02\t100\n"
        + "1A\ts3\tA\tT\t0.20\t100\n"
        + "1B\ts4\tA\tT\t0.30\t100\n"
        + "1B\ts5\tA\tT\t0.04\t100\n",
        encoding="utf-8",
    )
    ld.write_text(
        LD_HEADER
        + "chr1A 1000 s1 chr1A 2000 s2 0.91 0.99\n"
        + "chr1A 2000 s2 chr1A 3000 s3 0.12 0.993\n"
        + "chr1B 500 s4 chr1B 1500000 s5 0.04 0.10\n",
        encoding="utf-8",
    )
    prefix = str(tmp_path / "all.variant.circos_ld_map")
    ana_circos_ld_map(
        pvar_files=[str(pvar)],
        afreq_files=[str(afreq)],
        ld_files=[str(ld)],
        output_prefix=prefix,
        bin_bp=1_000_000,
        chrom_lengths={"1A": 5_000_000, "1B": 4_000_000},
    )
    sites = load_sites_from_pvar_afreq(str(pvar), str(afreq))
    assert set(sites["chrom"]) == {"1A", "1B"}
    assert int((sites["maf"] < 0.05).sum()) == 3
    blocks = pd.read_csv(f"{prefix}.blocks.info.tsv", sep="\t")
    density = pd.read_csv(f"{prefix}.density.info.tsv", sep="\t")
    assert {"1A", "1B"} <= set(blocks["chrom"])
    assert (blocks["end_pos"] >= blocks["start_pos"]).all()
    assert {"n_rare", "n_common"} <= set(density.columns)
    assert Path(f"{prefix}.png").is_file()
    assert Path(f"{prefix}.png").stat().st_size > 1000


def test_ld_blocks_match_when_pair_orientation_is_reversed():
    snps = pd.DataFrame(
        {
            "chrom": ["1A", "1A", "1A"],
            "id": ["s1", "s2", "s3"],
            "pos": [100, 200, 300],
        }
    )
    forward = pd.DataFrame(
        {
            "id_a": ["s1", "s2"],
            "id_b": ["s2", "s3"],
            "r2": [0.90, 0.10],
            "dprime": [0.99, 0.50],
        }
    )
    reversed_pairs = forward.rename(columns={"id_a": "id_b", "id_b": "id_a"})
    fwd_blocks = ld_blocks_from_pairs(snps, forward, r2_threshold=0.8)
    rev_blocks = ld_blocks_from_pairs(snps, reversed_pairs, r2_threshold=0.8)
    assert list(zip(fwd_blocks["start_id"], fwd_blocks["end_id"])) == list(
        zip(rev_blocks["start_id"], rev_blocks["end_id"])
    )


def test_ana_circos_ld_map_parallel_matches_serial(tmp_path: Path):
    def write_sg(sg: str, chrom: str, rows: list[tuple]):
        pvar = tmp_path / f"{sg}_test.plink2.pvar"
        afreq = tmp_path / f"{sg}.info.afreq"
        ld = tmp_path / f"{sg}.info.ld_block.ld"
        pvar.write_text(
            PVAR_HEADER + "".join(f"chr{chrom}\t{pos}\t{vid}\tA\tT\n" for pos, vid, _maf, _r2, _dp in rows),
            encoding="utf-8",
        )
        afreq.write_text(
            AFREQ_HEADER
            + "".join(f"{chrom}\t{vid}\tA\tT\t{maf}\t100\n" for _pos, vid, maf, _r2, _dp in rows),
            encoding="utf-8",
        )
        ld_lines = [LD_HEADER]
        for i in range(len(rows) - 1):
            pos_a, id_a, _m1, r2, dp = rows[i]
            pos_b, id_b, _m2, _r2b, _dpb = rows[i + 1]
            ld_lines.append(f"chr{chrom} {pos_a} {id_a} chr{chrom} {pos_b} {id_b} {r2} {dp}\n")
        ld.write_text("".join(ld_lines), encoding="utf-8")
        return str(pvar), str(afreq), str(ld)

    a_files = write_sg(
        "A",
        "1A",
        [(1000, "s1", 0.01, 0.91, 0.99), (2000, "s2", 0.02, 0.12, 0.993), (3000, "s3", 0.20, 0.0, 0.0)],
    )
    b_files = write_sg(
        "B",
        "1B",
        [(500, "s4", 0.30, 0.04, 0.10), (1_500_000, "s5", 0.04, 0.0, 0.0)],
    )
    pvars, afreqs, lds = zip(a_files, b_files)
    kwargs = dict(
        pvar_files=list(pvars),
        afreq_files=list(afreqs),
        ld_files=list(lds),
        bin_bp=1_000_000,
        chrom_lengths={"1A": 5_000_000, "1B": 4_000_000},
    )
    serial_prefix = str(tmp_path / "serial.variant.circos_ld_map")
    parallel_prefix = str(tmp_path / "parallel.variant.circos_ld_map")
    ana_circos_ld_map(**kwargs, output_prefix=serial_prefix, n_jobs=1)
    ana_circos_ld_map(**kwargs, output_prefix=parallel_prefix, n_jobs=2)
    serial = pd.read_csv(f"{serial_prefix}.blocks.info.tsv", sep="\t")
    parallel = pd.read_csv(f"{parallel_prefix}.blocks.info.tsv", sep="\t")
    serial = serial.sort_values(["chrom", "start_pos", "end_pos"]).reset_index(drop=True)
    parallel = parallel.sort_values(["chrom", "start_pos", "end_pos"]).reset_index(drop=True)
    pd.testing.assert_frame_equal(serial, parallel, check_dtype=False)


def test_circos_track_labels_omit_maf_thresholds():
    for label in (RARE_TRACK_LABEL, COMMON_TRACK_LABEL, LD_BLOCK_TRACK_LABEL):
        assert "MAF" not in label
        assert "(" not in label


def test_ana_circos_ld_map_density_and_blocks_use_separate_pvars(tmp_path: Path):
    dens_pvar = tmp_path / "thin.pvar"
    dens_afreq = tmp_path / "thin.afreq"
    ld_pvar = tmp_path / "common.pvar"
    ld = tmp_path / "common.ld"
    dens_pvar.write_text(
        PVAR_HEADER
        + "chr1A\t1000\trare1\tA\tT\n"
        + "chr1A\t2000\trare2\tA\tT\n"
        + "chr1A\t3000000\tcom1\tA\tT\n",
        encoding="utf-8",
    )
    dens_afreq.write_text(
        AFREQ_HEADER
        + "1A\trare1\tA\tT\t0.01\t100\n"
        + "1A\trare2\tA\tT\t0.02\t100\n"
        + "1A\tcom1\tA\tT\t0.20\t100\n",
        encoding="utf-8",
    )
    ld_pvar.write_text(
        PVAR_HEADER
        + "chr1A\t10000\tc1\tA\tT\n"
        + "chr1A\t20000\tc2\tA\tT\n"
        + "chr1A\t400000\tc3\tA\tT\n",
        encoding="utf-8",
    )
    ld.write_text(
        LD_HEADER
        + "chr1A 10000 c1 chr1A 20000 c2 0.91 0.99\n"
        + "chr1A 20000 c2 chr1A 400000 c3 0.05 0.20\n",
        encoding="utf-8",
    )
    prefix = str(tmp_path / "hybrid.variant.circos_ld_map")
    ana_circos_ld_map(
        pvar_files=[str(dens_pvar)],
        afreq_files=[str(dens_afreq)],
        ld_files=[str(ld)],
        ld_pvar_files=[str(ld_pvar)],
        output_prefix=prefix,
        bin_bp=1_000_000,
        chrom_lengths={"1A": 5_000_000},
        n_jobs=1,
    )
    density = pd.read_csv(f"{prefix}.density.info.tsv", sep="\t")
    blocks = pd.read_csv(f"{prefix}.blocks.info.tsv", sep="\t")
    assert int(density["n_rare"].sum()) == 2
    assert int(density["n_common"].sum()) == 1
    assert set(blocks["start_id"]).isdisjoint({"rare1", "rare2"})
    assert "c1" in set(blocks["start_id"]) | set(blocks["end_id"])


def test_plot_circos_ld_map_from_tables_writes_png(tmp_path: Path):
    density_path = tmp_path / "all.variant.circos_ld_map.density.info.tsv"
    blocks_path = tmp_path / "all.variant.circos_ld_map.blocks.info.tsv"
    pd.DataFrame(
        {
            "chrom": ["1A", "1A"],
            "bin_start": [0, 1_000_000],
            "bin_end": [1_000_000, 2_000_000],
            "bin_mid": [500_000.0, 1_500_000.0],
            "n_rare": [10, 2],
            "n_common": [1, 8],
            "n_sites": [11, 10],
        }
    ).to_csv(density_path, sep="\t", index=False)
    pd.DataFrame(
        {
            "chrom": ["1A", "1A"],
            "start_id": ["s1", "s3"],
            "end_id": ["s2", "s3"],
            "start_pos": [1000, 80_000],
            "end_pos": [50_000, 80_000],
            "n_snps": [2, 1],
        }
    ).to_csv(blocks_path, sep="\t", index=False)
    prefix = str(tmp_path / "from_tables.variant.circos_ld_map")
    plot_circos_ld_map_from_tables(
        str(density_path),
        str(blocks_path),
        prefix,
        chrom_lengths={"1A": 5_000_000},
    )
    assert Path(f"{prefix}.png").is_file()
    assert Path(f"{prefix}.png").stat().st_size > 1000
    assert Path(f"{prefix}.blocks_1mb.info.tsv").is_file()


def test_ld_block_coverage_intervals_merge_1mb_bins_at_half_coverage():
    blocks = pd.DataFrame(
        {
            "chrom": ["1A", "1A", "1A", "1A", "1A"],
            "start_id": ["a", "b", "c", "d", "e"],
            "end_id": ["a2", "b2", "c2", "d2", "e2"],
            "start_pos": [0, 1_000_000, 2_000_000, 3_500_000, 4_000_000],
            "end_pos": [599_999, 1_599_999, 2_099_999, 3_500_000, 4_699_999],
            "n_snps": [10, 10, 3, 1, 8],
        }
    )
    out = ld_block_coverage_intervals(
        blocks,
        chrom_lengths={"1A": 5_000_000},
        bin_bp=1_000_000,
        coverage_frac=0.5,
    )
    assert list(zip(out["start_pos"].astype(int), out["end_pos"].astype(int))) == [
        (0, 2_000_000),
        (4_000_000, 5_000_000),
    ]
    assert list(out["n_bins"].astype(int)) == [2, 1]
