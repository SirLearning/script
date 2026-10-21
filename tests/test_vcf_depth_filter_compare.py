"""Tests for MiDR/MaDR VCF genotype-mask comparison."""

from pathlib import Path

from genetics.genomics.variant.vcf_depth_filter_compare import (
    compare_filter_vcfs,
    is_het_gt,
    is_missing_gt,
    parse_ad,
    site_depth,
)


def _write_tbm(path: Path) -> Path:
    path.write_text(
        "Taxa\tCoverage-Of-All-Bams\tBams\n"
        "s1\t2.0\t/tmp/s1.bam\n"
        "s2\t10.0\t/tmp/s2.bam\n",
        encoding="utf-8",
    )
    return path


def _write_vcf(path: Path, samples: list[str], records: list[str]) -> Path:
    header = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples)
        + "\n"
    )
    path.write_text(header + "".join(records), encoding="utf-8")
    return path


def test_gt_helpers():
    assert is_missing_gt("./.")
    assert is_missing_gt("./.:0,0")
    assert not is_missing_gt("0/1:3,1")
    assert is_het_gt("0/1:3,1")
    assert not is_het_gt("1/1:0,4")
    assert parse_ad("0/1:3,1:0,1,4", ["GT", "AD", "GL"]) == [3, 1]
    assert site_depth([3, 1]) == 4


def test_midr_madr_cell_classes(tmp_path: Path):
    tbm = _write_tbm(tmp_path / "tbm.txt")
    # s1 cov=2 → MiDR 0.2 is 0.4, so DP=1 is called; MaDR 3 is 6.
    # s2 cov=10 → MiDR 2, MaDR 30.
    unf = _write_vcf(
        tmp_path / "u.vcf",
        ["s1", "s2"],
        [
            "1\t10\t.\tA\tC\t.\t.\tDP=41;NZ=2;AD=21,20;HT=1;MAF=0.5\tGT:AD\t0/1:1,0\t0/0:40,0\n",
            "1\t20\t.\tG\tT\t.\t.\tDP=42;NZ=2;AD=0,42;HT=0;MAF=0.5\tGT:AD\t1/1:0,2\t1/1:0,40\n",
        ],
    )
    flt = _write_vcf(
        tmp_path / "f.vcf",
        ["s1", "s2"],
        [
            "1\t10\t.\tA\tC\t.\t.\tDP=1;NZ=1;AD=1,0;HT=1;MAF=0.5\tGT:AD\t0/1:1,0\t./.\n",
            "1\t20\t.\tG\tT\t.\t.\tDP=2;NZ=1;AD=0,2;HT=0;MAF=1.0\tGT:AD\t1/1:0,2\t./.\n",
        ],
    )
    result = compare_filter_vcfs(unf, flt, tbm, midr=0.2, madr=3.0, max_sites=10)
    assert result.samples_equal
    assert result.n_sites_matched == 2
    assert result.n_sites_only_unfiltered == 0
    assert result.n_sites_only_filtered == 0
    assert result.cells.both_called == 2
    assert result.cells.both_called_gt_mismatch == 0
    assert result.cells.both_called_ad_mismatch == 0
    assert result.cells.unfiltered_called_filtered_missing == 2
    assert result.cells.unfiltered_missing_filtered_called == 0
    assert result.cells.madr_explained == 2
    assert result.cells.midr_explained == 0
    assert result.cells.called_to_missing_unexplained == 0
    assert result.cells.leak_still_called == 0
    assert result.cells.het_to_hom == 0
    assert result.cells.het_unfiltered == 1
    assert result.cells.het_filtered == 1
