"""Selfing approximation: rewrite unphased VCF GT '/' to '|'."""

from pathlib import Path
import gzip

import pytest

from genetics.genomics.variant.selfing_gt import rewrite_gt_slash_to_pipe


VCF_BODY = """\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##INFO=<ID=XX,Number=1,Type=String,Description="Has a slash">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1A	100	1-100	A	T	.	.	XX=a/b	GT:AD	0/0:2,0	0/1:1,1	./.
chr1A	200	1-200	C	G	.	.	.	GT	1/1	0|1	.
"""


def test_slash_to_pipe_rewrites_gt_only(tmp_path: Path):
    src = tmp_path / "in.vcf"
    dst = tmp_path / "out.vcf"
    src.write_text(VCF_BODY, encoding="utf-8")

    n_sites = rewrite_gt_slash_to_pipe(src, dst)
    assert n_sites == 2

    text = dst.read_text(encoding="utf-8")
    assert "XX=a/b" in text
    assert "0|0:2,0" in text
    assert "0|1:1,1" in text
    assert ".|." in text
    assert "1|1" in text
    assert "0|1" in text
    assert "0/0" not in text
    assert "0/1" not in text
    assert "./." not in text


def test_slash_to_pipe_roundtrip_gzip(tmp_path: Path):
    src = tmp_path / "in.vcf.gz"
    dst = tmp_path / "out.vcf.gz"
    with gzip.open(src, "wt", encoding="utf-8") as handle:
        handle.write(VCF_BODY)

    rewrite_gt_slash_to_pipe(src, dst)
    with gzip.open(dst, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert "0|0:2,0" in text
    assert "XX=a/b" in text


def test_ne_from_two_diploids_one_het(tmp_path: Path):
    from genetics.genomics.variant.singer_arg import estimate_ne_from_vcf

    vcf = tmp_path / "tiny.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
                "1\t50\ta\tA\tT\t.\t.\t.\tGT\t0|0\t0|1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    # Two diploids, allele counts 3 ref / 1 alt → freqs 0.75, 0.25 → π_site = 0.375
    # window length 100 → diversity 0.00375 → Ne = 0.00375 / 4 / 1e-8
    ne = estimate_ne_from_vcf(vcf, start=0, end=100, mut_rate=1e-8)
    assert ne == pytest.approx(0.00375 / 4e-8)
