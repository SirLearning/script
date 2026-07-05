"""Unit tests for nucleotide diversity (pi) from PLINK2 gcount columns."""

import pytest

from genetics.genomics.plink.nucleotide_diversity import (
    is_biallelic_snp,
    merge_thin_common_pi_by_ref,
    site_pi_from_gcount,
)


def test_is_biallelic_snp_accepts_substitution():
    assert is_biallelic_snp('A', 'G') is True


def test_is_biallelic_snp_rejects_symbolic_indel():
    assert is_biallelic_snp('C', 'D') is False
    assert is_biallelic_snp('T', 'I') is False


def test_is_biallelic_snp_rejects_multichar():
    assert is_biallelic_snp('AT', 'A') is False


def test_site_pi_fixed_het():
    # N=2 diploid, both het -> x=2 alt alleles, 2N=4 -> pi = 2/(4*3)*2*2 = 8/12
    pi = site_pi_from_gcount(hom_ref=0, het=2, hom_alt=0)
    assert pi == pytest.approx(8 / 12)


def test_site_pi_fixed_hom_alt():
    # N=2, both hom alt -> x=4, 2N=4 -> pi = 2/(4*3)*4*0 = 0
    pi = site_pi_from_gcount(hom_ref=0, het=0, hom_alt=2)
    assert pi == pytest.approx(0.0)


def test_site_pi_insufficient_samples():
    import math
    assert math.isnan(site_pi_from_gcount(hom_ref=1, het=0, hom_alt=0))


def test_merge_thin_common_pi_by_ref(tmp_path):
    thin = tmp_path / 'thin.chr_pi.by_ref.tsv'
    common = tmp_path / 'common.chr_pi.by_ref.tsv'
    thin.write_text(
        'ref_name\tlength_bp\tn_snp_sites\tpi_sum_sites\tpi\tpi_per_mb\n'
        'chr1A\t100\t10\t0.01\t0.0001\t100\n'
    )
    common.write_text(
        'ref_name\tlength_bp\tn_snp_sites\tpi_sum_sites\tpi\tpi_per_mb\n'
        'chr1A\t100\t8\t0.008\t0.00008\t80\n'
    )
    merged = merge_thin_common_pi_by_ref(thin, common)
    assert merged.loc[0, 'common_fraction'] == pytest.approx(0.8)
