"""Chr-count readers must agree on arm-numbered and ref-chr named pvars."""

import pandas as pd
import pytest

from genetics.genomics.plink.results_io import (
    iter_process_pvar_positions,
    summarize_chr_variant_counts,
)
from infra.wheat.ref_v1 import pvar_chrom_to_plink_chr_id


def _write_pvar(path, rows):
    lines = ['#CHROM\tPOS\tID\tREF\tALT']
    lines += [f'{chrom}\t{pos}\t{chrom}-{pos}\tA\tG' for chrom, pos in rows]
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def test_pvar_chrom_to_plink_chr_id_maps_names_to_first_arm():
    assert pvar_chrom_to_plink_chr_id('13') == '13'
    assert pvar_chrom_to_plink_chr_id('chr3A') == '13'
    assert pvar_chrom_to_plink_chr_id('chrUn') == '0'
    assert pvar_chrom_to_plink_chr_id('Mit') == '43'
    assert pvar_chrom_to_plink_chr_id('Chl') == '44'
    assert pvar_chrom_to_plink_chr_id('not_a_chrom') is None


def test_named_pvar_counts_match_arm_numbered_pvar(tmp_path):
    """Renaming #CHROM must not zero out the per-chromosome counts."""
    numeric_dir = tmp_path / 'numeric'
    named_dir = tmp_path / 'named'
    numeric_dir.mkdir()
    named_dir.mkdir()

    _write_pvar(numeric_dir / 'A_test.plink2.pvar', [('13', 100), ('13', 200), ('1', 50)])
    _write_pvar(named_dir / 'A_test.plink2.pvar', [('chr3A', 100), ('chr3A', 200), ('chr1A', 50)])

    numeric = summarize_chr_variant_counts(numeric_dir, str(tmp_path / 'num'))
    named = summarize_chr_variant_counts(named_dir, str(tmp_path / 'nam'))

    per_chr_numeric = pd.read_csv(f'{tmp_path / "num"}.chr_variant_counts.tsv', sep='\t')
    per_chr_named = pd.read_csv(f'{tmp_path / "nam"}.chr_variant_counts.tsv', sep='\t')

    named_rows = per_chr_named[per_chr_named['plink_chr'] != 'total']
    assert named_rows['n_variants'].sum() == 3
    assert named_rows.loc[named_rows['plink_chr'] == '1', 'n_variants'].item() == 1
    assert named_rows.loc[named_rows['plink_chr'] == '13', 'n_variants'].item() == 2
    pd.testing.assert_series_equal(
        per_chr_numeric['n_variants'], per_chr_named['n_variants']
    )
    assert numeric is not None and named is not None


def test_iter_process_pvar_positions_reads_named_chroms(tmp_path):
    _write_pvar(tmp_path / 'A_test.plink2.pvar', [('chr3A', 100), ('Mit', 7)])
    assert sorted(iter_process_pvar_positions(tmp_path)) == [(13, 100), (43, 7)]


def test_popdep_rejects_merged_pvar(tmp_path):
    """Arm-local frozen grids cannot be joined against merged coordinates."""
    from genetics.genomics.variant.popdep import annotate_variants_popdep_from_pvar

    pvar = tmp_path / 'A.pvar'
    _write_pvar(pvar, [('chr3A', 100)])
    with pytest.raises(ValueError, match='merged ref-chr names'):
        annotate_variants_popdep_from_pvar(
            str(pvar), str(tmp_path), str(tmp_path / 'out.tsv')
        )
