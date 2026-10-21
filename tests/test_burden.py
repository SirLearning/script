"""Gene burden collapsing and association."""

import numpy as np
import pandas as pd
import pytest

from genetics.gwas.burden import (
    COLLAPSE_CARRIER,
    COLLAPSE_COUNT,
    align_burden_inputs,
    collapse_gene_scores,
    collapse_gene_scores_from_raw,
    collapse_gene_scores_from_traw,
    load_gene_map,
    load_plink_raw_genotypes,
    run_burden_association,
    run_burden_test,
)


def _write_raw(path, sample_ids, columns, dosages):
    header = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'] + list(columns)
    lines = ['\t'.join(header)]
    for i, sample in enumerate(sample_ids):
        row = ['0', sample, '0', '0', '0', '-9'] + [str(v) for v in dosages[i]]
        lines.append('\t'.join(row))
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def _write_traw(path, sample_ids, variant_ids, counted_alleles, alt_alleles, dosages):
    header = ['CHR', 'SNP', '(C)M', 'POS', 'COUNTED', 'ALT'] + [f'0_{s}' for s in sample_ids]
    lines = ['\t'.join(header)]
    for variant_id, counted, alt, row in zip(variant_ids, counted_alleles, alt_alleles, dosages):
        chrom, pos = variant_id.split('-', 1)
        lines.append('\t'.join([chrom, variant_id, '0', pos, counted, alt] + [str(v) for v in row]))
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def _write_variant_sets(path, rows, set_columns=('set_plof_rare_wap',)):
    header = ['ID', 'Gene'] + list(set_columns)
    lines = ['\t'.join(header)]
    for variant_id, gene, *flags in rows:
        lines.append('\t'.join([variant_id, gene] + [str(f) for f in flags]))
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def _gene_info(n_variants=1, missing=0.0, gene='GENE1'):
    return pd.DataFrame(
        {'n_variants': [n_variants], 'mean_missing_rate': [missing]}, index=[gene]
    )


def _write_pheno(path, values, trait='PH_M_cm_mean'):
    lines = [f'#FID\tIID\t{trait}']
    for sample, value in values.items():
        lines.append(f'0\t{sample}\t{value}')
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def test_load_gene_map_keeps_only_set_members_with_a_gene(tmp_path):
    table = tmp_path / 'A.variant_sets.tsv'
    _write_variant_sets(table, [
        ('1-100', 'GENE1', 'True'),
        ('1-200', 'GENE1', 'True'),
        ('1-300', 'GENE2', 'False'),   # not in the set
        ('1-400', '-', 'True'),        # VEP assigned no gene
    ])
    gene_map = load_gene_map(table, 'plof_rare_wap')
    assert list(gene_map['ID']) == ['1-100', '1-200']
    assert set(gene_map['gene']) == {'GENE1'}


def test_load_gene_map_accepts_prefixed_and_bare_set_name(tmp_path):
    table = tmp_path / 'A.variant_sets.tsv'
    _write_variant_sets(table, [('1-100', 'GENE1', 'True')])
    assert len(load_gene_map(table, 'plof_rare_wap')) == 1
    assert len(load_gene_map(table, 'set_plof_rare_wap')) == 1


def test_load_gene_map_rejects_unknown_set(tmp_path):
    table = tmp_path / 'A.variant_sets.tsv'
    _write_variant_sets(table, [('1-100', 'GENE1', 'True')])
    with pytest.raises(ValueError, match='no column set_nope'):
        load_gene_map(table, 'nope')


def test_collapse_count_sums_minor_alleles_across_a_gene(tmp_path):
    raw = tmp_path / 'A.raw'
    _write_raw(
        raw,
        ['S1', 'S2', 'S3', 'S4'],
        ['1-100_A', '1-200_G'],
        [[0, 0], [1, 0], [1, 1], [2, 0]],
    )
    gene_map = pd.DataFrame({'ID': ['1-100', '1-200'], 'gene': ['GENE1', 'GENE1']})
    genotypes = load_plink_raw_genotypes(raw)
    scores, gene_info = collapse_gene_scores(genotypes, gene_map, COLLAPSE_COUNT)
    assert list(scores['GENE1']) == [0.0, 1.0, 2.0, 2.0]
    assert gene_info.loc['GENE1', 'n_variants'] == 2


def test_collapse_carrier_reduces_to_indicator(tmp_path):
    raw = tmp_path / 'A.raw'
    _write_raw(
        raw,
        ['S1', 'S2', 'S3', 'S4'],
        ['1-100_A', '1-200_G'],
        [[0, 0], [1, 0], [1, 1], [2, 0]],
    )
    gene_map = pd.DataFrame({'ID': ['1-100', '1-200'], 'gene': ['GENE1', 'GENE1']})
    genotypes = load_plink_raw_genotypes(raw)
    scores, _ = collapse_gene_scores(genotypes, gene_map, COLLAPSE_CARRIER)
    assert list(scores['GENE1']) == [0.0, 1.0, 1.0, 1.0]


def test_collapse_flips_when_counted_allele_is_major(tmp_path):
    """--export A counts ALT, so a major ALT must be flipped to keep scores minor-allele."""
    raw = tmp_path / 'A.raw'
    _write_raw(
        raw,
        ['S1', 'S2', 'S3', 'S4'],
        ['1-100_A'],
        [[2], [2], [2], [0]],
    )
    gene_map = pd.DataFrame({'ID': ['1-100'], 'gene': ['GENE1']})
    genotypes = load_plink_raw_genotypes(raw)
    scores, _ = collapse_gene_scores(genotypes, gene_map, COLLAPSE_COUNT)
    assert list(scores['GENE1']) == [0.0, 0.0, 0.0, 2.0]


def test_missing_call_on_a_flipped_variant_is_not_a_carrier(tmp_path):
    """Zero-filling before the flip would turn every missing call into a homozygote.

    This is the real shape of a rare pLoF site: ALT is the major allele, so the
    column is flipped, and roughly a tenth of the calls are missing.
    """
    n = 100
    dosages = [[2]] * 88 + [['NA']] * 10 + [[1]] * 2
    sample_ids = [f'S{i:03d}' for i in range(n)]
    raw = tmp_path / 'A.raw'
    _write_raw(raw, sample_ids, ['1-100_A'], dosages)

    gene_map = pd.DataFrame({'ID': ['1-100'], 'gene': ['GENE1']})
    genotypes = load_plink_raw_genotypes(raw)
    scores, gene_info = collapse_gene_scores(genotypes, gene_map, COLLAPSE_COUNT)

    # Only the two heterozygotes carry the minor allele; the 10 missing calls must
    # not be counted, and must certainly not be counted as homozygotes.
    assert int((scores['GENE1'] > 0).sum()) == 2
    assert scores['GENE1'].max() == 1.0
    assert gene_info.loc['GENE1', 'mean_missing_rate'] == pytest.approx(0.10)


def test_collapse_treats_missing_call_as_non_carrier(tmp_path):
    raw = tmp_path / 'A.raw'
    _write_raw(
        raw,
        ['S1', 'S2', 'S3'],
        ['1-100_A'],
        [['NA'], [1], [0]],
    )
    gene_map = pd.DataFrame({'ID': ['1-100'], 'gene': ['GENE1']})
    genotypes = load_plink_raw_genotypes(raw)
    scores, _ = collapse_gene_scores(genotypes, gene_map, COLLAPSE_COUNT)
    assert list(scores['GENE1']) == [0.0, 1.0, 0.0]


def test_collapse_from_raw_matches_in_memory_on_wide_and_flipped_sites(tmp_path):
    """Streaming collapse must not require the full sample-by-variant matrix in RAM."""
    raw = tmp_path / 'A.raw'
    _write_raw(
        raw,
        ['S1', 'S2', 'S3', 'S4'],
        ['1-100_A', '1-200_G', '1-300_T', '1-400_C'],
        [
            [2, 0, 'NA', 1],
            [2, 1, 0, 0],
            [2, 0, 1, 1],
            [0, 1, 0, 0],
        ],
    )
    gene_map = pd.DataFrame({
        'ID': ['1-100', '1-200', '1-300', '1-400'],
        'gene': ['GENE1', 'GENE1', 'GENE2', 'GENE2'],
    })
    genotypes = load_plink_raw_genotypes(raw)
    for rule in (COLLAPSE_COUNT, COLLAPSE_CARRIER):
        expected_scores, expected_info = collapse_gene_scores(genotypes, gene_map, rule)
        got_scores, got_info = collapse_gene_scores_from_raw(raw, gene_map, collapse=rule)
        pd.testing.assert_frame_equal(
            got_scores[expected_scores.columns], expected_scores, check_dtype=False
        )
        pd.testing.assert_frame_equal(got_info, expected_info, check_dtype=False)


def test_collapse_from_traw_matches_in_memory_on_wide_and_flipped_sites(tmp_path):
    """Variant-major .traw must match the sample-major in-memory collapse."""
    raw = tmp_path / 'A.raw'
    traw = tmp_path / 'A.traw'
    _write_raw(
        raw,
        ['S1', 'S2', 'S3', 'S4'],
        ['1-100_A', '1-200_G', '1-300_T', '1-400_C'],
        [
            [2, 0, 'NA', 1],
            [2, 1, 0, 0],
            [2, 0, 1, 1],
            [0, 1, 0, 0],
        ],
    )
    _write_traw(
        traw,
        ['S1', 'S2', 'S3', 'S4'],
        ['1-100', '1-200', '1-300', '1-400'],
        ['A', 'G', 'T', 'C'],
        ['G', 'A', 'A', 'T'],
        [
            [2, 2, 2, 0],
            [0, 1, 0, 1],
            ['NA', 0, 1, 0],
            [1, 0, 1, 0],
        ],
    )
    gene_map = pd.DataFrame({
        'ID': ['1-100', '1-200', '1-300', '1-400'],
        'gene': ['GENE1', 'GENE1', 'GENE2', 'GENE2'],
    })
    genotypes = load_plink_raw_genotypes(raw)
    for rule in (COLLAPSE_COUNT, COLLAPSE_CARRIER):
        expected_scores, expected_info = collapse_gene_scores(genotypes, gene_map, rule)
        got_scores, got_info = collapse_gene_scores_from_traw(traw, gene_map, collapse=rule)
        pd.testing.assert_frame_equal(
            got_scores[expected_scores.columns], expected_scores, check_dtype=False
        )
        pd.testing.assert_frame_equal(got_info, expected_info, check_dtype=False)


def test_collapse_keeps_genes_apart(tmp_path):
    raw = tmp_path / 'A.raw'
    _write_raw(
        raw,
        ['S1', 'S2'],
        ['1-100_A', '1-200_G'],
        [[1, 0], [0, 1]],
    )
    gene_map = pd.DataFrame({'ID': ['1-100', '1-200'], 'gene': ['GENE1', 'GENE2']})
    genotypes = load_plink_raw_genotypes(raw)
    scores, gene_info = collapse_gene_scores(genotypes, gene_map, COLLAPSE_COUNT)
    assert list(scores['GENE1']) == [1.0, 0.0]
    assert list(scores['GENE2']) == [0.0, 1.0]
    assert gene_info.loc['GENE1', 'n_variants'] == 1 and gene_info.loc['GENE2', 'n_variants'] == 1


def test_burden_beta_and_p_match_statsmodels_ols():
    """The Frisch-Waugh shortcut must reproduce a full OLS fit with covariates."""
    statsmodels = pytest.importorskip('statsmodels.api')
    rng = np.random.default_rng(20260903)
    n = 200
    samples = pd.Index([f'S{i:03d}' for i in range(n)])
    score = rng.poisson(0.4, n).astype(float)
    pc1 = rng.normal(size=n)
    trait = 3.0 + 1.7 * score - 0.8 * pc1 + rng.normal(scale=2.0, size=n)

    scores = pd.DataFrame({'GENE1': score}, index=samples)
    covariates = pd.DataFrame({'PC1': pc1}, index=samples)
    inputs = align_burden_inputs(
        scores, _gene_info(1), pd.Series(trait, index=samples), covariates
    )
    got = run_burden_association(inputs).iloc[0]

    design = statsmodels.add_constant(np.column_stack([score, pc1]))
    fit = statsmodels.OLS(trait, design).fit()

    assert got['beta'] == pytest.approx(fit.params[1], rel=1e-9)
    assert got['se'] == pytest.approx(fit.bse[1], rel=1e-9)
    assert got['p_value'] == pytest.approx(fit.pvalues[1], rel=1e-9)
    assert got['n_samples'] == n


def test_burden_reports_untestable_gene_without_crashing():
    """A gene nobody carries has no contrast; it must be reported, not divided by zero."""
    samples = pd.Index([f'S{i}' for i in range(10)])
    scores = pd.DataFrame({'GENE1': np.zeros(10)}, index=samples)
    trait = pd.Series(np.arange(10, dtype=float), index=samples)
    inputs = align_burden_inputs(scores, _gene_info(2), trait)
    got = run_burden_association(inputs).iloc[0]
    assert got['n_carriers'] == 0
    assert np.isnan(got['beta']) and np.isnan(got['p_value'])
    assert got['n_variants'] == 2


def test_align_intersects_samples_across_tables():
    scores = pd.DataFrame({'GENE1': [1.0, 2.0, 3.0, 4.0]}, index=pd.Index(['S1', 'S2', 'S3', 'S4']))
    trait = pd.Series([1.0, 2.0, 3.0, 4.0], index=pd.Index(['S2', 'S3', 'S4', 'S5']))
    covar = pd.DataFrame({'PC1': [0.1, 0.2, 0.3]}, index=pd.Index(['S2', 'S3', 'S4']))
    inputs = align_burden_inputs(scores, _gene_info(1), trait, covar)
    assert list(inputs.samples) == ['S2', 'S3', 'S4']


def test_align_rejects_a_panel_too_small_to_fit():
    scores = pd.DataFrame({'GENE1': [1.0, 2.0]}, index=pd.Index(['S1', 'S2']))
    trait = pd.Series([1.0, 2.0], index=pd.Index(['S1', 'S2']))
    with pytest.raises(ValueError, match='Only 2 samples'):
        align_burden_inputs(scores, _gene_info(1), trait)


def test_run_burden_test_writes_full_precision_p_values(tmp_path):
    """save_df_to_tsv rounds to 4 decimals by default, which would flatten p-values."""
    rng = np.random.default_rng(7)
    n = 120
    sample_ids = [f'S{i:03d}' for i in range(n)]
    dosage = rng.integers(0, 2, size=(n, 1))
    trait = 100.0 + 25.0 * dosage[:, 0] + rng.normal(scale=1.0, size=n)

    raw = tmp_path / 'A.raw'
    _write_raw(raw, sample_ids, ['1-100_A'], dosage.tolist())
    sets = tmp_path / 'A.variant_sets.tsv'
    _write_variant_sets(sets, [('1-100', 'GENE1', 'True')])
    pheno = tmp_path / 'pheno.tsv'
    _write_pheno(pheno, dict(zip(sample_ids, trait)))

    results = run_burden_test(
        str(raw), str(sets), 'plof_rare_wap', str(pheno),
        'PH_M_cm_mean', str(tmp_path / 'A.plof_rare_wap'),
    )
    assert len(results) == 1
    # A 25 cm effect on 120 samples is unmissable, so the p-value must survive
    # rounding as something far below the 1e-4 the default format would clamp to.
    written = pd.read_csv(tmp_path / 'A.plof_rare_wap.burden.tsv', sep='\t')
    assert written.loc[0, 'p_value'] < 1e-20
    assert written.loc[0, 'beta'] == pytest.approx(25.0, abs=1.0)
