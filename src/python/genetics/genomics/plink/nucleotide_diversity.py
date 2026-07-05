"""Nucleotide diversity (pi, Tajima theta_pi) from PLINK2 pfile artefacts.

PLINK2 does not compute pi directly. This module reads merged subgenome pvar +
``*.info.gcount`` under a process directory, keeps biallelic SNPs only (excludes
indel symbolic alleles D/I and multi-character alleles), and applies the
pairwise nucleotide-difference estimator per assayed site:

    pi_site = (2 / (2N(2N-1))) * x * (2N - x)

where N is the number of called diploid samples at the site and x is the alt-allele
count among 2N chromosomes (x = 2 * HOM_ALT + HET).

Regional pi divides the sum of pi_site values by the reference length L in bp
(chromosome, genome bin, or whole genome). Unassayed positions contribute zero,
consistent with VCFtools ``--window-pi`` on SNP-only input.
"""

from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.plink.results_io import (
    _ref_name_at_global_bp,
    _ref_name_sort_key,
)
from genetics.genomics.variant.variant_utils import load_df_from_plink_gcount
from infra.utils.io import load_df_from_tsv, save_df_to_tsv
from infra.wheat.ref_v1 import (
    all_ref_v1_plink_chr_ids,
    get_ref_v1_chr_name,
    get_ref_v1_genome_segment_layout,
    get_ref_v1_plink_chr_length,
)


def is_biallelic_snp(ref: str, alt: str) -> bool:
    """Return True for single-nucleotide biallelic SNPs (excludes D/I indels)."""
    alt_primary = str(alt).split(',')[0]
    if alt_primary in ('D', 'I'):
        return False
    if len(ref) != 1 or len(alt_primary) != 1:
        return False
    valid = frozenset('ACGTacgt')
    return ref in valid and alt_primary in valid


def site_pi_from_gcount(hom_ref: int, het: int, hom_alt: int) -> float:
    """
    Tajima per-site pi from PLINK2 ``--geno-counts`` columns.

    Uses 2N chromosomes (N diploid called samples). Returns NaN when N < 2.
    """
    n = int(hom_ref) + int(het) + int(hom_alt)
    if n < 2:
        return float('nan')
    chromosomes = 2 * n
    alt_alleles = 2 * int(hom_alt) + int(het)
    return (2.0 / (chromosomes * (chromosomes - 1))) * alt_alleles * (chromosomes - alt_alleles)


def _subgenome_from_merged_pvar(pvar_path: Path) -> str:
    stem = pvar_path.name.replace('.plink2.pvar', '')
    return stem.split('_', 1)[0]


def _list_merged_pvar_gcount_pairs(process_dir: Path) -> list[tuple[Path, Path]]:
    pairs: list[tuple[Path, Path]] = []
    for pvar in sorted(process_dir.glob('*_test.plink2.pvar')):
        sub = _subgenome_from_merged_pvar(pvar)
        gcount = process_dir / 'variant' / f'{sub}.info.gcount'
        if not gcount.exists():
            raise FileNotFoundError(f'Missing gcount for {pvar.name}: {gcount}')
        pairs.append((pvar, gcount))
    if not pairs:
        raise FileNotFoundError(
            f'No *_test.plink2.pvar + variant/*.info.gcount under {process_dir}'
        )
    return pairs


def _load_snp_site_table(pvar_path: Path, gcount_path: Path) -> pd.DataFrame:
    """Join pvar REF/ALT with gcount rows on variant ID; keep SNPs with pi_site."""
    pvar_rows = []
    with open(pvar_path, encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            plink_chr = int(parts[0])
            pos = int(parts[1])
            ref, alt = parts[3], parts[4]
            if not is_biallelic_snp(ref, alt):
                continue
            pvar_rows.append({
                'ID': parts[2],
                'plink_chr': plink_chr,
                'pos': pos,
                'ref_name': get_ref_v1_chr_name(plink_chr),
            })

    if not pvar_rows:
        return pd.DataFrame(columns=['ID', 'plink_chr', 'pos', 'ref_name', 'pi_site'])

    pvar_df = pd.DataFrame(pvar_rows)
    gcount = load_df_from_plink_gcount(str(gcount_path))
    required = {'ID', 'HOM_REF_CT', 'HET_REF_ALT_CTS', 'TWO_ALT_GENO_CTS'}
    missing = required - set(gcount.columns)
    if missing:
        raise ValueError(f'{gcount_path} missing columns: {sorted(missing)}')

    merged = pvar_df.merge(
        gcount[list(required)],
        on='ID',
        how='inner',
        validate='one_to_one',
    )
    merged['pi_site'] = merged.apply(
        lambda row: site_pi_from_gcount(
            row['HOM_REF_CT'],
            row['HET_REF_ALT_CTS'],
            row['TWO_ALT_GENO_CTS'],
        ),
        axis=1,
    )
    return merged.dropna(subset=['pi_site']).reset_index(drop=True)


def _load_all_snp_site_pi_table(process_dir: str | Path) -> pd.DataFrame:
    """Load SNP-only per-site pi for all merged subgenome pfiles under process_dir."""
    frames = []
    for pvar, gcount in _list_merged_pvar_gcount_pairs(Path(process_dir)):
        frames.append(_load_snp_site_table(pvar, gcount))
    if not frames:
        return pd.DataFrame(columns=['plink_chr', 'pos', 'ref_name', 'pi_site'])
    return pd.concat(frames, ignore_index=True)


@lru_cache(maxsize=4)
def _cached_site_pi_table(process_dir: str) -> pd.DataFrame:
    return _load_all_snp_site_pi_table(process_dir)


def iter_snp_site_pi(process_dir: str | Path):
    """Yield (plink_chr, pos, ref_name, pi_site) for every assayed biallelic SNP."""
    table = _cached_site_pi_table(str(process_dir))
    for row in table.itertuples(index=False):
        yield int(row.plink_chr), int(row.pos), str(row.ref_name), float(row.pi_site)


def get_ref_v1_ref_name_length_bp(ref_name: str) -> int:
    """Sum PLINK segment lengths that map to a vmap4 ref name (e.g. chr1A)."""
    total = 0
    for cid in all_ref_v1_plink_chr_ids():
        if get_ref_v1_chr_name(cid) == ref_name:
            total += get_ref_v1_plink_chr_length(cid)
    return total


def _ref_name_lengths() -> dict[str, int]:
    lengths: dict[str, int] = {}
    for cid in all_ref_v1_plink_chr_ids():
        name = get_ref_v1_chr_name(cid)
        lengths[name] = lengths.get(name, 0) + get_ref_v1_plink_chr_length(cid)
    return lengths


def summarize_pi_by_ref(process_dir: str | Path, output_prefix: str) -> pd.DataFrame:
    """
    Compute pi per vmap4 ref chromosome and write ``{output_prefix}.chr_pi.by_ref.tsv``.

    pi_ref = sum(pi_site over SNPs on ref) / L_ref_bp
    """
    ref_lengths = _ref_name_lengths()
    table = _cached_site_pi_table(str(process_dir))
    pi_sum: dict[str, float] = defaultdict(float)
    site_counts: dict[str, int] = defaultdict(int)
    for row in table.itertuples(index=False):
        pi_sum[str(row.ref_name)] += float(row.pi_site)
        site_counts[str(row.ref_name)] += 1

    rows = []
    for ref_name in sorted(ref_lengths, key=_ref_name_sort_key):
        length_bp = ref_lengths[ref_name]
        total_pi = pi_sum.get(ref_name, 0.0)
        rows.append({
            'ref_name': ref_name,
            'length_bp': length_bp,
            'n_snp_sites': site_counts.get(ref_name, 0),
            'pi_sum_sites': total_pi,
            'pi': total_pi / length_bp if length_bp > 0 else 0.0,
            'pi_per_mb': (total_pi / length_bp * 1_000_000.0) if length_bp > 0 else 0.0,
        })

    df = pd.DataFrame(rows)
    total_length = sum(ref_lengths.values())
    total_pi_sum = df['pi_sum_sites'].sum()
    total_pi = total_pi_sum / total_length if total_length > 0 else 0.0
    total_row = {
        'ref_name': 'total',
        'length_bp': total_length,
        'n_snp_sites': int(df['n_snp_sites'].sum()),
        'pi_sum_sites': total_pi_sum,
        'pi': total_pi,
        'pi_per_mb': total_pi * 1_000_000.0,
    }
    df = pd.concat([df, pd.DataFrame([total_row])], ignore_index=True)
    save_df_to_tsv(df, f'{output_prefix}.chr_pi.by_ref.tsv', float_format='%.8g')
    return df


def compute_genome_bin_pi(process_dir: str | Path, bin_size_bp: int = 5_000_000) -> pd.DataFrame:
    """
    Bin SNP pi contributions into fixed-width windows along the vmap4 genome.

    Returns bin_index, ref_name, bin_start_mb, bin_end_mb, pi_sum_sites, pi,
    and pi_per_mb (= pi * 1e6, parallel to variant density_per_mb).
    """
    segments, total_bp = get_ref_v1_genome_segment_layout()
    starts = {seg['plink_chr']: seg['start_bp'] for seg in segments}
    n_bins = (total_bp + bin_size_bp - 1) // bin_size_bp
    pi_sums = np.zeros(n_bins, dtype=np.float64)

    table = _cached_site_pi_table(str(process_dir))
    for row in table.itertuples(index=False):
        global_bp = starts[int(row.plink_chr)] + int(row.pos) - 1
        if global_bp < 0 or global_bp >= total_bp:
            continue
        pi_sums[global_bp // bin_size_bp] += float(row.pi_site)

    bin_indices = np.arange(n_bins)
    bin_starts_mb = bin_indices * bin_size_bp / 1e6
    ref_names = [_ref_name_at_global_bp(segments, int(i * bin_size_bp)) for i in bin_indices]
    pi = pi_sums / float(bin_size_bp)
    return pd.DataFrame({
        'bin_index': bin_indices,
        'ref_name': ref_names,
        'bin_start_mb': bin_starts_mb,
        'bin_end_mb': bin_starts_mb + bin_size_bp / 1e6,
        'pi_sum_sites': pi_sums,
        'pi': pi,
        'pi_per_mb': pi * 1_000_000.0,
    })


def load_ordered_chr_pi_by_ref(tsv_path: str | Path) -> pd.DataFrame:
    """Load ``*.chr_pi.by_ref.tsv`` excluding the total row, in display order."""
    df = load_df_from_tsv(tsv_path)
    if df is None or df.empty:
        raise ValueError(f'Empty pi by ref table: {tsv_path}')
    df = df.loc[df['ref_name'] != 'total'].copy()
    if 'pi_per_mb' not in df.columns:
        df['pi_per_mb'] = pd.to_numeric(df['pi'], errors='coerce').fillna(0.0) * 1_000_000.0
    df['pi'] = pd.to_numeric(df['pi'], errors='coerce').fillna(0.0)
    df['pi_per_mb'] = pd.to_numeric(df['pi_per_mb'], errors='coerce').fillna(0.0)
    return df.sort_values(
        by='ref_name',
        key=lambda s: s.map(_ref_name_sort_key),
    ).reset_index(drop=True)


def merge_thin_common_pi_by_ref(thin_by_ref_tsv: str | Path, common_by_ref_tsv: str | Path) -> pd.DataFrame:
    """Merge per-ref pi tables and compute common/thin pi ratio (uses pi_per_mb)."""
    thin = load_ordered_chr_pi_by_ref(thin_by_ref_tsv).rename(
        columns={'pi': 'pi_thin', 'pi_per_mb': 'test_thin'},
    )
    common = load_ordered_chr_pi_by_ref(common_by_ref_tsv).rename(
        columns={'pi': 'pi_common', 'pi_per_mb': 'test_common_thin'},
    )
    merged = thin[['ref_name', 'test_thin', 'pi_thin']].merge(
        common[['ref_name', 'test_common_thin', 'pi_common']],
        on='ref_name',
        how='outer',
    ).fillna(0.0)
    merged['test_thin'] = merged['test_thin'].astype(float)
    merged['test_common_thin'] = merged['test_common_thin'].astype(float)
    merged = merged.sort_values(
        by='ref_name',
        key=lambda s: s.map(_ref_name_sort_key),
    ).reset_index(drop=True)
    merged['common_fraction'] = np.where(
        merged['test_thin'] > 0,
        merged['test_common_thin'] / merged['test_thin'],
        0.0,
    )
    return merged[['ref_name', 'test_thin', 'test_common_thin', 'pi_thin', 'pi_common', 'common_fraction']]


def build_genome_pi_compare_table(thin_bins: pd.DataFrame, common_bins: pd.DataFrame) -> pd.DataFrame:
    """Align two genome pi bin tables on bin_index (uses pi_per_mb columns)."""
    table = thin_bins[
        ['bin_index', 'ref_name', 'bin_start_mb', 'bin_end_mb', 'pi_per_mb']
    ].rename(columns={'pi_per_mb': 'test_thin'}).copy()
    common_by_bin = common_bins.set_index('bin_index')['pi_per_mb']
    table['test_common_thin'] = table['bin_index'].map(common_by_bin).fillna(0.0)
    return table.sort_values('bin_index').reset_index(drop=True)


def plot_thin_common_pi_compare(
    thin_process_dir: str | Path,
    common_process_dir: str | Path,
    output_prefix: str,
    bin_size_bp: int = 5_000_000,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute pi for both mods and write info TSVs + two comparison plots.

    Mirrors ``plot_thin_common_chr_variant_compare`` deliverables:
      - ``{prefix}.thin_common_pi_compare.info.tsv``
      - ``{prefix}.pi.genome_mb_density.info.tsv``
      - ``{prefix}.pi.chr_distribution.line.png``
      - ``{prefix}.pi.genome_mb_density.line.png`` (stacked panels)
      - ``{prefix}.pi.common_fraction.bar.png``
    """
    from genetics.genomics.plink.results_io import plot_thin_common_by_ref_distribution_line
    from infra.utils.graph import plot_bar_chart, plot_genome_binned_density_panels

    thin_pi_path = f'{output_prefix}.test_thin.chr_pi.by_ref.tsv'
    common_pi_path = f'{output_prefix}.test_common_thin.chr_pi.by_ref.tsv'
    summarize_pi_by_ref(thin_process_dir, f'{output_prefix}.test_thin')
    summarize_pi_by_ref(common_process_dir, f'{output_prefix}.test_common_thin')

    merged = merge_thin_common_pi_by_ref(thin_pi_path, common_pi_path)
    save_df_to_tsv(merged, f'{output_prefix}.thin_common_pi_compare.info.tsv')

    plot_thin_common_by_ref_distribution_line(
        merged,
        filename=f'{output_prefix}.pi.chr_distribution.line.png',
        title='Nucleotide diversity pi by chromosome (test_thin vs test_common_thin, SNP-only)',
        y_label='pi per Mb',
    )

    thin_bins = compute_genome_bin_pi(thin_process_dir, bin_size_bp=bin_size_bp)
    common_bins = compute_genome_bin_pi(common_process_dir, bin_size_bp=bin_size_bp)
    density = build_genome_pi_compare_table(thin_bins, common_bins)
    save_df_to_tsv(density, f'{output_prefix}.pi.genome_mb_density.info.tsv', float_format='%.4f')

    plot_genome_binned_density_panels(
        density,
        x_col='bin_index',
        ref_name_col='ref_name',
        panel_specs=[
            {'y_col': 'test_thin', 'title': 'test_thin', 'color': '#4c72b0'},
            {'y_col': 'test_common_thin', 'title': 'test_common_thin', 'color': '#c44e52'},
        ],
        filename=f'{output_prefix}.pi.genome_mb_density.line.png',
        suptitle='Genome-wide nucleotide diversity pi (5 Mb bins, SNP-only, vmap4 order)',
        y_label='pi per Mb (sum site pi / bin width)',
        bin_size_mb=bin_size_bp // 1_000_000,
    )

    plot_bar_chart(
        merged['ref_name'].tolist(),
        merged['common_fraction'].tolist(),
        title='Common-thin pi fraction of test_thin by chromosome',
        ylabel='Fraction (test_common_thin pi / test_thin pi)',
        filename=f'{output_prefix}.pi.common_fraction.bar.png',
        ylim=(0.0, None),
        figure_size=(14, 5),
        rotate_xlabels=45,
    )

    return merged, thin_bins, density
