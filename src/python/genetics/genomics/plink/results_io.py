"""Load PLINK / PLINK2 / GCTA tabular outputs (format parsing only; no association math)."""

from collections import Counter
from pathlib import Path

import pandas as pd
import numpy as np

from infra.utils.io import load_df_generic, save_df_to_tsv
from infra.wheat.ref_v1 import all_ref_v1_plink_chr_ids, get_ref_v1_chr_name, get_ref_v1_genome_segment_layout


def load_plink_eigenvec(path):
    """Load PLINK 1.9 or 2 ``.eigenvec`` (headerless FID IID PC… or ``#FID`` header)."""

    with open(path, encoding='utf-8') as handle:
        first_line = handle.readline()
    if first_line.lstrip().startswith('#'):
        df = load_df_generic(path)
    else:
        df = pd.read_csv(path, sep=r'\s+', header=None, comment='#')
        if df.shape[1] < 3:
            raise ValueError(f'Expected FID IID PC… columns in eigenvec: {path}')
        df.columns = ['FID', 'Sample'] + [f'PC{i + 1}' for i in range(df.shape[1] - 2)]
    if df is None or df.empty:
        raise ValueError(f'Empty eigenvec: {path}')
    if '#FID' in df.columns:
        df = df.rename(columns={'#FID': 'FID'})
    if '#IID' in df.columns and 'Sample' not in df.columns:
        df = df.rename(columns={'#IID': 'Sample'})
    elif 'IID' in df.columns and 'Sample' not in df.columns:
        df = df.rename(columns={'IID': 'Sample'})
    return df


def load_plink2_eigenvec(path):
    """Backward-compatible alias for :func:`load_plink_eigenvec`."""

    return load_plink_eigenvec(path)


def load_plink2_eigenval(path):
    vals = []
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals.append(float(line.split()[0]))
    if not vals:
        raise ValueError(f'No eigenvalues in {path}')
    total = float(sum(vals))
    ratios = [v / total if total > 0 else 0.0 for v in vals]
    return pd.DataFrame({
        'PC': [f'PC{i + 1}' for i in range(len(vals))],
        'Eigenvalue': vals,
        'ExplainedVarianceRatio': ratios,
    })


def load_plink2_glm(path):
    df = load_df_generic(path)
    if df is None or df.empty:
        raise ValueError(f'Empty GLM result: {path}')
    rename = {'#CHROM': 'CHR', 'POS': 'BP', 'ID': 'SNP', 'P': 'PVALUE', 'OBS_CT': 'N'}
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    if 'PVALUE' not in df.columns and 'P' in df.columns:
        df = df.rename(columns={'P': 'PVALUE'})
    df['PVALUE'] = pd.to_numeric(df['PVALUE'], errors='coerce')
    df = df.dropna(subset=['PVALUE']).sort_values('PVALUE').reset_index(drop=True)
    df['Index'] = np.arange(1, len(df) + 1)
    df['LOG10P'] = -np.log10(df['PVALUE'].clip(lower=1e-300))
    return df


def load_gcta_gwas(path):
    df = load_df_generic(path)
    if df is None or df.empty:
        raise ValueError(f'Empty GCTA GWAS result: {path}')
    if 'P' in df.columns and 'PVALUE' not in df.columns:
        df = df.rename(columns={'P': 'PVALUE'})
    df['PVALUE'] = pd.to_numeric(df['PVALUE'], errors='coerce')
    df = df.dropna(subset=['PVALUE']).sort_values('PVALUE').reset_index(drop=True)
    df['Index'] = np.arange(1, len(df) + 1)
    df['LOG10P'] = -np.log10(df['PVALUE'].clip(lower=1e-300))
    return df


def summarize_tagsnp_from_prune(prune_in_path, output_prefix, max_tags=1000):
    tags = []
    with open(prune_in_path, encoding='utf-8') as handle:
        for line in handle:
            vid = line.strip()
            if vid:
                tags.append(vid)
            if len(tags) >= int(max_tags):
                break
    out = pd.DataFrame({'TagSNP': tags})
    save_df_to_tsv(out, f'{output_prefix}.tagsnp.tsv')
    return out


def _count_variants_in_pvar(pvar_path):
    counts = Counter()
    with open(pvar_path, encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            chrom = line.split('\t', 1)[0]
            if chrom:
                counts[chrom] += 1
    return counts


def _plink_chr_from_thin_pvar_name(pvar_name):
    """Parse PLINK chr id from chrNNN.common.thin.pvar basename."""
    stem = Path(pvar_name).name.replace('.common.thin.pvar', '')
    num = stem.replace('chr', '').lstrip('0')
    return num if num else '0'


def _collect_chr_variant_counts(process_dir):
    root = Path(process_dir)
    thin_pvars = sorted(root.glob('chr*.common.thin.pvar'))
    if thin_pvars:
        counter = Counter()
        for pvar in thin_pvars:
            plink_chr = _plink_chr_from_thin_pvar_name(pvar.name)
            counter[plink_chr] = sum(
                1 for line in open(pvar, encoding='utf-8') if not line.startswith('#')
            )
        return counter, 'per_chr_thin_pvar'

    merged_pvars = sorted(root.glob('*_test.plink2.pvar'))
    if not merged_pvars:
        merged_pvars = sorted(root.glob('*.plink2.pvar'))
    if not merged_pvars:
        raise FileNotFoundError(
            f'No chr*.common.thin.pvar or *_test.plink2.pvar under {process_dir}'
        )

    counter = Counter()
    for pvar in merged_pvars:
        counter.update(_count_variants_in_pvar(pvar))
    return counter, 'merged_subgenome_pvar'


def _ref_name_sort_key(ref_name):
    """chrUn first, chr1A–chr7D in name order, Mit/Chl last (before total)."""
    if ref_name == 'chrUn':
        return (0, '')
    if ref_name == 'Mit':
        return (2, 0)
    if ref_name == 'Chl':
        return (2, 1)
    return (1, ref_name)


def _plink_chr_row_sort_key(plink_chr):
    """PLINK chr 0 first, 1–42 numeric, Mit (43) and Chl (44) last (before total)."""
    if plink_chr == 'total':
        return (3, 0)
    if plink_chr == '0':
        return (0, 0)
    if plink_chr == '43':
        return (2, 0)
    if plink_chr == '44':
        return (2, 1)
    return (1, int(plink_chr))


def _append_total_row(df, count_col, label_col, total_label='total'):
    total = int(df[count_col].sum())
    total_row = {label_col: total_label, count_col: total}
    for col in df.columns:
        if col not in total_row:
            total_row[col] = ''
    return pd.concat([df, pd.DataFrame([total_row])], ignore_index=True)


def load_ordered_chr_variant_by_ref(tsv_path):
    """Load ``*.chr_variant_counts.by_ref.tsv`` excluding the total row, in display order."""
    df = load_df_generic(tsv_path)
    if df is None or df.empty:
        raise ValueError(f'Empty chr variant counts: {tsv_path}')
    df = df.loc[df['ref_name'] != 'total'].copy()
    df['n_variants'] = pd.to_numeric(df['n_variants'], errors='coerce').fillna(0).astype(int)
    return df.sort_values(
        by='ref_name',
        key=lambda s: s.map(_ref_name_sort_key),
    ).reset_index(drop=True)


def summarize_chr_variant_counts(process_dir, output_prefix):
    """
    Count variant sites per PLINK chromosome under a process directory.

    Prefers per-chromosome ``chr*.common.thin.pvar`` when present; otherwise
    aggregates ``*_test.plink2.pvar`` (or ``*.plink2.pvar``) merged subgenome files.
    Maps PLINK chr ids to vmap4 ref names via ``get_ref_v1_chr_name`` (infra_ref_v1.nf parity).

    Writes:
      - ``{output_prefix}.chr_variant_counts.tsv`` (plink_chr, ref_name, n_variants)
      - ``{output_prefix}.chr_variant_counts.by_ref.tsv`` (ref_name, n_variants)
    """
    counter, source = _collect_chr_variant_counts(process_dir)

    per_chr_rows = []
    for plink_chr in all_ref_v1_plink_chr_ids():
        per_chr_rows.append({
            'plink_chr': str(plink_chr),
            'ref_name': get_ref_v1_chr_name(plink_chr),
            'n_variants': int(counter.get(str(plink_chr), 0)),
            'source': source,
        })
    per_chr = pd.DataFrame(per_chr_rows)
    per_chr = per_chr.sort_values(
        by='plink_chr',
        key=lambda s: s.map(_plink_chr_row_sort_key),
    ).reset_index(drop=True)
    per_chr = _append_total_row(per_chr, 'n_variants', 'plink_chr')
    per_chr.loc[per_chr['plink_chr'] == 'total', 'ref_name'] = 'total'
    save_df_to_tsv(per_chr, f'{output_prefix}.chr_variant_counts.tsv')

    by_ref = (
        per_chr.loc[per_chr['plink_chr'] != 'total']
        .groupby('ref_name', as_index=False)['n_variants']
        .sum()
    )
    by_ref = by_ref.sort_values(
        by='ref_name',
        key=lambda s: s.map(_ref_name_sort_key),
    ).reset_index(drop=True)
    by_ref = _append_total_row(by_ref, 'n_variants', 'ref_name')
    save_df_to_tsv(by_ref, f'{output_prefix}.chr_variant_counts.by_ref.tsv')
    return per_chr, by_ref


def iter_process_pvar_positions(process_dir):
    """Yield (plink_chr_id, pos) for every variant row under a process directory."""
    root = Path(process_dir)
    thin_pvars = sorted(root.glob('chr*.common.thin.pvar'))
    if thin_pvars:
        for pvar in thin_pvars:
            plink_chr = int(_plink_chr_from_thin_pvar_name(pvar.name))
            with open(pvar, encoding='utf-8') as handle:
                for line in handle:
                    if line.startswith('#'):
                        continue
                    pos = int(line.split('\t', 3)[1])
                    yield plink_chr, pos
        return

    merged_pvars = sorted(root.glob('*_test.plink2.pvar'))
    if not merged_pvars:
        merged_pvars = sorted(root.glob('*.plink2.pvar'))
    if not merged_pvars:
        raise FileNotFoundError(
            f'No chr*.common.thin.pvar or *_test.plink2.pvar under {process_dir}'
        )
    for pvar in merged_pvars:
        with open(pvar, encoding='utf-8') as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                parts = line.split('\t', 3)
                yield int(parts[0]), int(parts[1])


def compute_genome_bin_counts(process_dir, bin_size_bp=5_000_000):
    """
    Bin variants into fixed-width windows along the ordered vmap4 genome (PLINK 0-44).

    Returns a DataFrame with bin_index, ref_name, bin_start_mb, bin_end_mb, n_variants,
    and density_per_mb (variants per Mb of bin width).
    """
    segments, total_bp = get_ref_v1_genome_segment_layout()
    starts = {seg['plink_chr']: seg['start_bp'] for seg in segments}
    n_bins = (total_bp + bin_size_bp - 1) // bin_size_bp
    counts = np.zeros(n_bins, dtype=np.int64)

    for plink_chr, pos in iter_process_pvar_positions(process_dir):
        global_bp = starts[plink_chr] + pos - 1
        if global_bp < 0 or global_bp >= total_bp:
            continue
        counts[global_bp // bin_size_bp] += 1

    bin_indices = np.arange(n_bins)
    bin_starts_mb = bin_indices * bin_size_bp / 1e6
    ref_names = [_ref_name_at_global_bp(segments, int(i * bin_size_bp)) for i in bin_indices]
    return pd.DataFrame({
        'bin_index': bin_indices,
        'ref_name': ref_names,
        'bin_start_mb': bin_starts_mb,
        'bin_end_mb': bin_starts_mb + bin_size_bp / 1e6,
        'n_variants': counts,
        'density_per_mb': counts.astype(float) * (1_000_000 / bin_size_bp),
    })


def _ref_name_at_global_bp(segments, bp):
    for seg in segments:
        if seg['start_bp'] <= bp < seg['end_bp']:
            return seg['ref_name']
    return segments[-1]['ref_name']


def compute_genome_mb_bin_counts(process_dir, bin_size_bp=5_000_000):
    """Backward-compatible alias for :func:`compute_genome_bin_counts`."""
    return compute_genome_bin_counts(process_dir, bin_size_bp=bin_size_bp)


def merge_thin_common_chr_by_ref(thin_by_ref_tsv, common_by_ref_tsv):
    """Merge per-ref chromosome counts and compute common/thin fraction."""
    thin = load_ordered_chr_variant_by_ref(thin_by_ref_tsv).rename(
        columns={'n_variants': 'test_thin'}
    )
    common = load_ordered_chr_variant_by_ref(common_by_ref_tsv).rename(
        columns={'n_variants': 'test_common_thin'}
    )
    merged = thin.merge(common, on='ref_name', how='outer').fillna(0)
    merged['test_thin'] = merged['test_thin'].astype(int)
    merged['test_common_thin'] = merged['test_common_thin'].astype(int)
    merged = merged.sort_values(
        by='ref_name',
        key=lambda s: s.map(_ref_name_sort_key),
    ).reset_index(drop=True)
    merged['common_fraction'] = np.where(
        merged['test_thin'] > 0,
        merged['test_common_thin'] / merged['test_thin'],
        0.0,
    )
    return merged


def build_genome_density_compare_table(thin_bins, common_bins):
    """Align two genome bin tables on bin_index for stacked density plotting."""
    table = thin_bins[
        ['bin_index', 'ref_name', 'bin_start_mb', 'bin_end_mb', 'density_per_mb']
    ].rename(columns={'density_per_mb': 'test_thin'}).copy()
    common_by_bin = common_bins.set_index('bin_index')['density_per_mb']
    table['test_common_thin'] = table['bin_index'].map(common_by_bin).fillna(0.0)
    return table.sort_values('bin_index').reset_index(drop=True)


def plot_thin_common_by_ref_distribution_line(
    merged,
    filename,
    title,
    y_label,
    figure_size=(14, 5),
):
    """Dual-line per-ref chromosome plot (test_thin vs test_common_thin)."""
    from infra.utils.graph import plot_multi_line_series

    plot_multi_line_series(
        merged,
        x_col='ref_name',
        y_specs=[
            {'y_col': 'test_thin', 'label': 'test_thin', 'color': '#4c72b0', 'linewidth': 2.0},
            {'y_col': 'test_common_thin', 'label': 'test_common_thin', 'color': '#c44e52', 'linewidth': 2.0},
        ],
        title=title,
        filename=filename,
        x_label='Chromosome',
        y_label=y_label,
        figure_size=figure_size,
        rotate_xlabels=45,
    )


def prepare_plink_phenotype_table(phenotype_file, trait, output_path):
    pheno = load_df_generic(phenotype_file)
    if pheno is None or pheno.empty:
        raise ValueError(f'Empty phenotype: {phenotype_file}')
    if 'Sample' not in pheno.columns:
        raise ValueError('Phenotype table requires Sample column')
    if trait not in pheno.columns:
        raise ValueError(f'Missing trait column: {trait}')
    out = pheno[['Sample', trait]].copy()
    out.columns = ['IID', trait]
    out.insert(0, '#FID', out['IID'].astype(str))
    save_df_to_tsv(out, output_path)
    return output_path
