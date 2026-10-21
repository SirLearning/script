"""Plot GWAS association results from PLINK2 --glm or GCTA."""

import numpy as np
import pandas as pd

from genetics.genomics.plink.results_io import load_gcta_gwas, load_plink2_glm
from infra.utils.graph import plot_gwas_qq, plot_manhattan
from infra.utils.io import save_df_to_tsv

GENOME_WIDE_P = 5e-8


def plot_gwas_results(association_path, output_prefix, source='plink2'):
    if source == 'gcta':
        out = load_gcta_gwas(association_path)
    else:
        out = load_plink2_glm(association_path)
    save_df_to_tsv(out, f'{output_prefix}.gwas.tsv')
    exp = -np.log10((np.arange(1, len(out) + 1) - 0.5) / (len(out) + 0.5))
    obs = -np.log10(np.sort(out['PVALUE'].values))
    plot_gwas_qq(exp, obs, title='GWAS QQ Plot', filename=f'{output_prefix}.qq.png')
    plot_manhattan(
        chromosomes=out['CHR'], positions=out['BP'], values=out['LOG10P'],
        title='GWAS Manhattan', filename=f'{output_prefix}.manhattan.png',
        thresholds={f'P = {GENOME_WIDE_P:g}': -np.log10(GENOME_WIDE_P)},
    )


def plot_kgwas_results(glm_path, output_prefix):
    plot_gwas_results(glm_path, output_prefix, source='plink2')


def plot_burden_results(burden_path, output_prefix, gff_path):
    """Gene-level QQ and Manhattan for a ``*.burden.tsv``.

    Points sit at the GFF gene midpoint. The dashed line is Bonferroni
    ``0.05 / n_tested`` — genome-wide 5e-8 is a single-variant threshold and
    does not apply to a few thousand genes.
    """
    from genetics.genomics.gene.expression import GENE_ID_COLUMN, load_gene_loci
    from infra.utils.io import load_df_generic

    burden = load_df_generic(str(burden_path))
    if burden is None or burden.empty:
        raise ValueError(f"Empty burden table: {burden_path}")
    if "gene" not in burden.columns or "p_value" not in burden.columns:
        raise ValueError(f"{burden_path} needs gene and p_value columns")

    tested = burden.copy()
    tested["p_value"] = pd.to_numeric(tested["p_value"], errors="coerce")
    tested = tested.dropna(subset=["p_value"]).reset_index(drop=True)
    if tested.empty:
        raise ValueError(f"{burden_path} has no finite p-values")
    tested["LOG10P"] = -np.log10(tested["p_value"].clip(lower=1e-300))

    n = len(tested)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / (n + 0.5))
    obs = -np.log10(np.sort(tested["p_value"].to_numpy()))
    plot_gwas_qq(exp, obs, title="Burden QQ Plot", filename=f"{output_prefix}.qq.png")

    loci = load_gene_loci(gff_path)
    placed = tested.merge(loci, left_on="gene", right_on=GENE_ID_COLUMN, how="left")
    missing = int(placed["pos"].isna().sum())
    if missing:
        print(f"[Info] {missing} of {n} genes have no GFF locus and are omitted from the Manhattan")
    placed = placed.dropna(subset=["chrom", "pos"]).reset_index(drop=True)
    if placed.empty:
        raise ValueError(f"No burden gene matched a locus in {gff_path}")

    bonferroni = 0.05 / n
    plot_manhattan(
        chromosomes=placed["chrom"],
        positions=placed["pos"],
        values=placed["LOG10P"],
        title="Burden Manhattan",
        filename=f"{output_prefix}.manhattan.png",
        thresholds={f"Bonferroni {bonferroni:.2e}": -np.log10(bonferroni)},
    )
    keep = [c for c in ("gene", "chrom", "pos", "n_variants", "n_carriers", "p_value", "LOG10P") if c in placed.columns]
    save_df_to_tsv(
        placed[keep],
        f"{output_prefix}.burden.plot.tsv",
        float_format="%.6g",
        na_rep="NA",
    )
    return placed
