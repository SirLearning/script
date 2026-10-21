"""Rare-variant gene burden association for the WAP panel.

Single-variant GLM has no power on variants carried by a handful of accessions,
so the burden track collapses every variant a gene carries in one ``variant_sets``
set into a single per-sample score and tests that score against the trait.

Inputs come from tracks that already exist:

* genotypes as a PLINK2 ``--export Av`` ``.traw`` (variant-major) or a
  ``--export A`` ``.raw``, restricted to the set's ``*.extract`` list and to
  ``assoc.keep``. ``run_burden_test`` streams either file and never holds the
  full sample-by-variant matrix.
* the gene assignment from ``{id}.variant_sets.tsv`` (VEP ``Gene``), which is the
  only place variant-to-gene mapping lives
* the same phenotype table and structure PCs the single-variant track uses

Genes are grouped per subgenome. A/B/D homoeologs are *not* collapsed into triads:
the subgenomes live in separate pfiles and no triad map is wired in yet.

Aggregation is deliberately unweighted. Frequency weighting (Madsen-Browning,
beta(1,25)) only pays off once a gene holds many variants spanning a range of
frequencies; on the current sets most genes hold one.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from infra.utils.io import load_df_generic, save_df_to_tsv

COLLAPSE_COUNT = "count"
COLLAPSE_CARRIER = "carrier"
COLLAPSE_RULES = (COLLAPSE_COUNT, COLLAPSE_CARRIER)

GENE_COLUMN_CANDIDATES = ("Gene", "gene_id")
MISSING_GENE_TOKENS = frozenset({"", "-", "NA", "."})

# PLINK2 .raw carries six sample columns before the genotype block.
RAW_LEADING_COLUMNS = ("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
RAW_MISSING_TOKENS = frozenset({"NA", "na", "NaN", "nan", ".", ""})
# Row chunks keep peak RAM to tens of samples × N variants, not the full matrix.
RAW_ROW_CHUNK = 64
# --export Av is variant-major: a few thousand sites × the WAP panel is tiny.
TRAW_LEADING_COLUMNS = ("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT")
TRAW_VARIANT_CHUNK = 8192

# p-values need far more precision than the repo-wide 4-decimal default.
BURDEN_FLOAT_FORMAT = "%.6g"

BURDEN_COLUMNS = (
    "gene",
    "n_variants",
    "n_carriers",
    "carrier_fraction",
    "mean_missing_rate",
    "beta",
    "se",
    "t_stat",
    "p_value",
    "n_samples",
)


@dataclass(frozen=True)
class BurdenInputs:
    """Sample-aligned matrices feeding the association step."""

    samples: pd.Index
    scores: pd.DataFrame  # gene x sample
    trait: np.ndarray
    covariates: np.ndarray  # n_samples x n_covariates, intercept included
    gene_info: pd.DataFrame  # n_variants, mean_missing_rate, indexed by gene


def _strip_raw_allele_suffix(column, known_ids):
    """Recover the variant ID from a ``.raw`` column named ``{ID}_{allele}``.

    Variant IDs here are ``chrom-pos`` and carry no underscore, but the suffix is
    stripped by lookup rather than by blind splitting so an ID that does contain
    one still resolves.
    """
    if column in known_ids:
        return column
    head, sep, _ = column.rpartition("_")
    if sep and head in known_ids:
        return head
    return None


def load_plink_raw_genotypes(raw_path, keep_ids=None):
    """Load a PLINK2 ``--export A`` matrix as a samples x variant-column frame.

    Returns the frame indexed by IID with the six leading sample columns dropped.
    """
    raw = pd.read_csv(raw_path, sep=r"\s+")
    if "IID" not in raw.columns:
        raise ValueError(f"{raw_path} has no IID column; is it a PLINK2 --export A file?")
    genotype_columns = [c for c in raw.columns if c not in RAW_LEADING_COLUMNS]
    if not genotype_columns:
        raise ValueError(f"{raw_path} carries no genotype columns")
    out = raw.set_index(raw["IID"].astype(str))[genotype_columns]
    out.index.name = "Sample"
    if keep_ids is not None:
        out = out.loc[out.index.isin(set(keep_ids))]
    return out


def _raw_header_columns(raw_path):
    header = pd.read_csv(raw_path, sep=r"\s+", nrows=0)
    if "IID" not in header.columns:
        raise ValueError(f"{raw_path} has no IID column; is it a PLINK2 --export A file?")
    genotype_columns = [c for c in header.columns if c not in RAW_LEADING_COLUMNS]
    if not genotype_columns:
        raise ValueError(f"{raw_path} carries no genotype columns")
    return genotype_columns


def _map_raw_columns_to_genes(genotype_columns, gene_map):
    known_ids = set(gene_map["ID"])
    keep_columns = []
    variant_ids = []
    for column in genotype_columns:
        variant_id = _strip_raw_allele_suffix(column, known_ids)
        if variant_id is not None:
            keep_columns.append(column)
            variant_ids.append(variant_id)
    if not keep_columns:
        raise ValueError("No .raw genotype column matched a variant in the set")
    gene_by_id = gene_map.drop_duplicates(subset=["ID"]).set_index("ID")["gene"]
    genes = gene_by_id.loc[variant_ids]
    unique_genes = sorted(pd.unique(genes.to_numpy()))
    gene_index = {gene: i for i, gene in enumerate(unique_genes)}
    col_gene_i = np.fromiter((gene_index[gene] for gene in genes), dtype=np.int32, count=len(genes))
    return keep_columns, variant_ids, genes.to_numpy(), unique_genes, col_gene_i


def _iter_raw_dosage_chunks(raw_path, keep_columns, chunksize=RAW_ROW_CHUNK):
    dtype = {column: np.float32 for column in keep_columns}
    dtype["IID"] = str
    yield from pd.read_csv(
        raw_path,
        sep=r"\s+",
        usecols=["IID"] + keep_columns,
        dtype=dtype,
        na_values=list(RAW_MISSING_TOKENS),
        chunksize=chunksize,
        low_memory=True,
    )


def collapse_gene_scores_from_raw(raw_path, gene_map, collapse=COLLAPSE_COUNT, chunksize=RAW_ROW_CHUNK):
    """Collapse a PLINK2 ``.raw`` file without holding the full dosage matrix.

    Two sequential scans: allele frequency (for the major-ALT flip) then the
    per-gene sums. Each scan only materialises ``chunksize`` samples at a time.
    """
    if collapse not in COLLAPSE_RULES:
        raise ValueError(f"collapse must be one of {COLLAPSE_RULES}, got {collapse!r}")

    keep_columns, _variant_ids, _genes, unique_genes, col_gene_i = _map_raw_columns_to_genes(
        _raw_header_columns(raw_path), gene_map
    )
    n_var = len(keep_columns)
    n_genes = len(unique_genes)
    alt_sum = np.zeros(n_var, dtype=np.float64)
    n_called = np.zeros(n_var, dtype=np.int32)
    iids = []
    for chunk in _iter_raw_dosage_chunks(raw_path, keep_columns, chunksize=chunksize):
        iids.extend(chunk["IID"].astype(str).tolist())
        dosage = chunk[keep_columns].to_numpy(dtype=np.float64, copy=False)
        called = ~np.isnan(dosage)
        alt_sum += np.nansum(dosage, axis=0)
        n_called += called.sum(axis=0).astype(np.int32)

    n_samples = len(iids)
    frequency = np.zeros(n_var, dtype=np.float64)
    np.divide(alt_sum, 2.0 * n_called, out=frequency, where=n_called > 0)
    flip = frequency > 0.5
    missing_rate = 1.0 - (n_called.astype(np.float64) / max(n_samples, 1))

    cols_by_gene = [[] for _ in range(n_genes)]
    for j, gene_i in enumerate(col_gene_i):
        cols_by_gene[int(gene_i)].append(j)
    cols_by_gene = [np.asarray(cols, dtype=np.int32) for cols in cols_by_gene]

    scores = np.zeros((n_samples, n_genes), dtype=np.float64)
    offset = 0
    for chunk in _iter_raw_dosage_chunks(raw_path, keep_columns, chunksize=chunksize):
        n_row = len(chunk)
        dosage = chunk[keep_columns].to_numpy(dtype=np.float64, copy=True)
        if flip.any():
            dosage[:, flip] = 2.0 - dosage[:, flip]
        np.nan_to_num(dosage, copy=False, nan=0.0)
        for gene_i, cols in enumerate(cols_by_gene):
            scores[offset:offset + n_row, gene_i] += dosage[:, cols].sum(axis=1)
        offset += n_row

    if collapse == COLLAPSE_CARRIER:
        scores = (scores > 0).astype(np.float64)

    n_variants = np.bincount(col_gene_i, minlength=n_genes).astype(np.int64)
    miss_sum = np.bincount(col_gene_i, weights=missing_rate, minlength=n_genes)
    mean_missing = np.zeros(n_genes, dtype=np.float64)
    np.divide(miss_sum, n_variants, out=mean_missing, where=n_variants > 0)

    score_df = pd.DataFrame(
        scores, index=pd.Index(iids, name="Sample"), columns=unique_genes
    )
    gene_info = pd.DataFrame(
        {"n_variants": n_variants.astype(int), "mean_missing_rate": mean_missing},
        index=unique_genes,
    )
    return score_df, gene_info


def _traw_sample_columns(header_columns):
    """Map PLINK2 ``{FID}_{IID}`` header names to sample IIDs."""
    sample_columns = [c for c in header_columns if c not in TRAW_LEADING_COLUMNS]
    if not sample_columns:
        raise ValueError("traw file carries no sample columns")
    iids = []
    for column in sample_columns:
        fid, sep, iid = column.partition("_")
        iids.append(iid if sep else column)
    return sample_columns, iids


def collapse_gene_scores_from_traw(traw_path, gene_map, collapse=COLLAPSE_COUNT, chunksize=TRAW_VARIANT_CHUNK):
    """Collapse a PLINK2 ``--export Av`` ``.traw`` without a second file scan.

    Each row is one variant, so allele frequency, the major-ALT flip, and the
    per-gene sum happen on that row before the next chunk is read.
    """
    if collapse not in COLLAPSE_RULES:
        raise ValueError(f"collapse must be one of {COLLAPSE_RULES}, got {collapse!r}")

    header = pd.read_csv(traw_path, sep=r"\s+", nrows=0)
    if "SNP" not in header.columns:
        raise ValueError(f"{traw_path} has no SNP column; is it a PLINK2 --export Av file?")
    sample_columns, iids = _traw_sample_columns(header.columns)
    n_samples = len(iids)
    known_ids = set(gene_map["ID"])
    gene_by_id = gene_map.drop_duplicates(subset=["ID"]).set_index("ID")["gene"]
    unique_genes = sorted(pd.unique(gene_by_id.to_numpy()))
    gene_index = {gene: i for i, gene in enumerate(unique_genes)}
    n_genes = len(unique_genes)

    scores = np.zeros((n_samples, n_genes), dtype=np.float64)
    n_variants = np.zeros(n_genes, dtype=np.int64)
    miss_sum = np.zeros(n_genes, dtype=np.float64)
    matched = 0
    dtype = {column: np.float32 for column in sample_columns}
    dtype["SNP"] = str
    for chunk in pd.read_csv(
        traw_path,
        sep=r"\s+",
        usecols=["SNP"] + sample_columns,
        dtype=dtype,
        na_values=list(RAW_MISSING_TOKENS),
        chunksize=chunksize,
        low_memory=True,
    ):
        dosage = chunk[sample_columns].to_numpy(dtype=np.float64, copy=True)
        called = ~np.isnan(dosage)
        n_called = called.sum(axis=1)
        frequency = np.zeros(len(chunk), dtype=np.float64)
        np.divide(np.nansum(dosage, axis=1), 2.0 * n_called, out=frequency, where=n_called > 0)
        flip = frequency > 0.5
        if flip.any():
            dosage[flip] = 2.0 - dosage[flip]
        np.nan_to_num(dosage, copy=False, nan=0.0)
        missing_rate = 1.0 - (n_called.astype(np.float64) / max(n_samples, 1))
        for i, snp in enumerate(chunk["SNP"].astype(str)):
            variant_id = snp if snp in known_ids else _strip_raw_allele_suffix(snp, known_ids)
            if variant_id is None:
                continue
            gene = gene_by_id.loc[variant_id]
            gene_i = gene_index[gene]
            scores[:, gene_i] += dosage[i]
            n_variants[gene_i] += 1
            miss_sum[gene_i] += missing_rate[i]
            matched += 1

    if matched == 0:
        raise ValueError("No .traw SNP column matched a variant in the set")

    keep = n_variants > 0
    unique_genes = [gene for gene, flag in zip(unique_genes, keep) if flag]
    scores = scores[:, keep]
    n_variants = n_variants[keep]
    miss_sum = miss_sum[keep]
    mean_missing = np.zeros(len(unique_genes), dtype=np.float64)
    np.divide(miss_sum, n_variants, out=mean_missing, where=n_variants > 0)

    if collapse == COLLAPSE_CARRIER:
        scores = (scores > 0).astype(np.float64)

    score_df = pd.DataFrame(
        scores, index=pd.Index(iids, name="Sample"), columns=unique_genes
    )
    gene_info = pd.DataFrame(
        {"n_variants": n_variants.astype(int), "mean_missing_rate": mean_missing},
        index=unique_genes,
    )
    return score_df, gene_info


def load_gene_map(variant_sets_table, set_name):
    """Map variant ID to gene for one ``variant_sets`` set column.

    Variants the set excludes, and variants VEP left without a gene, are dropped:
    a burden score is only defined for variants that belong to a gene.
    """
    # load_df_generic sniffs the suffix with str.endswith, so it needs a real string.
    sets = load_df_generic(str(variant_sets_table))
    if sets is None or sets.empty:
        raise ValueError(f"Empty variant_sets table: {variant_sets_table}")
    set_column = set_name if set_name.startswith("set_") else f"set_{set_name}"
    if set_column not in sets.columns:
        available = sorted(c for c in sets.columns if c.startswith("set_"))
        raise ValueError(
            f"{variant_sets_table} has no column {set_column}; available: {', '.join(available)}"
        )
    gene_column = next((c for c in GENE_COLUMN_CANDIDATES if c in sets.columns), None)
    if gene_column is None:
        raise ValueError(f"{variant_sets_table} has no Gene column; burden needs gene grouping")
    if "ID" not in sets.columns:
        raise ValueError(f"{variant_sets_table} has no ID column")

    in_set = sets[set_column].astype(str).str.lower().isin({"true", "1"})
    out = sets.loc[in_set, ["ID", gene_column]].rename(columns={gene_column: "gene"})
    out["ID"] = out["ID"].astype(str)
    out["gene"] = out["gene"].astype(str).str.strip()
    out = out.loc[~out["gene"].isin(MISSING_GENE_TOKENS)]
    return out.drop_duplicates(subset=["ID", "gene"]).reset_index(drop=True)


def collapse_gene_scores(genotypes, gene_map, collapse=COLLAPSE_COUNT):
    """Collapse per-variant dosages into one score per gene per sample.

    Dosages are flipped where the counted allele is the major one, so the score
    always counts minor alleles. Missing calls count as non-carrier, which is the
    conservative reading for a rare variant.

    ``count`` sums minor alleles across the gene; ``carrier`` reduces that to a
    0/1 indicator.
    """
    if collapse not in COLLAPSE_RULES:
        raise ValueError(f"collapse must be one of {COLLAPSE_RULES}, got {collapse!r}")

    known_ids = set(gene_map["ID"])
    column_to_id = {}
    for column in genotypes.columns:
        variant_id = _strip_raw_allele_suffix(column, known_ids)
        if variant_id is not None:
            column_to_id[column] = variant_id
    if not column_to_id:
        raise ValueError("No .raw genotype column matched a variant in the set")

    dosage = genotypes[list(column_to_id)].apply(pd.to_numeric, errors="coerce")
    # --export A counts the ALT allele, which is not always the minor one. Order
    # matters twice here: the frequency must ignore missing calls, and the flip
    # must run before they are zero-filled. Flipping a zero-filled missing call
    # yields 2 - 0 = 2, inventing a homozygous carrier of the rare allele.
    frequency = dosage.mean(axis=0) / 2.0
    flip = frequency > 0.5
    if flip.any():
        dosage.loc[:, flip] = 2.0 - dosage.loc[:, flip]
    missing_rate = dosage.isna().mean(axis=0)
    dosage = dosage.fillna(0.0)

    id_by_column = pd.Series(column_to_id)
    gene_by_id = gene_map.set_index("ID")["gene"]
    gene_by_column = id_by_column.map(gene_by_id)

    scores = dosage.T.groupby(gene_by_column.values).sum().T
    if collapse == COLLAPSE_CARRIER:
        scores = (scores > 0).astype(float)
    gene_info = pd.DataFrame({
        "n_variants": id_by_column.groupby(gene_by_column.values).size(),
        # Missing calls score as non-carrier, so a gene built from poorly called
        # variants is biased toward the null. Carry the rate so that is visible.
        "mean_missing_rate": missing_rate.groupby(gene_by_column.values).mean(),
    }).reindex(scores.columns)
    gene_info["n_variants"] = gene_info["n_variants"].astype(int)
    return scores, gene_info


def _load_trait(phenotype_table, trait):
    pheno = load_df_generic(str(phenotype_table))
    if pheno is None or pheno.empty:
        raise ValueError(f"Empty phenotype table: {phenotype_table}")
    pheno = pheno.rename(columns={"#FID": "FID", "#IID": "IID"})
    if "IID" not in pheno.columns:
        raise ValueError(f"{phenotype_table} has no IID column")
    if trait not in pheno.columns:
        available = ", ".join(c for c in pheno.columns if c not in {"FID", "IID"})
        raise ValueError(f"{phenotype_table} has no trait {trait}; available: {available}")
    out = pheno[["IID", trait]].copy()
    out["IID"] = out["IID"].astype(str)
    out[trait] = pd.to_numeric(out[trait], errors="coerce")
    return out.dropna(subset=[trait]).set_index("IID")[trait]


def _load_covariates(covar_table):
    if covar_table is None:
        return None
    covar = load_df_generic(str(covar_table))
    if covar is None or covar.empty:
        raise ValueError(f"Empty covariate table: {covar_table}")
    covar = covar.rename(columns={"#FID": "FID", "#IID": "IID"})
    id_column = "IID" if "IID" in covar.columns else None
    if id_column is None:
        raise ValueError(f"{covar_table} has no IID column")
    value_columns = [c for c in covar.columns if c not in {"FID", "IID"}]
    if not value_columns:
        raise ValueError(f"{covar_table} carries no covariate columns")
    out = covar[[id_column] + value_columns].copy()
    out[id_column] = out[id_column].astype(str)
    out = out.set_index(id_column)
    return out.apply(pd.to_numeric, errors="coerce").dropna()


def align_burden_inputs(scores, gene_info, trait_series, covariates=None):
    """Intersect genotype, phenotype and covariate samples into aligned matrices."""
    samples = scores.index.intersection(trait_series.index)
    if covariates is not None:
        samples = samples.intersection(covariates.index)
    if len(samples) < 3:
        raise ValueError(
            f"Only {len(samples)} samples shared by genotypes, phenotype and covariates"
        )
    samples = pd.Index(sorted(samples))

    design = [np.ones((len(samples), 1))]
    if covariates is not None:
        design.append(covariates.loc[samples].to_numpy(dtype=float))
    return BurdenInputs(
        samples=samples,
        scores=scores.loc[samples].T,
        trait=trait_series.loc[samples].to_numpy(dtype=float),
        covariates=np.hstack(design),
        gene_info=gene_info,
    )


def run_burden_association(inputs):
    """Regress the trait on each gene score, adjusting for the covariates.

    Covariates are projected out once and each gene is then a simple regression on
    the residuals (Frisch-Waugh-Lovell), which is exact for the full model and
    keeps the per-gene cost to a couple of dot products.
    """
    q, _ = np.linalg.qr(inputs.covariates)
    n_samples = len(inputs.samples)
    df_resid = n_samples - q.shape[1] - 1
    if df_resid <= 0:
        raise ValueError(f"No residual degrees of freedom: {n_samples} samples")

    y = inputs.trait
    y_resid = y - q @ (q.T @ y)
    y_ss = float(y_resid @ y_resid)

    rows = []
    for gene, score in inputs.scores.iterrows():
        g = score.to_numpy(dtype=float)
        n_carriers = int((g > 0).sum())
        g_resid = g - q @ (q.T @ g)
        g_ss = float(g_resid @ g_resid)
        info = inputs.gene_info.loc[gene] if gene in inputs.gene_info.index else {}
        row = {
            "gene": gene,
            "n_variants": int(info.get("n_variants", 0)),
            "n_carriers": n_carriers,
            "carrier_fraction": n_carriers / n_samples,
            "mean_missing_rate": float(info.get("mean_missing_rate", np.nan)),
            "n_samples": n_samples,
        }
        # A gene whose score is constant once covariates are removed carries no
        # contrast; report it rather than dividing by zero.
        if g_ss <= 0 or n_carriers == 0:
            row.update({"beta": np.nan, "se": np.nan, "t_stat": np.nan, "p_value": np.nan})
            rows.append(row)
            continue
        beta = float(g_resid @ y_resid) / g_ss
        rss = y_ss - beta * beta * g_ss
        se = float(np.sqrt(max(rss, 0.0) / df_resid / g_ss))
        t_stat = beta / se if se > 0 else np.nan
        p_value = float(2.0 * stats.t.sf(abs(t_stat), df_resid)) if se > 0 else np.nan
        row.update({"beta": beta, "se": se, "t_stat": t_stat, "p_value": p_value})
        rows.append(row)

    out = pd.DataFrame(rows, columns=list(BURDEN_COLUMNS))
    return out.sort_values("p_value", na_position="last").reset_index(drop=True)


def run_burden_test(
    raw_path,
    variant_sets_table,
    set_name,
    phenotype_table,
    trait,
    output_prefix,
    covar_table=None,
    collapse=COLLAPSE_COUNT,
):
    """Gene burden association for one subgenome and one variant set.

    Writes ``{output_prefix}.burden.tsv`` and returns the result frame.
    """
    gene_map = load_gene_map(variant_sets_table, set_name)
    print(f"[Info] {set_name}: {len(gene_map)} gene-assigned variants in the set")
    if gene_map.empty:
        raise ValueError(f"Set {set_name} has no gene-assigned variants; nothing to collapse")

    collapse_fn = (
        collapse_gene_scores_from_traw
        if str(raw_path).endswith(".traw")
        else collapse_gene_scores_from_raw
    )
    scores, gene_info = collapse_fn(raw_path, gene_map, collapse=collapse)
    print(
        f"[Info] collapsed {int(gene_info['n_variants'].sum())} variants into "
        f"{scores.shape[1]} genes ({collapse}) across {scores.shape[0]} samples; "
        f"mean {gene_info['n_variants'].mean():.2f} variants per gene, "
        f"mean missing rate {gene_info['mean_missing_rate'].mean():.3f}"
    )

    trait_series = _load_trait(phenotype_table, trait)
    covariates = _load_covariates(covar_table)
    inputs = align_burden_inputs(scores, gene_info, trait_series, covariates)
    print(
        f"[Info] testing {inputs.scores.shape[0]} genes on {len(inputs.samples)} samples "
        f"with {inputs.covariates.shape[1] - 1} covariates"
    )

    results = run_burden_association(inputs)
    tested = results["p_value"].notna().sum()
    print(f"[Info] {tested} genes testable; {int((results['p_value'] < 0.05).sum())} at p<0.05")
    save_df_to_tsv(
        results, f"{output_prefix}.burden.tsv", float_format=BURDEN_FLOAT_FORMAT, na_rep="NA"
    )
    return results
