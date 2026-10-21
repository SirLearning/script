"""Gene expression atlases normalized onto IWGSC RefSeq v1.1 gene IDs.

Two TPM sources sit under ``00data/12cs_tpm`` and they do not agree on gene ID
version:

``camp/``
    ``TraesCS1D02G000100`` — already v1.1, joins to the VEP GFF directly.
``ramirez2018_refseq1/``
    ``TraesCS1A01G000100`` — v1.0, plus 158k ``…LC`` low-confidence genes.

The v1.1 GFF carries ``previous_id=`` on every gene row, which is the official
v1.0 to v1.1 map, so no string surgery is needed. Do **not** normalize by
stripping the ``LC`` suffix: HC and LC genes are numbered independently in
v1.0, so ``TraesCS1A01G000100`` and ``TraesCS1A01G000100LC`` are different
genes and stripping merges them. The GFF is HC-only, so LC genes carry no
annotated variants and are dropped.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from infra.utils.graph import plot_categorical_bar_chart, plot_distribution_with_stats
from infra.utils.io import load_df_generic, save_df_to_tsv, save_thresholds

GENE_ID_COLUMN = "gene_id"
LEGACY_ID_COLUMN = "previous_gene_id"

# Non-sample columns seen in the TPM matrices.
NON_SAMPLE_COLUMNS = frozenset({"Geneid", "gene", "Length", "length"})

DEFAULT_DETECT_MIN = 0.5


def load_gene_id_map(gff_path: str | Path) -> pd.DataFrame:
    """Read gene rows from a v1.1 GFF3 into a ``gene_id`` / ``previous_gene_id`` table."""

    path = Path(gff_path)
    opener = gzip.open if path.suffix == ".gz" else open
    rows: list[dict[str, str]] = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            gene_id = ""
            previous = ""
            for attribute in fields[8].split(";"):
                if attribute.startswith("ID="):
                    gene_id = attribute[3:]
                elif attribute.startswith("previous_id="):
                    previous = attribute[12:]
            if gene_id:
                rows.append({GENE_ID_COLUMN: gene_id, LEGACY_ID_COLUMN: previous})
    if not rows:
        raise ValueError(f"No gene rows found in {gff_path}")
    return pd.DataFrame(rows)


def load_gene_loci(gff_path: str | Path) -> pd.DataFrame:
    """Gene midpoints on whole-chromosome labels (``chr1A``, ``Mit``).

    Official HC GFF already uses those names and linear coordinates. The Lulab
    GFF uses arm-segment ids; those rows are shifted with
    ``plink_segment_to_ref_position`` so a gene-level Manhattan lands on the
    same axis as a ``*_refchr`` pfile.
    """
    from infra.wheat.ref_v1 import (
        get_ref_v1_chr_name,
        is_ref_v1_chr_name,
        plink_segment_to_ref_position,
        pvar_chrom_to_plink_chr_id,
    )

    path = Path(gff_path)
    opener = gzip.open if path.suffix == ".gz" else open
    rows: list[dict] = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            gene_id = ""
            for attribute in fields[8].split(";"):
                if attribute.startswith("ID="):
                    gene_id = attribute[3:]
                    if gene_id.startswith("gene:"):
                        gene_id = gene_id[5:]
                    break
            if not gene_id:
                continue
            chrom_raw = fields[0].strip()
            start = int(fields[3])
            end = int(fields[4])
            if is_ref_v1_chr_name(chrom_raw):
                chrom = chrom_raw
                mid = (start + end) // 2
            else:
                plink_chr = pvar_chrom_to_plink_chr_id(chrom_raw)
                if plink_chr is None:
                    continue
                chrom = get_ref_v1_chr_name(plink_chr)
                _, start_lin = plink_segment_to_ref_position(plink_chr, start)
                _, end_lin = plink_segment_to_ref_position(plink_chr, end)
                mid = (start_lin + end_lin) // 2
            rows.append({
                GENE_ID_COLUMN: gene_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "pos": mid,
            })
    if not rows:
        raise ValueError(f"No gene rows found in {gff_path}")
    return pd.DataFrame(rows).drop_duplicates(subset=[GENE_ID_COLUMN], keep="first")


def load_tpm_matrix(path: str | Path) -> pd.DataFrame:
    """Read a TPM matrix; first column is the gene, the rest are samples."""

    table = load_df_generic(str(path))
    if table is None or table.empty:
        raise ValueError(f"Empty TPM matrix: {path}")
    gene_column = table.columns[0]
    sample_columns = [
        column
        for column in table.columns[1:]
        if column not in NON_SAMPLE_COLUMNS and not str(column).startswith("Unnamed")
    ]
    if not sample_columns:
        raise ValueError(f"TPM matrix has no sample columns: {path}")
    matrix = table[[gene_column, *sample_columns]].rename(columns={gene_column: GENE_ID_COLUMN})
    matrix[GENE_ID_COLUMN] = matrix[GENE_ID_COLUMN].astype(str)
    for column in sample_columns:
        matrix[column] = pd.to_numeric(matrix[column], errors="coerce").astype("float32")
    return matrix.drop_duplicates(subset=[GENE_ID_COLUMN], keep="first")


def normalize_gene_ids(matrix: pd.DataFrame, id_map: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Map a matrix onto v1.1 gene IDs, whichever version it arrived in.

    Returns the mapped matrix and a small report of how the rows resolved.
    """

    current = set(id_map[GENE_ID_COLUMN])
    legacy = id_map[id_map[LEGACY_ID_COLUMN].astype(bool)]
    legacy_to_current = dict(zip(legacy[LEGACY_ID_COLUMN], legacy[GENE_ID_COLUMN]))

    ids = matrix[GENE_ID_COLUMN]
    n_direct = int(ids.isin(current).sum())
    n_legacy = int(ids.isin(legacy_to_current).sum())
    n_low_confidence = int(ids.str.endswith("LC").sum())

    mapped = ids.where(ids.isin(current), ids.map(legacy_to_current))
    out = matrix.assign(**{GENE_ID_COLUMN: mapped})
    out = out[out[GENE_ID_COLUMN].notna()].drop_duplicates(subset=[GENE_ID_COLUMN], keep="first")
    report = {
        "n_input_genes": int(len(matrix)),
        "n_direct_v11": n_direct,
        "n_mapped_from_v10": n_legacy,
        "n_low_confidence_dropped": n_low_confidence,
        "n_mapped_genes": int(len(out)),
    }
    return out.reset_index(drop=True), report


def load_sample_groups(
    metadata_path: str | Path,
    *,
    sample_col: str,
    group_col: str,
) -> pd.Series:
    """Read a sample-to-group mapping (e.g. run accession to high-level tissue)."""

    table = load_df_generic(str(metadata_path))
    for column in (sample_col, group_col):
        if column not in table.columns:
            raise ValueError(f"metadata missing column {column!r}; has {list(table.columns)[:8]}")
    groups = table[[sample_col, group_col]].dropna().drop_duplicates(subset=[sample_col])
    return pd.Series(
        groups[group_col].astype(str).to_numpy(),
        index=groups[sample_col].astype(str),
        name=group_col,
    )


def tissue_specificity_tau(values: np.ndarray) -> np.ndarray:
    """Tau over the group means: 0 is uniform expression, 1 is one-group specific."""

    peak = np.nanmax(values, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        scaled = values / peak[:, None]
    n_groups = values.shape[1]
    if n_groups < 2:
        return np.full(values.shape[0], np.nan, dtype="float32")
    tau = np.nansum(1.0 - scaled, axis=1) / (n_groups - 1)
    return np.where(peak > 0, tau, np.nan).astype("float32")


def summarize_gene_expression(
    matrix: pd.DataFrame,
    *,
    sample_groups: pd.Series | None = None,
    detect_min: float = DEFAULT_DETECT_MIN,
) -> pd.DataFrame:
    """Per-gene expression summary, with group means and tau when groups are known."""

    sample_columns = [column for column in matrix.columns if column != GENE_ID_COLUMN]
    values = matrix[sample_columns].to_numpy(dtype="float32")
    summary = pd.DataFrame({GENE_ID_COLUMN: matrix[GENE_ID_COLUMN].to_numpy()})
    summary["max_tpm"] = np.nanmax(values, axis=1)
    summary["mean_tpm"] = np.nanmean(values, axis=1)
    summary["median_tpm"] = np.nanmedian(values, axis=1)
    summary["n_samples"] = len(sample_columns)
    summary["n_samples_expressed"] = (values >= detect_min).sum(axis=1)
    summary["expressed"] = summary["n_samples_expressed"] > 0

    if sample_groups is None:
        summary["breadth"] = np.nan
        summary["tau"] = np.nan
        return summary

    shared = [column for column in sample_columns if column in sample_groups.index]
    if not shared:
        raise ValueError("No TPM sample column matched the metadata sample ids")
    grouped = (
        matrix[shared]
        .T.groupby(sample_groups.reindex(shared).to_numpy())
        .mean()
        .T.astype("float32")
    )
    group_values = grouped.to_numpy(dtype="float32")
    summary["breadth"] = (group_values >= detect_min).sum(axis=1)
    summary["n_groups"] = group_values.shape[1]
    summary["tau"] = tissue_specificity_tau(group_values)
    for column in grouped.columns:
        summary[f"tpm_{column}"] = grouped[column].to_numpy()
    return summary


def plot_gene_expression(summary: pd.DataFrame, output_prefix: str | Path) -> None:
    """Peak-expression distribution, plus expressed-gene counts per group."""

    logged = summary.assign(log10_max_tpm=np.log10(summary["max_tpm"].clip(lower=0) + 1.0))
    plot_distribution_with_stats(
        data=logged,
        col="log10_max_tpm",
        title="Peak gene expression across samples",
        filename=f"{output_prefix}.gene_expr.max_tpm.png",
        x_label="log10(max TPM + 1)",
        median_val=float(logged["log10_max_tpm"].median()),
    )

    group_columns = [column for column in summary.columns if column.startswith("tpm_")]
    if not group_columns:
        return
    labels = [column[4:] for column in group_columns]
    counts = [float((summary[column] >= DEFAULT_DETECT_MIN).sum()) for column in group_columns]
    order = np.argsort(counts)[::-1]
    plot_categorical_bar_chart(
        categories=[labels[i] for i in order],
        values=[counts[i] for i in order],
        title="Genes expressed per group",
        ylabel="Genes",
        filename=f"{output_prefix}.gene_expr.by_group.png",
        rotate_xlabels=45,
    )


def run_gene_expression(
    tpm_table: str | Path,
    gff_path: str | Path,
    output_prefix: str | Path,
    *,
    metadata_path: str | Path | None = None,
    sample_col: str = "run_accession",
    group_col: str = "High level tissue",
    detect_min: float = DEFAULT_DETECT_MIN,
) -> pd.DataFrame:
    """Normalize one TPM matrix to v1.1 gene IDs and write the per-gene summary."""

    id_map = load_gene_id_map(gff_path)
    matrix, report = normalize_gene_ids(load_tpm_matrix(tpm_table), id_map)

    sample_groups = None
    if metadata_path:
        sample_groups = load_sample_groups(
            metadata_path,
            sample_col=sample_col,
            group_col=group_col,
        )
    summary = summarize_gene_expression(
        matrix,
        sample_groups=sample_groups,
        detect_min=detect_min,
    )

    save_df_to_tsv(summary, f"{output_prefix}.gene_expr.tsv")
    stats: dict[str, int | float] = dict(report)
    stats["detect_min"] = detect_min
    stats["n_expressed_genes"] = int(summary["expressed"].sum())
    stats["median_max_tpm"] = float(summary["max_tpm"].median())
    if summary["tau"].notna().any():
        stats["median_tau"] = float(summary["tau"].median())
        stats["n_group_specific"] = int((summary["tau"] >= 0.85).sum())
    save_thresholds(stats, f"{output_prefix}.gene_expr.th.tsv")
    plot_gene_expression(summary, output_prefix)
    print(
        f"gene_expr: input={report['n_input_genes']} mapped={report['n_mapped_genes']} "
        f"(direct={report['n_direct_v11']} via_v10={report['n_mapped_from_v10']} "
        f"LC={report['n_low_confidence_dropped']}) expressed={stats['n_expressed_genes']}"
    )
    return summary
