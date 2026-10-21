"""WAP phenotype matrix aligned to the association keep list.

Turns a raw phenotype table into a PLINK2 ``--pheno`` file whose sample order
and membership follow ``assoc.keep`` (WAP only, per the cohort contract in
``doc/project_knowledge/domain/sample_cohorts.yaml``).

This module deliberately performs no mixed-model adjustment: BLUE/BLUP across
field replicates and environments is a separate decision and belongs upstream
of this alignment step. Rows are taken as given; duplicated sample IDs are
averaged only when ``aggregate_duplicates`` is set.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from infra.utils.graph import plot_bar_chart
from infra.utils.io import load_df_generic, save_df_to_tsv, save_thresholds

FID_COLUMN = "#FID"
IID_COLUMN = "IID"
DEFAULT_ID_COLUMN = "IID"
DEFAULT_FID = "0"
PLINK_MISSING = "NA"


def load_keep_samples(keep_file: str | Path) -> pd.DataFrame:
    """Read a PLINK2 keep file (FID IID, no header) preserving its order."""

    raw = pd.read_csv(keep_file, sep=r"\s+", header=None, dtype=str)
    if raw.empty:
        raise ValueError(f"Empty keep file: {keep_file}")
    if raw.shape[1] == 1:
        return pd.DataFrame({FID_COLUMN: DEFAULT_FID, IID_COLUMN: raw.iloc[:, 0].astype(str)})
    return pd.DataFrame(
        {
            FID_COLUMN: raw.iloc[:, 0].astype(str),
            IID_COLUMN: raw.iloc[:, 1].astype(str),
        }
    )


def load_phenotype_table(
    phenotype_file: str | Path,
    *,
    id_col: str = DEFAULT_ID_COLUMN,
) -> pd.DataFrame:
    """Load a wide phenotype table and normalise its sample identifier column."""

    table = load_df_generic(str(phenotype_file))
    if table is None or table.empty:
        raise ValueError(f"Empty phenotype table: {phenotype_file}")
    if id_col not in table.columns:
        raise ValueError(f"Phenotype table missing id column {id_col}: {phenotype_file}")
    out = table.copy()
    out[id_col] = out[id_col].astype(str)
    if id_col != IID_COLUMN:
        out = out.rename(columns={id_col: IID_COLUMN})
    return out


def select_traits(phenotype: pd.DataFrame, traits: str | None = None) -> list[str]:
    """Trait columns to carry through; ``None`` keeps every non-ID column."""

    available = [col for col in phenotype.columns if col not in {IID_COLUMN, FID_COLUMN}]
    if traits is None or not str(traits).strip():
        return available
    wanted = [part.strip() for part in str(traits).split(",") if part.strip()]
    missing = [name for name in wanted if name not in available]
    if missing:
        raise ValueError(f"Phenotype table missing trait(s): {', '.join(missing)}")
    return wanted


def align_phenotype_to_keep(
    phenotype: pd.DataFrame,
    keep: pd.DataFrame,
    traits: list[str],
    *,
    aggregate_duplicates: bool = False,
) -> pd.DataFrame:
    """Left-join traits onto the keep list, preserving keep order."""

    work = phenotype[[IID_COLUMN] + traits].copy()
    if aggregate_duplicates:
        work = work.groupby(IID_COLUMN, as_index=False, sort=False).mean(numeric_only=True)
    else:
        work = work.drop_duplicates(subset=[IID_COLUMN], keep="first")
    merged = keep.merge(work, on=IID_COLUMN, how="left")
    for trait in traits:
        merged[trait] = pd.to_numeric(merged[trait], errors="coerce")
    return merged


def summarize_phenotype_matrix(matrix: pd.DataFrame, traits: list[str]) -> pd.DataFrame:
    n_samples = int(len(matrix))
    rows: list[dict[str, object]] = []
    for trait in traits:
        values = pd.to_numeric(matrix[trait], errors="coerce")
        observed = int(values.notna().sum())
        rows.append(
            {
                "trait": trait,
                "n_samples": n_samples,
                "n_observed": observed,
                "missing_rate": (n_samples - observed) / n_samples if n_samples else pd.NA,
                "mean": float(values.mean()) if observed else pd.NA,
                "sd": float(values.std()) if observed > 1 else pd.NA,
                "min": float(values.min()) if observed else pd.NA,
                "max": float(values.max()) if observed else pd.NA,
            }
        )
    return pd.DataFrame(rows)


def plot_phenotype_coverage(summary: pd.DataFrame, output_prefix: str | Path) -> None:
    if summary.empty:
        return
    names = summary["trait"].astype(str).tolist()
    values = [float(v) for v in summary["n_observed"].tolist()]
    plot_bar_chart(
        names=names,
        values=values,
        title="WAP phenotype coverage by trait",
        ylabel="Samples with a value",
        filename=f"{output_prefix}.pheno.coverage.png",
        ylim=(0.0, (max(values) * 1.15) if values and max(values) > 0 else 1.05),
        rotate_xlabels=45,
    )


def run_wap_phenotype_matrix(
    phenotype_file: str | Path,
    keep_file: str | Path,
    output_prefix: str | Path,
    *,
    id_col: str = DEFAULT_ID_COLUMN,
    traits: str | None = None,
    aggregate_duplicates: bool = False,
) -> pd.DataFrame:
    """Write a PLINK2 pheno file, a per-trait QC table, thresholds, and a plot."""

    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    keep = load_keep_samples(keep_file)
    phenotype = load_phenotype_table(phenotype_file, id_col=id_col)
    trait_names = select_traits(phenotype, traits)
    if not trait_names:
        raise ValueError("No trait columns left after selection")

    matrix = align_phenotype_to_keep(
        phenotype, keep, trait_names, aggregate_duplicates=aggregate_duplicates
    )
    # PLINK2 splits on whitespace, so blanks must be written as its NA code.
    save_df_to_tsv(matrix, f"{prefix}.pheno.tsv", na_rep=PLINK_MISSING)
    summary = summarize_phenotype_matrix(matrix, trait_names)
    save_df_to_tsv(summary, f"{prefix}.pheno.qc.tsv")
    stats: dict[str, int | float] = {
        "n_keep_samples": int(len(keep)),
        "n_traits": len(trait_names),
    }
    for row in summary.itertuples(index=False):
        stats[f"n_observed_{row.trait}"] = int(row.n_observed)
    save_thresholds(stats, f"{prefix}.pheno.th.tsv")
    plot_phenotype_coverage(summary, prefix)
    print(
        f"phenotype_wap: keep={stats['n_keep_samples']} traits={stats['n_traits']} "
        f"({', '.join(trait_names)})"
    )
    return matrix
