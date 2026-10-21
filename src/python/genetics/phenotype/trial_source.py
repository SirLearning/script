"""Field-trial plot tables into project phenotype sources.

Input is the DBone export shape, one row per field plot::

    sampleName,trialYear,trait,PH_1,PH_2,...,PH_10

``sampleName`` is ``Z<yy>-<plot_id>`` for an association accession and
``Z<yy>-<plot_id>-<check_name>`` for a repeated check cultivar.

Two facts drive the split:

* A plot slot used by a check in any year is a check slot in every year. Only
  the 2019 rows carry the ``-Jimai22`` label; 2018 and 2020 record the same
  slots unlabelled. Those slots have no WAP genotype, so keeping them would
  inject check-cultivar measurements into the association panel.
* Genotype IDs come in two forms, ``P0001`` and ``P0030A`` (registry
  ``domain/sample_cohorts.yaml``, ``id_prefix: "P / PA"``). The numeric parts
  do not collide, so stripping a trailing ``A`` gives an unambiguous join key.

Replicate values are averaged per accession and year. That is a plain mean, not
a BLUE/BLUP: no block, spatial, or year-effect model is fitted here.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from infra.utils.graph import plot_bar_chart
from infra.utils.io import load_df_generic, save_df_to_tsv, save_thresholds

ID_COLUMNS = ("sampleName", "trialYear", "trait")
TRIAL_NAME_RE = re.compile(r"^Z(?P<yy>\d{2})-(?P<plot>P\d+)(?:-(?P<check>.+))?$")
GENOTYPE_SUFFIX_RE = re.compile(r"A$")
MEAN_TRAIT_SUFFIX = "mean"


def parse_trial_sample_name(sample_name: str) -> tuple[str, str, str]:
    """Split ``Z19-P0025-Jimai22`` into (year code, plot id, check name)."""

    match = TRIAL_NAME_RE.match(str(sample_name).strip())
    if match is None:
        raise ValueError(f"Unparseable trial sample name: {sample_name}")
    return match.group("yy"), match.group("plot"), match.group("check") or ""


def normalize_genotype_id(sample_id: str) -> str:
    """Join key shared by phenotype plot ids and genotype IIDs (``P0030A`` -> ``P0030``)."""

    return GENOTYPE_SUFFIX_RE.sub("", str(sample_id).strip())


def load_trial_plot_csv(csv_file: str | Path) -> pd.DataFrame:
    """Melt a wide plot table into one row per (plot, replicate) with a value."""

    raw = load_df_generic(str(csv_file))
    if raw is None or raw.empty:
        raise ValueError(f"Empty trial table: {csv_file}")
    missing = [col for col in ID_COLUMNS if col not in raw.columns]
    if missing:
        raise ValueError(f"Trial table missing column(s) {', '.join(missing)}: {csv_file}")
    rep_columns = [col for col in raw.columns if col not in ID_COLUMNS]
    if not rep_columns:
        raise ValueError(f"Trial table has no replicate columns: {csv_file}")

    tidy = raw.melt(
        id_vars=list(ID_COLUMNS),
        value_vars=rep_columns,
        var_name="rep",
        value_name="value",
    )
    tidy["value"] = pd.to_numeric(tidy["value"], errors="coerce")
    tidy = tidy[tidy["value"].notna()].copy()

    parsed = [parse_trial_sample_name(name) for name in tidy["sampleName"]]
    tidy["year_code"] = [item[0] for item in parsed]
    tidy["plot_id"] = [item[1] for item in parsed]
    tidy["check_name"] = [item[2] for item in parsed]
    tidy["trial_year"] = pd.to_numeric(tidy["trialYear"], errors="coerce").astype("Int64")
    return tidy.reset_index(drop=True)


def check_plot_slots(tidy: pd.DataFrame) -> set[str]:
    """Plot slots labelled with a check cultivar in any year."""

    labelled = tidy.loc[tidy["check_name"].astype(str) != "", "plot_id"]
    return set(labelled.astype(str))


def split_check_rows(tidy: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Separate association-accession rows from check-cultivar rows.

    Unlabelled rows sitting in a labelled slot inherit that slot's check name,
    so the check source keeps every year of the same cultivar.
    """

    labelled = tidy.loc[tidy["check_name"].astype(str) != ""]
    slot_names = (
        labelled.drop_duplicates(subset=["plot_id"]).set_index("plot_id")["check_name"].to_dict()
    )
    is_check = tidy["plot_id"].astype(str).isin(slot_names)
    checks = tidy[is_check].copy()
    checks["check_name"] = checks["plot_id"].map(slot_names)
    return tidy[~is_check].copy(), checks


def build_genotype_id_map(genotype_ids: list[str]) -> dict[str, str]:
    """Map a normalised plot id to the genotype IID that carries it."""

    mapping: dict[str, str] = {}
    for genotype_id in genotype_ids:
        key = normalize_genotype_id(genotype_id)
        existing = mapping.get(key)
        if existing is not None and existing != genotype_id:
            raise ValueError(f"Ambiguous genotype ids for plot {key}: {existing} and {genotype_id}")
        mapping[key] = str(genotype_id)
    return mapping


def load_genotype_ids(keep_file: str | Path) -> list[str]:
    """Read genotype IIDs from a PLINK2 keep file (FID IID) or a one-column list."""

    raw = pd.read_csv(keep_file, sep=r"\s+", header=None, dtype=str)
    if raw.empty:
        raise ValueError(f"Empty keep file: {keep_file}")
    column = 0 if raw.shape[1] == 1 else 1
    return raw.iloc[:, column].astype(str).tolist()


def attach_genotype_ids(accessions: pd.DataFrame, id_map: dict[str, str]) -> pd.DataFrame:
    """Add the genotype ``IID`` column; rows without a genotype get NA."""

    out = accessions.copy()
    out["IID"] = [id_map.get(normalize_genotype_id(plot)) for plot in out["plot_id"]]
    return out


def pivot_plot_means(matched: pd.DataFrame) -> pd.DataFrame:
    """Per-accession replicate means: one column per trial year plus an overall mean."""

    if matched.empty:
        raise ValueError("No genotyped accession rows to pivot")
    per_year = matched.groupby(["IID", "trait", "trial_year"], as_index=False)["value"].mean()
    per_year["column"] = [
        f"{row.trait}_{int(row.trial_year)}" for row in per_year.itertuples(index=False)
    ]
    wide = per_year.pivot(index="IID", columns="column", values="value")

    overall = matched.groupby(["IID", "trait"], as_index=False)["value"].mean()
    overall["column"] = [f"{trait}_{MEAN_TRAIT_SUFFIX}" for trait in overall["trait"]]
    wide = wide.join(overall.pivot(index="IID", columns="column", values="value"))
    wide.columns.name = None
    return wide.reset_index().sort_values("IID").reset_index(drop=True)


def summarize_trial_source(
    tidy: pd.DataFrame,
    accessions: pd.DataFrame,
    checks: pd.DataFrame,
    matched: pd.DataFrame,
) -> dict[str, int]:
    return {
        "n_plot_values": int(len(tidy)),
        "n_accession_values": int(len(accessions)),
        "n_check_values": int(len(checks)),
        "n_check_slots": len(check_plot_slots(tidy)),
        "n_accession_plots": int(accessions["plot_id"].nunique()),
        "n_genotyped_accessions": int(matched["IID"].nunique()),
        "n_ungenotyped_plots": int(
            accessions.loc[~accessions["plot_id"].isin(matched["plot_id"]), "plot_id"].nunique()
        ),
        "n_trial_years": int(tidy["trial_year"].nunique()),
    }


def plot_trial_coverage(wide: pd.DataFrame, output_prefix: str | Path) -> None:
    traits = [col for col in wide.columns if col != "IID"]
    values = [float(wide[col].notna().sum()) for col in traits]
    plot_bar_chart(
        names=traits,
        values=values,
        title="Genotyped accessions per trial year",
        ylabel="Accessions with a value",
        filename=f"{output_prefix}.trial.coverage.png",
        ylim=(0.0, (max(values) * 1.15) if values and max(values) > 0 else 1.05),
        rotate_xlabels=30,
    )


def run_trial_source_extract(
    csv_file: str | Path,
    keep_file: str | Path,
    output_prefix: str | Path,
    *,
    check_prefix: str | Path | None = None,
) -> dict[str, int]:
    """Write the accession source, the check-cultivar source, and a wide trait table."""

    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    tidy = load_trial_plot_csv(csv_file)
    accessions, checks = split_check_rows(tidy)
    id_map = build_genotype_id_map(load_genotype_ids(keep_file))
    accessions = attach_genotype_ids(accessions, id_map)
    matched = accessions[accessions["IID"].notna()].copy()

    columns = ["IID", "plot_id", "trial_year", "trait", "rep", "value"]
    save_df_to_tsv(matched[columns], f"{prefix}.plot.tsv")
    wide = pivot_plot_means(matched)
    # PLINK2 splits on whitespace, so blanks must be written as its NA code.
    save_df_to_tsv(wide, f"{prefix}.pheno.tsv", na_rep="NA")

    check_out = Path(check_prefix) if check_prefix else prefix.with_name(f"{prefix.name}.check")
    check_out.parent.mkdir(parents=True, exist_ok=True)
    check_columns = ["plot_id", "check_name", "trial_year", "trait", "rep", "value"]
    save_df_to_tsv(checks[check_columns], f"{check_out}.plot.tsv")

    stats = summarize_trial_source(tidy, accessions, checks, matched)
    save_thresholds(stats, f"{prefix}.src.th.tsv")
    plot_trial_coverage(wide, prefix)
    print(
        f"trial_source: plots={stats['n_plot_values']} accessions={stats['n_genotyped_accessions']} "
        f"checks={stats['n_check_values']} unmatched_plots={stats['n_ungenotyped_plots']}"
    )
    return stats
