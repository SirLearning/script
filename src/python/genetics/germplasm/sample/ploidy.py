"""Genome-formula ploidy: previous sample_group records, corrected by breachr.

``sample_group`` is the ploidy assumed at alignment / mosdepth comparison
(AA, SS, DD, AABB, AABBDD). Full-cohort breachr segment calls then drop a
claimed subgenome only when every homoeologue is ``subgenome_absent`` or
``subgenome_reference_mismatch``. Flags on subgenomes the previous label
never claimed (AB D-null, A-genome B/D-null) are expected, not a correction.
Focal CNV and chromosome aneuploidy do not change the formula.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from genetics.germplasm.sample.cohort_manifest import load_sample_group_table
from genetics.germplasm.sample.mosdepth_aneuploidy import normalize_sample_id
from infra.utils.io import load_df_from_tsv, save_df_to_tsv

GROUP_TO_PLOIDY = {
    "A": "AA",
    "S": "SS",
    "D": "DD",
    "AB": "AABB",
    "ABD": "AABBDD",
    "WAP": "AABBDD",
    "Watkins": "AABBDD",
    "Nature": "AABBDD",
    "HZNU": "AABBDD",
}

PLOIDY_SUBGENOMES = {
    "AA": frozenset("A"),
    "SS": frozenset("B"),
    "DD": frozenset("D"),
    "AABB": frozenset("AB"),
    "AABBDD": frozenset("ABD"),
}

PLOIDY_FROM_SUBGENOMES = {sgs: ploidy for ploidy, sgs in PLOIDY_SUBGENOMES.items()}

SUBGENOME_PLOIDIES = {
    "A": ("AA", "AABB", "AABBDD"),
    "B": ("SS", "AABB", "AABBDD"),
    "D": ("DD", "AABBDD"),
}

PLOIDY_ORDER = ("AA", "SS", "DD", "AABB", "AABBDD")

PLOIDY_PALETTE = {
    "AA": "#4c78a8",
    "SS": "#f58518",
    "DD": "#54a24b",
    "AABB": "#e45756",
    "AABBDD": "#72b7b2",
}

SUBGENOMES = ("A", "B", "D")

PLOIDY_CALL_TYPES = frozenset({"subgenome_absent", "subgenome_reference_mismatch"})
MIN_CHR_FOR_SUBGENOME_FLAG = 7


def nominal_ploidy_from_group(sample_group: str) -> str:
    """Map a sample_group label to the recorded genome formula."""
    key = str(sample_group).strip()
    if key not in GROUP_TO_PLOIDY:
        raise ValueError(f"Unknown sample_group for ploidy: {sample_group!r}")
    return GROUP_TO_PLOIDY[key]


def formula_from_subgenomes(present: frozenset[str]) -> str | None:
    """Map an observed present-subgenome set to AA/SS/DD/AABB/AABBDD, or None."""
    return PLOIDY_FROM_SUBGENOMES.get(frozenset(present))


def subgenome_ploidies(subgenome: str) -> tuple[str, ...]:
    """Ploidy classes that carry this subgenome in the panel."""
    key = str(subgenome).strip().upper()
    if key not in SUBGENOME_PLOIDIES:
        raise ValueError(f"Unknown subgenome: {subgenome!r}")
    return SUBGENOME_PLOIDIES[key]


def _join_key(sample: str) -> str:
    return normalize_sample_id(str(sample))


def load_breachr_segments(path: str | Path) -> pd.DataFrame:
    """Load a concatenated breachr segment table (full cohort or a test slice)."""
    table = load_df_from_tsv(str(path))
    if table is None or table.empty:
        raise ValueError(f"Could not load breachr segments: {path}")
    return table


def subgenome_flags_from_breachr(
    segments: pd.DataFrame,
    *,
    min_chr: int = MIN_CHR_FOR_SUBGENOME_FLAG,
) -> pd.DataFrame:
    """One row per sample with A/B/D call_type when all homoeologues share a ploidy call."""
    empty = pd.DataFrame(columns=["Sample", "A", "B", "D"])
    if segments is None or segments.empty:
        return empty
    if "call_type" not in segments.columns or "logical_chr" not in segments.columns:
        return empty
    sample_col = "sample" if "sample" in segments.columns else "Sample"
    work = segments[segments["call_type"].isin(PLOIDY_CALL_TYPES)].copy()
    if work.empty:
        return empty
    work["Sample"] = work[sample_col].astype(str)
    work["sg"] = work["logical_chr"].astype(str).str[-1]
    work = work[work["sg"].isin(SUBGENOMES)]
    counts = (
        work.groupby(["Sample", "sg", "call_type"], sort=False)
        .size()
        .reset_index(name="n_chr")
    )
    counts = counts[counts["n_chr"] >= int(min_chr)]
    if counts.empty:
        return empty
    flags = (
        counts.groupby(["Sample", "sg"], sort=False)["call_type"]
        .agg(lambda s: str(s.iloc[0]))
        .unstack("sg")
    )
    for sg in SUBGENOMES:
        if sg not in flags.columns:
            flags[sg] = pd.NA
    return flags.reset_index()[["Sample", "A", "B", "D"]]


def _decide_ploidy(nominal: str, drop: dict[str, str]) -> tuple[str, str, frozenset[str]]:
    nominal_sgs = PLOIDY_SUBGENOMES[nominal]
    drop_sgs = frozenset(drop)
    claimed_drop = nominal_sgs & drop_sgs
    present = nominal_sgs - claimed_drop
    if not claimed_drop:
        return nominal, "agrees", present
    observed = formula_from_subgenomes(present)
    if observed is None:
        return nominal, "breachr_unmapped", present
    return observed, "breachr_ploidy_mismatch", present


def assign_corrected_ploidy(
    groups: pd.DataFrame,
    breachr_segments: pd.DataFrame,
) -> pd.DataFrame:
    """Left-join sample_group records to breachr subgenome flags and correct formula."""
    if groups is None or groups.empty:
        raise ValueError("sample_group table is empty")
    work = groups.copy()
    if "Sample" not in work.columns or "Group" not in work.columns:
        raise ValueError("sample_group table needs Sample and Group columns")
    work["Sample"] = work["Sample"].astype(str)
    work["Group"] = work["Group"].astype(str)
    work["Sample_norm"] = work["Sample"].map(_join_key)
    work["nominal_ploidy"] = work["Group"].map(nominal_ploidy_from_group)

    flags = subgenome_flags_from_breachr(breachr_segments)
    if flags.empty:
        flags = pd.DataFrame(columns=["Sample", "A", "B", "D"])
    flags = flags.copy()
    flags["Sample_norm"] = flags["Sample"].map(_join_key)
    flags = flags.rename(columns={"Sample": "sample_breachr"})
    flags = flags.drop_duplicates(subset=["Sample_norm"], keep="first")
    work = work.merge(flags, on="Sample_norm", how="left").reset_index(drop=True)

    decisions = []
    records = work.to_dict(orient="records")
    for rec in records:
        drop: dict[str, str] = {}
        for sg in SUBGENOMES:
            value = rec.get(sg)
            if value is None or (isinstance(value, float) and pd.isna(value)):
                continue
            if pd.isna(value):
                continue
            text = str(value).strip()
            if text and text.lower() != "nan":
                drop[sg] = text
        ploidy, reason, present = _decide_ploidy(str(rec["nominal_ploidy"]), drop)
        decisions.append(
            {
                "ploidy": ploidy,
                "reason": reason,
                "present_subgenomes": "".join(sg for sg in SUBGENOMES if sg in present),
                "breachr_drop": "".join(sg for sg in SUBGENOMES if sg in drop),
                "breachr_A": drop.get("A", ""),
                "breachr_B": drop.get("B", ""),
                "breachr_D": drop.get("D", ""),
            }
        )
    decided = pd.DataFrame(decisions)
    out = pd.concat([work.reset_index(drop=True), decided], axis=1)
    keep = [
        "Sample",
        "Group",
        "nominal_ploidy",
        "ploidy",
        "reason",
        "present_subgenomes",
        "breachr_drop",
        "breachr_A",
        "breachr_B",
        "breachr_D",
        "sample_breachr",
    ]
    for col in keep:
        if col not in out.columns:
            out[col] = pd.NA
    return out[keep]


def build_corrected_ploidy_table(
    group_file: str | Path,
    breachr_segments_path: str | Path,
    output_path: str | Path | None = None,
) -> pd.DataFrame:
    """Load sample_group + breachr segments and optionally write the ploidy TSV."""
    groups = load_sample_group_table(group_file)
    segments = load_breachr_segments(breachr_segments_path)
    table = assign_corrected_ploidy(groups, segments)
    if output_path is not None:
        save_df_to_tsv(table, str(output_path))
    return table
