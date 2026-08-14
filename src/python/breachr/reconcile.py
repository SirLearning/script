"""Reconcile mosdepth aneuploidy flags with breachr MAPQ-depth segment calls."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from breachr.constants import SEGMENT_OUTPUT_COLUMNS

_ABNORMAL_ENTRY_RE = re.compile(r"(\w+):([\w+]+)\(([0-9.]+)\)")

MOSDEPTH_GAIN_LABELS = frozenset({"trisomic", "tetrasomic+"})
MOSDEPTH_LOSS_LABELS = frozenset({"monosomic", "null"})

RECONCILE_OUTPUT_COLUMNS = (
    "sample",
    "logical_chr",
    "mosdepth_label",
    "mosdepth_rel_depth",
    "reconciliation",
    "breachr_call_type",
    "breachr_length_bp",
    "breachr_n_bins",
    "breachr_median_depth_rel",
    "breachr_confidence",
    "breachr_evidence",
    "notes",
)


@dataclass(frozen=True)
class ReconcileConfig:
    """Thresholds for mosdepth vs breachr comparison."""

    min_cnv_length_bp: int = 2_000_000
    min_cnv_bins: int = 2
    min_cnv_confidence: float = 0.05


def parse_abnormal_chr_list(text: str | float | None) -> list[dict[str, object]]:
    """Parse ``abnormal_chr_list`` field from mosdepth flagged table."""
    if text is None or (isinstance(text, float) and pd.isna(text)):
        return []
    entries: list[dict[str, object]] = []
    for match in _ABNORMAL_ENTRY_RE.finditer(str(text)):
        entries.append(
            {
                "logical_chr": match.group(1),
                "mosdepth_label": match.group(2),
                "mosdepth_rel_depth": float(match.group(3)),
            },
        )
    return entries


def _large_focal_cnv(chr_segments: pd.DataFrame, config: ReconcileConfig) -> pd.DataFrame:
    focal = chr_segments[chr_segments["call_type"].isin(["deletion", "duplication"])].copy()
    if focal.empty:
        return focal
    return focal[
        (focal["length_bp"] >= config.min_cnv_length_bp)
        & (focal["n_bins"] >= config.min_cnv_bins)
        & (focal["confidence"] >= config.min_cnv_confidence)
    ]


def reconcile_flagged_chr(
    sample: str,
    logical_chr: str,
    mosdepth_label: str,
    mosdepth_rel_depth: float,
    chr_segments: pd.DataFrame,
    config: ReconcileConfig | None = None,
) -> dict[str, object]:
    """Compare one mosdepth-flagged chromosome with breachr calls on that chr."""
    config = config or ReconcileConfig()
    base = {
        "sample": sample,
        "logical_chr": logical_chr,
        "mosdepth_label": mosdepth_label,
        "mosdepth_rel_depth": mosdepth_rel_depth,
        "breachr_call_type": "",
        "breachr_length_bp": 0,
        "breachr_n_bins": 0,
        "breachr_median_depth_rel": float("nan"),
        "breachr_confidence": 0.0,
        "breachr_evidence": "",
        "notes": "",
    }

    whole_chr = chr_segments[
        chr_segments["call_type"].isin(["aneuploidy_chr_gain", "aneuploidy_chr_loss"])
    ]
    subgenome_call = chr_segments[
        chr_segments["call_type"].isin(["subgenome_reference_mismatch", "subgenome_absent"])
    ]
    if not subgenome_call.empty:
        row = subgenome_call.iloc[0]
        return {
            **base,
            "reconciliation": "reference_ploidy_mismatch_not_aneuploidy",
            "breachr_call_type": str(row["call_type"]),
            "breachr_length_bp": int(row["length_bp"]),
            "breachr_n_bins": int(row["n_bins"]),
            "breachr_median_depth_rel": float(row["median_depth_rel"]),
            "breachr_confidence": float(row["confidence"]),
            "breachr_evidence": str(row["evidence"]),
            "notes": "subgenome absent or reference ploidy mismatch",
        }
    if not whole_chr.empty:
        row = whole_chr.iloc[0]
        reconciliation = "consistent_whole_chr_aneuploidy"
        notes = "breachr whole-chromosome aneuploidy call"
        if mosdepth_label in MOSDEPTH_GAIN_LABELS and row["call_type"] == "aneuploidy_chr_loss":
            notes = "mosdepth gain vs breachr loss — direction mismatch"
            reconciliation = "direction_mismatch"
        elif mosdepth_label in MOSDEPTH_LOSS_LABELS and row["call_type"] == "aneuploidy_chr_gain":
            notes = "mosdepth loss vs breachr gain — direction mismatch"
            reconciliation = "direction_mismatch"
        return {
            **base,
            "reconciliation": reconciliation,
            "breachr_call_type": row["call_type"],
            "breachr_length_bp": int(row["length_bp"]),
            "breachr_n_bins": int(row["n_bins"]),
            "breachr_median_depth_rel": float(row["median_depth_rel"]),
            "breachr_confidence": float(row["confidence"]),
            "breachr_evidence": str(row["evidence"]),
            "notes": notes,
        }

    focal_cnv = _large_focal_cnv(chr_segments, config)
    if not focal_cnv.empty:
        row = focal_cnv.sort_values(["length_bp", "confidence"], ascending=False).iloc[0]
        gain_flag = mosdepth_label in MOSDEPTH_GAIN_LABELS
        loss_flag = mosdepth_label in MOSDEPTH_LOSS_LABELS
        cnv_type = str(row["call_type"])

        if gain_flag and cnv_type == "duplication":
            reconciliation = "likely_cnv_not_aneuploidy"
            notes = "mosdepth trisomic/tetrasomic flag explained by focal duplication"
        elif loss_flag and cnv_type == "deletion":
            reconciliation = "likely_cnv_not_aneuploidy"
            notes = "mosdepth monosomic/null flag explained by focal deletion"
        elif cnv_type in {"subgenome_reference_mismatch", "subgenome_absent"}:
            reconciliation = "reference_ploidy_mismatch_not_aneuploidy"
            notes = f"subgenome-level {cnv_type}, not focal chr aneuploidy"
        elif gain_flag and cnv_type == "deletion":
            reconciliation = "focal_cnv_direction_mismatch"
            notes = "mosdepth gain flag but breachr finds large deletion"
        elif loss_flag and cnv_type == "duplication":
            reconciliation = "focal_cnv_direction_mismatch"
            notes = "mosdepth loss flag but breachr finds large duplication"
        else:
            reconciliation = "focal_cnv_ambiguous"
            notes = f"large focal {cnv_type} without matching mosdepth direction"

        return {
            **base,
            "reconciliation": reconciliation,
            "breachr_call_type": cnv_type,
            "breachr_length_bp": int(row["length_bp"]),
            "breachr_n_bins": int(row["n_bins"]),
            "breachr_median_depth_rel": float(row["median_depth_rel"]),
            "breachr_confidence": float(row["confidence"]),
            "breachr_evidence": str(row["evidence"]),
            "notes": notes,
        }

    ih = chr_segments[chr_segments["call_type"] == "inter_homoeologue"]
    if not ih.empty and mosdepth_label in MOSDEPTH_GAIN_LABELS:
        row = ih.sort_values("length_bp", ascending=False).iloc[0]
        return {
            **base,
            "reconciliation": "possible_inter_homoeologue_not_aneuploidy",
            "breachr_call_type": "inter_homoeologue",
            "breachr_length_bp": int(row["length_bp"]),
            "breachr_n_bins": int(row["n_bins"]),
            "breachr_median_depth_rel": float(row["median_depth_rel"]),
            "breachr_confidence": float(row["confidence"]),
            "breachr_evidence": str(row["evidence"]),
            "notes": "disomic depth with elevated low-MQ — shadow mapping not trisomy",
        }

    if mosdepth_label == "null" and logical_chr.endswith("D"):
        return {
            **base,
            "reconciliation": "expected_null_d_subgenome",
            "notes": "D-null in tetraploid/absent-D context — not treated as misclassification",
        }

    return {
        **base,
        "reconciliation": "no_breachr_call",
        "notes": "no whole-chr aneuploidy or large focal CNV on this chromosome",
    }


def load_segment_tables(segment_paths: list[Path]) -> pd.DataFrame:
    """Load and concatenate per-sample breachr segment TSVs."""
    frames: list[pd.DataFrame] = []
    for path in segment_paths:
        if not path.exists() or path.stat().st_size == 0:
            continue
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            continue
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=list(SEGMENT_OUTPUT_COLUMNS))
    return pd.concat(frames, ignore_index=True)


def reconcile_flagged_samples(
    flagged_df: pd.DataFrame,
    segments_df: pd.DataFrame,
    samples_with_mapq: set[str] | None = None,
    config: ReconcileConfig | None = None,
) -> pd.DataFrame:
    """Build per-(sample, chr) reconciliation table for aneuploidy-flagged rows."""
    config = config or ReconcileConfig()
    samples_with_mapq = samples_with_mapq or set()

    flagged = flagged_df.copy()
    if "aneuploid_flag" in flagged.columns:
        flagged = flagged[
            flagged["aneuploid_flag"].astype(str).str.lower().isin({"true", "1", "yes"})
        ]

    records: list[dict[str, object]] = []
    seg_by_sample: dict[str, pd.DataFrame] = {}
    if not segments_df.empty:
        for sample, grp in segments_df.groupby("sample", sort=False):
            seg_by_sample[str(sample)] = grp

    for row in flagged.itertuples(index=False):
        sample = str(row.sample)
        if sample not in samples_with_mapq:
            for entry in parse_abnormal_chr_list(getattr(row, "abnormal_chr_list", "")):
                records.append(
                    {
                        "sample": sample,
                        "logical_chr": entry["logical_chr"],
                        "mosdepth_label": entry["mosdepth_label"],
                        "mosdepth_rel_depth": entry["mosdepth_rel_depth"],
                        "reconciliation": "missing_mapq",
                        "breachr_call_type": "",
                        "breachr_length_bp": 0,
                        "breachr_n_bins": 0,
                        "breachr_median_depth_rel": float("nan"),
                        "breachr_confidence": 0.0,
                        "breachr_evidence": "",
                        "notes": "no published mapq_depth_bin.tsv for this sample",
                    },
                )
            continue

        sample_segments = seg_by_sample.get(sample, pd.DataFrame(columns=list(SEGMENT_OUTPUT_COLUMNS)))
        for entry in parse_abnormal_chr_list(getattr(row, "abnormal_chr_list", "")):
            logical_chr = str(entry["logical_chr"])
            chr_segments = sample_segments[sample_segments["logical_chr"] == logical_chr]
            records.append(
                reconcile_flagged_chr(
                    sample=sample,
                    logical_chr=logical_chr,
                    mosdepth_label=str(entry["mosdepth_label"]),
                    mosdepth_rel_depth=float(entry["mosdepth_rel_depth"]),
                    chr_segments=chr_segments,
                    config=config,
                ),
            )

    if not records:
        return pd.DataFrame(columns=list(RECONCILE_OUTPUT_COLUMNS))
    out = pd.DataFrame(records)
    return out[list(RECONCILE_OUTPUT_COLUMNS)]


def summarize_reconciliation(reconcile_df: pd.DataFrame) -> pd.DataFrame:
    """Count reconciliation categories overall and for likely misclassifications."""
    if reconcile_df.empty:
        return pd.DataFrame(
            columns=["scope", "category", "n_rows", "n_samples", "n_chromosomes"],
        )

    rows: list[dict[str, object]] = []

    def _add(scope: str, subset: pd.DataFrame, category: str) -> None:
        rows.append(
            {
                "scope": scope,
                "category": category,
                "n_rows": int(len(subset)),
                "n_samples": int(subset["sample"].nunique()) if not subset.empty else 0,
                "n_chromosomes": int(subset[["sample", "logical_chr"]].drop_duplicates().shape[0])
                if not subset.empty
                else 0,
            },
        )

    for category, grp in reconcile_df.groupby("reconciliation", sort=True):
        _add("all_flagged_chr", grp, str(category))

    mis = reconcile_df[
        reconcile_df["reconciliation"].isin(
            [
                "likely_cnv_not_aneuploidy",
                "possible_inter_homoeologue_not_aneuploidy",
            ],
        )
    ]
    _add("cnv_misclassification_candidates", mis, "all_candidates")
    for category, grp in mis.groupby("reconciliation", sort=True):
        _add("cnv_misclassification_candidates", grp, str(category))

    return pd.DataFrame(rows)


def write_reconcile_outputs(
    segment_paths: list[Path],
    flagged_path: Path,
    out_segments_all: Path,
    out_misclassified: Path,
    out_summary: Path,
    config: ReconcileConfig | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Merge segment calls and write reconciliation artefacts."""
    config = config or ReconcileConfig()
    flagged_df = pd.read_csv(flagged_path, sep="\t")
    segments_df = load_segment_tables(segment_paths)
    samples_with_mapq = {p.name.replace(".breachr_segments.tsv", "") for p in segment_paths}

    segments_df.to_csv(out_segments_all, sep="\t", index=False)

    reconcile_df = reconcile_flagged_samples(
        flagged_df=flagged_df,
        segments_df=segments_df,
        samples_with_mapq=samples_with_mapq,
        config=config,
    )
    reconcile_df.to_csv(out_misclassified.parent / "breachr_reconcile_all.tsv", sep="\t", index=False)

    misclassified = reconcile_df[
        reconcile_df["reconciliation"].isin(
            [
                "likely_cnv_not_aneuploidy",
                "possible_inter_homoeologue_not_aneuploidy",
            ],
        )
    ].copy()
    misclassified.to_csv(out_misclassified, sep="\t", index=False)

    summary_df = summarize_reconciliation(reconcile_df)
    summary_df.to_csv(out_summary, sep="\t", index=False)

    return segments_df, reconcile_df, summary_df


def resolve_mapq_sample_aliases(sample: str) -> list[str]:
    """Return candidate sample IDs for locating a published MAPQ TSV."""
    sample = sample.strip()
    candidates = [sample]
    if sample.endswith("_deduped"):
        candidates.append(sample[: -len("_deduped")])
    return candidates
