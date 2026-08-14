"""Main entry: call segments from a MAPQ-depth bin table."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from breachr.classify import chr_aneuploidy_call, label_bins, segment_confidence
from breachr.config import CallConfig
from breachr.constants import LOGICAL_AUTOSOMES, SEGMENT_OUTPUT_COLUMNS
from breachr.dosage import chr_dosage_profile, is_uniform_monosomic, is_uniform_trisomic
from breachr.features import (
    add_depth_relatives,
    add_mapq_fractions,
    add_smoothed_depth,
    compute_subgenome_baselines,
    sample_mq_low_frac_reference,
    subgenome_from_logical_chr,
)
from breachr.io import load_mapq_depth_bin_tsv
from breachr.partition import partition_chromosome, segments_to_records
from breachr.subgenome import SubgenomeState, detect_subgenome_states


def _whole_chr_record(
    sample: str,
    chr_df: pd.DataFrame,
    call_type: str,
    confidence: float,
    evidence: str,
) -> dict[str, object]:
    depth_med = float(pd.to_numeric(chr_df["depth_rel"], errors="coerce").median())
    mq_med = float(pd.to_numeric(chr_df["mq_low_frac"], errors="coerce").median())
    return {
        "sample": sample,
        "logical_chr": str(chr_df["logical_chr"].iloc[0]),
        "bin_start": int(chr_df["bin_start"].min()),
        "bin_end": int(chr_df["bin_end"].max()),
        "length_bp": int(chr_df["bin_end"].max() - chr_df["bin_start"].min() + 1),
        "n_bins": int(len(chr_df)),
        "median_depth_rel": depth_med,
        "median_mq_low_frac": mq_med,
        "call_type": call_type,
        "confidence": confidence,
        "evidence": evidence,
    }


def _subgenome_record(
    sample: str,
    chr_df: pd.DataFrame,
    state: SubgenomeState,
) -> dict[str, object]:
    return _whole_chr_record(
        sample=sample,
        chr_df=chr_df,
        call_type=state.state,
        confidence=0.85,
        evidence=state.evidence,
    )


def _consolidate_segment_records(
    records: list[dict[str, object]],
    chr_df: pd.DataFrame,
    sample: str,
    mq_reference: float,
    config: CallConfig,
) -> list[dict[str, object]]:
    """Upgrade fragmented deletion/dup runs to whole-chr aneuploidy when uniform."""
    if not records:
        return records
    profile = chr_dosage_profile(chr_df, config)
    if profile is None:
        return records

    del_records = [r for r in records if r["call_type"] == "deletion"]
    dup_records = [r for r in records if r["call_type"] == "duplication"]

    if len(del_records) >= config.consolidate_deletion_segments_min and is_uniform_monosomic(profile, config):
        call_type, conf, evidence = chr_aneuploidy_call(chr_df, config)
        if call_type is None:
            call_type = "aneuploidy_chr_loss"
            conf = segment_confidence(
                call_type,
                profile.median_rel,
                float(pd.to_numeric(chr_df["mq_low_frac"], errors="coerce").median()),
                mq_reference,
                config,
            )
            evidence = (
                f"consolidated_from_{len(del_records)}_deletion_segments;"
                f"uniform_monosomic;median_rel={profile.median_rel:.3f}"
            )
        return [_whole_chr_record(sample, chr_df, call_type, conf, evidence)]

    if len(dup_records) >= config.consolidate_deletion_segments_min and is_uniform_trisomic(profile, config):
        call_type, conf, evidence = chr_aneuploidy_call(chr_df, config)
        if call_type is None:
            call_type = "aneuploidy_chr_gain"
            conf = segment_confidence(
                call_type,
                profile.median_rel,
                float(pd.to_numeric(chr_df["mq_low_frac"], errors="coerce").median()),
                mq_reference,
                config,
            )
            evidence = (
                f"consolidated_from_{len(dup_records)}_duplication_segments;"
                f"uniform_trisomic;median_rel={profile.median_rel:.3f}"
            )
        return [_whole_chr_record(sample, chr_df, call_type, conf, evidence)]

    return records


def call_mapq_depth_bins(
    df: pd.DataFrame,
    sample: str,
    config: CallConfig | None = None,
) -> pd.DataFrame:
    """Call CNV / inter-homoeologue / aneuploidy segments from 1 Mb MAPQ-depth bins."""
    config = config or CallConfig()
    work = add_mapq_fractions(df)
    baselines = compute_subgenome_baselines(work, config)
    work = add_depth_relatives(work, baselines)
    work = add_smoothed_depth(work, config)
    mq_reference = sample_mq_low_frac_reference(work, config)
    work = label_bins(work, mq_reference, config)
    subgenome_states = detect_subgenome_states(work, config)

    records: list[dict[str, object]] = []
    for logical_chr in LOGICAL_AUTOSOMES:
        chr_df = work[work["logical_chr"] == logical_chr]
        if chr_df.empty:
            continue

        sg = subgenome_from_logical_chr(logical_chr)
        if sg in subgenome_states:
            records.append(_subgenome_record(sample, chr_df, subgenome_states[sg]))
            continue

        call_type, conf, evidence = chr_aneuploidy_call(chr_df, config)
        if call_type is not None:
            records.append(_whole_chr_record(sample, chr_df, call_type, conf, evidence))
            continue

        chr_records: list[dict[str, object]] = []
        for seg in partition_chromosome(chr_df, config):
            chr_records.extend(
                segments_to_records([seg], chr_df, sample, mq_reference, config),
            )
        records.extend(_consolidate_segment_records(chr_records, chr_df, sample, mq_reference, config))

    if not records:
        return pd.DataFrame(columns=list(SEGMENT_OUTPUT_COLUMNS))
    out = pd.DataFrame(records)
    return out[list(SEGMENT_OUTPUT_COLUMNS)].sort_values(
        ["logical_chr", "bin_start"],
    ).reset_index(drop=True)


def call_mapq_depth_bin_tsv(
    tsv_path: str | Path,
    sample: str | None = None,
    config: CallConfig | None = None,
) -> pd.DataFrame:
    """Load a published TSV and return segment calls."""
    path = Path(tsv_path)
    sample = sample or path.name.replace(".mapq_depth_bin.tsv", "")
    df = load_mapq_depth_bin_tsv(path)
    return call_mapq_depth_bins(df, sample=sample, config=config)
