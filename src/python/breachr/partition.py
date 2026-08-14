"""Segmentation along chromosomes (CNVpytor-style partition on 1 Mb bins)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from breachr.config import CallConfig


@dataclass(frozen=True)
class RawSegment:
    """Contiguous run of bins sharing a non-normal preliminary label."""

    logical_chr: str
    start_bin_idx: int
    end_bin_idx: int
    bin_label: str


def _merge_same_label_runs(chr_df: pd.DataFrame, config: CallConfig) -> list[RawSegment]:
    """Merge consecutive bins with identical ``bin_label`` (excluding normal)."""
    segments: list[RawSegment] = []
    if chr_df.empty:
        return segments

    logical_chr = str(chr_df["logical_chr"].iloc[0])
    current_label: str | None = None
    start_idx: int | None = None
    prev_bin_idx: int | None = None

    for row in chr_df.itertuples(index=False):
        label = str(row.bin_label)
        bin_idx = int(row.bin_idx)
        if label in {"normal", "low_support"}:
            if current_label is not None and start_idx is not None and prev_bin_idx is not None:
                segments.append(
                    RawSegment(logical_chr, start_idx, prev_bin_idx, current_label),
                )
            current_label = None
            start_idx = None
            prev_bin_idx = None
            continue

        if label == current_label:
            prev_bin_idx = bin_idx
            continue

        if current_label is not None and start_idx is not None and prev_bin_idx is not None:
            segments.append(
                RawSegment(logical_chr, start_idx, prev_bin_idx, current_label),
            )
        current_label = label
        start_idx = bin_idx
        prev_bin_idx = bin_idx

    if current_label is not None and start_idx is not None and prev_bin_idx is not None:
        segments.append(RawSegment(logical_chr, start_idx, prev_bin_idx, current_label))

    return [s for s in segments if (s.end_bin_idx - s.start_bin_idx + 1) >= config.min_segment_bins]


def merge_adjacent_segments(segments: list[RawSegment], *, max_gap_bins: int = 0) -> list[RawSegment]:
    """Merge consecutive segments with the same label on one chromosome."""
    if not segments:
        return []
    ordered = sorted(segments, key=lambda s: s.start_bin_idx)
    merged: list[RawSegment] = [ordered[0]]
    for seg in ordered[1:]:
        prev = merged[-1]
        gap = seg.start_bin_idx - prev.end_bin_idx - 1
        if seg.bin_label == prev.bin_label and gap <= max_gap_bins:
            merged[-1] = RawSegment(prev.logical_chr, prev.start_bin_idx, seg.end_bin_idx, prev.bin_label)
        else:
            merged.append(seg)
    return merged


def partition_chromosome(chr_df: pd.DataFrame, config: CallConfig) -> list[RawSegment]:
    """Partition one chromosome into MAPQ-depth segments.

    Uses per-bin labels from joint classification, then optionally splits long
    runs where smoothed depth jumps (mean-shift inspired breakpoint scan).
    """
    labeled = chr_df.sort_values("bin_idx").reset_index(drop=True)
    from breachr.dosage import chr_dosage_profile, is_partial_cnv_pattern

    profile = chr_dosage_profile(labeled, config)
    partial_cnv = profile is not None and is_partial_cnv_pattern(profile, config)
    base_segments = merge_adjacent_segments(
        _merge_same_label_runs(labeled, config),
        max_gap_bins=2 if partial_cnv else 0,
    )
    if not base_segments or config.partition_bandwidth_bins <= 1:
        return base_segments

    if partial_cnv:
        return base_segments

    # Refine only within homogeneous label runs that span many bins.
    refined: list[RawSegment] = []
    for seg in base_segments:
        sub = labeled[
            (labeled["bin_idx"] >= seg.start_bin_idx)
            & (labeled["bin_idx"] <= seg.end_bin_idx)
        ]
        if len(sub) < config.min_segment_bins * 2:
            refined.append(seg)
            continue

        depth = pd.to_numeric(sub["depth_rel_smooth"], errors="coerce").to_numpy(dtype=float)
        split_at = _depth_breakpoints(depth, config.partition_bandwidth_bins)
        if not split_at:
            refined.append(seg)
            continue

        bin_indices = sub["bin_idx"].astype(int).to_numpy()
        edges = [0] + split_at + [len(sub)]
        for i in range(len(edges) - 1):
            lo = edges[i]
            hi = edges[i + 1] - 1
            if hi - lo + 1 < config.min_segment_bins:
                continue
            refined.append(
                RawSegment(
                    seg.logical_chr,
                    int(bin_indices[lo]),
                    int(bin_indices[hi]),
                    seg.bin_label,
                ),
            )
    return merge_adjacent_segments(refined)


def _depth_breakpoints(depth: np.ndarray, bandwidth: int) -> list[int]:
    """Find internal indices where smoothed depth shifts beyond local bandwidth."""
    if len(depth) < 2 * bandwidth + 1:
        return []
    breakpoints: list[int] = []
    for i in range(bandwidth, len(depth) - bandwidth):
        left = depth[i - bandwidth : i]
        right = depth[i : i + bandwidth]
        if not np.isfinite(left).any() or not np.isfinite(right).any():
            continue
        left_med = float(np.nanmedian(left))
        right_med = float(np.nanmedian(right))
        if abs(left_med - right_med) >= 0.35:
            breakpoints.append(i)
    return breakpoints


def segments_to_records(
    segments: list[RawSegment],
    chr_df: pd.DataFrame,
    sample: str,
    mq_reference: float,
    config: CallConfig,
) -> list[dict[str, object]]:
    """Summarise raw segments into call records."""
    from breachr.classify import segment_confidence

    records: list[dict[str, object]] = []
    for seg in segments:
        sub = chr_df[
            (chr_df["bin_idx"] >= seg.start_bin_idx)
            & (chr_df["bin_idx"] <= seg.end_bin_idx)
        ]
        if sub.empty:
            continue
        bin_start = int(sub["bin_start"].min())
        bin_end = int(sub["bin_end"].max())
        median_depth = float(pd.to_numeric(sub["depth_rel"], errors="coerce").median())
        median_mq = float(pd.to_numeric(sub["mq_low_frac"], errors="coerce").median())
        evidence = ";".join(sorted(set(sub["bin_evidence"].astype(str))))
        call_type = seg.bin_label
        conf = segment_confidence(
            call_type,
            median_depth,
            median_mq,
            mq_reference,
            config,
        )
        records.append(
            {
                "sample": sample,
                "logical_chr": seg.logical_chr,
                "bin_start": bin_start,
                "bin_end": bin_end,
                "length_bp": bin_end - bin_start + 1,
                "n_bins": len(sub),
                "median_depth_rel": median_depth,
                "median_mq_low_frac": median_mq,
                "call_type": call_type,
                "confidence": conf,
                "evidence": evidence,
            },
        )
    return records
