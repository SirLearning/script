"""Shared constants for breachr MAPQ-depth segment calling."""

from __future__ import annotations

MQ_BIN_LABELS = (
    "mq_0_10",
    "mq_10_20",
    "mq_20_30",
    "mq_30_40",
    "mq_40_50",
    "mq_50_60",
)

MQ_LOW_LABELS = ("mq_0_10", "mq_10_20")
MQ_MID_LABELS = ("mq_20_30", "mq_30_40")
MQ_HIGH_LABELS = ("mq_40_50", "mq_50_60")

MAPQ_DEPTH_BIN_COLUMNS = (
    "logical_chr",
    "bin_idx",
    "bin_start",
    "bin_end",
    *MQ_BIN_LABELS,
    "total_reads",
    "bin_mean_depth",
)

# Wheat logical autosomes chr1A–chr7D (same order as mosdepth aneuploidy).
LOGICAL_AUTOSOMES = tuple(
    f"chr{i}{genome}" for i in range(1, 8) for genome in ("A", "B", "D")
)

SUBGENOME_AUTOSOMES = {
    sg: tuple(f"chr{i}{sg}" for i in range(1, 8)) for sg in ("A", "B", "D")
}

CALL_TYPES = (
    "normal",
    "deletion",
    "duplication",
    "inter_homoeologue",
    "aneuploidy_chr_loss",
    "aneuploidy_chr_gain",
    "subgenome_absent",
    "subgenome_reference_mismatch",
    "low_support",
)

SEGMENT_OUTPUT_COLUMNS = (
    "sample",
    "logical_chr",
    "bin_start",
    "bin_end",
    "length_bp",
    "n_bins",
    "median_depth_rel",
    "median_mq_low_frac",
    "call_type",
    "confidence",
    "evidence",
)
