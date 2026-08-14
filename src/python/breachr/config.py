"""Thresholds and tunables for MAPQ-depth joint calling."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CallConfig:
    """MAPQ- and depth-aware calling thresholds.

    Depth bounds align loosely with mosdepth ``CN_CLASS_BOUNDS`` but operate on
    per-1Mb bins. MAPQ gates encode the idea that misplaced reads (deletions,
    inter-homeologous shadow, some aneuploidy bins) carry excess low-MQ mass.
    """

    # Relative depth (bin_mean_depth / subgenome baseline).
    deletion_depth_rel_max: float = 0.55
    duplication_depth_rel_min: float = 1.40
    disomic_depth_min: float = 0.75
    disomic_depth_max: float = 1.25
    aneuploidy_chr_loss_median_max: float = 0.60
    aneuploidy_chr_gain_median_min: float = 1.40

    # MAPQ: low MAPQ fraction = (mq_0_10 + mq_10_20) / total_reads.
    mq_low_frac_delta: float = 0.12
    mq_low_frac_absolute_min: float = 0.30
    duplication_mq_low_frac_max: float = 0.28

    min_reads_per_bin: int = 50
    min_segment_bins: int = 2
    smooth_window_bins: int = 3
    chr_aneuploidy_bin_fraction: float = 0.65
    partition_bandwidth_bins: int = 3

    # Uniform whole-chromosome aneuploidy vs focal CNV.
    min_chr_bins_for_aneuploidy: int = 20
    chr_disomic_frac_max_for_aneuploidy: float = 0.12
    chr_partial_cnv_disomic_frac_min: float = 0.15
    chr_partial_cnv_altered_frac_min: float = 0.20
    chr_partial_cnv_disomic_run_frac_min: float = 0.08
    uniform_monosomic_median_min: float = 0.42
    uniform_monosomic_median_max: float = 0.68
    uniform_monosomic_loss_frac_min: float = 0.35
    uniform_trisomic_median_min: float = 1.32
    uniform_trisomic_median_max: float = 1.80
    uniform_trisomic_gain_frac_min: float = 0.35
    uniform_aneuploidy_iqr_max: float = 0.20
    uniform_aneuploidy_cv_max: float = 0.35
    consolidate_deletion_segments_min: int = 4

    # Subgenome absent / reference ploidy mismatch (e.g. AB in ABD ref).
    subgenome_depth_ratio_mismatch_max: float = 0.65
    subgenome_mq_low_mismatch_min: float = 0.28
    subgenome_absent_reads_median_max: int = 10

    # Confidence scaling.
    depth_confidence_scale: float = 0.35
    mq_confidence_scale: float = 0.25
