"""Per-bin and segment classification using joint MAPQ + depth evidence."""

from __future__ import annotations

import numpy as np
import pandas as pd

from breachr.config import CallConfig
from breachr.dosage import (
    ChrDosageProfile,
    chr_dosage_profile,
    is_partial_cnv_pattern,
    is_uniform_monosomic,
    is_uniform_trisomic,
)


def mq_is_elevated(mq_low_frac: float, mq_reference: float, config: CallConfig) -> bool:
    """True when low-MAPQ fraction exceeds sample reference plus margin."""
    return mq_low_frac >= max(
        mq_reference + config.mq_low_frac_delta,
        config.mq_low_frac_absolute_min,
    )


def classify_bin_row(
    depth_rel: float,
    mq_low_frac: float,
    total_reads: float,
    mq_reference: float,
    config: CallConfig,
) -> tuple[str, str]:
    """Return preliminary bin label and short evidence string."""
    if not np.isfinite(depth_rel) or total_reads < config.min_reads_per_bin:
        return "low_support", "insufficient_reads_or_depth"

    mq_high = mq_is_elevated(mq_low_frac, mq_reference, config)
    mq_ok = mq_low_frac <= config.duplication_mq_low_frac_max

    if depth_rel <= config.deletion_depth_rel_max:
        if mq_high:
            return (
                "deletion",
                "low_depth_with_elevated_low_mq_misplaced_reads",
            )
        return "deletion", "low_depth_sparse_signal"

    if depth_rel >= config.duplication_depth_rel_min:
        if mq_ok:
            return "duplication", "high_depth_with_consistent_mq"
        return "duplication", "high_depth_with_noisy_mq"

    if config.disomic_depth_min <= depth_rel <= config.disomic_depth_max and mq_high:
        return (
            "inter_homoeologue",
            "disomic_depth_with_elevated_low_mq_shadow_mapping",
        )

    return "normal", "disomic_depth_with_normal_mq"


def label_bins(df: pd.DataFrame, mq_reference: float, config: CallConfig) -> pd.DataFrame:
    """Attach ``bin_label`` and ``bin_evidence`` columns."""
    out = df.copy()
    labels = []
    evidence = []
    for row in out.itertuples(index=False):
        label, ev = classify_bin_row(
            depth_rel=float(row.depth_rel) if pd.notna(row.depth_rel) else float("nan"),
            mq_low_frac=float(row.mq_low_frac),
            total_reads=float(row.total_reads),
            mq_reference=mq_reference,
            config=config,
        )
        labels.append(label)
        evidence.append(ev)
    out["bin_label"] = labels
    out["bin_evidence"] = evidence
    return out


def chr_aneuploidy_call(
    chr_df: pd.DataFrame,
    config: CallConfig,
) -> tuple[str | None, float, str]:
    """Detect whole-chromosome dosage shift; reject partial-CNV patterns."""
    profile = chr_dosage_profile(chr_df, config)
    if profile is None:
        return None, 0.0, ""

    mq_med = float(pd.to_numeric(chr_df["mq_low_frac"], errors="coerce").median())

    if is_partial_cnv_pattern(profile, config):
        return None, 0.0, (
            f"partial_cnv_pattern;disomic_frac={profile.disomic_frac:.2f};"
            f"loss_frac={profile.loss_frac:.2f};gain_frac={profile.gain_frac:.2f}"
        )

    if is_uniform_monosomic(profile, config):
        conf = min(
            1.0,
            (config.uniform_monosomic_median_max - profile.median_rel) / config.depth_confidence_scale,
        )
        ev = (
            f"uniform_monosomic;chr_median_depth_rel={profile.median_rel:.3f};"
            f"loss_bin_fraction={profile.loss_frac:.2f};disomic_frac={profile.disomic_frac:.2f};"
            f"iqr={profile.iqr_rel:.3f};mq_low_frac_median={mq_med:.3f}"
        )
        if mq_med >= config.mq_low_frac_absolute_min:
            ev += ";elevated_low_mq"
        return "aneuploidy_chr_loss", conf, ev

    if is_uniform_trisomic(profile, config):
        conf = min(
            1.0,
            (profile.median_rel - config.uniform_trisomic_median_min) / config.depth_confidence_scale,
        )
        ev = (
            f"uniform_trisomic;chr_median_depth_rel={profile.median_rel:.3f};"
            f"gain_bin_fraction={profile.gain_frac:.2f};disomic_frac={profile.disomic_frac:.2f};"
            f"iqr={profile.iqr_rel:.3f};mq_low_frac_median={mq_med:.3f}"
        )
        return "aneuploidy_chr_gain", conf, ev

    depth = pd.to_numeric(chr_df["depth_rel"], errors="coerce")
    valid = depth.dropna()
    if valid.empty:
        return None, 0.0, ""

    median_rel = float(valid.median())
    n_bins = len(valid)
    loss_bins = int((valid <= config.deletion_depth_rel_max).sum())
    gain_bins = int((valid >= config.duplication_depth_rel_min).sum())
    loss_frac = loss_bins / n_bins
    gain_frac = gain_bins / n_bins

    if (
        loss_frac >= config.chr_aneuploidy_bin_fraction
        and median_rel <= config.aneuploidy_chr_loss_median_max
        and profile.disomic_frac <= config.chr_disomic_frac_max_for_aneuploidy
    ):
        conf = min(1.0, (config.aneuploidy_chr_loss_median_max - median_rel) / config.depth_confidence_scale)
        ev = (
            f"chr_median_depth_rel={median_rel:.3f};loss_bin_fraction={loss_frac:.2f};"
            f"mq_low_frac_median={mq_med:.3f}"
        )
        if mq_med >= config.mq_low_frac_absolute_min:
            ev += ";elevated_low_mq_supports_mis_mapping_on_reduced_dose"
        return "aneuploidy_chr_loss", conf, ev

    if (
        gain_frac >= config.chr_aneuploidy_bin_fraction
        and median_rel >= config.aneuploidy_chr_gain_median_min
        and profile.disomic_frac <= config.chr_disomic_frac_max_for_aneuploidy
    ):
        conf = min(1.0, (median_rel - config.aneuploidy_chr_gain_median_min) / config.depth_confidence_scale)
        ev = (
            f"chr_median_depth_rel={median_rel:.3f};gain_bin_fraction={gain_frac:.2f};"
            f"mq_low_frac_median={mq_med:.3f};extra_copy_reads_map_with_typical_mq"
        )
        return "aneuploidy_chr_gain", conf, ev

    return None, 0.0, ""


def segment_confidence(
    call_type: str,
    median_depth_rel: float,
    median_mq_low_frac: float,
    mq_reference: float,
    config: CallConfig,
) -> float:
    """Heuristic confidence in [0, 1] from depth and MAPQ separation."""
    if call_type in {"normal", "low_support"}:
        return 0.0

    depth_margin = 0.0
    if call_type == "deletion":
        depth_margin = max(0.0, config.deletion_depth_rel_max - median_depth_rel)
    elif call_type == "duplication":
        depth_margin = max(0.0, median_depth_rel - config.duplication_depth_rel_min)
    elif call_type == "inter_homoeologue":
        depth_margin = max(
            0.0,
            1.0 - abs(median_depth_rel - 1.0),
        ) * 0.5
    elif call_type.startswith("aneuploidy") or call_type.startswith("subgenome"):
        depth_margin = abs(median_depth_rel - 1.0)

    mq_margin = max(0.0, median_mq_low_frac - mq_reference)
    score = depth_margin / config.depth_confidence_scale + mq_margin / config.mq_confidence_scale
    return float(min(1.0, max(0.05, score)))
