"""Plot mosdepth vs breachr reconciliation summaries."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from infra.utils.graph import (
    plot_categorical_bar_chart,
    plot_stacked_histogram,
    plot_stacked_proportion_bar,
)

# Short labels for stacked / bar charts.
_RECONCILE_LABELS = {
    "likely_cnv_not_aneuploidy": "CNV not aneuploidy",
    "possible_inter_homoeologue_not_aneuploidy": "Inter-homoeologue",
    "consistent_whole_chr_aneuploidy": "Consistent aneuploidy",
    "reference_ploidy_mismatch_not_aneuploidy": "Reference ploidy mismatch",
    "focal_cnv_direction_mismatch": "CNV direction mismatch",
    "focal_cnv_ambiguous": "Focal CNV ambiguous",
    "no_breachr_call": "No breachr call",
    "expected_null_d_subgenome": "Expected D-null",
    "direction_mismatch": "Direction mismatch",
    "missing_mapq": "Missing MAPQ",
}

_RECONCILE_COLORS = {
    "likely_cnv_not_aneuploidy": "#d95f02",
    "possible_inter_homoeologue_not_aneuploidy": "#e6ab02",
    "consistent_whole_chr_aneuploidy": "#1b9e77",
    "reference_ploidy_mismatch_not_aneuploidy": "#a6761d",
    "focal_cnv_direction_mismatch": "#7570b3",
    "focal_cnv_ambiguous": "#66a61e",
    "no_breachr_call": "#999999",
    "expected_null_d_subgenome": "#cccccc",
    "direction_mismatch": "#e7298a",
    "missing_mapq": "#dddddd",
}

_BREACHR_CALL_LABELS = {
    "aneuploidy_chr_loss": "Chr loss (monosomy)",
    "aneuploidy_chr_gain": "Chr gain (trisomy)",
    "subgenome_reference_mismatch": "Subgenome ref mismatch",
    "subgenome_absent": "Subgenome absent",
    "deletion": "Focal deletion",
    "duplication": "Focal duplication",
    "inter_homoeologue": "Inter-homoeologue",
}

_BREACHR_CALL_COLORS = {
    "aneuploidy_chr_loss": "#1b9e77",
    "aneuploidy_chr_gain": "#7570b3",
    "subgenome_reference_mismatch": "#a6761d",
    "subgenome_absent": "#666666",
    "deletion": "#d95f02",
    "duplication": "#e7298a",
    "inter_homoeologue": "#e6ab02",
}

_COHORT_STACK_CATEGORIES = [
    "likely_cnv_not_aneuploidy",
    "consistent_whole_chr_aneuploidy",
    "reference_ploidy_mismatch_not_aneuploidy",
    "focal_cnv_direction_mismatch",
    "no_breachr_call",
]

_CNV_MISCLASS = {
    "likely_cnv_not_aneuploidy",
    "possible_inter_homoeologue_not_aneuploidy",
}


def summarize_reconcile_by_sample(reconcile_df: pd.DataFrame) -> pd.DataFrame:
    """One row per flagged sample with chr-flag counts by reconciliation category."""
    if reconcile_df.empty:
        return pd.DataFrame()

    rows: list[dict[str, object]] = []
    for sample, grp in reconcile_df.groupby("sample", sort=True):
        counts = grp["reconciliation"].value_counts()
        n_cnv_misclass = int(counts.reindex(list(_CNV_MISCLASS), fill_value=0).sum())
        n_consistent = int(counts.get("consistent_whole_chr_aneuploidy", 0))
        n_ref_mismatch = int(counts.get("reference_ploidy_mismatch_not_aneuploidy", 0))
        n_total = int(len(grp))
        if n_ref_mismatch == n_total and n_total > 0:
            sample_class = "reference_ploidy_mismatch_only"
        elif n_cnv_misclass > 0 and n_consistent > 0:
            sample_class = "mixed_cnv_and_aneuploidy"
        elif n_cnv_misclass > 0:
            sample_class = "cnv_misclassification"
        elif n_consistent > 0:
            sample_class = "consistent_aneuploidy_only"
        elif n_ref_mismatch > 0:
            sample_class = "mixed_reference_and_other"
        else:
            sample_class = "other_or_ambiguous"

        rows.append(
            {
                "sample": sample,
                "n_chr_flags": n_total,
                "n_cnv_misclassified_chr": n_cnv_misclass,
                "n_consistent_aneuploidy_chr": n_consistent,
                "n_reference_ploidy_mismatch_chr": n_ref_mismatch,
                "n_other_chr": n_total - n_cnv_misclass - n_consistent - n_ref_mismatch,
                "sample_class": sample_class,
            },
        )
    return pd.DataFrame(rows)


def plot_breachr_cnv_length_distribution(
    segments_all_path: str | Path,
    plot_dir: str | Path,
    *,
    max_length_mb: float = 120.0,
    bins: int = 60,
) -> list[Path]:
    """Histogram of focal CNV segment lengths (deletion / duplication) across all samples."""
    segments_all_path = Path(segments_all_path)
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    segments = pd.read_csv(segments_all_path, sep="\t")
    focal = segments[segments["call_type"].isin(["deletion", "duplication"])].copy()
    if focal.empty:
        return []

    focal["length_mb"] = pd.to_numeric(focal["length_bp"], errors="coerce") / 1e6
    focal = focal[focal["length_mb"].notna()]
    focal["call_type"] = focal["call_type"].map(
        {"deletion": "Deletion", "duplication": "Duplication"},
    )

    n_segments = len(focal)
    n_samples = focal["sample"].nunique()
    median_mb = float(focal["length_mb"].median())
    mean_mb = float(focal["length_mb"].mean())
    title_base = (
        f"Focal CNV length distribution (n={n_segments} segments, "
        f"{n_samples} samples)"
    )
    common = dict(
        data=focal,
        x_col="length_mb",
        hue_col="call_type",
        hue_order=["Deletion", "Duplication"],
        palette_map={"Deletion": "#d95f02", "Duplication": "#e7298a"},
        xlabel="Segment length (Mb)",
        ylabel="Number of CNV segments",
        mean_val=mean_mb,
        median_val=median_mb,
        xlim=(0.0, max_length_mb),
        legend_title="CNV type",
    )

    written: list[Path] = []
    out_log = plot_dir / "cnv_length_distribution.png"
    plot_stacked_histogram(
        **common,
        title=title_base,
        filename=str(out_log),
        log_scale=True,
    )
    written.append(out_log)

    out_linear = plot_dir / "cnv_length_distribution_linear.png"
    plot_stacked_histogram(
        **common,
        title=title_base,
        filename=str(out_linear),
        log_scale=False,
    )
    written.append(out_linear)

    return written


def plot_breachr_reconcile(
    reconcile_all_path: str | Path,
    flagged_path: str | Path,
    plot_dir: str | Path,
    *,
    segments_all_path: str | Path | None = None,
    sample_summary_path: str | Path | None = None,
) -> list[Path]:
    """Write reconciliation summary plots for flagged aneuploidy benchmark."""
    reconcile_all_path = Path(reconcile_all_path)
    flagged_path = Path(flagged_path)
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    reconcile = pd.read_csv(reconcile_all_path, sep="\t")
    flagged = pd.read_csv(flagged_path, sep="\t")
    flagged = flagged[flagged["aneuploid_flag"].astype(str).str.lower().isin({"true", "1", "yes"})]
    cohort_map = flagged.set_index("sample")["cohort"].astype(str).to_dict()
    reconcile = reconcile.copy()
    reconcile["cohort"] = reconcile["sample"].map(cohort_map).fillna("unknown")

    sample_summary = summarize_reconcile_by_sample(reconcile)
    flagged_samples = flagged["sample"].astype(str).unique()
    missing_samples = set(flagged_samples) - set(sample_summary["sample"].astype(str))
    if missing_samples:
        extra = pd.DataFrame(
            [
                {
                    "sample": s,
                    "n_chr_flags": 0,
                    "n_cnv_misclassified_chr": 0,
                    "n_consistent_aneuploidy_chr": 0,
                    "n_reference_ploidy_mismatch_chr": 0,
                    "n_other_chr": 0,
                    "sample_class": "other_or_ambiguous",
                }
                for s in sorted(missing_samples)
            ],
        )
        sample_summary = pd.concat([sample_summary, extra], ignore_index=True)
    if sample_summary_path is not None:
        sample_summary.to_csv(sample_summary_path, sep="\t", index=False)

    written: list[Path] = []

    # 1) Chr-level reconciliation counts.
    chr_counts = reconcile["reconciliation"].value_counts()
    categories = [str(c) for c in chr_counts.index]
    labels = [_RECONCILE_LABELS.get(c, c) for c in categories]
    colors = [_RECONCILE_COLORS.get(c, "steelblue") for c in categories]
    out1 = plot_dir / "reconciliation_chr_level_counts.png"
    plot_categorical_bar_chart(
        labels,
        [float(chr_counts[c]) for c in categories],
        title="Flagged chr outcomes (mosdepth vs breachr)",
        ylabel="Number of chr flags",
        filename=str(out1),
        colors=colors,
        ylim=(0, max(chr_counts.max() * 1.15, 1)),
        rotate_xlabels=35,
        figure_size=(12, 6),
    )
    written.append(out1)

    # 2) CNV misclassification by logical chromosome.
    mis = reconcile[reconcile["reconciliation"].isin(_CNV_MISCLASS)]
    if not mis.empty:
        by_chr = mis.groupby("logical_chr").size().sort_values(ascending=False)
        out2 = plot_dir / "cnv_misclassified_by_chromosome.png"
        plot_categorical_bar_chart(
            [str(c) for c in by_chr.index],
            [float(v) for v in by_chr.values],
            title="CNV misclassified as aneuploidy — by chromosome",
            ylabel="Number of chr flags",
            filename=str(out2),
            colors=["#d95f02"] * len(by_chr),
            ylim=(0, max(by_chr.max() * 1.15, 1)),
            rotate_xlabels=45,
            figure_size=(12, 6),
        )
        written.append(out2)

    # 3) Sample-level class counts (261 flagged samples).
    if not sample_summary.empty:
        cls_counts = sample_summary["sample_class"].value_counts()
        cls_labels = {
            "cnv_misclassification": "CNV misclassification",
            "consistent_aneuploidy_only": "Consistent aneuploidy only",
            "mixed_cnv_and_aneuploidy": "Mixed CNV + aneuploidy",
            "reference_ploidy_mismatch_only": "Reference ploidy mismatch only",
            "mixed_reference_and_other": "Mixed reference + other",
            "other_or_ambiguous": "Other / ambiguous",
        }
        cls_colors = {
            "cnv_misclassification": "#d95f02",
            "consistent_aneuploidy_only": "#1b9e77",
            "mixed_cnv_and_aneuploidy": "#7570b3",
            "reference_ploidy_mismatch_only": "#a6761d",
            "mixed_reference_and_other": "#666666",
            "other_or_ambiguous": "#999999",
        }
        cats = [str(c) for c in cls_counts.index]
        out3 = plot_dir / "sample_level_outcome_counts.png"
        plot_categorical_bar_chart(
            [cls_labels.get(c, c) for c in cats],
            [float(cls_counts[c]) for c in cats],
            title=f"Flagged samples by primary breachr outcome (n={len(sample_summary)})",
            ylabel="Number of samples",
            filename=str(out3),
            colors=[cls_colors.get(c, "steelblue") for c in cats],
            ylim=(0, max(cls_counts.max() * 1.15, 1)),
            rotate_xlabels=25,
            figure_size=(10, 6),
        )
        written.append(out3)

    # 4) Cohort stacked proportions (chr-level).
    cohort_order = sorted(reconcile["cohort"].unique())
    cohort_rows: dict[str, dict[str, float]] = {}
    for cohort in cohort_order:
        sub = reconcile[reconcile["cohort"] == cohort]
        total = len(sub)
        if total == 0:
            continue
        counts = sub["reconciliation"].value_counts()
        other = total - counts.reindex(_COHORT_STACK_CATEGORIES, fill_value=0).sum()
        cohort_rows[cohort] = {
            _RECONCILE_LABELS.get(cat, cat): counts.get(cat, 0) / total
            for cat in _COHORT_STACK_CATEGORIES
        }
        cohort_rows[cohort]["Other"] = other / total
    if cohort_rows:
        prop_df = pd.DataFrame(cohort_rows).T
        out4 = plot_dir / "reconciliation_by_cohort_stacked.png"
        plot_stacked_proportion_bar(
            prop_df,
            title="Chr-flag reconciliation fraction by cohort",
            filename=str(out4),
            x_label="Cohort",
            y_label="Fraction of chr flags",
            stack_col_title="Outcome",
            figure_size=(12, 6),
        )
        written.append(out4)

    # 5) Samples with ≥1 CNV-misclassified chr — chr count distribution.
    cnv_samples = sample_summary[sample_summary["n_cnv_misclassified_chr"] > 0]
    if not cnv_samples.empty:
        hist = cnv_samples["n_cnv_misclassified_chr"].value_counts().sort_index()
        out5 = plot_dir / "cnv_misclassified_chr_per_sample.png"
        plot_categorical_bar_chart(
            [f"{int(k)} chr" for k in hist.index],
            [float(v) for v in hist.values],
            title="CNV-misclassified samples — chr flags per sample",
            ylabel="Number of samples",
            filename=str(out5),
            colors=["#d95f02"] * len(hist),
            ylim=(0, max(hist.max() * 1.15, 1)),
            rotate_xlabels=0,
            figure_size=(8, 5),
        )
        written.append(out5)

    # 6) Breachr call type on flagged chr rows (primary call per reconciliation row).
    call_mask = reconcile["breachr_call_type"].astype(str).str.len() > 0
    call_subset = reconcile.loc[call_mask].copy()
    if not call_subset.empty:
        call_counts = call_subset["breachr_call_type"].value_counts()
        call_cats = [str(c) for c in call_counts.index]
        out6 = plot_dir / "breachr_call_type_on_flagged_chr.png"
        plot_categorical_bar_chart(
            [_BREACHR_CALL_LABELS.get(c, c) for c in call_cats],
            [float(call_counts[c]) for c in call_cats],
            title="Breachr primary call type on mosdepth-flagged chromosomes",
            ylabel="Number of chr flags",
            filename=str(out6),
            colors=[_BREACHR_CALL_COLORS.get(c, "steelblue") for c in call_cats],
            ylim=(0, max(call_counts.max() * 1.15, 1)),
            rotate_xlabels=35,
            figure_size=(12, 6),
        )
        written.append(out6)

    segments_path = Path(segments_all_path) if segments_all_path is not None else (
        reconcile_all_path.parent / "breachr_segments_all.tsv"
    )
    if segments_path.exists():
        written.extend(
            plot_breachr_cnv_length_distribution(
                segments_all_path=segments_path,
                plot_dir=plot_dir,
            ),
        )

    return written
