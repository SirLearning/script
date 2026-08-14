"""Genome-grid MAPQ depth cap plots with breachr annotations for flagged samples."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from breachr.io import load_mapq_depth_bin_tsv
from breachr.reconcile import parse_abnormal_chr_list
from genetics.germplasm.sample.bam_mapq_depth_bin import plot_mapq_depth_cap_genome_grid
from infra.utils.graph import TICK_FONT_SIZE

_RECONCILE_SHORT = {
    "likely_cnv_not_aneuploidy": "CNV not aneuploidy",
    "possible_inter_homoeologue_not_aneuploidy": "Inter-homoeologue",
    "consistent_whole_chr_aneuploidy": "Consistent aneuploidy",
    "focal_cnv_direction_mismatch": "CNV dir. mismatch",
    "focal_cnv_ambiguous": "Focal CNV ambiguous",
    "no_breachr_call": "No breachr call",
    "expected_null_d_subgenome": "Expected D-null",
    "direction_mismatch": "Direction mismatch",
    "missing_mapq": "Missing MAPQ",
}

_SEGMENT_COLORS = {
    "deletion": "#2166ac",
    "duplication": "#d95f02",
    "inter_homoeologue": "#e6ab02",
    "aneuploidy_chr_loss": "#1b7837",
    "aneuploidy_chr_gain": "#5aae61",
    "subgenome_reference_mismatch": "#984ea3",
    "subgenome_absent": "#999999",
}


def _primary_segment_label(chr_segments: pd.DataFrame) -> str:
    if chr_segments.empty:
        return ""
    focal = chr_segments[chr_segments["call_type"].isin(_SEGMENT_COLORS)].copy()
    if focal.empty:
        return ""
    row = focal.sort_values(["length_bp", "confidence"], ascending=False).iloc[0]
    length_mb = int(row["length_bp"]) / 1_000_000
    return (
        f"{row['call_type']} {length_mb:.1f} Mb "
        f"(depth_rel={float(row['median_depth_rel']):.2f})"
    )


def build_chr_annotations(
    sample: str,
    flagged_df: pd.DataFrame,
    reconcile_df: pd.DataFrame,
    segments_df: pd.DataFrame,
) -> dict[str, dict[str, str]]:
    """Build per-chromosome annotation payloads for genome cap plots."""
    flagged_row = flagged_df[flagged_df["sample"].astype(str) == sample]
    abnormal: dict[str, dict[str, object]] = {}
    if not flagged_row.empty:
        for entry in parse_abnormal_chr_list(flagged_row.iloc[0]["abnormal_chr_list"]):
            abnormal[str(entry["logical_chr"])] = entry

    sample_reconcile = reconcile_df[reconcile_df["sample"].astype(str) == sample]
    sample_segments = segments_df[segments_df["sample"].astype(str) == sample]

    annotations: dict[str, dict[str, str]] = {}
    for logical_chr in abnormal:
        mosdepth_label = str(abnormal[logical_chr]["mosdepth_label"])
        mosdepth_rel = float(abnormal[logical_chr]["mosdepth_rel_depth"])
        rec_rows = sample_reconcile[sample_reconcile["logical_chr"] == logical_chr]
        chr_segments = sample_segments[sample_segments["logical_chr"] == logical_chr]

        reconciliation = ""
        breachr_call = ""
        if not rec_rows.empty:
            rec = rec_rows.iloc[0]
            reconciliation = _RECONCILE_SHORT.get(str(rec["reconciliation"]), str(rec["reconciliation"]))
            if str(rec["breachr_call_type"]):
                length_mb = int(rec["breachr_length_bp"]) / 1_000_000 if int(rec["breachr_length_bp"]) > 0 else 0.0
                breachr_call = str(rec["breachr_call_type"])
                if length_mb > 0:
                    breachr_call += f" {length_mb:.1f} Mb"
        else:
            breachr_call = _primary_segment_label(chr_segments)

        annotations[logical_chr] = {
            "mosdepth_title": f"{logical_chr}:{mosdepth_label}({mosdepth_rel:.2f})",
            "reconciliation": reconciliation,
            "breachr_call": breachr_call,
        }
    return annotations


def _overlay_breachr_segments(
    ax: plt.Axes,
    chr_segments: pd.DataFrame,
    y_cap: float,
) -> None:
    for row in chr_segments.itertuples(index=False):
        call_type = str(row.call_type)
        if call_type not in _SEGMENT_COLORS:
            continue
        start_mb = float(row.bin_start) / 1_000_000
        end_mb = float(row.bin_end) / 1_000_000
        ax.axvspan(
            start_mb,
            end_mb,
            color=_SEGMENT_COLORS[call_type],
            alpha=0.12,
            linewidth=0,
        )


def plot_breachr_flagged_mapq_depth_cap(
    mapq_tsv_path: str | Path,
    sample: str,
    flagged_path: str | Path,
    reconcile_all_path: str | Path,
    segments_path: str | Path,
    output_png: str | Path,
    *,
    baseline_depth: float | None = None,
    y_cap_multiplier: float = 3.0,
) -> Path:
    """Plot genome MAPQ depth cap grid with mosdepth + breachr chr annotations."""
    mapq_tsv_path = Path(mapq_tsv_path)
    output_png = Path(output_png)
    flagged_df = pd.read_csv(flagged_path, sep="\t")
    reconcile_df = pd.read_csv(reconcile_all_path, sep="\t")
    segments_df = pd.read_csv(segments_path, sep="\t")

    if baseline_depth is None:
        flagged_row = flagged_df[flagged_df["sample"].astype(str) == sample]
        if not flagged_row.empty and pd.notna(flagged_row.iloc[0].get("baseline_depth")):
            baseline_depth = float(flagged_row.iloc[0]["baseline_depth"])

    df = load_mapq_depth_bin_tsv(mapq_tsv_path)
    annotations = build_chr_annotations(sample, flagged_df, reconcile_df, segments_df)

    def _panel_hook(logical_chr: str, ax: plt.Axes, chr_df: pd.DataFrame, y_cap: float) -> None:
        if logical_chr not in annotations:
            return
        ann = annotations[logical_chr]
        ax.set_title(ann["mosdepth_title"], fontsize=TICK_FONT_SIZE, color="#cc0000", pad=14)
        lines = []
        if ann["reconciliation"]:
            lines.append(f"breachr: {ann['reconciliation']}")
        if ann["breachr_call"]:
            lines.append(str(ann["breachr_call"]))
        if lines:
            ax.text(
                0.02,
                0.98,
                "\n".join(lines),
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=max(6, TICK_FONT_SIZE - 3),
                color="#333333",
                bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none", "pad": 2},
            )
        chr_segments = segments_df[
            (segments_df["sample"].astype(str) == sample)
            & (segments_df["logical_chr"] == logical_chr)
        ]
        _overlay_breachr_segments(ax, chr_segments, y_cap)

    return plot_mapq_depth_cap_genome_grid(
        df,
        sample=sample,
        output_png=output_png,
        baseline_depth=baseline_depth,
        y_cap_multiplier=y_cap_multiplier,
        panel_hook=_panel_hook,
    )


def plot_all_flagged_mapq_depth_caps(
    flagged_path: str | Path,
    mapq_dir: str | Path,
    segments_dir: str | Path,
    reconcile_all_path: str | Path,
    plot_dir: str | Path,
    *,
    resolve_mapq_sample=None,
) -> list[Path]:
    """Generate annotated cap plots for all aneuploidy-flagged samples."""
    flagged_path = Path(flagged_path)
    mapq_dir = Path(mapq_dir)
    segments_dir = Path(segments_dir)
    plot_dir = Path(plot_dir)
    flagged_df = pd.read_csv(flagged_path, sep="\t")
    flagged_df = flagged_df[flagged_df["aneuploid_flag"].astype(str).str.lower().isin({"true", "1", "yes"})]

    resolve = resolve_mapq_sample or (lambda sample: sample)
    written: list[Path] = []
    for sample in flagged_df["sample"].astype(str):
        mapq_id = str(resolve(sample))
        mapq_tsv = mapq_dir / f"{mapq_id}.mapq_depth_bin.tsv"
        segments_tsv = segments_dir / f"{sample}.breachr_segments.tsv"
        if not mapq_tsv.exists() or not segments_tsv.exists():
            continue
        out_dir = plot_dir / sample
        out_dir.mkdir(parents=True, exist_ok=True)
        out_png = out_dir / f"{sample}.mapq_depth_cap.genome.png"
        written.append(
            plot_breachr_flagged_mapq_depth_cap(
                mapq_tsv_path=mapq_tsv,
                sample=sample,
                flagged_path=flagged_path,
                reconcile_all_path=reconcile_all_path,
                segments_path=segments_tsv,
                output_png=out_png,
            ),
        )
    return written
