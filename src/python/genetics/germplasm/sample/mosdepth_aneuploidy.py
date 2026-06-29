"""Detect chromosome-level aneuploidy from mosdepth ``*.summary.txt`` files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.germplasm.sample.anno import anno_group
from infra.utils.graph import (
    LEGEND_FONT_SIZE,
    TITLE_FONT_SIZE,
    TICK_FONT_SIZE,
    X_LABEL_FONT_SIZE,
    Y_LABEL_FONT_SIZE,
    plot_categorical_bar_chart,
    plot_correlation_with_regression,
    plot_stacked_distribution,
    sample_group_hue_config,
)
from infra.utils.io import load_df_generic
from infra.wheat import ref_v1

# Logical chromosomes chr1A–chr7D (21 wheat autosomes in ref v1 naming).
LOGICAL_AUTOSOMES = tuple(
    f"chr{i}{genome}" for i in range(1, 8) for genome in ("A", "B", "D")
)

SUBGENOME_AUTOSOMES = {
    sg: tuple(f"chr{i}{sg}" for i in range(1, 8)) for sg in ("A", "B", "D")
}

CHR_DEPTH_COLUMNS = (
    "logical_chr",
    "subgenome",
    "depth",
    "length",
    "sample",
    "cohort",
    "expected",
    "baseline",
    "rel_depth",
    "cn_class",
)

# Cohort folder -> subgenomes expected at disomic dosage (relative to cohort baseline).
COHORT_EXPECTED_SUBGENOMES = {
    "01A": ("A",),
    "02AB": ("A", "B"),
    "03ABD": ("A", "B", "D"),
    "04D": ("D",),
    "05HZNU": ("A", "B", "D"),
    "06Nature": ("A", "B", "D"),
    "07S": ("A", "B", "D"),
    "08WAP": ("A", "B", "D"),
    "09Watkins": ("A", "B", "D"),
}

# Relative depth vs cohort-autosome median: copy-number style labels.
CN_CLASS_BOUNDS = (
    ("null", 0.0, 0.20),
    ("monosomic", 0.20, 0.60),
    ("disomic", 0.60, 1.40),
    ("trisomic", 1.40, 1.85),
    ("tetrasomic+", 1.85, float("inf")),
)

CN_CLASS_BG_COLORS = {
    "null": "#d9d9d9",
    "monosomic": "#9ecae1",
    "disomic": "#b6d7a8",
    "trisomic": "#f9cb9c",
    "tetrasomic+": "#ea9999",
}

UNEXPECTED_SUBGENOME_RATIO = 0.35  # off-subgenome median / expected median


@dataclass(frozen=True)
class ChrDepthRecord:
    sample: str
    cohort: str
    logical_chr: str
    subgenome: str
    depth: float
    baseline: float
    rel_depth: float
    cn_class: str


def _logical_chr_segment_ids(logical_chr: str) -> list[int]:
    return [
        cid
        for cid, name in ref_v1._REF_V1_PLINK_TO_NAME.items()
        if name == logical_chr
    ]


def expected_logical_chromosomes(subgenomes: tuple[str, ...]) -> list[str]:
    genomes = set(subgenomes)
    return [c for c in LOGICAL_AUTOSOMES if c[-1] in genomes]


def parse_logical_chromosome_depths(summary_path: str | Path) -> pd.DataFrame:
    """Length-weighted mean depth per logical chromosome (chr1A, chr2B, …)."""
    df = load_df_generic(str(summary_path))
    if df is None or df.empty:
        return pd.DataFrame(columns=["logical_chr", "subgenome", "depth", "length"])

    df = df.copy()
    df["chrom"] = df["chrom"].astype(str)
    df = df[df["chrom"] != "total"]
    df["chrom"] = df["chrom"].astype(int)

    rows = []
    for logical in LOGICAL_AUTOSOMES:
        seg_ids = _logical_chr_segment_ids(logical)
        sub = df[df["chrom"].isin(seg_ids)]
        total_len = float(sub["length"].sum())
        if total_len <= 0:
            depth = float("nan")
        else:
            depth = float(sub["bases"].sum() / total_len)
        rows.append(
            {
                "logical_chr": logical,
                "subgenome": logical[-1],
                "depth": depth,
                "length": total_len,
            }
        )
    return pd.DataFrame(rows)


def classify_relative_depth(rel: float) -> str:
    if not np.isfinite(rel):
        return "missing"
    for label, lo, hi in CN_CLASS_BOUNDS:
        if lo <= rel < hi:
            return label
    return "missing"


def analyze_sample_summary(
    summary_path: str | Path,
    sample: str,
    cohort: str,
    expected_subgenomes: tuple[str, ...] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Return (per-chromosome table, sample-level summary) for one mosdepth summary.

    Baseline = median depth across expected logical autosomes for the cohort.
    """
    expected_subgenomes = expected_subgenomes or COHORT_EXPECTED_SUBGENOMES.get(
        cohort, ("A", "B", "D")
    )
    expected_chrs = expected_logical_chromosomes(expected_subgenomes)
    chr_df = parse_logical_chromosome_depths(summary_path)
    if chr_df.empty:
        return chr_df, pd.DataFrame()

    exp = chr_df[chr_df["logical_chr"].isin(expected_chrs)].copy()
    baseline = float(exp["depth"].median())
    if not np.isfinite(baseline) or baseline <= 0:
        baseline = float("nan")

    chr_df["sample"] = sample
    chr_df["cohort"] = cohort
    chr_df["expected"] = chr_df["logical_chr"].isin(expected_chrs)
    chr_df["baseline"] = baseline
    chr_df["rel_depth"] = chr_df["depth"] / baseline
    chr_df["cn_class"] = chr_df["rel_depth"].map(classify_relative_depth)

    off = chr_df[~chr_df["expected"] & chr_df["subgenome"].isin({"A", "B", "D"})]
    off_medians = off.groupby("subgenome")["rel_depth"].median() if not off.empty else pd.Series(dtype=float)

    abnormal = chr_df[
        chr_df["expected"] & ~chr_df["cn_class"].isin(("disomic",))
    ]
    sample_row = {
        "sample": sample,
        "cohort": cohort,
        "summary_path": str(summary_path),
        "baseline_depth": baseline,
        "n_expected_chr": len(expected_chrs),
        "n_abnormal_chr": int((chr_df["expected"] & ~chr_df["cn_class"].eq("disomic")).sum()),
        "abnormal_chr_list": ",".join(
            f"{r.logical_chr}:{r.cn_class}({r.rel_depth:.2f})"
            for r in abnormal.itertuples()
        ),
        "max_unexpected_subgenome": float(off_medians.max()) if len(off_medians) else 0.0,
        "unexpected_subgenome_flag": bool(
            len(off_medians) and float(off_medians.max()) >= UNEXPECTED_SUBGENOME_RATIO
        ),
        "aneuploid_flag": bool(
            (chr_df["expected"] & ~chr_df["cn_class"].eq("disomic")).any()
            or (len(off_medians) and float(off_medians.max()) >= UNEXPECTED_SUBGENOME_RATIO)
        ),
    }
    for sg in ("A", "B", "D"):
        sample_row[f"unexpected_{sg}_rel_median"] = float(off_medians.get(sg, np.nan))

    return chr_df, pd.DataFrame([sample_row])


def scan_cohort_depth_dir(
    cohort_dir: str | Path,
    cohort: str | None = None,
    expected_subgenomes: tuple[str, ...] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Scan all ``*.mosdepth.summary.txt`` files in a cohort folder."""
    cohort_dir = Path(cohort_dir)
    cohort = cohort or cohort_dir.name
    expected_subgenomes = expected_subgenomes or COHORT_EXPECTED_SUBGENOMES.get(
        cohort, ("A", "B", "D")
    )

    chr_parts: list[pd.DataFrame] = []
    sample_parts: list[pd.DataFrame] = []
    for path in sorted(cohort_dir.glob("*.mosdepth.summary.txt")):
        sample = normalize_sample_from_summary_stem(path.name)
        chr_df, sample_df = analyze_sample_summary(
            path, sample=sample, cohort=cohort, expected_subgenomes=expected_subgenomes
        )
        if not chr_df.empty:
            chr_parts.append(chr_df)
        if not sample_df.empty:
            sample_parts.append(sample_df)

    chr_all = pd.concat(chr_parts, ignore_index=True) if chr_parts else pd.DataFrame()
    sample_all = pd.concat(sample_parts, ignore_index=True) if sample_parts else pd.DataFrame()
    return chr_all, sample_all


def normalize_sample_id(sample: str) -> str:
    """Canonical taxa id for joins (handles ``*_deduped.bam`` mosdepth summary stems)."""
    s = str(sample)
    suffixes = (
        ".mosdepth.summary.txt",
        ".rmdup.bam",
        "_deduped.bam",
        ".bam",
        "_deduped",
        ".rmdup",
    )
    changed = True
    while changed:
        changed = False
        for suffix in suffixes:
            if s.endswith(suffix):
                s = s[: -len(suffix)]
                changed = True
                break
    return s


def normalize_sample_from_summary_stem(stem: str) -> str:
    """Derive sample id from a mosdepth summary filename stem."""
    return normalize_sample_id(stem)


def aggregate_aneuploidy_tables(
    chr_paths: list[str | Path],
    sample_paths: list[str | Path],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Merge per-sample TSVs into cohort-wide tables and a flagged subset."""
    chr_parts = [pd.read_csv(p, sep="\t") for p in chr_paths if Path(p).exists()]
    sample_parts = [pd.read_csv(p, sep="\t") for p in sample_paths if Path(p).exists()]
    chr_all = pd.concat(chr_parts, ignore_index=True) if chr_parts else pd.DataFrame()
    sample_all = pd.concat(sample_parts, ignore_index=True) if sample_parts else pd.DataFrame()
    if sample_all.empty:
        flagged = pd.DataFrame()
    else:
        flagged = sample_all[sample_all["aneuploid_flag"]].sort_values(
            ["cohort", "n_abnormal_chr"], ascending=[True, False]
        )
    return chr_all, sample_all, flagged


def split_chr_depth_by_chromosome(chr_all: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """One table per logical chromosome; rows are samples, columns match CHR_DEPTH_COLUMNS."""
    if chr_all.empty:
        return {}
    out: dict[str, pd.DataFrame] = {}
    for logical_chr in LOGICAL_AUTOSOMES:
        sub = chr_all[chr_all["logical_chr"] == logical_chr].copy()
        if sub.empty:
            continue
        sub = sub.sort_values(["cohort", "sample"]).reset_index(drop=True)
        for col in CHR_DEPTH_COLUMNS:
            if col not in sub.columns:
                sub[col] = np.nan
        out[logical_chr] = sub[list(CHR_DEPTH_COLUMNS)]
    return out


def attach_sample_group(
    chr_all: pd.DataFrame,
    group_file: str | Path | None,
) -> pd.DataFrame:
    """Add ``Group`` from ``sample_group.txt`` (via ``anno_group``)."""
    if chr_all.empty or not group_file or not Path(group_file).exists():
        out = chr_all.copy()
        if "Group" not in out.columns:
            out["Group"] = "Unknown"
        return out
    tmp = chr_all.rename(columns={"sample": "Sample"})
    tmp["Sample"] = tmp["Sample"].map(normalize_sample_id)
    annotated = anno_group(tmp, group_file=str(group_file), save_tsv=False)
    if annotated is None:
        out = chr_all.copy()
        out["Group"] = "Unknown"
        return out
    if "sample" not in annotated.columns and "Sample" in annotated.columns:
        annotated = annotated.rename(columns={"Sample": "sample"})
    return annotated


def write_per_chromosome_depth_tables(
    chr_all: pd.DataFrame,
    out_dir: str | Path,
) -> list[Path]:
    """Write ``{logical_chr}.depth.tsv`` — one file per chromosome, one row per sample."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for logical_chr, sub in split_chr_depth_by_chromosome(chr_all).items():
        path = out_dir / f"{logical_chr}.depth.tsv"
        sub.to_csv(path, sep="\t", index=False)
        written.append(path)
    return written


def cn_class_background_intervals(xlim: tuple[float, float] | None = (0.0, 3.0)) -> list[dict]:
    """Background spans for rel_depth CN classes (clipped to ``xlim``)."""
    x_hi = xlim[1] if xlim is not None else 3.0
    intervals = []
    for label, lo, hi in CN_CLASS_BOUNDS:
        intervals.append(
            {
                "label": label,
                "lo": lo,
                "hi": hi if np.isfinite(hi) else x_hi,
                "color": CN_CLASS_BG_COLORS.get(label, "#eeeeee"),
                "alpha": 0.4,
            }
        )
    return intervals


def _evaluable_on_chromosome(sub: pd.DataFrame) -> pd.DataFrame:
    """Expected chromosome rows with finite relative depth."""
    rel = pd.to_numeric(sub["rel_depth"], errors="coerce")
    return sub[sub["expected"] & np.isfinite(rel)].copy()


def _abnormal_cn_mask(sub: pd.DataFrame) -> pd.Series:
    """True when expected chromosome copy class is not disomic."""
    return ~sub["cn_class"].isin(("disomic", "missing"))


def plot_per_chromosome_flagged_group_charts(
    chr_all: pd.DataFrame,
    plot_dir: str | Path,
    group_file: str | Path | None = None,
) -> list[Path]:
    """
    Per chromosome: (1) sample_group composition among flagged samples;
    (2) flagged fraction within each sample_group.
    """
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)
    if chr_all.empty:
        return []

    plot_df = attach_sample_group(chr_all, group_file)
    written: list[Path] = []

    for logical_chr in LOGICAL_AUTOSOMES:
        sub = _evaluable_on_chromosome(plot_df[plot_df["logical_chr"] == logical_chr])
        if sub.empty:
            continue

        sub = sub.copy()
        sub["_abnormal"] = _abnormal_cn_mask(sub)
        flagged = sub[sub["_abnormal"]]
        n_flag = len(flagged)

        hue_order, palette_map = sample_group_hue_config(sub["Group"])

        if n_flag > 0:
            flag_counts = flagged["Group"].value_counts()
            groups_present = [g for g in hue_order if g in flag_counts.index]
            groups_present.extend(sorted(g for g in flag_counts.index if g not in groups_present))
            count_vals = [int(flag_counts.get(g, 0)) for g in groups_present]
            count_colors = [palette_map.get(g, (0.5, 0.5, 0.5)) for g in groups_present]
            count_path = plot_dir / f"{logical_chr}.flagged_group.count.png"
            plot_categorical_bar_chart(
                groups_present,
                count_vals,
                title=f"{logical_chr} flagged samples by Group (n={n_flag})",
                ylabel="Flagged sample count",
                filename=str(count_path),
                colors=count_colors,
                ylim=(0, max(count_vals) * 1.15 if count_vals else 1),
                value_formatter=lambda v: f"{int(v)}",
            )
            written.append(count_path)

        rate_rows = []
        for g in hue_order:
            gsub = sub[sub["Group"] == g]
            if gsub.empty:
                continue
            rate_rows.append(
                {"Group": g, "rate_pct": 100.0 * float(gsub["_abnormal"].mean()), "n": len(gsub)}
            )
        extras = sorted(set(sub["Group"]) - set(hue_order))
        for g in extras:
            gsub = sub[sub["Group"] == g]
            rate_rows.append(
                {"Group": g, "rate_pct": 100.0 * float(gsub["_abnormal"].mean()), "n": len(gsub)}
            )
        if not rate_rows:
            continue
        rate_df = pd.DataFrame(rate_rows)
        rate_groups = rate_df["Group"].tolist()
        rate_vals = rate_df["rate_pct"].tolist()
        rate_colors = [palette_map.get(g, (0.5, 0.5, 0.5)) for g in rate_groups]
        rate_path = plot_dir / f"{logical_chr}.flagged_rate_by_group.png"
        plot_categorical_bar_chart(
            rate_groups,
            rate_vals,
            title=f"{logical_chr} flagged rate by Group",
            ylabel="Flagged samples (%)",
            filename=str(rate_path),
            colors=rate_colors,
            ylim=(0, min(100.0, max(rate_vals) * 1.2 + 1) if rate_vals else 100),
            value_formatter=lambda v: f"{v:.1f}%",
        )
        written.append(rate_path)

    return written


def plot_per_chromosome_rel_depth_distributions(
    chr_all: pd.DataFrame,
    plot_dir: str | Path,
    group_file: str | Path | None = None,
    *,
    bins: int = 80,
    xlim: tuple[float, float] | None = (0.0, 3.0),
) -> list[Path]:
    """Stacked rel_depth histograms per logical chromosome, colored by sample ``Group``."""
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)
    if chr_all.empty:
        return []

    plot_df = attach_sample_group(chr_all, group_file)
    written: list[Path] = []
    bg_intervals = cn_class_background_intervals(xlim)

    for logical_chr in LOGICAL_AUTOSOMES:
        sub = _evaluable_on_chromosome(plot_df[plot_df["logical_chr"] == logical_chr])
        if sub.empty:
            continue
        n_samples = sub["sample"].nunique()
        title_base = f"{logical_chr} relative depth (n={n_samples})"
        plot_specs = (
            (f"{logical_chr}.rel_depth.dist.png", False),
            (f"{logical_chr}.rel_depth.dist.logy.png", True),
        )
        for filename, log_scale in plot_specs:
            path = plot_dir / filename
            plot_stacked_distribution(
                sub,
                col="rel_depth",
                group_col="Group",
                title=title_base,
                filename=str(path),
                x_label="Relative depth (vs cohort baseline)",
                y_label="Sample count",
                bins=bins,
                xlim=xlim,
                auto_xlim=xlim is None,
                log_scale=log_scale,
                background_intervals=bg_intervals,
            )
            written.append(path)
    return written


def plot_per_chromosome_rel_depth_vs_baseline(
    chr_all: pd.DataFrame,
    plot_dir: str | Path,
    group_file: str | Path | None = None,
    *,
    ylim: tuple[float, float] | None = (0.0, 3.0),
) -> list[Path]:
    """
    Per logical chromosome: scatter of genome-wide mean depth (baseline) vs relative depth.

    One point per sample where the chromosome is expected and depths are finite.
    Points are colored by sample ``Group`` when ``group_file`` is provided.
    """
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)
    if chr_all.empty:
        return []

    plot_df = attach_sample_group(chr_all, group_file)
    written: list[Path] = []
    for logical_chr in LOGICAL_AUTOSOMES:
        sub = _evaluable_on_chromosome(plot_df[plot_df["logical_chr"] == logical_chr])
        if sub.empty:
            continue
        sub = sub.copy()
        sub["baseline"] = pd.to_numeric(sub["baseline"], errors="coerce")
        sub["rel_depth"] = pd.to_numeric(sub["rel_depth"], errors="coerce")
        sub = sub.dropna(subset=["baseline", "rel_depth"])
        if len(sub) < 2:
            continue

        n_samples = sub["sample"].nunique()
        path = plot_dir / f"{logical_chr}.rel_depth.vs_baseline.png"
        plot_correlation_with_regression(
            sub,
            x_col="baseline",
            y_col="rel_depth",
            title=f"{logical_chr} relative depth vs genome-wide mean depth (n={n_samples})",
            filename=str(path),
            x_label="Genome-wide mean depth (coverage)",
            y_label="Relative depth (chr depth / genome-wide mean)",
            ylim=ylim,
            group_col="Group",
        )
        written.append(path)
    return written


def plot_flagged_sample_chr_rel_depth_profiles(
    chr_all: pd.DataFrame,
    flagged_samples: pd.DataFrame | str | Path,
    output_path: str | Path,
    group_file: str | Path | None = None,
    info_path: str | Path | None = None,
    *,
    ylim: tuple[float, float] | None = (0.0, 3.0),
    sample_line_alpha: float = 0.22,
    include_subgenome_panels: bool = True,
) -> list[Path]:
    """
    Line plot of rel_depth across logical chromosomes for flagged samples only.

    Each flagged sample is one line; color follows ``Group`` (same palette as other
    aneuploidy plots). Writes the genome-wide panel plus optional A/B/D subgenome
    panels (chr1A–chr7A, chr1B–chr7B, chr1D–chr7D).
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if isinstance(flagged_samples, (str, Path)):
        flagged_df = pd.read_csv(flagged_samples, sep="\t")
    else:
        flagged_df = flagged_samples.copy()
    if flagged_df.empty:
        print("[Warning] No flagged samples; skipping chr rel_depth profile plot")
        return []

    flagged_ids = set(flagged_df["sample"].map(normalize_sample_id))

    plot_df = attach_sample_group(chr_all, group_file)
    plot_df = plot_df.copy()
    plot_df["sample_norm"] = plot_df["sample"].map(normalize_sample_id)
    sub = _evaluable_on_chromosome(plot_df)
    sub = sub[sub["sample_norm"].isin(flagged_ids)]
    if sub.empty:
        print("[Warning] No evaluable flagged chr rows; skipping profile plot")
        return []

    if info_path is not None:
        info_path = Path(info_path)
        info_path.parent.mkdir(parents=True, exist_ok=True)
        sub.sort_values(["Group", "sample_norm", "logical_chr"]).to_csv(
            info_path, sep="\t", index=False,
        )

    written: list[Path] = []
    panel_specs: list[tuple[tuple[str, ...], Path, str]] = [
        (
            LOGICAL_AUTOSOMES,
            output_path,
            f"Flagged samples: relative depth by chromosome (n={sub['sample_norm'].nunique()})",
        ),
    ]
    if include_subgenome_panels:
        stem = output_path.stem
        suffix = output_path.suffix
        parent = output_path.parent
        for sg in ("A", "B", "D"):
            chrs = SUBGENOME_AUTOSOMES[sg]
            sg_sub = sub[sub["logical_chr"].isin(chrs)]
            if sg_sub.empty:
                continue
            panel_specs.append(
                (
                    chrs,
                    parent / f"{stem}.sub{sg}{suffix}",
                    (
                        f"Flagged samples: {sg} subgenome relative depth "
                        f"(n={sg_sub['sample_norm'].nunique()})"
                    ),
                )
            )

    for logical_chrs, png_path, title in panel_specs:
        sub_panel = sub[sub["logical_chr"].isin(logical_chrs)]
        if sub_panel.empty:
            continue
        _save_flagged_chr_rel_depth_profile_figure(
            sub_panel,
            logical_chrs,
            png_path,
            title=title,
            ylim=ylim,
            sample_line_alpha=sample_line_alpha,
        )
        written.append(png_path)
    return written


def _save_flagged_chr_rel_depth_profile_figure(
    sub: pd.DataFrame,
    logical_chrs: tuple[str, ...],
    output_path: Path,
    *,
    title: str,
    ylim: tuple[float, float] | None,
    sample_line_alpha: float,
) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    hue_order, palette_map = sample_group_hue_config(sub["Group"])
    group_counts = sub.groupby("Group")["sample_norm"].nunique().to_dict()

    fig_w = max(8.0, 1.1 * len(logical_chrs))
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(fig_w, 6))

    for sample_id, sdf in sub.groupby("sample_norm", sort=False):
        group = str(sdf["Group"].iloc[0])
        by_chr = sdf.set_index("logical_chr")["rel_depth"]
        ys = [float(by_chr.get(c, np.nan)) for c in logical_chrs]
        ax.plot(
            list(logical_chrs),
            ys,
            color=palette_map.get(group, (0.5, 0.5, 0.5)),
            alpha=sample_line_alpha,
            linewidth=0.9,
            zorder=1,
        )

    for group in hue_order:
        if group not in group_counts:
            continue
        ax.plot(
            [],
            [],
            color=palette_map.get(group, (0.5, 0.5, 0.5)),
            linewidth=2.0,
            label=f"{group} (n={group_counts[group]})",
        )
    ax.axhline(1.0, color="#666666", linestyle="--", linewidth=1.2, label="disomic (1.0)", zorder=0)

    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.set_xlabel("Chromosome", fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel("Relative depth (vs cohort baseline)", fontsize=Y_LABEL_FONT_SIZE)
    ax.set_xticks(list(range(len(logical_chrs))))
    ax.set_xticklabels(list(logical_chrs), rotation=45, ha="right", fontsize=TICK_FONT_SIZE)
    ax.tick_params(axis="y", labelsize=TICK_FONT_SIZE)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, -0.22),
        ncol=min(5, len(hue_order) + 1),
        fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved to {output_path}")


def plot_all_per_chromosome_figures(
    chr_all: pd.DataFrame,
    plot_dir: str | Path,
    group_file: str | Path | None = None,
    *,
    bins: int = 80,
    xlim: tuple[float, float] | None = (0.0, 3.0),
) -> list[Path]:
    """Rel_depth distributions (linear + log-y), baseline regression scatter, and flagged group charts."""
    plot_dir = Path(plot_dir)
    paths: list[Path] = []
    paths.extend(
        plot_per_chromosome_rel_depth_distributions(
            chr_all, plot_dir, group_file=group_file, bins=bins, xlim=xlim,
        )
    )
    paths.extend(
        plot_per_chromosome_rel_depth_vs_baseline(
            chr_all, plot_dir, group_file=group_file, ylim=xlim,
        )
    )
    paths.extend(
        plot_per_chromosome_flagged_group_charts(chr_all, plot_dir, group_file=group_file)
    )
    return paths


def publish_aneuploidy_collect_outputs(
    chr_paths: list[str | Path],
    sample_paths: list[str | Path],
    out_dir: str | Path = ".",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[Path]]:
    """Merge per-sample tables; write summary and per-chromosome depth TSVs under ``out_dir``."""
    out_dir = Path(out_dir)
    chr_all, sample_all, flagged = aggregate_aneuploidy_tables(chr_paths, sample_paths)

    sample_all.to_csv(out_dir / "sample_aneuploidy_summary.tsv", sep="\t", index=False)
    flagged.to_csv(out_dir / "flagged_samples.tsv", sep="\t", index=False)
    chr_table_paths = write_per_chromosome_depth_tables(chr_all, out_dir / "chr_depth")
    return chr_all, sample_all, flagged, chr_table_paths


def load_chr_depth_from_dir(chr_depth_dir: str | Path) -> pd.DataFrame:
    """Load flat ``*.depth.tsv`` files from a single chr_depth directory."""
    chr_depth_dir = Path(chr_depth_dir)
    paths = sorted(chr_depth_dir.glob("*.depth.tsv"))
    if not paths:
        return pd.DataFrame(columns=list(CHR_DEPTH_COLUMNS))
    parts = [pd.read_csv(p, sep="\t") for p in paths]
    return pd.concat(parts, ignore_index=True)


def replot_rel_depth_distributions_from_chr_depth(
    chr_depth_dir: str | Path,
    plot_dir: str | Path,
    group_file: str | Path | None = None,
    flagged_samples_path: str | Path | None = None,
    *,
    bins: int = 80,
    xlim: tuple[float, float] | None = (0.0, 3.0),
) -> list[Path]:
    """Replot all per-chromosome figures from existing ``*.depth.tsv`` (no collect merge)."""
    chr_all = load_chr_depth_from_dir(chr_depth_dir)
    if chr_all.empty:
        raise FileNotFoundError(f"No *.depth.tsv under {chr_depth_dir}")
    plot_dir = Path(plot_dir)
    paths: list[Path] = list(
        plot_all_per_chromosome_figures(
            chr_all,
            plot_dir,
            group_file=group_file,
            bins=bins,
            xlim=xlim,
        )
    )
    if flagged_samples_path and Path(flagged_samples_path).exists():
        profile_png = plot_dir / "flagged_sample_group.rel_depth.chr_distribution.line.png"
        info_dir = plot_dir.parent / "info" if plot_dir.name == "plots" else plot_dir
        profile_info = info_dir / "flagged_sample_chr_rel_depth.profile.tsv"
        paths.extend(
            plot_flagged_sample_chr_rel_depth_profiles(
                chr_all,
                flagged_samples_path,
                profile_png,
                group_file=group_file,
                info_path=profile_info,
                ylim=xlim,
            )
        )
    return paths


def scan_depth_root(
    depth_root: str | Path,
    cohorts: tuple[str, ...] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    depth_root = Path(depth_root)
    cohorts = cohorts or tuple(COHORT_EXPECTED_SUBGENOMES.keys())
    chr_parts, sample_parts = [], []
    for cohort in cohorts:
        cohort_dir = depth_root / cohort
        if not cohort_dir.is_dir():
            continue
        c_df, s_df = scan_cohort_depth_dir(cohort_dir, cohort=cohort)
        if not c_df.empty:
            chr_parts.append(c_df)
        if not s_df.empty:
            sample_parts.append(s_df)
    chr_all = pd.concat(chr_parts, ignore_index=True) if chr_parts else pd.DataFrame()
    sample_all = pd.concat(sample_parts, ignore_index=True) if sample_parts else pd.DataFrame()
    return chr_all, sample_all
