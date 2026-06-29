"""Plot mosdepth summary mean depth and per-chromosome depth density for BAM benchmarks."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from infra.utils.graph import (
    LEGEND_FONT_SIZE,
    TICK_FONT_SIZE,
    TITLE_FONT_SIZE,
    X_LABEL_FONT_SIZE,
    Y_LABEL_FONT_SIZE,
)
from infra.utils.io import load_df_generic
from infra.wheat import ref_v1

DEFAULT_INPUT_DIR = Path("/data/home/tusr1/01projects/vmap4/00data/04depth/03ABD")
DEFAULT_OUTPUT_DIR = Path(
    "/data1/dazheng_tusr1/vmap4.VCF.v1/benchmark/popdep_bench/plots/dep_dist"
)

CS_BAM_SAMPLES = (
    "CS_mp_2018_8X",
    "CS_sg_2014_3X",
    "CS_sg_2017_60X",
)

ABD_DEPTH_TARGETS = (1, 3, 5, 8, 10)

DEFAULT_PALETTE = (
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
)

CS_SAMPLE_COLORS = {
    "CS_mp_2018_8X": "#1f77b4",
    "CS_sg_2014_3X": "#ff7f0e",
    "CS_sg_2017_60X": "#2ca02c",
}

SEGMENT_CHR_IDS = tuple(i for i in ref_v1.get_chr_ids() if 1 <= i <= 42)


def _segment_label(chr_id: int) -> str:
    chrom = ref_v1.get_chromosome(chr_id)
    part = 1 if chr_id % 2 == 1 else 2
    return f"{chrom}-{part}"


def _summary_path(input_dir: Path, sample: str) -> Path:
    return input_dir / f"{sample}.bam.mosdepth.summary.txt"


def _dist_path(input_dir: Path, sample: str) -> Path:
    return input_dir / f"{sample}.bam.mosdepth.global.dist.txt"


def _sample_colors(samples: tuple[str, ...]) -> dict[str, str]:
    return {sample: DEFAULT_PALETTE[i % len(DEFAULT_PALETTE)] for i, sample in enumerate(samples)}


def load_abd_total_mean_depths(input_dir: Path) -> list[tuple[str, float]]:
    rows: list[tuple[str, float]] = []
    for path in sorted(input_dir.glob("ABD_*.bam.mosdepth.summary.txt")):
        df = load_df_generic(str(path))
        total = df.loc[df["chrom"].astype(str) == "total"]
        if total.empty:
            continue
        taxa = path.name.replace(".bam.mosdepth.summary.txt", "")
        rows.append((taxa, float(total.iloc[0]["mean"])))
    rows.sort(key=lambda item: item[1])
    return rows


def select_abd_taxa_by_target_depth(
    input_dir: Path,
    targets: tuple[int, ...] = ABD_DEPTH_TARGETS,
) -> tuple[tuple[str, ...], pd.DataFrame]:
    """Pick one ABD taxon per target genome-wide mean depth (greedy, no reuse)."""
    depths = load_abd_total_mean_depths(input_dir)
    if not depths:
        raise FileNotFoundError(f"No ABD_*.bam.mosdepth.summary.txt under {input_dir}")

    used: set[str] = set()
    picks: list[dict[str, object]] = []
    for target in targets:
        candidates = [(taxa, mean) for taxa, mean in depths if taxa not in used]
        taxa, mean = min(candidates, key=lambda item: abs(item[1] - target))
        used.add(taxa)
        picks.append(
            {
                "target_depth_x": target,
                "taxa": taxa,
                "genome_mean_depth": round(mean, 3),
                "abs_error_x": round(abs(mean - target), 3),
            }
        )

    selection = pd.DataFrame(picks)
    samples = tuple(selection["taxa"].astype(str))
    return samples, selection


def load_mosdepth_summary(input_dir: Path, sample: str) -> pd.DataFrame:
    path = _summary_path(input_dir, sample)
    df = load_df_generic(str(path))
    df["chrom"] = df["chrom"].astype(str)
    df = df[df["chrom"].str.isdigit()].copy()
    df["chr_id"] = df["chrom"].astype(int)
    df = df[df["chr_id"].isin(SEGMENT_CHR_IDS)].copy()
    df["segment"] = df["chr_id"].map(_segment_label)
    df["sample"] = sample
    return df.sort_values("chr_id").reset_index(drop=True)


def load_all_summaries(input_dir: Path, samples: tuple[str, ...]) -> pd.DataFrame:
    frames = [load_mosdepth_summary(input_dir, sample) for sample in samples]
    return pd.concat(frames, ignore_index=True)


def build_mean_depth_wide(summary: pd.DataFrame, samples: tuple[str, ...]) -> pd.DataFrame:
    order = [_segment_label(cid) for cid in SEGMENT_CHR_IDS]
    wide = summary.pivot(index="segment", columns="sample", values="mean")
    wide = wide.reindex(order)
    wide = wide.reset_index().rename(columns={"segment": "segment_label"})
    wide["segment_idx"] = np.arange(len(wide))
    return wide


def plot_mean_depth_by_segment(
    wide: pd.DataFrame,
    samples: tuple[str, ...],
    output_dir: Path,
    *,
    sample_colors: dict[str, str],
    title: str,
    stem: str,
) -> Path:
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(16, 6))
    x = np.arange(len(wide))
    ax.set_xticks(x)
    ax.set_xticklabels(wide["segment_label"], rotation=90, fontsize=8)
    for sample in samples:
        ax.plot(
            x,
            wide[sample],
            color=sample_colors[sample],
            linewidth=1.8,
            label=sample,
        )
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.set_xlabel("Genome segment", fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel("Mean depth", fontsize=Y_LABEL_FONT_SIZE)
    ax.tick_params(axis="y", labelsize=TICK_FONT_SIZE)
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, -0.22),
        fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )
    out = output_dir / f"{stem}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot to {out}")
    return out


def _cumulative_to_density(rows: list[tuple[int, float, float]]) -> tuple[np.ndarray, np.ndarray]:
    """Convert mosdepth global.dist cumulative rows to depth / density arrays."""
    if not rows:
        return np.array([]), np.array([])

    rows_sorted = sorted(rows, key=lambda r: r[1], reverse=True)
    depths: list[float] = []
    densities: list[float] = []
    prev_cum = 0.0
    for _, depth, cum in rows_sorted:
        density = float(cum) - prev_cum
        if density > 0:
            depths.append(float(depth))
            densities.append(density)
        prev_cum = float(cum)
    return np.asarray(depths), np.asarray(densities)


def load_mosdepth_density(
    input_dir: Path,
    sample: str,
    chr_ids: tuple[int, ...] = SEGMENT_CHR_IDS,
) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    path = _dist_path(input_dir, sample)
    wanted = set(chr_ids)
    by_chr: dict[int, list[tuple[int, float, float]]] = {cid: [] for cid in chr_ids}

    with open(path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3 or not parts[0].isdigit():
                continue
            cid = int(parts[0])
            if cid not in wanted:
                continue
            by_chr[cid].append((cid, float(parts[1]), float(parts[2])))

    return {cid: _cumulative_to_density(by_chr[cid]) for cid in chr_ids}


def _depth_xlim(
    depth: np.ndarray,
    density: np.ndarray,
    mean_depth: float | None,
    max_density_frac: float = 0.999,
) -> tuple[float, float]:
    if depth.size == 0:
        return 0.0, 1.0
    cum = np.cumsum(density)
    if cum[-1] <= 0:
        return 0.0, float(depth.max())
    norm = cum / cum[-1]
    keep = norm <= max_density_frac
    hi = float(depth[keep].max()) if keep.any() else float(depth.max())
    if mean_depth is not None and mean_depth > 0:
        hi = max(hi, mean_depth * 3.0)
    return 0.0, hi


def plot_chr_depth_density(
    densities: dict[str, dict[int, tuple[np.ndarray, np.ndarray]]],
    summary: pd.DataFrame,
    samples: tuple[str, ...],
    output_dir: Path,
    *,
    sample_colors: dict[str, str],
    chr_subdir: str = "by_chromosome",
) -> list[Path]:
    chr_dir = output_dir / chr_subdir
    chr_dir.mkdir(parents=True, exist_ok=True)
    mean_lookup = summary.set_index(["sample", "chr_id"])["mean"].to_dict()
    written: list[Path] = []

    for chr_id in SEGMENT_CHR_IDS:
        chrom = ref_v1.get_chromosome(chr_id)
        part = 1 if chr_id % 2 == 1 else 2
        title = f"Depth density — {chrom} segment {part}"
        fname = chr_dir / f"{chrom}_seg{part}_depth_density.png"

        sns.set_style("white")
        fig, ax = plt.subplots(figsize=(10, 6))
        x_max = 0.0
        for sample in samples:
            depth, density = densities[sample][chr_id]
            if depth.size == 0:
                continue
            mean_depth = mean_lookup.get((sample, chr_id))
            _, hi = _depth_xlim(depth, density, mean_depth)
            x_max = max(x_max, hi)
            ax.plot(
                depth,
                density,
                color=sample_colors[sample],
                linewidth=1.5,
                label=sample,
            )

        ax.set_xlim(0, x_max if x_max > 0 else 1.0)
        ax.set_ylim(0, None)
        ax.set_title(title, fontsize=TITLE_FONT_SIZE)
        ax.set_xlabel("Depth", fontsize=X_LABEL_FONT_SIZE)
        ax.set_ylabel("Fraction of bases", fontsize=Y_LABEL_FONT_SIZE)
        ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
        ax.legend(
            loc="upper right",
            fontsize=LEGEND_FONT_SIZE,
            frameon=False,
        )
        fig.savefig(fname, dpi=300, bbox_inches="tight")
        plt.close(fig)
        written.append(fname)

    print(f"Saved {len(written)} per-segment depth density plots under {chr_dir}")
    return written


def plot_sample_density_panels(
    densities: dict[str, dict[int, tuple[np.ndarray, np.ndarray]]],
    summary: pd.DataFrame,
    samples: tuple[str, ...],
    output_dir: Path,
    *,
    sample_colors: dict[str, str],
) -> list[Path]:
    """One faceted figure per sample (6×7 grid of segments)."""
    mean_lookup = summary.set_index(["sample", "chr_id"])["mean"].to_dict()
    written: list[Path] = []
    n_cols = 7
    n_rows = 6

    for sample in samples:
        sns.set_style("white")
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(24, 18), sharex=False, sharey=False)
        fig.suptitle(
            f"Per-segment depth density — {sample}",
            fontsize=TITLE_FONT_SIZE + 2,
            y=0.995,
        )

        for idx, chr_id in enumerate(SEGMENT_CHR_IDS):
            row, col = divmod(idx, n_cols)
            ax = axes[row, col]
            depth, density = densities[sample][chr_id]
            label = _segment_label(chr_id)
            if depth.size:
                mean_depth = mean_lookup.get((sample, chr_id))
                _, hi = _depth_xlim(depth, density, mean_depth)
                mask = depth <= hi
                ax.plot(
                    depth[mask],
                    density[mask],
                    color=sample_colors[sample],
                    linewidth=0.9,
                )
                ax.set_xlim(0, hi)
            ax.set_title(label, fontsize=9)
            ax.tick_params(axis="both", labelsize=7)

        for ax in axes.flat[len(SEGMENT_CHR_IDS) :]:
            ax.set_visible(False)

        for ax in axes[-1, :]:
            ax.set_xlabel("Depth", fontsize=8)
        for ax in axes[:, 0]:
            ax.set_ylabel("Fraction", fontsize=8)

        out = output_dir / f"{sample}_depth_density_panels.png"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        written.append(out)
        print(f"Saved plot to {out}")

    return written


def run_mosdepth_bam_depth_plots(
    input_dir: str | Path = DEFAULT_INPUT_DIR,
    output_dir: str | Path = DEFAULT_OUTPUT_DIR,
    samples: tuple[str, ...] = CS_BAM_SAMPLES,
    *,
    sample_colors: dict[str, str] | None = None,
    title: str = "Mean sequencing depth by wheat segment (CS reference BAMs)",
    stem: str = "cs_mean_depth_by_segment",
    chr_subdir: str = "by_chromosome",
) -> pd.DataFrame:
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    colors = sample_colors or _sample_colors(samples)

    summary = load_all_summaries(input_dir, samples)
    wide = build_mean_depth_wide(summary, samples)
    plot_mean_depth_by_segment(
        wide,
        samples,
        output_dir,
        sample_colors=colors,
        title=title,
        stem=stem,
    )

    densities = {sample: load_mosdepth_density(input_dir, sample) for sample in samples}
    plot_chr_depth_density(
        densities,
        summary,
        samples,
        output_dir,
        sample_colors=colors,
        chr_subdir=chr_subdir,
    )
    plot_sample_density_panels(
        densities,
        summary,
        samples,
        output_dir,
        sample_colors=colors,
    )

    summary_out = output_dir / f"{stem}.tsv"
    wide.to_csv(summary_out, sep="\t", index=False)
    print(f"Wrote summary table to {summary_out}")
    return wide


def run_abd_depth_tier_plots(
    input_dir: str | Path = DEFAULT_INPUT_DIR,
    output_dir: str | Path = DEFAULT_OUTPUT_DIR / "abd_tiers",
    targets: tuple[int, ...] = ABD_DEPTH_TARGETS,
) -> pd.DataFrame:
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    samples, selection = select_abd_taxa_by_target_depth(input_dir, targets)
    selection.to_csv(output_dir / "abd_depth_tier_selection.tsv", sep="\t", index=False)
    print("Selected ABD taxa by target genome mean depth:")
    for _, row in selection.iterrows():
        print(
            f"  ~{int(row['target_depth_x'])}x -> {row['taxa']} "
            f"(genome mean {row['genome_mean_depth']}x)"
        )

    run_mosdepth_bam_depth_plots(
        input_dir,
        output_dir,
        samples,
        title="Mean sequencing depth by wheat segment (ABD depth tiers)",
        stem="abd_mean_depth_by_segment",
        chr_subdir="by_chromosome",
    )
    return selection


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot mosdepth mean depth and per-chromosome depth density for BAMs",
    )
    parser.add_argument(
        "--input-dir",
        default=str(DEFAULT_INPUT_DIR),
        help="Directory containing *.bam.mosdepth.{summary,global.dist}.txt",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory for PNG outputs",
    )
    parser.add_argument(
        "--cohort",
        choices=("cs", "abd_tiers", "both"),
        default="cs",
        help="Sample set to plot (default: cs reference BAMs)",
    )
    args = parser.parse_args()

    if args.cohort in ("cs", "both"):
        run_mosdepth_bam_depth_plots(
            args.input_dir,
            args.output_dir,
            CS_BAM_SAMPLES,
            sample_colors=CS_SAMPLE_COLORS,
        )
    if args.cohort in ("abd_tiers", "both"):
        run_abd_depth_tier_plots(args.input_dir, Path(args.output_dir) / "abd_tiers")


if __name__ == "__main__":
    main()
