"""Per-bin MAPQ and depth profiles from BAM alignments on vmap4 ref v1."""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from genetics.germplasm.sample.mosdepth_aneuploidy import LOGICAL_AUTOSOMES
from infra.utils.io import load_df_generic
from infra.wheat import ref_v1

SAMTOOLS_VIEW_FLAGS = 1796  # unmapped, secondary, QC-fail, duplicate
MQ_BIN_LABELS = (
    "mq_0_10",
    "mq_10_20",
    "mq_20_30",
    "mq_30_40",
    "mq_40_50",
    "mq_50_60",
)
DEFAULT_BIN_SIZE_BP = 1_000_000
TSV_COLUMNS = (
    "logical_chr",
    "bin_idx",
    "bin_start",
    "bin_end",
    *MQ_BIN_LABELS,
    "total_reads",
    "bin_mean_depth",
)


@dataclass
class _BinAccumulator:
    mq_counts: np.ndarray = field(default_factory=lambda: np.zeros(6, dtype=np.int64))
    ref_bases: float = 0.0
    read_count: int = 0


def _cigar_ref_len(cigar: str) -> int:
    total = 0
    for size, op in re.findall(r"(\d+)([=XMDN])", cigar):
        if op in "=XMDN":
            total += int(size)
    return total


def _mq_bin_index(mapq: int) -> int:
    if mapq < 0:
        return -1
    if mapq >= 60:
        return 5
    return min(mapq // 10, 5)


def _logical_position(plink_chr: int, pos: int) -> tuple[str, int] | None:
    meta = ref_v1._refSegmentIntervalMap.get(int(plink_chr))
    if meta is None:
        return None
    logical = meta["ref_name"]
    if logical not in LOGICAL_AUTOSOMES:
        return None
    logical_pos = int(meta["ref_start"]) + int(pos)
    return logical, logical_pos


def logical_chromosome_lengths() -> dict[str, int]:
    lengths: dict[str, int] = {}
    for logical in LOGICAL_AUTOSOMES:
        chrom = logical.replace("chr", "", 1)
        intervals = ref_v1._ref_segment_intervals_for_chromosome(chrom)
        lengths[logical] = max(meta["ref_end"] for _, meta in intervals)
    return lengths


def _init_bin_store(
    bin_size_bp: int,
) -> tuple[dict[tuple[str, int], _BinAccumulator], dict[str, int]]:
    chr_lengths = logical_chromosome_lengths()
    store: dict[tuple[str, int], _BinAccumulator] = {}
    for logical, length in chr_lengths.items():
        n_bins = int(np.ceil(length / bin_size_bp))
        for bin_idx in range(n_bins):
            store[(logical, bin_idx)] = _BinAccumulator()
    return store, chr_lengths


def scan_bam_mapq_depth_bins(
    bam_path: str | Path,
    bin_size_bp: int = DEFAULT_BIN_SIZE_BP,
    samtools: str = "samtools",
    threads: int = 4,
    progress_every: int = 5_000_000,
) -> pd.DataFrame:
    """Scan primary alignments and aggregate MAPQ/read depth into fixed-width bins."""
    bam_path = Path(bam_path)
    if not bam_path.exists():
        raise FileNotFoundError(bam_path)

    store, chr_lengths = _init_bin_store(bin_size_bp)
    cmd = [
        samtools,
        "view",
        "-@",
        str(int(threads)),
        "-F",
        str(SAMTOOLS_VIEW_FLAGS),
        str(bam_path),
    ]
    print(f"[Info] scan_bam_mapq_depth_bins: {' '.join(cmd)}")

    n_align = 0
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:
            continue
        rname, pos_str, mapq_str, cigar = parts[2], parts[3], parts[4], parts[5]
        mapped = _logical_position(int(rname), int(pos_str))
        if mapped is None:
            continue
        logical, logical_pos = mapped
        if logical_pos < 1:
            continue
        bin_idx = (logical_pos - 1) // int(bin_size_bp)
        key = (logical, int(bin_idx))
        if key not in store:
            continue

        mapq = int(mapq_str)
        mq_idx = _mq_bin_index(mapq)
        if mq_idx < 0:
            continue
        ref_len = max(_cigar_ref_len(cigar), 1)

        acc = store[key]
        acc.mq_counts[mq_idx] += 1
        acc.read_count += 1
        acc.ref_bases += ref_len
        n_align += 1
        if progress_every and n_align % progress_every == 0:
            print(f"[Info] processed {n_align:,} alignments...")

    stderr = proc.stderr.read() if proc.stderr is not None else ""
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"samtools view failed ({rc}): {stderr.strip()}")
    print(f"[Info] finished BAM pass: {n_align:,} alignments")

    rows = []
    for logical in LOGICAL_AUTOSOMES:
        length = chr_lengths[logical]
        n_bins = int(np.ceil(length / bin_size_bp))
        for bin_idx in range(n_bins):
            acc = store[(logical, bin_idx)]
            bin_start = bin_idx * bin_size_bp + 1
            bin_end = min((bin_idx + 1) * bin_size_bp, length)
            span = max(bin_end - bin_start + 1, 1)
            rows.append(
                {
                    "logical_chr": logical,
                    "bin_idx": bin_idx,
                    "bin_start": bin_start,
                    "bin_end": bin_end,
                    **{
                        label: int(acc.mq_counts[i])
                        for i, label in enumerate(MQ_BIN_LABELS)
                    },
                    "total_reads": int(acc.read_count),
                    "bin_mean_depth": acc.ref_bases / span,
                }
            )
    return pd.DataFrame(rows, columns=list(TSV_COLUMNS))


def write_mapq_depth_bin_tsv(df: pd.DataFrame, output_tsv: str | Path) -> Path:
    out = Path(output_tsv)
    df.to_csv(out, sep="\t", index=False)
    return out


def run_sample_mapq_depth_bins(
    sample: str,
    bam_path: str | Path,
    output_tsv: str | Path,
    plot_dir: str | Path | None = None,
    bin_size_bp: int = DEFAULT_BIN_SIZE_BP,
    samtools: str = "samtools",
    threads: int = 4,
    make_plots: bool = False,
) -> tuple[Path, list[Path]]:
    """Scan one BAM, write TSV, optionally emit QC plots."""
    df = scan_bam_mapq_depth_bins(
        bam_path,
        bin_size_bp=bin_size_bp,
        samtools=samtools,
        threads=threads,
    )
    tsv_path = write_mapq_depth_bin_tsv(df, output_tsv)
    plot_paths: list[Path] = []
    if make_plots and plot_dir is not None:
        plot_paths = plot_mapq_depth_profiles(
            df,
            sample=sample,
            plot_dir=plot_dir,
            bin_size_bp=bin_size_bp,
        )
    return tsv_path, plot_paths


def run_flagged_sample_mapq_depth_check(
    sample: str,
    bam_path: str | Path,
    output_tsv: str | Path,
    plot_dir: str | Path = ".",
    bin_size_bp: int = DEFAULT_BIN_SIZE_BP,
    samtools: str = "samtools",
    threads: int = 4,
    make_plots: bool = True,
) -> tuple[Path, list[Path]]:
    """Backward-compatible wrapper used by earlier mapq_check runs."""
    return run_sample_mapq_depth_bins(
        sample=sample,
        bam_path=bam_path,
        output_tsv=output_tsv,
        plot_dir=plot_dir,
        bin_size_bp=bin_size_bp,
        samtools=samtools,
        threads=threads,
        make_plots=make_plots,
    )


def _subgenome_frame(df: pd.DataFrame, subgenome: str) -> pd.DataFrame:
    return df[df["logical_chr"].str.endswith(subgenome)].copy()


def plot_mapq_depth_profiles(
    df: pd.DataFrame,
    sample: str,
    plot_dir: str | Path,
    bin_size_bp: int = DEFAULT_BIN_SIZE_BP,
) -> list[Path]:
    """Emit genome/subgenome depth, MAPQ fraction, capped depth, and 5 Mb heatmap."""
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    def _save(fig: plt.Figure, name: str) -> None:
        path = plot_dir / f"{sample}.{name}.png"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        written.append(path)
        print(f"[Info] wrote {path}")

    cap_path = plot_mapq_depth_cap_genome_grid(
        df,
        sample=sample,
        output_png=plot_dir / f"{sample}.mapq_depth_cap.genome.png",
        bin_size_bp=bin_size_bp,
    )
    written.append(cap_path)

    for tag, sub_df in [
        ("subA", _subgenome_frame(df, "A")),
        ("subB", _subgenome_frame(df, "B")),
        ("subD", _subgenome_frame(df, "D")),
    ]:
        if sub_df.empty:
            continue
        sub_df = sub_df.sort_values(["logical_chr", "bin_idx"]).reset_index(drop=True)
        x = np.arange(len(sub_df))
        depth = sub_df["bin_mean_depth"].astype(float)

        fig, ax = plt.subplots(figsize=(12, 3))
        ax.plot(x, depth, linewidth=0.6, color="#1f77b4")
        ax.set_title(f"{sample} bin mean depth ({tag})")
        ax.set_ylabel("Depth")
        ax.set_xlabel("Bin index")
        _save(fig, f"mapq_depth.{tag}")

        cap = np.clip(depth, 0, np.nanpercentile(depth[np.isfinite(depth)], 99))
        fig, ax = plt.subplots(figsize=(12, 3))
        ax.plot(x, cap, linewidth=0.6, color="#d62728")
        ax.set_title(f"{sample} bin mean depth capped ({tag})")
        ax.set_ylabel("Depth")
        ax.set_xlabel("Bin index")
        _save(fig, f"mapq_depth_cap.{tag}")

        mq = sub_df[list(MQ_BIN_LABELS)].astype(float)
        frac = mq.div(mq.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
        fig, ax = plt.subplots(figsize=(12, 3))
        ax.stackplot(
            x,
            [frac[col] for col in MQ_BIN_LABELS],
            labels=MQ_BIN_LABELS,
            alpha=0.85,
        )
        ax.set_ylim(0, 1)
        ax.set_title(f"{sample} MAPQ fraction ({tag})")
        ax.set_xlabel("Bin index")
        ax.legend(loc="upper right", fontsize=6, ncol=3)
        _save(fig, f"mapq_frac.{tag}")

    heatmap_bin = max(int(bin_size_bp), 5_000_000)
    if heatmap_bin != bin_size_bp:
        agg = (
            df.assign(
                heat_bin=(df["bin_idx"] * bin_size_bp) // heatmap_bin,
            )
            .groupby(["logical_chr", "heat_bin"], as_index=False)["bin_mean_depth"]
            .mean()
        )
    else:
        agg = df[["logical_chr", "bin_idx", "bin_mean_depth"]].rename(
            columns={"bin_idx": "heat_bin"},
        )
    pivot = agg.pivot(index="logical_chr", columns="heat_bin", values="bin_mean_depth")
    pivot = pivot.reindex(list(LOGICAL_AUTOSOMES))
    path = plot_dir / f"{sample}.mapq_depth.heatmap5mb.png"
    fig, ax = plt.subplots(figsize=(14, 6))
    data = pivot.to_numpy(dtype=float)
    im = ax.imshow(data, aspect="auto", interpolation="nearest", cmap="YlOrRd")
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(list(pivot.index), fontsize=6)
    ax.set_xlabel("5 Mb bin index")
    ax.set_title(f"{sample} bin mean depth (5 Mb)")
    fig.colorbar(im, ax=ax, label="Depth")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    written.append(path)
    print(f"[Info] wrote {path}")
    return written


MQ_STACK_COLORS = (
    "#8b0000",
    "#ff7f0e",
    "#ffd700",
    "#cfe2f3",
    "#6fa8dc",
    "#1f4e79",
)


def plot_mapq_depth_cap_genome_grid(
    df: pd.DataFrame,
    sample: str,
    output_png: str | Path,
    *,
    baseline_depth: float | None = None,
    y_cap_multiplier: float = 3.0,
    bin_size_bp: int = DEFAULT_BIN_SIZE_BP,
    panel_hook=None,
) -> Path:
    """7x3 genome grid of MAPQ-decomposed capped depth (mapq_check style)."""
    from infra.utils.graph import LEGEND_FONT_SIZE, TICK_FONT_SIZE, TITLE_FONT_SIZE, Y_LABEL_FONT_SIZE

    output_png = Path(output_png)
    output_png.parent.mkdir(parents=True, exist_ok=True)

    work = df.sort_values(["logical_chr", "bin_idx"]).reset_index(drop=True)
    finite_depth = pd.to_numeric(work["bin_mean_depth"], errors="coerce")
    finite_depth = finite_depth[np.isfinite(finite_depth)]
    if baseline_depth is not None and np.isfinite(baseline_depth) and baseline_depth > 0:
        y_cap = float(baseline_depth) * float(y_cap_multiplier)
        mean_cov = float(baseline_depth)
    else:
        y_cap = float(np.nanpercentile(finite_depth, 99)) if not finite_depth.empty else 1.0
        mean_cov = float(finite_depth.median()) if not finite_depth.empty else float("nan")

    fig, axes = plt.subplots(7, 3, figsize=(14, 16), sharey=True)
    fig.suptitle(
        f"{sample} genome (A/B/D) (y≤{y_cap:.1f})",
        fontsize=TITLE_FONT_SIZE,
        y=0.995,
    )

    for row_idx in range(7):
        for col_idx, subgenome in enumerate(("A", "B", "D")):
            logical_chr = f"chr{row_idx + 1}{subgenome}"
            ax = axes[row_idx, col_idx]
            chr_df = work[work["logical_chr"] == logical_chr].copy()
            if chr_df.empty:
                ax.set_xticks([])
                ax.set_yticks([])
                ax.text(0.5, 0.5, logical_chr, transform=ax.transAxes, ha="center", va="center")
                continue

            x_mb = chr_df["bin_start"].astype(float) / 1_000_000
            depth = pd.to_numeric(chr_df["bin_mean_depth"], errors="coerce").fillna(0.0)
            depth_cap = np.clip(depth, 0, y_cap)
            totals = pd.to_numeric(chr_df["total_reads"], errors="coerce").replace(0, np.nan)
            mq_depths = []
            for label in MQ_BIN_LABELS:
                mq_counts = pd.to_numeric(chr_df[label], errors="coerce").fillna(0.0)
                mq_depths.append((depth_cap * (mq_counts / totals)).fillna(0.0).to_numpy())

            ax.stackplot(
                x_mb,
                mq_depths,
                colors=MQ_STACK_COLORS,
                labels=list(MQ_BIN_LABELS),
                linewidth=0,
            )
            ax.set_ylim(0, y_cap)
            ax.set_xlabel("Position (Mb)", fontsize=TICK_FONT_SIZE - 1)
            if col_idx == 0:
                ax.set_ylabel("Depth", fontsize=Y_LABEL_FONT_SIZE)
            ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE - 2)
            ax.set_title(logical_chr, fontsize=TICK_FONT_SIZE, color="#222222")
            if panel_hook is not None:
                panel_hook(logical_chr, ax, chr_df, y_cap)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="lower center",
            ncol=3,
            fontsize=LEGEND_FONT_SIZE - 1,
            frameon=False,
            bbox_to_anchor=(0.5, 0.01),
        )
    if np.isfinite(mean_cov):
        fig.text(0.98, 0.01, f"Mean coverage: {mean_cov:.2f}x", ha="right", fontsize=LEGEND_FONT_SIZE - 1)

    fig.tight_layout(rect=(0, 0.04, 1, 0.98))
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[Info] wrote {output_png}")
    return output_png


def replot_mapq_depth_from_tsv(
    tsv_path: str | Path,
    sample: str | None = None,
    plot_dir: str | Path = ".",
    bin_size_bp: int = DEFAULT_BIN_SIZE_BP,
) -> list[Path]:
    """Rebuild PNGs from a published mapq_depth_bin TSV."""
    tsv_path = Path(tsv_path)
    sample = sample or tsv_path.name.replace(".mapq_depth_bin.tsv", "")
    df = load_df_generic(str(tsv_path))
    if df is None or df.empty:
        raise ValueError(f"Empty or missing TSV: {tsv_path}")
    return plot_mapq_depth_profiles(
        df,
        sample=sample,
        plot_dir=plot_dir,
        bin_size_bp=bin_size_bp,
    )


def write_taxa_bam_index_tsv(taxa_bam_dir: str | Path, output_tsv: str | Path) -> Path:
    """Build a flat sample→BAM table from taxaBamMap files (legacy helper)."""
    taxa_bam_dir = Path(taxa_bam_dir)
    rows: list[dict[str, str]] = []
    for path in sorted(taxa_bam_dir.glob("*.taxaBamMap.txt")):
        if path.name.startswith("all.") or path.name.startswith("final."):
            continue
        for line in path.read_text().splitlines():
            if not line.strip() or line.startswith("Taxa"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            rows.append({"sample": parts[0], "bam_path": parts[2]})
    out = Path(output_tsv)
    pd.DataFrame(rows).drop_duplicates(subset=["sample"]).to_csv(out, sep="\t", index=False)
    return out
