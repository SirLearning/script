"""Population depth reference from cohort ``*.mapq_depth_bin.tsv`` files (1 Mb bins)."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import polars as pl

from breachr.config import CallConfig
from breachr.constants import LOGICAL_AUTOSOMES, SUBGENOME_AUTOSOMES
from breachr.io import load_mapq_depth_bin_tsv

# MAPQ >= 30 uses the top three histogram bins (30–39, 40–49, 50–60).
MQ_GE30_LABELS = ("mq_30_40", "mq_40_50", "mq_50_60")

POPDEP_BIN_COLUMNS = (
    "logical_chr",
    "bin_idx",
    "bin_start",
    "bin_end",
    "n_samples",
    "depth_mean",
    "depth_sd",
    "rel_mean",
    "rel_sd",
    "n_cnv",
    "cnv_frac",
)

SAMPLE_REL_COLUMNS = (
    "sample",
    "logical_chr",
    "bin_idx",
    "bin_start",
    "bin_end",
    "bin_depth",
    "depth_rel",
    "is_cnv",
)

BuildEngine = Literal["polars", "streaming", "pandas"]
DepthMode = Literal["mq30", "all"]


def normalize_depth_mode(depth_mode: str | DepthMode) -> DepthMode:
    mode = str(depth_mode).strip().lower()
    if mode in {"all", "none", "unfiltered", "raw", "bin_mean"}:
        return "all"
    if mode in {"mq30", "mq30+", "ge30"}:
        return "mq30"
    raise ValueError(f"Unsupported depth_mode: {depth_mode!r} (use 'all' or 'mq30')")


def sample_bin_depth(df: pd.DataFrame, *, depth_mode: DepthMode = "mq30") -> pd.Series:
    """
    Per-bin absolute depth for popdep aggregation.

    ``all``: published ``bin_mean_depth`` (primary alignments, no MAPQ filter).
    ``mq30``: approximate MAPQ>=30 depth from the MAPQ histogram (legacy).
    """
    mode = normalize_depth_mode(depth_mode)
    depth = pd.to_numeric(df["bin_mean_depth"], errors="coerce").astype(float)
    if mode == "all":
        return depth.replace([np.inf, -np.inf], np.nan)
    return estimate_mq30_bin_depth(df)


def flag_bin_cnv_from_depth_rel(
    depth_rel: pd.Series | np.ndarray,
    *,
    config: CallConfig | None = None,
) -> np.ndarray:
    """True when relative depth indicates deletion or duplication (breachr bounds)."""
    cfg = config or CallConfig()
    rel = pd.to_numeric(depth_rel, errors="coerce").to_numpy(dtype=float)
    return np.isfinite(rel) & (
        (rel <= cfg.deletion_depth_rel_max) | (rel >= cfg.duplication_depth_rel_min)
    )


def estimate_mq30_bin_depth(df: pd.DataFrame) -> pd.Series:
    """
    Estimate per-bin depth from reads with MAPQ >= 30.

    Published TSVs store ``bin_mean_depth`` over all primary alignments and a
    MAPQ histogram only (not depth split by MAPQ). We approximate MQ30+ depth as:

        depth_mq30 = bin_mean_depth * (reads_mq30+ / total_reads)

    when ``total_reads > 0``.
    """
    reads = pd.to_numeric(df["total_reads"], errors="coerce").astype(float)
    depth = pd.to_numeric(df["bin_mean_depth"], errors="coerce").astype(float)
    mq30_reads = df[list(MQ_GE30_LABELS)].sum(axis=1).astype(float)
    frac = mq30_reads / reads.replace(0, np.nan)
    out = depth * frac
    return out.replace([np.inf, -np.inf], np.nan)


def _subgenome_from_logical_chr(logical_chr: str) -> str:
    return logical_chr[-1] if len(logical_chr) >= 2 else ""


def sample_subgenome_baselines(
    df: pd.DataFrame,
    *,
    min_reads: int = 1,
    depth_mode: DepthMode = "mq30",
) -> dict[str, float]:
    """Median bin depth per subgenome for one sample (positive bins only)."""
    work = df.copy()
    work["bin_depth"] = sample_bin_depth(work, depth_mode=depth_mode)
    work["subgenome"] = work["logical_chr"].map(_subgenome_from_logical_chr)
    reads = pd.to_numeric(work["total_reads"], errors="coerce").fillna(0)
    supported = work.loc[reads >= min_reads].copy()

    baselines: dict[str, float] = {}
    for sg, chroms in SUBGENOME_AUTOSOMES.items():
        sub = supported[supported["logical_chr"].isin(chroms)]
        depths = pd.to_numeric(sub["bin_depth"], errors="coerce")
        depths = depths[(depths > 0) & np.isfinite(depths)]
        baselines[sg] = float(depths.median()) if not depths.empty else float("nan")
    return baselines


def sample_subgenome_mq30_baselines(
    df: pd.DataFrame,
    *,
    min_reads: int = 1,
) -> dict[str, float]:
    """Backward-compatible alias for MQ30 subgenome baselines."""
    return sample_subgenome_baselines(df, min_reads=min_reads, depth_mode="mq30")


def add_sample_relative_depth(
    df: pd.DataFrame,
    baselines: dict[str, float] | None = None,
    *,
    min_reads: int = 1,
    config: CallConfig | None = None,
    depth_mode: DepthMode = "mq30",
) -> pd.DataFrame:
    """Add ``bin_depth``, ``depth_rel``, and ``is_cnv`` (bin depth / subgenome median)."""
    out = df.copy()
    out["bin_depth"] = sample_bin_depth(out, depth_mode=depth_mode)
    out["subgenome"] = out["logical_chr"].map(_subgenome_from_logical_chr)
    if baselines is None:
        baselines = sample_subgenome_baselines(out, min_reads=min_reads, depth_mode=depth_mode)
    base = out["subgenome"].map(baselines).astype(float)
    out["depth_rel"] = out["bin_depth"] / base.replace(0, np.nan)
    out["depth_rel"] = out["depth_rel"].replace([np.inf, -np.inf], np.nan)
    reads = pd.to_numeric(out["total_reads"], errors="coerce").fillna(0)
    out.loc[reads < min_reads, ["bin_depth", "depth_rel"]] = np.nan
    out["is_cnv"] = flag_bin_cnv_from_depth_rel(out["depth_rel"], config=config)
    out.loc[reads < min_reads, "is_cnv"] = False
    return out


def add_sample_mq30_relative_depth(
    df: pd.DataFrame,
    baselines: dict[str, float] | None = None,
    *,
    min_reads: int = 1,
    config: CallConfig | None = None,
) -> pd.DataFrame:
    """Backward-compatible alias using MQ30+ depth approximation."""
    return add_sample_relative_depth(
        df,
        baselines,
        min_reads=min_reads,
        config=config,
        depth_mode="mq30",
    )


def _load_one_sample_tsv(path: str, min_reads: int, depth_mode: DepthMode) -> pd.DataFrame:
    sample = Path(path).name.replace(".mapq_depth_bin.tsv", "")
    df = load_mapq_depth_bin_tsv(path)
    df = add_sample_relative_depth(df, min_reads=min_reads, depth_mode=depth_mode)
    df["sample"] = sample
    return df[list(SAMPLE_REL_COLUMNS)]


def _bin_depth_expr_polars(*, depth_mode: DepthMode) -> pl.Expr:
    reads = pl.col("total_reads").cast(pl.Float64)
    depth = pl.col("bin_mean_depth").cast(pl.Float64)
    if normalize_depth_mode(depth_mode) == "all":
        return depth
    mq30_reads = pl.sum_horizontal([pl.col(c) for c in MQ_GE30_LABELS]).cast(pl.Float64)
    return pl.when(reads > 0).then(depth * mq30_reads / reads).otherwise(None)


def _prepare_sample_polars(
    path: str,
    *,
    min_reads: int,
    config: CallConfig,
    depth_mode: DepthMode,
) -> pl.DataFrame:
    """Load one sample TSV with Polars and attach bin depth / rel / CNV flags."""
    sample = Path(path).name.replace(".mapq_depth_bin.tsv", "")
    df = pl.read_csv(path, separator="\t")
    reads = pl.col("total_reads").cast(pl.Float64)
    bin_depth = _bin_depth_expr_polars(depth_mode=depth_mode)

    df = df.with_columns(
        bin_depth=bin_depth,
        subgenome=pl.col("logical_chr").cast(pl.Utf8).str.slice(-1),
    )

    supported = df.filter(reads >= min_reads)
    baselines: dict[str, float | None] = {}
    for sg, chroms in SUBGENOME_AUTOSOMES.items():
        sub = supported.filter(pl.col("logical_chr").is_in(chroms)).filter(
            pl.col("bin_depth") > 0
        )
        baselines[sg] = (
            float(sub.select(pl.col("bin_depth").median()).item())
            if sub.height
            else None
        )

    baseline_expr = (
        pl.when(pl.col("subgenome") == "A")
        .then(pl.lit(baselines["A"], dtype=pl.Float64))
        .when(pl.col("subgenome") == "B")
        .then(pl.lit(baselines["B"], dtype=pl.Float64))
        .when(pl.col("subgenome") == "D")
        .then(pl.lit(baselines["D"], dtype=pl.Float64))
        .otherwise(None)
    )

    depth_rel = pl.when((reads >= min_reads) & (baseline_expr > 0)).then(
        pl.col("bin_depth") / baseline_expr
    ).otherwise(None)

    df = df.with_columns(depth_rel=depth_rel).with_columns(
        is_cnv=(
            pl.col("depth_rel").is_finite()
            & (
                (pl.col("depth_rel") <= config.deletion_depth_rel_max)
                | (pl.col("depth_rel") >= config.duplication_depth_rel_min)
            )
        )
    )

    return df.select(
        pl.lit(sample).alias("sample"),
        pl.col("logical_chr").cast(pl.Utf8),
        pl.col("bin_idx").cast(pl.Int64),
        pl.col("bin_start").cast(pl.Int64),
        pl.col("bin_end").cast(pl.Int64),
        pl.col("bin_depth"),
        pl.col("depth_rel"),
        pl.col("is_cnv"),
    )


def _aggregate_long_polars(long: pl.DataFrame, *, min_samples: int) -> pl.DataFrame:
    key = ["logical_chr", "bin_idx", "bin_start", "bin_end"]
    return (
        long.group_by(key, maintain_order=True)
        .agg(
            pl.col("bin_depth").mean().alias("depth_mean"),
            pl.col("bin_depth").std().alias("depth_sd"),
            pl.col("depth_rel").mean().alias("rel_mean"),
            pl.col("depth_rel").std().alias("rel_sd"),
            pl.col("bin_depth").count().alias("n_samples"),
            pl.col("is_cnv").sum().alias("n_cnv"),
            pl.col("is_cnv").mean().alias("cnv_frac"),
        )
        .filter(pl.col("n_samples") >= min_samples)
        .sort(["logical_chr", "bin_idx"])
    )


def build_mapq_bin_popdep_polars(
    cohort_long: pl.DataFrame,
    *,
    min_samples: int = 3,
) -> pd.DataFrame:
    """Aggregate a Polars long table to the population grid."""
    grid = _aggregate_long_polars(cohort_long, min_samples=min_samples)
    return (
        grid.with_columns(
            pl.col("n_samples").cast(pl.Int64),
            pl.col("n_cnv").cast(pl.Int64),
        )
        .to_pandas()[list(POPDEP_BIN_COLUMNS)]
    )


def load_cohort_mapq_bin_tables_polars(
    mapq_dir: str | Path,
    *,
    pattern: str = "*.mapq_depth_bin.tsv",
    max_workers: int = 16,
    min_reads: int = 1,
    config: CallConfig | None = None,
    depth_mode: DepthMode = "mq30",
) -> pl.DataFrame:
    """Parallel Polars load of all sample TSVs."""
    cfg = config or CallConfig()
    mode = normalize_depth_mode(depth_mode)
    mapq_dir = Path(mapq_dir)
    paths = sorted(mapq_dir.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No files matching {pattern} under {mapq_dir}")

    workers = max(1, min(int(max_workers), len(paths)))
    print(
        f"[Info] load_cohort_mapq_bin_tables_polars: {len(paths)} files, "
        f"workers={workers}, depth_mode={mode}"
    )
    frames: list[pl.DataFrame] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(
                _prepare_sample_polars,
                str(p),
                min_reads=min_reads,
                config=cfg,
                depth_mode=mode,
            ): p
            for p in paths
        }
        for i, fut in enumerate(as_completed(futures), start=1):
            frames.append(fut.result())
            if i % 500 == 0 or i == len(paths):
                print(f"[Info] loaded {i}/{len(paths)} sample TSVs (polars)")

    return pl.concat(frames, how="vertical_relaxed").sort(
        ["logical_chr", "bin_idx", "sample"]
    )


def load_cohort_mapq_bin_tables(
    mapq_dir: str | Path,
    *,
    pattern: str = "*.mapq_depth_bin.tsv",
    max_workers: int = 16,
    min_reads: int = 1,
    depth_mode: DepthMode = "mq30",
) -> pd.DataFrame:
    """Load all sample TSVs in parallel; return long table with bin depth + rel."""
    mode = normalize_depth_mode(depth_mode)
    mapq_dir = Path(mapq_dir)
    paths = sorted(mapq_dir.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No files matching {pattern} under {mapq_dir}")

    frames: list[pd.DataFrame] = []
    workers = max(1, min(int(max_workers), len(paths)))
    print(
        f"[Info] load_cohort_mapq_bin_tables: {len(paths)} files, "
        f"workers={workers}, depth_mode={mode}"
    )
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(_load_one_sample_tsv, str(p), min_reads, mode): p for p in paths
        }
        for i, fut in enumerate(as_completed(futures), start=1):
            frames.append(fut.result())
            if i % 500 == 0 or i == len(paths):
                print(f"[Info] loaded {i}/{len(paths)} sample TSVs")

    out = pd.concat(frames, ignore_index=True)
    return out.sort_values(["logical_chr", "bin_idx", "sample"]).reset_index(drop=True)


@dataclass
class _BinRunningStats:
    """Online mean / sample variance (Welford) for one 1 Mb bin."""

    bin_start: int
    bin_end: int
    n: int = 0
    n_cnv: int = 0
    depth_mean: float = 0.0
    depth_m2: float = 0.0
    rel_mean: float = 0.0
    rel_m2: float = 0.0

    def update(self, bin_depth: float, depth_rel: float, is_cnv: bool) -> None:
        if not (np.isfinite(bin_depth) and np.isfinite(depth_rel)):
            return
        self.n += 1
        if is_cnv:
            self.n_cnv += 1
        self.depth_mean, self.depth_m2 = _welford_update(
            self.n, self.depth_mean, self.depth_m2, float(bin_depth)
        )
        self.rel_mean, self.rel_m2 = _welford_update(
            self.n, self.rel_mean, self.rel_m2, float(depth_rel)
        )


def _welford_update(
    n: int, mean: float, m2: float, value: float
) -> tuple[float, float]:
    delta = value - mean
    mean += delta / n
    delta2 = value - mean
    m2 += delta * delta2
    return mean, m2


def _running_stats_to_frame(
    stats: dict[tuple[str, int], _BinRunningStats],
    *,
    min_samples: int,
) -> pd.DataFrame:
    rows: list[dict] = []
    for (logical_chr, bin_idx), agg in sorted(
        stats.items(), key=lambda item: (item[0][0], item[0][1])
    ):
        if agg.n < int(min_samples):
            continue
        depth_sd = float(np.sqrt(agg.depth_m2 / (agg.n - 1))) if agg.n > 1 else float("nan")
        rel_sd = float(np.sqrt(agg.rel_m2 / (agg.n - 1))) if agg.n > 1 else float("nan")
        rows.append(
            {
                "logical_chr": logical_chr,
                "bin_idx": bin_idx,
                "bin_start": agg.bin_start,
                "bin_end": agg.bin_end,
                "n_samples": agg.n,
                "depth_mean": agg.depth_mean,
                "depth_sd": depth_sd,
                "rel_mean": agg.rel_mean,
                "rel_sd": rel_sd,
                "n_cnv": agg.n_cnv,
                "cnv_frac": agg.n_cnv / agg.n if agg.n else float("nan"),
            }
        )
    return pd.DataFrame(rows, columns=list(POPDEP_BIN_COLUMNS))


def build_mapq_bin_popdep(
    cohort_long: pd.DataFrame,
    *,
    min_samples: int = 3,
) -> pd.DataFrame:
    """
    Aggregate cohort long table to one row per 1 Mb bin.

    Columns (TIGER-compatible naming shortened):
    - ``depth_mean`` / ``depth_sd``: mean and sample SD of per-sample bin depth
    - ``rel_mean`` / ``rel_sd``: mean and sample SD of per-sample ``depth_rel``
    - ``n_cnv`` / ``cnv_frac``: cohort bins flagged deletion/duplication by depth
    """
    key = ["logical_chr", "bin_idx", "bin_start", "bin_end"]
    grouped = cohort_long.groupby(key, sort=False, as_index=False)

    depth_mean = grouped["bin_depth"].mean()
    depth_sd = grouped["bin_depth"].std(ddof=1)
    rel_mean = grouped["depth_rel"].mean()
    rel_sd = grouped["depth_rel"].std(ddof=1)
    n_samples = grouped["bin_depth"].count()
    n_cnv = grouped["is_cnv"].sum()

    out = depth_mean[key].copy()
    out["n_samples"] = n_samples["bin_depth"].astype(int)
    out["depth_mean"] = depth_mean["bin_depth"]
    out["depth_sd"] = depth_sd["bin_depth"]
    out["rel_mean"] = rel_mean["depth_rel"]
    out["rel_sd"] = rel_sd["depth_rel"]
    out["n_cnv"] = n_cnv["is_cnv"].astype(int)
    out["cnv_frac"] = out["n_cnv"] / out["n_samples"]
    out = out.loc[out["n_samples"] >= int(min_samples)].copy()
    out = out.sort_values(["logical_chr", "bin_idx"]).reset_index(drop=True)
    return out[list(POPDEP_BIN_COLUMNS)]


def build_mapq_bin_popdep_from_dir(
    mapq_dir: str | Path,
    *,
    max_workers: int = 16,
    min_reads: int = 1,
    min_samples: int = 3,
    keep_long: bool = True,
    engine: BuildEngine = "polars",
    config: CallConfig | None = None,
    depth_mode: DepthMode = "mq30",
) -> tuple[pd.DataFrame, pd.DataFrame | pl.DataFrame | None]:
    """
    Load cohort TSVs and return (population grid, optional per-sample long table).

    ``depth_mode="all"`` uses published ``bin_mean_depth`` (no MAPQ filter).
    ``depth_mode="mq30"`` uses the MAPQ>=30 histogram approximation (legacy).
    """
    mode = normalize_depth_mode(depth_mode)
    if engine == "polars":
        long_pl = load_cohort_mapq_bin_tables_polars(
            mapq_dir,
            max_workers=max_workers,
            min_reads=min_reads,
            config=config,
            depth_mode=mode,
        )
        grid = build_mapq_bin_popdep_polars(long_pl, min_samples=min_samples)
        if keep_long:
            return grid, long_pl.to_pandas()
        return grid, None

    mapq_dir = Path(mapq_dir)
    paths = sorted(mapq_dir.glob("*.mapq_depth_bin.tsv"))
    if not paths:
        raise FileNotFoundError(f"No *.mapq_depth_bin.tsv under {mapq_dir}")

    if engine == "pandas":
        long = load_cohort_mapq_bin_tables(
            mapq_dir,
            max_workers=max_workers,
            min_reads=min_reads,
            depth_mode=mode,
        )
        grid = build_mapq_bin_popdep(long, min_samples=min_samples)
        return grid, long if keep_long else None

    workers = max(1, min(int(max_workers), len(paths)))
    print(
        f"[Info] build_mapq_bin_popdep_from_dir(streaming): {len(paths)} files, "
        f"workers={workers}, keep_long={keep_long}, depth_mode={mode}"
    )

    running: dict[tuple[str, int], _BinRunningStats] = {}
    long_frames: list[pd.DataFrame] | None = [] if keep_long else None

    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(_load_one_sample_tsv, str(p), min_reads, mode): p for p in paths
        }
        for i, fut in enumerate(as_completed(futures), start=1):
            sample_long = fut.result()
            if long_frames is not None:
                long_frames.append(sample_long)
            for row in sample_long.itertuples(index=False):
                key = (str(row.logical_chr), int(row.bin_idx))
                agg = running.get(key)
                if agg is None:
                    agg = _BinRunningStats(
                        bin_start=int(row.bin_start),
                        bin_end=int(row.bin_end),
                    )
                    running[key] = agg
                agg.update(float(row.bin_depth), float(row.depth_rel), bool(row.is_cnv))
            if i % 500 == 0 or i == len(paths):
                print(f"[Info] aggregated {i}/{len(paths)} sample TSVs")

    grid = _running_stats_to_frame(running, min_samples=min_samples)
    long = (
        pd.concat(long_frames, ignore_index=True).sort_values(
            ["logical_chr", "bin_idx", "sample"]
        ).reset_index(drop=True)
        if long_frames is not None
        else None
    )
    return grid, long


def write_mapq_bin_popdep(
    grid: pd.DataFrame,
    output_dir: str | Path,
    *,
    prefix: str = "mapq30_1mb",
    per_chromosome: bool = True,
) -> list[Path]:
    """Write combined and optional per-chromosome popdep TSVs."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    combined = output_dir / f"{prefix}.popdep_bins.tsv"
    grid.to_csv(combined, sep="\t", index=False)
    written.append(combined)
    print(f"[Info] wrote {combined} ({len(grid)} bins)")

    if per_chromosome:
        for logical_chr in LOGICAL_AUTOSOMES:
            sub = grid[grid["logical_chr"] == logical_chr]
            if sub.empty:
                continue
            path = output_dir / f"{prefix}.{logical_chr}.popdep_bins.tsv"
            sub.to_csv(path, sep="\t", index=False)
            written.append(path)
    return written


def write_sample_mq30_rel_depth(
    cohort_long: pd.DataFrame,
    output_path: str | Path,
) -> Path:
    """Write per-sample MQ30 relative depth (long format)."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cohort_long.to_csv(output_path, sep="\t", index=False)
    print(f"[Info] wrote {output_path} ({len(cohort_long)} rows)")
    return output_path


def run_mapq_bin_popdep_cohort(
    mapq_dir: str | Path,
    output_dir: str | Path,
    *,
    prefix: str = "mapq30_1mb",
    max_workers: int = 32,
    min_reads: int = 1,
    min_samples: int = 3,
    engine: BuildEngine = "polars",
    plot: bool = True,
    depth_mode: DepthMode = "all",
) -> tuple[pd.DataFrame, list[Path], list[Path]]:
    """
    Build cohort 1 Mb popdep grid, write TSVs, and optional QC scatters.

    Intended for Nextflow ``mapq_bin_popdep_build`` (work dir → ``publishDir``).
    """
    from genetics.germplasm.sample.mapq_bin_popdep_plot import (
        plot_mapq_bin_popdep_qc_scatters,
    )

    mode = normalize_depth_mode(depth_mode)
    output_dir = Path(output_dir)
    grid, _long = build_mapq_bin_popdep_from_dir(
        mapq_dir,
        max_workers=max_workers,
        min_reads=min_reads,
        min_samples=min_samples,
        keep_long=False,
        engine=engine,
        depth_mode=mode,
    )
    tsv_paths = write_mapq_bin_popdep(grid, output_dir, prefix=prefix)
    plot_paths: list[Path] = []
    if plot:
        plot_paths = plot_mapq_bin_popdep_qc_scatters(
            grid,
            output_dir / prefix,
            depth_mode=mode,
        )
    print(
        f"[Info] run_mapq_bin_popdep_cohort: depth_mode={mode}, bins={len(grid)}, "
        f"tsvs={len(tsv_paths)}, plots={len(plot_paths)}"
    )
    return grid, tsv_paths, plot_paths


def run_mapq_bin_popdep_replot(
    grid_tsv: str | Path,
    output_dir: str | Path,
    *,
    prefix: str = "mapq30_1mb",
    depth_mode: DepthMode = "all",
) -> list[Path]:
    """Replot QC scatters from an existing combined popdep grid TSV."""
    from genetics.germplasm.sample.mapq_bin_popdep_plot import (
        plot_mapq_bin_popdep_qc_scatters,
    )

    grid = pd.read_csv(grid_tsv, sep="\t")
    plot_paths = plot_mapq_bin_popdep_qc_scatters(
        grid,
        Path(output_dir) / prefix,
        depth_mode=depth_mode,
    )
    print(f"[Info] run_mapq_bin_popdep_replot: plots={len(plot_paths)}")
    return plot_paths
