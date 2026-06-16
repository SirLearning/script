import fcntl
import gzip
import os
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.variant.variant_utils import load_df_from_plink_variant
from infra.utils.graph import (
    LEGEND_FONT_SIZE,
    TICK_FONT_SIZE,
    TITLE_FONT_SIZE,
    X_LABEL_FONT_SIZE,
    Y_LABEL_FONT_SIZE,
    plot_distribution_with_stats,
    plot_regression_comparison,
)
from infra.utils.io import load_df_from_space_sep_no_header, load_df_from_tsv, save_df_to_tsv

# Frozen vmap4 test_plink site-MQ reference (partial abstract_mq_50_bams one-off build).
DEFAULT_MQ_DIR = "/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/abstract_mq_50_bams"

# Process-step variant MQ table: {subgenome}.mq.info.tsv (PLINK-style #CHROM header).
VARIANT_MQ_INFO_SUFFIX = ".mq.info.tsv"


def variant_mq_info_filename(subgenome_id: str) -> str:
    """Basename for per-subgenome variant MQ annotation under process/<mod>/variant/."""
    return f"{subgenome_id}{VARIANT_MQ_INFO_SUFFIX}"


def site_mq_ref_path(chr_id, mq_dir=None):
    """Path to padded reference grid (*.site_mq.ref.txt.gz preferred)."""
    root = Path(mq_dir or DEFAULT_MQ_DIR) / "reference"
    gz = root / f"{chr_id}.site_mq.ref.txt.gz"
    if gz.exists():
        return gz
    return root / f"{chr_id}.site_mq.ref.txt"


def site_mq_ref_bgz_path(chr_id, mq_dir=None):
    """BGZF copy of site_mq.ref for tabix random access (built once from *.txt.gz)."""
    return Path(mq_dir or DEFAULT_MQ_DIR) / "reference" / f"{chr_id}.site_mq.ref.txt.bgz"


def site_mq_ref_tbi_path(chr_id, mq_dir=None):
    """Tabix index for the BGZF site_mq.ref sidecar."""
    return Path(f"{site_mq_ref_bgz_path(chr_id, mq_dir)}.tbi")


def site_mq_calls_path(chr_id, mq_dir=None):
    """Path to mpileup calls table (*.site_mq.calls.tsv.gz)."""
    return Path(mq_dir or DEFAULT_MQ_DIR) / "reference" / f"{chr_id}.site_mq.calls.tsv.gz"


def load_mq_data(filepath, keep_missing=False):
    """
    Load site MQ (3 columns: Chrom, Position, MQ) from padded reference grid.
    """
    print(f"[Info] Loading MQ reference grid: {filepath}")
    col_names = ["Chrom", "Position", "MQ"]
    df = load_df_from_space_sep_no_header(filepath, col_names=col_names)
    if df is None:
        return None
    df["MQ"] = pd.to_numeric(df["MQ"], errors="coerce")
    if not keep_missing:
        df = df.dropna(subset=["MQ"])
    return df


def load_mq_calls(filepath, keep_missing=False):
    """
    Load per-site mpileup MQ table (5 columns: Chrom, Position, REF, ALT, MQ).
    MQ is float mean MAPQ from bcftools mpileup INFO/I16 (one row per mpileup site).
    """
    print(f"[Info] Loading MQ calls: {filepath}")
    col_names = ["Chrom", "Position", "REF", "ALT", "MQ"]
    df = load_df_from_space_sep_no_header(filepath, col_names=col_names)
    if df is None:
        return None
    df["MQ"] = pd.to_numeric(df["MQ"], errors="coerce")
    if not keep_missing:
        df = df.dropna(subset=["MQ"])
    return df


def load_variant_mq_info(filepath):
    """Load subgenome variant MQ annotation ({subgenome}.mq.info.tsv from process step)."""
    print(f"[Info] Loading variant MQ info: {filepath}")
    df = load_df_from_tsv(filepath)
    if df is None:
        return None
    df.columns = [c.replace("#", "") if isinstance(c, str) else c for c in df.columns]
    df["MQ"] = pd.to_numeric(df["MQ"], errors="coerce")
    return df


def save_variant_mq_info(df: pd.DataFrame, output_path: str) -> None:
    """Write variant MQ table with PLINK-style ``#CHROM`` header via ``save_df_to_tsv``."""
    out = df.copy()
    if "CHROM" in out.columns and "#CHROM" not in out.columns:
        out = out.rename(columns={"CHROM": "#CHROM"})
    cols = ["#CHROM", "ID", "POS", "MQ"]
    out = out[[c for c in cols if c in out.columns]]
    save_df_to_tsv(out, output_path)


def ensure_site_mq_ref_tabix(chr_id: str, mq_dir: str | None = None) -> Path:
    """
    Ensure BGZF + tabix index exist for one chromosome padded site_mq.ref grid.

    Standard ``gzip`` refs are converted once to ``*.site_mq.ref.txt.bgz`` (+ ``.tbi``).
    Tabix requires BGZF (not plain gzip); original ``*.txt.gz`` is left unchanged.
    """
    bgz = site_mq_ref_bgz_path(chr_id, mq_dir)
    tbi = site_mq_ref_tbi_path(chr_id, mq_dir)
    if bgz.exists() and tbi.exists():
        return bgz

    ref_path = site_mq_ref_path(chr_id, mq_dir)
    if not Path(ref_path).exists():
        raise FileNotFoundError(f"MQ reference missing: {ref_path}")

    bgz.parent.mkdir(parents=True, exist_ok=True)
    lock_path = bgz.with_suffix(".bgz.lock")
    with open(lock_path, "a", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        if bgz.exists() and tbi.exists():
            return bgz

        print(f"[Info] Building BGZF + tabix for {chr_id} from {ref_path} (one-time)...")
        ref_str = str(ref_path)
        if ref_str.endswith(".gz"):
            decompress = subprocess.Popen(
                ["gunzip", "-c", ref_str],
                stdout=subprocess.PIPE,
            )
            with open(bgz, "wb") as out_handle:
                bgzip_proc = subprocess.run(
                    ["bgzip", "-c"],
                    stdin=decompress.stdout,
                    stdout=out_handle,
                    check=True,
                )
            decompress.wait()
            if decompress.returncode != 0:
                raise RuntimeError(f"gunzip failed for {ref_path} (exit {decompress.returncode})")
            if bgzip_proc.returncode != 0:
                raise RuntimeError(f"bgzip failed for {chr_id}")
        else:
            subprocess.run(["bgzip", "-c", ref_str], stdout=open(bgz, "wb"), check=True)

        subprocess.run(
            ["tabix", "-s", "1", "-b", "2", "-e", "2", str(bgz)],
            check=True,
        )
        print(f"[Info] Tabix ready: {bgz} (+ .tbi)")
        return bgz


def _tabix_mq_lookup(chrom: int, positions: list[int], bgz_path: Path) -> dict[int, float | None]:
    """Batch tabix fetch for sorted variant positions on one chromosome."""
    if not positions:
        return {}

    result: dict[int, float | None] = {pos: None for pos in positions}
    regions_file = tempfile.NamedTemporaryFile(
        mode="w",
        suffix=f".chr{chrom}.regions.tsv",
        delete=False,
        encoding="utf-8",
    )
    chrom_s = str(chrom)
    try:
        for pos in positions:
            regions_file.write(f"{chrom_s}\t{pos}\t{pos}\n")
        regions_file.close()

        proc = subprocess.run(
            ["tabix", "-R", regions_file.name, str(bgz_path)],
            check=False,
            capture_output=True,
            text=True,
        )
        if proc.returncode not in (0, 1):
            raise RuntimeError(
                f"tabix failed for {bgz_path}: {proc.stderr.strip() or proc.stdout.strip()}"
            )
        for line in proc.stdout.splitlines():
            if not line.strip():
                continue
            _chrom, pos_s, mq_s = line.rstrip().split("\t", 2)
            pos = int(pos_s)
            if mq_s != ".":
                result[pos] = float(mq_s)
    finally:
        os.unlink(regions_file.name)

    return result


def _stream_mq_lookup(chrom: int, positions: list[int], mq_dir: str) -> dict[int, float | None]:
    """Fallback merge-join when tabix sidecar is unavailable."""
    ref_path = site_mq_ref_path(f"chr{chrom:03d}", mq_dir)
    if not Path(ref_path).exists():
        print(f"[Warning] MQ reference missing for chr{chrom}: {ref_path}")
        return {pos: None for pos in positions}

    targets = iter(positions)
    target = next(targets, None)
    result: dict[int, float | None] = {}

    open_fn = gzip.open if str(ref_path).endswith(".gz") else open
    with open_fn(ref_path, "rt") as handle:
        for line in handle:
            if target is None:
                break
            _chrom, pos_s, mq_s = line.rstrip().split("\t", 2)
            pos = int(pos_s)
            while target is not None and target < pos:
                result[target] = None
                target = next(targets, None)
            if target == pos:
                if mq_s == ".":
                    result[target] = None
                else:
                    result[target] = float(mq_s)
                target = next(targets, None)

    while target is not None:
        result[target] = None
        target = next(targets, None)
    return result


def _lookup_mq_for_chromosome(
    chrom: int,
    items: list[tuple[int, str]],
    mq_dir: str,
    use_tabix: bool,
) -> list[dict]:
    """Worker: tabix (or stream) lookup for one chromosome; returns row dicts."""
    positions = [pos for pos, _vid in items]
    chr_id = f"chr{chrom:03d}"
    if use_tabix:
        try:
            bgz = ensure_site_mq_ref_tabix(chr_id, mq_dir)
            mq_map = _tabix_mq_lookup(chrom, positions, bgz)
        except (FileNotFoundError, RuntimeError, subprocess.CalledProcessError) as exc:
            print(f"[Warning] Tabix lookup failed for {chr_id} ({exc}); falling back to stream scan.")
            mq_map = _stream_mq_lookup(chrom, positions, mq_dir)
    else:
        mq_map = _stream_mq_lookup(chrom, positions, mq_dir)

    rows = []
    for pos, variant_id in items:
        mq_val = mq_map.get(pos)
        rows.append(
            {
                "CHROM": chrom,
                "ID": variant_id,
                "POS": pos,
                "MQ": np.nan if mq_val is None else mq_val,
            }
        )
    return rows


def _annotate_chromosome_worker(args: tuple) -> tuple[int, list[dict]]:
    chrom, items, mq_dir, use_tabix = args
    rows = _lookup_mq_for_chromosome(chrom, items, mq_dir, use_tabix)
    return chrom, rows


def annotate_variants_mq_from_pvar(
    pvar_path: str,
    mq_dir: str,
    output_path: str,
    max_workers: int | None = None,
    use_tabix: bool = True,
) -> None:
    """
    Annotate merged subgenome PLINK2 variants with mean MAPQ from frozen site_mq.ref grids.

    Uses BGZF+tabix sidecars on ``*.site_mq.ref.txt.gz`` when available (built once per chr).
    Chromosomes are processed in parallel. Writes ``{subgenome}.mq.info.tsv`` with columns
    ``#CHROM, ID, POS, MQ`` (NaN when reference MQ is missing).
    """
    print(f"[Info] Annotating MQ from pvar={pvar_path} mq_dir={mq_dir} tabix={use_tabix}")
    by_chrom: dict[int, list[tuple[int, str]]] = defaultdict(list)
    with open(pvar_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 4)
            chrom = int(parts[0])
            pos = int(parts[1])
            variant_id = parts[2]
            by_chrom[chrom].append((pos, variant_id))

    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda item: item[0])

    chromosomes = sorted(by_chrom.keys())
    n_workers = max_workers or min(len(chromosomes), os.cpu_count() or 1)
    n_workers = max(1, n_workers)

    rows_by_chrom: dict[int, list[dict]] = {}
    if n_workers == 1 or len(chromosomes) == 1:
        for chrom in chromosomes:
            _, rows = _annotate_chromosome_worker(
                (chrom, by_chrom[chrom], mq_dir, use_tabix)
            )
            rows_by_chrom[chrom] = rows
    else:
        print(f"[Info] Parallel MQ lookup across {len(chromosomes)} chromosomes (workers={n_workers})")
        tasks = [
            (chrom, by_chrom[chrom], mq_dir, use_tabix) for chrom in chromosomes
        ]
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {pool.submit(_annotate_chromosome_worker, task): task[0] for task in tasks}
            for future in as_completed(futures):
                chrom, rows = future.result()
                rows_by_chrom[chrom] = rows

    rows: list[dict] = []
    for chrom in chromosomes:
        rows.extend(rows_by_chrom[chrom])

    out_df = pd.DataFrame(rows)
    n_numeric = int(out_df["MQ"].notna().sum())
    save_variant_mq_info(out_df, output_path)
    print(
        f"[Info] Wrote {len(out_df)} variant MQ rows ({n_numeric} numeric) to {output_path}"
    )


_MQ_BIN5_WIDTH = 5
_MQ_BIN5_SAMPLE_PER_BIN = 100


def _mq_bin5_label(mq: float) -> str:
    """MQ bins of width 5: [0,5), [5,10), [10,15), ..."""
    lo = int(np.floor(float(mq) / _MQ_BIN5_WIDTH) * _MQ_BIN5_WIDTH)
    hi = lo + _MQ_BIN5_WIDTH
    return f"{lo}-{hi}"


def _mq_bin5_center(bin_label: str) -> float:
    lo_s, hi_s = bin_label.split("-", 1)
    return (float(lo_s) + float(hi_s)) / 2.0


def _stratified_sample_mq_bin5(
    plot_df: pd.DataFrame,
    per_bin: int = _MQ_BIN5_SAMPLE_PER_BIN,
    random_seed: int = 42,
) -> pd.DataFrame:
    """Up to ``per_bin`` random variants per 5-wide MQ bin (all if bin is smaller)."""
    work = plot_df.dropna(subset=["MQ", "F_MISS"]).copy()
    if work.empty:
        return work
    work = work.assign(_mq_bin5=work["MQ"].map(_mq_bin5_label))
    parts: list[pd.DataFrame] = []
    for _, grp in work.groupby("_mq_bin5", sort=False):
        if len(grp) > per_bin:
            parts.append(grp.sample(n=per_bin, random_state=random_seed))
        else:
            parts.append(grp)
    return pd.concat(parts, ignore_index=True).drop(columns="_mq_bin5")


def _mq_miss_mean_by_mq_bin5(df: pd.DataFrame) -> pd.DataFrame:
    """Mean F_MISS per 5-MQ bin for overlay curve on regression plots."""
    work = df.dropna(subset=["MQ", "F_MISS"]).copy()
    if work.empty:
        return pd.DataFrame(columns=["MQ", "mean_F_MISS"])
    work["_mq_bin"] = work["MQ"].map(_mq_bin5_label)
    curve = work.groupby("_mq_bin", as_index=False)["F_MISS"].mean().rename(columns={"F_MISS": "mean_F_MISS"})
    curve["MQ"] = curve["_mq_bin"].map(_mq_bin5_center)
    return curve.sort_values("MQ").reset_index(drop=True)[["MQ", "mean_F_MISS"]]


def _plot_mq_miss_regression(
    plot_df: pd.DataFrame,
    output_path: str,
    title: str,
    y_lim: tuple[float, float] | None = (0.0, 1.0),
    stats_df: pd.DataFrame | None = None,
) -> None:
    """MQ vs F_MISS scatter with OLS/Huber regression and 5-MQ binned mean curve."""
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.linear_model import HuberRegressor, LinearRegression

    overlay_df = stats_df if stats_df is not None else plot_df
    local_valid = plot_df[["MQ", "F_MISS"]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(local_valid) < 10:
        print(f"Not enough valid data for {output_path}")
        return

    X = local_valid["MQ"].values.reshape(-1, 1)
    y = local_valid["F_MISS"].values

    ols = LinearRegression()
    ols.fit(X, y)
    ols_score = ols.score(X, y)
    ols_eq = f"y = {ols.coef_[0]:.4f}x + {ols.intercept_:.4f}"

    huber = HuberRegressor(epsilon=1.35)
    huber.fit(X, y)
    huber_score = huber.score(X, y)
    huber_eq = f"y = {huber.coef_[0]:.4f}x + {huber.intercept_:.4f}"

    sns.set_style("white")
    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        data=local_valid,
        x="MQ",
        y="F_MISS",
        alpha=0.3,
        s=15,
        color="#1f77b4",
        edgecolor="none",
        label="Variants",
    )

    x_min, x_max = local_valid["MQ"].min(), local_valid["MQ"].max()
    x_span = x_max - x_min or 1.0
    x_line_min = x_min - 0.05 * x_span
    x_line_max = x_max + 0.05 * x_span
    x_range = np.linspace(x_line_min, x_line_max, 100).reshape(-1, 1)
    plt.plot(
        x_range,
        ols.predict(x_range),
        color="blue",
        linewidth=2,
        linestyle="--",
        label=f"OLS: {ols_eq}, $R^2$={ols_score:.3f}",
    )
    plt.plot(
        x_range,
        huber.predict(x_range),
        color="green",
        linewidth=2,
        label=f"Huber: {huber_eq}, $R^2$={huber_score:.3f}",
    )

    mean_curve = _mq_miss_mean_by_mq_bin5(overlay_df)
    if len(mean_curve) >= 2:
        plt.plot(
            mean_curve["MQ"],
            mean_curve["mean_F_MISS"],
            color="darkorange",
            linewidth=2,
            linestyle="-",
            label=f"Mean F_MISS ({_MQ_BIN5_WIDTH}-MQ bins)",
        )

    plt.title(title, fontsize=TITLE_FONT_SIZE)
    plt.xlabel("Mean MAPQ", fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel("Missing rate (F_MISS)", fontsize=Y_LABEL_FONT_SIZE)
    plt.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
    plt.xlim(x_line_min, x_line_max)
    if y_lim is not None:
        plt.ylim(y_lim)
    plt.legend(loc="upper left", bbox_to_anchor=(0.0, -0.15), fontsize=LEGEND_FONT_SIZE, frameon=False)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved plot to {output_path}")
    plt.close()


def _load_mq_miss_merged(mq_path: str, vmiss_path: str) -> pd.DataFrame | None:
    df_mq = load_variant_mq_info(mq_path)
    df_miss = load_df_from_plink_variant(vmiss_path)
    if df_mq is None or df_miss is None:
        return None
    merged = pd.merge(
        df_miss[["ID", "F_MISS"]],
        df_mq[["ID", "MQ"]],
        on="ID",
        how="inner",
    ).dropna(subset=["MQ", "F_MISS"])
    return merged if not merged.empty else None


def _resolve_mq_output_paths(output_prefix: str, info_path: str | None = None) -> dict[str, str]:
    """
    Derive MQ distribution and info.tsv paths from the mq_miss plot prefix.

    Example: ``A.variant.mq_miss`` → dist ``A.variant.mq.dist.png``,
    info ``A.variant.mq_miss.info.tsv`` (or under sibling ``info/`` when prefix is in ``plots/``).
    """
    prefix_path = Path(output_prefix)
    stem = prefix_path.name
    if stem.endswith(".mq_miss"):
        mq_dist_stem = stem[: -len(".mq_miss")] + ".mq"
    else:
        mq_dist_stem = f"{stem}.mq"

    if prefix_path.parent and str(prefix_path.parent) not in ("", "."):
        mq_dist_prefix = str(prefix_path.parent / mq_dist_stem)
    else:
        mq_dist_prefix = mq_dist_stem

    if info_path:
        info_tsv = info_path
    elif prefix_path.parent.name == "plots" and prefix_path.parent.parent.is_dir():
        info_tsv = str(prefix_path.parent.parent / "info" / f"{stem}.info.tsv")
    else:
        info_tsv = f"{output_prefix}.info.tsv"

    return {
        "mq_miss": str(prefix_path),
        "mq_dist": mq_dist_prefix,
        "info": info_tsv,
    }


def ana_mq_missing_reg(
    mq_path: str,
    vmiss_path: str,
    output_prefix: str,
    info_path: str | None = None,
    regression_max_points: int = 50000,
    random_seed: int = 42,
) -> None:
    """MQ distribution, MQ vs F_MISS regression, log-missing-only plot, and 5-MQ bin sample plot."""
    paths = _resolve_mq_output_paths(output_prefix, info_path)
    merged = _load_mq_miss_merged(mq_path, vmiss_path)
    if merged is None:
        print("[Error] Could not load MQ or vmiss inputs, or no overlapping variants.")
        return

    n_variants = len(merged)
    mq_numeric = merged["MQ"].dropna()
    stats: dict[str, float | int] = {
        "Total_Variants": n_variants,
        "Pearson_r_MQ_F_MISS": float(merged["MQ"].corr(merged["F_MISS"])),
        "Mean_MQ": float(mq_numeric.mean()),
        "Median_MQ": float(mq_numeric.median()),
        "Mean_F_MISS": float(merged["F_MISS"].mean()),
        "Median_F_MISS": float(merged["F_MISS"].median()),
    }

    # MQ distribution (numeric MQ only): {subgenome}.variant.mq.dist.png
    mq_dist_df = merged.dropna(subset=["MQ"])
    if not mq_dist_df.empty:
        plot_distribution_with_stats(
            data=mq_dist_df,
            col="MQ",
            title="Distribution of variant mean MAPQ",
            filename=f"{paths['mq_dist']}.dist.png",
            mean_val=float(mq_dist_df["MQ"].mean()),
            median_val=float(mq_dist_df["MQ"].median()),
            std_val=float(mq_dist_df["MQ"].std()),
            x_label="Mean MAPQ",
            y_label="Count of variants",
            bins=50,
        )

    miss_prefix = paths["mq_miss"]
    # Full-range regression (subsampled when very large)
    if n_variants > regression_max_points:
        merged_plot = merged.sample(regression_max_points, random_state=random_seed)
    else:
        merged_plot = merged
    stats["Regression_Points"] = len(merged_plot)
    _plot_mq_miss_regression(
        merged_plot,
        f"{miss_prefix}.reg.png",
        title="Variant MQ vs missing rate",
        y_lim=(0.0, 1.0),
        stats_df=merged,
    )

    # MQ vs log(F_MISS): log missing rate only (F_MISS in (0,1) => log(F_MISS) < 0)
    log_miss_df = merged[(merged["F_MISS"] > 0) & (merged["F_MISS"] < 1)].copy()
    if len(log_miss_df) >= 10:
        if len(log_miss_df) > regression_max_points:
            log_miss_df = log_miss_df.sample(regression_max_points, random_state=random_seed)
        log_miss_df["Log_F_MISS"] = np.log(log_miss_df["F_MISS"])
        plot_regression_comparison(
            log_miss_df,
            "MQ",
            "Log_F_MISS",
            x_label="Mean MAPQ",
            y_label="log missing rate",
            filename=f"{miss_prefix}.reg.log_miss.png",
            y_lim=(None, 0.0),
            title="Variant MQ vs log missing rate",
        )
        stats["Regression_Points_Log_Miss"] = len(log_miss_df)

    # Stratified 5-MQ bin sample (up to 100 variants per bin), MAC bin50s-style
    sample_df = _stratified_sample_mq_bin5(merged, random_seed=random_seed)
    n_sampled = len(sample_df)
    n_bins = int(merged["MQ"].map(_mq_bin5_label).nunique()) if n_variants else 0
    stats["N_mq_bin5_bins"] = n_bins
    stats["N_sampled_bin5s"] = n_sampled
    if n_sampled >= 10:
        _plot_mq_miss_regression(
            sample_df,
            f"{miss_prefix}.reg.bin5s.png",
            title=(
                f"MQ vs missing (5-MQ bin sample, "
                f"n<={_MQ_BIN5_SAMPLE_PER_BIN}/bin)"
            ),
            y_lim=(0.0, 1.0),
            stats_df=merged,
        )
        stats["Pearson_r_sampled_bin5s"] = (
            float(sample_df["MQ"].corr(sample_df["F_MISS"])) if n_sampled > 1 else np.nan
        )

    save_df_to_tsv(pd.DataFrame([stats]), paths["info"])
