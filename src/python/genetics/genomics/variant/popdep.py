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
from scipy.optimize import curve_fit
from scipy.stats import chi2

from genetics.genomics.variant.mq import load_mq_data
from infra.wheat.ref_v1 import is_ref_v1_chr_name
from infra.utils.graph import (
    DEFAULT_AXIS_PADDING_FRACTION,
    axis_limits_with_adaptive_upper,
    padded_axis_limits,
    plot_bar_chart,
    plot_distribution_with_stats,
    plot_mean_variance_fit,
    plot_qq_residuals,
    plot_regression_comparison,
    plot_scatter_with_outliers,
    plot_scatter_with_thresholds,
)
from infra.utils.io import load_df_from_space_sep, load_df_from_tsv, save_df_to_tsv

from .variant_utils import load_df_from_plink_variant

# ==================================================================================
# Frozen reference paths and variant annotation (mirror site MQ workflow)
# ==================================================================================

def popdep_chr_ref_bgz_path(popdep_dir: str, chrom: int) -> Path:
    """BGZF popdepth grid for tabix random access."""
    return Path(popdep_dir) / "variant" / f"chr{int(chrom):03d}.popdep.txt.bgz"


def popdep_chr_ref_tbi_path(popdep_dir: str, chrom: int) -> Path:
    """Tabix index for the BGZF popdepth grid."""
    return Path(f"{popdep_chr_ref_bgz_path(popdep_dir, chrom)}.tbi")


def popdep_chr_ref_path(popdep_dir: str, chrom: int) -> str:
    """Path to frozen per-chromosome TIGER popdepth grid under ``{popdep_dir}/variant/``."""
    variant = Path(popdep_dir) / "variant"
    bgz = variant / f"chr{int(chrom):03d}.popdep.txt.bgz"
    if bgz.exists():
        return str(bgz)
    txt = variant / f"chr{int(chrom):03d}.popdep.txt"
    if txt.exists():
        return str(txt)
    gz = variant / f"chr{int(chrom):03d}.popdep.txt.gz"
    if gz.exists():
        return str(gz)
    return str(txt)


def ensure_popdep_ref_tabix(chr_id: str, popdep_dir: str) -> Path:
    """
    Ensure BGZF + tabix index exist for one chromosome popdepth grid.

    New builds publish ``*.popdep.txt.bgz`` directly. Legacy plain ``*.popdep.txt``
    or ``*.popdep.txt.gz`` refs are converted once to BGZF (+ ``.tbi``).
    """
    chrom = int(chr_id.replace("chr", ""))
    bgz = popdep_chr_ref_bgz_path(popdep_dir, chrom)
    tbi = popdep_chr_ref_tbi_path(popdep_dir, chrom)
    if bgz.exists() and tbi.exists():
        return bgz

    ref_path = Path(popdep_chr_ref_path(popdep_dir, chrom))
    if ref_path.suffix == ".bgz" and ref_path.exists():
        subprocess.run(
            ["tabix", "-S", "1", "-s", "1", "-b", "2", "-e", "2", "-f", str(ref_path)],
            check=True,
        )
        return ref_path
    if not ref_path.exists():
        raise FileNotFoundError(f"Popdep reference missing: {ref_path}")

    bgz.parent.mkdir(parents=True, exist_ok=True)
    lock_path = bgz.with_suffix(".bgz.lock")
    with open(lock_path, "a", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        if bgz.exists() and tbi.exists():
            return bgz

        print(f"[Info] Building BGZF + tabix for {chr_id} from {ref_path} (one-time)...")
        ref_str = str(ref_path)
        if ref_str.endswith(".bgz"):
            return ref_path
        if ref_str.endswith(".gz"):
            decompress = subprocess.Popen(
                ["gunzip", "-c", ref_str],
                stdout=subprocess.PIPE,
                text=True,
            )
            with open(bgz, "wb") as out_handle:
                bgzip_proc = subprocess.Popen(
                    ["bgzip", "-c"],
                    stdin=subprocess.PIPE,
                    stdout=out_handle,
                )
                assert decompress.stdout is not None
                header = decompress.stdout.readline()
                bgzip_proc.stdin.write(f"Chrom\t{header}".encode("utf-8"))
                for line in decompress.stdout:
                    bgzip_proc.stdin.write(f"{chrom}\t".encode("utf-8") + line.encode("utf-8"))
                bgzip_proc.stdin.close()
                bgzip_proc.wait()
            decompress.wait()
            if decompress.returncode != 0:
                raise RuntimeError(f"gunzip failed for {ref_path} (exit {decompress.returncode})")
            if bgzip_proc.returncode != 0:
                raise RuntimeError(f"bgzip failed for {chr_id}")
        else:
            with open(ref_str, encoding="utf-8") as src_handle, open(bgz, "wb") as out_handle:
                bgzip_proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=out_handle)
                header = src_handle.readline()
                bgzip_proc.stdin.write(f"Chrom\t{header}".encode("utf-8"))
                for line in src_handle:
                    bgzip_proc.stdin.write(f"{chrom}\t".encode("utf-8") + line.encode("utf-8"))
                bgzip_proc.stdin.close()
                bgzip_proc.wait()
                if bgzip_proc.returncode != 0:
                    raise RuntimeError(f"bgzip failed for {chr_id}")

        subprocess.run(
            ["tabix", "-S", "1", "-s", "1", "-b", "2", "-e", "2", str(bgz)],
            check=True,
        )
        print(f"[Info] Tabix ready: {bgz} (+ .tbi)")
        return bgz


def _popdep_grid_column_map(header_cols: list[str]) -> dict[str, int | None]:
    """Map TIGER / BGZF popdepth header names to column indices."""
    clean = [c.replace("#", "") if isinstance(c, str) else c for c in header_cols]

    def _idx(primary: str, aliases: list[str]) -> int | None:
        for name in [primary, *aliases]:
            if name in clean:
                return clean.index(name)
        return None

    return {
        "pos": _idx("Position", []),
        "mean": _idx("Depth_Mean", []),
        "sd": _idx("Depth_SD", []),
        "rel_mean": _idx("RelativeDepth_Mean", ["RelDepth_Mean"]),
        "rel_sd": _idx("RelativeDepth_SD", ["RelDepth_SD"]),
    }


def _parse_popdep_grid_row(
    parts: list[str],
    col_map: dict[str, int | None],
) -> tuple[int, float, float, float, float] | None:
    """Parse one popdepth grid row into position and absolute/relative depth stats."""
    pos_i = col_map.get("pos")
    mean_i = col_map.get("mean")
    sd_i = col_map.get("sd")
    if pos_i is None or mean_i is None or sd_i is None:
        return None
    if len(parts) <= max(pos_i, mean_i, sd_i):
        return None

    rel_mean = np.nan
    rel_sd = np.nan
    rel_mean_i = col_map.get("rel_mean")
    rel_sd_i = col_map.get("rel_sd")
    if rel_mean_i is not None and rel_sd_i is not None and len(parts) > rel_sd_i:
        rel_mean = float(parts[rel_mean_i])
        rel_sd = float(parts[rel_sd_i])

    return (
        int(parts[pos_i]),
        float(parts[mean_i]),
        float(parts[sd_i]),
        rel_mean,
        rel_sd,
    )


def _default_bgz_column_map() -> dict[str, int | None]:
    """Default indices after ``Chrom`` was prepended to the TIGER header."""
    return {
        "pos": 1,
        "mean": 2,
        "sd": 3,
        "rel_mean": 4,
        "rel_sd": 5,
    }


def _tabix_popdep_lookup(
    chrom: int,
    positions: list[int],
    bgz_path: Path,
) -> dict[int, tuple[float, float, float, float]]:
    """Batch tabix fetch: (Depth_Mean, Depth_SD, RelativeDepth_Mean, RelativeDepth_SD)."""
    if not positions:
        return {}

    col_map = _default_bgz_column_map()
    result: dict[int, tuple[float, float, float, float]] = {}
    regions_file = tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".popdep.regions.tsv",
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
            parsed = _parse_popdep_grid_row(line.rstrip().split("\t"), col_map)
            if parsed is None:
                continue
            pos, mean, sd, rel_mean, rel_sd = parsed
            result[pos] = (mean, sd, rel_mean, rel_sd)
    finally:
        os.unlink(regions_file.name)

    return result


def _stream_popdep_lookup(
    popdep_path: str,
    wanted_positions: list[int],
) -> dict[int, tuple[float, float, float, float]]:
    """Fallback merge-join when tabix sidecar is unavailable."""
    if not wanted_positions:
        return {}
    wanted_set = set(wanted_positions)
    max_pos = max(wanted_positions)
    result: dict[int, tuple[float, float, float, float]] = {}
    if not os.path.isfile(popdep_path):
        print(f"[Warn] Popdep reference missing: {popdep_path}")
        return result

    open_fn = gzip.open if popdep_path.endswith((".gz", ".bgz")) else open
    with open_fn(popdep_path, "rt", encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            return result
        header_cols = header.rstrip("\n").split("\t")
        col_map = _popdep_grid_column_map(header_cols)
        for line in handle:
            parsed = _parse_popdep_grid_row(line.rstrip("\n").split("\t"), col_map)
            if parsed is None:
                continue
            pos, mean, sd, rel_mean, rel_sd = parsed
            if pos in wanted_set:
                result[pos] = (mean, sd, rel_mean, rel_sd)
            if pos >= max_pos and len(result) == len(wanted_set):
                break
    return result


def save_variant_popdep_info(df: pd.DataFrame, output_path: str) -> None:
    """Write variant popdepth table with PLINK-style ``#CHROM`` header."""
    out = df.copy()
    if "CHROM" in out.columns and "#CHROM" not in out.columns:
        out = out.rename(columns={"CHROM": "#CHROM"})
    cols = [
        "#CHROM",
        "ID",
        "POS",
        "Depth_Mean",
        "Depth_SD",
        "Depth_CV",
        "RelativeDepth_Mean",
        "RelativeDepth_SD",
        "RelativeDepth_CV",
    ]
    out = out[[c for c in cols if c in out.columns]]
    save_df_to_tsv(out, output_path)


def _lookup_popdep_for_chromosome(
    chrom: int,
    items: list[tuple[int, str]],
    popdep_dir: str,
    use_tabix: bool,
) -> list[dict]:
    positions = [pos for pos, _variant_id in items]
    if use_tabix:
        try:
            bgz = ensure_popdep_ref_tabix(f"chr{chrom:03d}", popdep_dir)
            popdep_map = _tabix_popdep_lookup(chrom, positions, bgz)
        except (FileNotFoundError, RuntimeError, subprocess.CalledProcessError) as exc:
            print(f"[Warning] Tabix lookup failed for chr{chrom:03d} ({exc}); falling back to stream scan.")
            popdep_map = _stream_popdep_lookup(popdep_chr_ref_path(popdep_dir, chrom), positions)
    else:
        popdep_map = _stream_popdep_lookup(popdep_chr_ref_path(popdep_dir, chrom), positions)
    rows: list[dict] = []
    for pos, variant_id in items:
        depth = popdep_map.get(pos)
        depth_mean = np.nan
        depth_sd = np.nan
        depth_cv = np.nan
        rel_mean = np.nan
        rel_sd = np.nan
        rel_cv = np.nan
        if depth is not None:
            depth_mean, depth_sd, rel_mean, rel_sd = depth
            if depth_mean > 0:
                depth_cv = depth_sd / depth_mean
            if rel_mean > 0:
                rel_cv = rel_sd / rel_mean
        rows.append(
            {
                "CHROM": chrom,
                "ID": variant_id,
                "POS": pos,
                "Depth_Mean": depth_mean,
                "Depth_SD": depth_sd,
                "Depth_CV": depth_cv,
                "RelativeDepth_Mean": rel_mean,
                "RelativeDepth_SD": rel_sd,
                "RelativeDepth_CV": rel_cv,
            }
        )
    return rows


def _annotate_popdep_chromosome_worker(args: tuple) -> tuple[int, list[dict]]:
    chrom, items, popdep_dir, use_tabix = args
    rows = _lookup_popdep_for_chromosome(chrom, items, popdep_dir, use_tabix)
    return chrom, rows


def annotate_variants_popdep_from_pvar(
    pvar_path: str,
    popdep_dir: str,
    output_path: str,
    max_workers: int | None = None,
    use_tabix: bool = True,
) -> None:
    """
    Annotate merged subgenome PLINK2 variants with population depth from frozen main_raw grids.

    Uses BGZF+tabix on ``*.popdep.txt.bgz`` when available (built at freeze time or once from legacy
    ``*.popdep.txt`` / ``*.popdep.txt.gz``). Writes ``{subgenome}.popdep.info.tsv`` with columns
    ``#CHROM, ID, POS, Depth_Mean, Depth_SD, Depth_CV, RelativeDepth_Mean, RelativeDepth_SD, RelativeDepth_CV``
    (NaN when reference depth is missing).
    """
    print(f"[Info] Annotating popdep from pvar={pvar_path} popdep_dir={popdep_dir} tabix={use_tabix}")
    by_chrom: dict[int, list[tuple[int, str]]] = defaultdict(list)
    with open(pvar_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 4)
            # Frozen popdep grids are per arm segment with arm-local positions, so a
            # merged ref-chr pvar cannot be joined against them.
            if is_ref_v1_chr_name(parts[0]):
                raise ValueError(
                    f"{pvar_path} carries merged ref-chr names (#CHROM={parts[0]}); "
                    "annotate popdep on the arm-numbered source pfile instead."
                )
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
            _, rows = _annotate_popdep_chromosome_worker(
                (chrom, by_chrom[chrom], popdep_dir, use_tabix)
            )
            rows_by_chrom[chrom] = rows
    else:
        print(
            f"[Info] Parallel popdep lookup across {len(chromosomes)} chromosomes (workers={n_workers})"
        )
        tasks = [(chrom, by_chrom[chrom], popdep_dir, use_tabix) for chrom in chromosomes]
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {
                pool.submit(_annotate_popdep_chromosome_worker, task): task[0]
                for task in tasks
            }
            for future in as_completed(futures):
                chrom, rows = future.result()
                rows_by_chrom[chrom] = rows

    rows: list[dict] = []
    for chrom in chromosomes:
        rows.extend(rows_by_chrom[chrom])

    out_df = pd.DataFrame(rows)
    n_numeric = int(out_df["Depth_Mean"].notna().sum())
    save_variant_popdep_info(out_df, output_path)
    print(
        f"[Info] Wrote {len(out_df)} variant popdep rows ({n_numeric} numeric) to {output_path}"
    )

# ==================================================================================
# Data Loading Helpers
# ==================================================================================

def load_popdep_data(filepath):
    """
    Loads population depth file.
    Expected columns: Position, Depth_Mean, Depth_SD (tab-separated; plain, gzip, or BGZF).
    """
    print(f"[Info] Loading PopData: {filepath}")
    path_str = str(filepath)
    try:
        if path_str.endswith((".gz", ".bgz")):
            with gzip.open(path_str, "rt", encoding="utf-8") as handle:
                df = pd.read_csv(handle, sep=r"\s+")
        else:
            df = load_df_from_space_sep(path_str)
    except Exception as exc:
        print(f"[Error] Failed to read popdep file {path_str}: {exc}")
        return None
    if df is None: return None
    
    # Validation
    req_cols = ['Position', 'Depth_Mean', 'Depth_SD']
    if not all(c in df.columns for c in req_cols):
        print(f"[Error] Missing columns in popdep file. Required: {req_cols}")
        return None

    # Derived Metrics
    df['Depth_Var'] = df['Depth_SD'] ** 2
    # Avoid division by zero
    df['Depth_CV'] = np.where(
        df['Depth_Mean'] > 0, 
        df['Depth_SD'] / df['Depth_Mean'], 
        np.nan
    )
    return df

def load_df_from_tsv_popdep_info(filepath: str) -> pd.DataFrame | None:
    """Load variant-level ``*.popdep.info.tsv`` and normalize to popdep analysis columns."""
    print(f"[Info] Loading variant popdep info: {filepath}")
    df = load_df_from_space_sep(filepath)
    if df is None:
        return None
    df.columns = [c.replace("#", "") if isinstance(c, str) else c for c in df.columns]
    if "POS" in df.columns and "Position" not in df.columns:
        df = df.rename(columns={"POS": "Position"})
    req_cols = ["Position", "Depth_Mean", "Depth_SD"]
    if not all(c in df.columns for c in req_cols):
        print(f"[Error] Missing columns in popdep info. Required: {req_cols}")
        return None
    if "Depth_Var" not in df.columns:
        df["Depth_Var"] = df["Depth_SD"] ** 2
    if "Depth_CV" not in df.columns:
        df["Depth_CV"] = np.where(
            df["Depth_Mean"] > 0,
            df["Depth_SD"] / df["Depth_Mean"],
            np.nan,
        )
    if "RelativeDepth_Mean" in df.columns and "RelativeDepth_SD" in df.columns:
        if "RelativeDepth_Var" not in df.columns:
            df["RelativeDepth_Var"] = df["RelativeDepth_SD"] ** 2
        if "RelativeDepth_CV" not in df.columns:
            df["RelativeDepth_CV"] = np.where(
                df["RelativeDepth_Mean"] > 0,
                df["RelativeDepth_SD"] / df["RelativeDepth_Mean"],
                np.nan,
            )
    return df

# ==================================================================================
# Analysis 1: Mean-Variance Metrics
# ==================================================================================

def nb_variance_function(x, phi):
    """Negative Binomial Variance: Var = Mean + phi * Mean^2"""
    return x + phi * (x ** 2)

def calculate_criteria(y_true, y_pred, k, n):
    """Calculates AIC, BIC based on RSS."""
    rss = np.sum((y_true - y_pred) ** 2)
    aic = n * np.log(rss / n) + 2 * k
    bic = n * np.log(rss / n) + k * np.log(n)
    return aic, bic, rss

def ana_popdep_curve_fit(
    popdep_file, 
    output_prefix, 
    sample_size=50000
):
    """
    Fits Poisson vs Negative Binomial models to Mean-Variance relationship.
    """
    df = load_popdep_data(popdep_file)
    if df is None: return

    # Filter valid
    df_valid = df[(df['Depth_Mean'] > 0) & (df['Depth_Var'] > 0)].dropna()
    print(f"Valid sites for curve fit: {len(df_valid)}")

    x_data = df_valid['Depth_Mean'].values
    y_data = df_valid['Depth_Var'].values
    n_samples = len(y_data)

    # 1. Poisson
    y_pred_poisson = x_data

    # 2. Negative Binomial
    print("Fitting Negative Binomial Model...")
    phi_est = 0
    y_pred_nb = y_pred_poisson
    try:
        popt, pcov = curve_fit(nb_variance_function, x_data, y_data, bounds=(0, np.inf))
        phi_est = popt[0]
        y_pred_nb = nb_variance_function(x_data, phi_est)
        print(f"Estimated Phi: {phi_est:.6f}")
    except Exception as e:
        print(f"Curve fitting failed: {e}")
        # Continue with best effort

    # 3. Report
    aic_p, bic_p, rss_p = calculate_criteria(y_data, y_pred_poisson, 0, n_samples)
    aic_nb, bic_nb, rss_nb = calculate_criteria(y_data, y_pred_nb, 1, n_samples)
    
    print(f"Poisson | RSS: {rss_p:.2e} | AIC: {aic_p:.2f}")
    print(f"NegBin  | RSS: {rss_nb:.2e} | AIC: {aic_nb:.2f}")

    # 4. Plot Mean-Variance
    x_range = np.linspace(df_valid['Depth_Mean'].min(), df_valid['Depth_Mean'].max(), 500)
    
    plot_mean_variance_fit(
        x_data=df_valid['Depth_Mean'],
        y_data=df_valid['Depth_Var'],
        x_line=x_range,
        y_line_poisson=x_range,
        y_line_nb=nb_variance_function(x_range, phi_est),
        phi_est=phi_est,
        filename=f"{output_prefix}.mean_loss_fit.png",
        sample_size=sample_size
    )

    # 5. QQ Plot
    resid_nb = y_data - y_pred_nb
    if np.std(resid_nb) > 0:
        resid_std = (resid_nb - np.mean(resid_nb)) / np.std(resid_nb)
        plot_qq_residuals(resid_std, f"{output_prefix}.qq_residuals.png", title="QQ Plot of Residuals (NB)")

# ==================================================================================
# Analysis 2: Missing Rate vs Depth
# ==================================================================================

def ana_popdep_missing_reg(
    popdep_file, 
    vmiss_file, 
    output_prefix, 
    log_scale=True
):
    """
    Analyzes relationship between Depth Metrics and Variant Missing Rate.

    ``popdep_file`` may be a full per-chromosome TIGER grid or a variant-level
    ``*.popdep.info.tsv`` sidecar produced by ``annotate_variants_popdep_from_pvar``.
    """
    # 1. Load
    if str(popdep_file).endswith(".popdep.info.tsv"):
        df_dep = load_df_from_tsv_popdep_info(popdep_file)
    else:
        df_dep = load_popdep_data(popdep_file)
    df_miss = load_df_from_plink_variant(vmiss_file) # Using renamed function
    if df_dep is None or df_miss is None: return

    # Merge
    merged = pd.merge(df_miss[['Position', 'F_MISS']], df_dep, on='Position', how='inner')
    if merged.empty:
        print("[Error] No overlapping positions found.")
        return

    # Regression Plots
    metrics = [
        ('Depth_Mean', 'Mean Depth', False),
        ('Depth_CV', 'Depth CV', False),
    ]
    
    if log_scale:
        for col in ['Depth_Mean', 'Depth_CV', 'F_MISS']:
            merged[f'Log_{col}'] = np.log(merged[col].replace(0, np.nan))
        metrics.extend([
            ('Log_Depth_Mean', 'Log Mean Depth', True),
            ('Log_Depth_CV', 'Log Depth CV', True)
        ])

    for x_col, x_label, is_log in metrics:
        y_col = 'Log_F_MISS' if is_log else 'F_MISS'
        y_label = 'Log Missing Rate' if is_log else 'Missing Rate'

        plot_df = merged
        y_lim = None
        x_lim = None
        if y_col == 'F_MISS':
            y_lim = (0.0, 1.0)
        elif y_col == 'Log_F_MISS':
            plot_df = merged[(merged['F_MISS'] > 0) & (merged['F_MISS'] < 1)].copy()
            if len(plot_df) < 10:
                print(f"[Warning] Not enough F_MISS in (0,1) for log-scale plot: {x_col}")
                continue
            y_lim = (None, 0.0)
        if x_col == 'Depth_Mean':
            x_lim = (0.0, None)

        plot_regression_comparison(
            plot_df, x_col, y_col,
            x_label=x_label, y_label=y_label,
            filename=f"{output_prefix}.reg_{y_col}_vs_{x_col}.png",
            title=f"{y_label} vs {x_label}",
            x_lim=x_lim,
            y_lim=y_lim,
        )

# ==================================================================================
# Analysis 3: Mahalanobis Distance for Variant Depth Distribution
# ==================================================================================

def ana_popdep_mahalanobis(
    popdep_file,
    output_prefix,
    threshold_p_value=0.99,
    relative_depth=False,
):
    """
    Mahalanobis outlier detection on log(mean) vs log(CV).

    When ``relative_depth=True``, uses RelativeDepth_Mean / RelativeDepth_CV columns.
    """
    metric_tag = "reldepth" if relative_depth else "depth"
    mean_col = "RelativeDepth_Mean" if relative_depth else "Depth_Mean"
    cv_col = "RelativeDepth_CV" if relative_depth else "Depth_CV"
    sd_col = "RelativeDepth_SD" if relative_depth else "Depth_SD"
    log_mean_col = f"Log_{mean_col}"
    log_cv_col = f"Log_{cv_col}"

    print(
        f"[Info] Running Mahalanobis ({metric_tag}) on: {popdep_file}"
    )
    if str(popdep_file).endswith(".popdep.info.tsv"):
        df = load_df_from_tsv_popdep_info(popdep_file)
    else:
        df = load_popdep_data(popdep_file)
    if df is None:
        return
    if relative_depth and mean_col not in df.columns:
        print(f"[Warning] Skipping relative-depth Mahalanobis; missing {mean_col}")
        return
    if "Depth_Var" not in df.columns:
        df["Depth_Var"] = df["Depth_SD"] ** 2
    if relative_depth and "RelativeDepth_Var" not in df.columns and sd_col in df.columns:
        df["RelativeDepth_Var"] = df[sd_col] ** 2

    df_valid = df[(df[mean_col] > 0) & (df[cv_col] > 0)].copy()
    if df_valid.empty:
        print(f"[Error] No valid data for Mahalanobis ({metric_tag}).")
        return

    df_valid[log_mean_col] = np.log(df_valid[mean_col])
    df_valid[log_cv_col] = np.log(df_valid[cv_col])

    data = df_valid[[log_mean_col, log_cv_col]].values
    
    # Fast vectorized Mahalanobis calculation
    mean_vec = np.mean(data, axis=0)
    cov_matrix = np.cov(data, rowvar=False)
    
    try:
        inv_cov_matrix = np.linalg.inv(cov_matrix)
    except np.linalg.LinAlgError:
        print("[Error] Covariance matrix is singular, cannot calculate Mahalanobis distance.")
        return

    diff = data - mean_vec
    # Vectorized computation: Mahalanobis squared distance
    mahalanobis_sq = np.sum(np.dot(diff, inv_cov_matrix) * diff, axis=1)
    
    df_valid['Mahalanobis_Sq'] = mahalanobis_sq
    df_valid['Mahalanobis'] = np.sqrt(mahalanobis_sq)

    # Threshold based on Chi-Square distribution (df=2)
    threshold = chi2.ppf(threshold_p_value, df=2)
    df_valid['Outlier'] = df_valid['Mahalanobis_Sq'] > threshold

    out_file = f"{output_prefix}.info.tsv"
    out_cols = [
        "Position",
        mean_col,
        sd_col,
        cv_col,
        "Mahalanobis",
        "Mahalanobis_Sq",
        "Outlier",
    ]
    save_df_to_tsv(df_valid[out_cols], out_file)

    metric_label = "Relative Depth" if relative_depth else "Depth"
    plot_file = f"{output_prefix}.mahalanobis.png"
    plot_scatter_with_outliers(
        data=df_valid,
        x_col=log_mean_col,
        y_col=log_cv_col,
        outlier_col="Outlier",
        title=(
            f"Mahalanobis Outlier Detection ({metric_label}, p={threshold_p_value})"
        ),
        filename=plot_file,
        xlabel=f"Log {metric_label} Mean",
        ylabel=f"Log {metric_label} CV",
        color_normal="blue",
        color_outlier="red",
        s=2,
    )


_POPDEP_QC_METRIC_SPECS = (
    {
        "tag": "depth",
        "mean": "Depth_Mean",
        "sd": "Depth_SD",
        "var": "Depth_Var",
        "cv": "Depth_CV",
        "label_mean": "Mean Depth",
        "label_sd": "Depth SD",
        "label_var": "Depth Variance",
        "label_cv": "Depth CV",
    },
    {
        "tag": "reldepth",
        "mean": "RelativeDepth_Mean",
        "sd": "RelativeDepth_SD",
        "var": "RelativeDepth_Var",
        "cv": "RelativeDepth_CV",
        "label_mean": "Relative Mean Depth",
        "label_sd": "Relative Depth SD",
        "label_var": "Relative Depth Variance",
        "label_cv": "Relative Depth CV",
    },
)


def _load_popdep_for_analysis(popdep_file: str) -> pd.DataFrame | None:
    if str(popdep_file).endswith(".popdep.info.tsv"):
        return load_df_from_tsv_popdep_info(popdep_file)
    return load_popdep_data(popdep_file)


def _add_log_metric_columns(df: pd.DataFrame, base_cols: list[str]) -> pd.DataFrame:
    for col in base_cols:
        if col in df.columns:
            df[f"Log_{col}"] = np.log(df[col].replace(0, np.nan))
    return df


def _nonnegative_depth_metric_floor(col: str, spec: dict) -> float | None:
    """Raw mean / SD / variance metrics are non-negative; log-scale columns are exempt."""
    if col.startswith("Log_"):
        return None
    if col in (spec["mean"], spec["sd"], spec["var"], spec["cv"]):
        return 0.0
    return None


def _scatter_axis_limits(
    x_col: str,
    y_col: str,
    spec: dict,
) -> tuple[tuple[float | None, float | None] | None, tuple[float | None, float | None] | None]:
    x_floor = _nonnegative_depth_metric_floor(x_col, spec)
    y_floor = _nonnegative_depth_metric_floor(y_col, spec)
    x_lim = (x_floor, None) if x_floor is not None else None
    y_lim = (y_floor, None) if y_floor is not None else None
    return x_lim, y_lim


def _distribution_xlim(series: pd.Series) -> tuple[float, float] | None:
    lims = padded_axis_limits(series)
    if lims is None:
        return (0.0, None)
    return 0.0, lims[1]


def _plot_popdep_missing_regressions(
    merged: pd.DataFrame,
    spec: dict,
    output_prefix: str,
) -> None:
    mean_col = spec["mean"]
    cv_col = spec["cv"]
    tag = spec["tag"]
    if mean_col not in merged.columns or cv_col not in merged.columns:
        print(f"[Warning] Skipping missing-reg plots for {tag}; columns absent.")
        return

    log_mean_col = f"Log_{mean_col}"
    log_cv_col = f"Log_{cv_col}"
    merged = _add_log_metric_columns(merged, [mean_col, cv_col, "F_MISS"])
    log_miss_col = "Log_F_MISS"
    merged[log_miss_col] = np.log(merged["F_MISS"].replace(0, np.nan))

    missing_specs = [
        (log_mean_col, "F_MISS", spec["label_mean"], "Missing Rate", (0.0, 1.0)),
        (mean_col, "F_MISS", spec["label_mean"], "Missing Rate", (0.0, 1.0)),
        (cv_col, "F_MISS", spec["label_cv"], "Missing Rate", (0.0, 1.0)),
        (log_mean_col, log_miss_col, f"Log {spec['label_mean']}", "Log Missing Rate", None),
        (log_cv_col, log_miss_col, f"Log {spec['label_cv']}", "Log Missing Rate", None),
    ]

    for x_col, y_col, x_label, y_label, y_lim in missing_specs:
        plot_df = merged
        x_lim = None
        x_floor = _nonnegative_depth_metric_floor(x_col, spec)
        if x_floor is not None:
            x_lim = (x_floor, None)
        if y_col == log_miss_col:
            plot_df = merged[(merged["F_MISS"] > 0) & (merged["F_MISS"] < 1)].copy()
            if len(plot_df) < 10:
                print(f"[Warning] Not enough F_MISS in (0,1) for {tag} {x_col}")
                continue
            y_lim = (None, 0.0)

        plot_regression_comparison(
            plot_df,
            x_col,
            y_col,
            x_label=x_label,
            y_label=y_label,
            filename=f"{output_prefix}.{tag}.reg_{y_col}_vs_{x_col}.png",
            title=f"{y_label} vs {x_label} ({tag})",
            x_lim=x_lim,
            y_lim=y_lim,
        )


def _adaptive_scatter_axis_limits(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    spec: dict,
    *,
    min_points: int = 10,
) -> tuple[tuple[float, float] | None, tuple[float, float] | None]:
    """Per-axis upper trim limits; floors at 0 for raw mean / SD / variance."""
    valid = df[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(valid) < min_points:
        return None, None

    x_lim = axis_limits_with_adaptive_upper(
        valid[x_col],
        min_points=min_points,
        lower_override=_nonnegative_depth_metric_floor(x_col, spec),
    )
    y_lim = axis_limits_with_adaptive_upper(
        valid[y_col],
        min_points=min_points,
        lower_override=_nonnegative_depth_metric_floor(y_col, spec),
    )
    return x_lim, y_lim


def _percentile_upper_scatter_limits(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    spec: dict,
    *,
    upper_pct: float = 99.9,
    min_points: int = 10,
    padding_fraction: float = DEFAULT_AXIS_PADDING_FRACTION,
) -> tuple[pd.DataFrame | None, tuple[float, float] | None, tuple[float, float] | None]:
    """
    Max-only percentile trim: drop rows above ``upper_pct`` on either axis.

    Lower axis bound is 0 for raw mean / SD / variance / CV; no lower-percentile
    row filter (unlike legacy trim99 which also excluded the bottom 1%).
    """
    valid = df[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(valid) < min_points:
        return None, None, None

    hi_q = upper_pct / 100.0
    x_hi = float(valid[x_col].quantile(hi_q))
    y_hi = float(valid[y_col].quantile(hi_q))
    keep = (valid[x_col] <= x_hi) & (valid[y_col] <= y_hi)
    trimmed = df.loc[valid.index[keep]]
    if len(trimmed) < min_points:
        return None, None, None

    def _axis_lim(series: pd.Series, col: str, hi_val: float) -> tuple[float, float]:
        lo_raw = float(series.min())
        floor = _nonnegative_depth_metric_floor(col, spec)
        span = hi_val - lo_raw
        if span <= 0:
            span = max(abs(lo_raw), abs(hi_val), 1.0) * 0.01
        if floor is not None:
            lo_axis = floor
        else:
            lo_axis = lo_raw - span * padding_fraction
        hi_axis = hi_val + span * padding_fraction
        return lo_axis, hi_axis

    valid_trim = trimmed[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    x_lim = _axis_lim(valid_trim[x_col], x_col, x_hi)
    y_lim = _axis_lim(valid_trim[y_col], y_col, y_hi)
    return trimmed, x_lim, y_lim


def _plot_percentile_upper_scatter(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    spec: dict,
    *,
    x_label: str,
    y_label: str,
    base_title: str,
    base_filename: str,
    upper_pct: float,
    suffix: str,
) -> None:
    """Write a max-only percentile scatter variant (axis floor 0, no lower row filter)."""
    trim_df, trim_x_lim, trim_y_lim = _percentile_upper_scatter_limits(
        df,
        x_col,
        y_col,
        spec,
        upper_pct=upper_pct,
    )
    if trim_df is None:
        return
    pct_label = (
        f"{int(upper_pct)}th" if upper_pct == int(upper_pct) else f"{upper_pct}th"
    )
    plot_regression_comparison(
        trim_df,
        x_col,
        y_col,
        x_label=x_label,
        y_label=y_label,
        filename=base_filename.replace(".png", f".{suffix}.png"),
        title=f"{base_title}, excluding above {pct_label} percentile on each axis",
        x_lim=trim_x_lim,
        y_lim=trim_y_lim,
    )


def _plot_popdep_depth_scatter(df: pd.DataFrame, spec: dict, output_prefix: str) -> None:
    mean_col = spec["mean"]
    sd_col = spec["sd"]
    var_col = spec["var"]
    cv_col = spec["cv"]
    tag = spec["tag"]
    if mean_col not in df.columns:
        print(f"[Warning] Skipping depth scatter for {tag}; {mean_col} absent.")
        return

    df = _add_log_metric_columns(df, [mean_col, sd_col, var_col, cv_col])
    scatter_specs = [
        (sd_col, mean_col, spec["label_sd"], spec["label_mean"]),
        (var_col, mean_col, spec["label_var"], spec["label_mean"]),
        (cv_col, mean_col, spec["label_cv"], spec["label_mean"]),
        (f"Log_{var_col}", f"Log_{mean_col}", f"Log {spec['label_var']}", f"Log {spec['label_mean']}"),
    ]

    for y_col, x_col, y_label, x_label in scatter_specs:
        if y_col not in df.columns or x_col not in df.columns:
            continue
        x_lim, y_lim = _scatter_axis_limits(x_col, y_col, spec)
        base_title = f"{y_label} vs {x_label} ({tag})"
        base_filename = f"{output_prefix}.{tag}.scatter_{y_col}_vs_{x_col}.png"

        plot_regression_comparison(
            df,
            x_col,
            y_col,
            x_label=x_label,
            y_label=y_label,
            filename=base_filename,
            title=base_title,
            x_lim=x_lim,
            y_lim=y_lim,
        )

        trim_x_lim, trim_y_lim = _adaptive_scatter_axis_limits(
            df,
            x_col,
            y_col,
            spec,
        )
        if trim_x_lim is not None or trim_y_lim is not None:
            plot_regression_comparison(
                df,
                x_col,
                y_col,
                x_label=x_label,
                y_label=y_label,
                filename=base_filename.replace(".png", ".trim.png"),
                title=f"{base_title}, adaptive upper-axis trim",
                x_lim=trim_x_lim if trim_x_lim is not None else x_lim,
                y_lim=trim_y_lim if trim_y_lim is not None else y_lim,
            )

        for upper_pct, suffix in ((99.0, "trim99"), (99.9, "trim999")):
            _plot_percentile_upper_scatter(
                df,
                x_col,
                y_col,
                spec,
                x_label=x_label,
                y_label=y_label,
                base_title=base_title,
                base_filename=base_filename,
                upper_pct=upper_pct,
                suffix=suffix,
            )


def _plot_popdep_depth_distributions(
    df: pd.DataFrame,
    spec: dict,
    output_prefix: str,
) -> None:
    mean_col = spec["mean"]
    sd_col = spec["sd"]
    var_col = spec["var"]
    tag = spec["tag"]
    if mean_col not in df.columns or var_col not in df.columns:
        print(f"[Warning] Skipping distributions for {tag}; columns absent.")
        return

    for col, label in (
        (mean_col, spec["label_mean"]),
        (sd_col, spec["label_sd"]),
        (var_col, spec["label_var"]),
    ):
        if col not in df.columns:
            continue
        series = df[col].replace([np.inf, -np.inf], np.nan).dropna()
        if series.empty:
            print(f"[Warning] No data for {tag} distribution: {col}")
            continue
        mean_val = float(series.mean())
        median_val = float(series.median())
        std_val = float(series.std())
        xlim = _distribution_xlim(series)

        plot_distribution_with_stats(
            data=df,
            col=col,
            title=f"{label} Distribution ({tag})",
            filename=f"{output_prefix}.{tag}.dist_{col}.png",
            mean_val=mean_val,
            median_val=median_val,
            std_val=std_val,
            x_label=label,
            log_scale=False,
            xlim=xlim,
        )
        plot_distribution_with_stats(
            data=df,
            col=col,
            title=f"{label} Distribution ({tag}, log count)",
            filename=f"{output_prefix}.{tag}.dist_{col}_logy.png",
            mean_val=mean_val,
            median_val=median_val,
            std_val=std_val,
            x_label=label,
            log_scale=True,
            xlim=xlim,
        )


def ana_popdep_qc_plots(popdep_file, vmiss_file, output_prefix):
    """
    Extended popdep QC plots: depth vs missing regressions, mean-vs-dispersion scatter,
    and mean/variance distributions — for both absolute and relative depth metrics.
    """
    df_dep = _load_popdep_for_analysis(popdep_file)
    df_miss = load_df_from_plink_variant(vmiss_file)
    if df_dep is None or df_miss is None:
        return

    merged = pd.merge(
        df_miss[["Position", "F_MISS"]],
        df_dep,
        on="Position",
        how="inner",
    )
    if merged.empty:
        print("[Error] No overlapping positions for popdep QC plots.")
        return

    for spec in _POPDEP_QC_METRIC_SPECS:
        _plot_popdep_missing_regressions(merged, spec, output_prefix)
        _plot_popdep_depth_scatter(df_dep, spec, output_prefix)
        _plot_popdep_depth_distributions(df_dep, spec, output_prefix)


# ==================================================================================
# Leave-one-out influence (compare PopDep grids vs baseline)
# ==================================================================================

POPDEP_LOO_METRIC_COLS = (
    "Depth_Mean",
    "Depth_SD",
    "RelativeDepth_Mean",
    "RelativeDepth_SD",
)

POPDEP_LOO_ABSOLUTE_METRIC_COLS = ("Depth_Mean", "Depth_SD")
POPDEP_LOO_PATTERN_METRIC_COLS = ("RelativeDepth_Mean", "RelativeDepth_SD")
PATTERN_INFLUENCE_PCT_WEIGHT = 0.01


def append_influence_components(stats: dict[str, float | int]) -> dict[str, float | int]:
    """
    Split LOO influence into absolute depth scale vs relative coverage pattern.

    ``influence_absolute`` = RMSE(Depth_Mean) + RMSE(Depth_SD)
    ``influence_pattern`` = RMSE(RelMean) + RMSE(RelSD) + 0.01 × pct sites rel>1%
    ``influence_score`` = sum of both (legacy combined index)
    """
    abs_rmse = sum(float(stats.get(f"rmse_{c}", 0.0) or 0.0) for c in POPDEP_LOO_ABSOLUTE_METRIC_COLS)
    pat_rmse = sum(float(stats.get(f"rmse_{c}", 0.0) or 0.0) for c in POPDEP_LOO_PATTERN_METRIC_COLS)
    pct = float(stats.get("pct_any_metric_rel_gt_1pct", 0.0) or 0.0)
    pattern = pat_rmse + PATTERN_INFLUENCE_PCT_WEIGHT * pct
    stats["influence_absolute"] = float(abs_rmse)
    stats["influence_pattern"] = float(pattern)
    stats["influence_score"] = float(abs_rmse + pattern)
    return stats


def ensure_influence_components_df(df: pd.DataFrame) -> pd.DataFrame:
    """Backfill split influence columns from RMSE fields when missing."""
    if df is None or df.empty:
        return df
    out = df.copy()
    rmse_cols = [f"rmse_{c}" for c in POPDEP_LOO_METRIC_COLS]
    if not all(c in out.columns for c in rmse_cols):
        return out
    if "influence_absolute" in out.columns and "influence_pattern" in out.columns:
        return out
    rows = []
    for _, row in out.iterrows():
        stats = append_influence_components(row.to_dict())
        rows.append(
            {
                "influence_absolute": stats["influence_absolute"],
                "influence_pattern": stats["influence_pattern"],
                "influence_score": stats["influence_score"],
            }
        )
    comp = pd.DataFrame(rows, index=out.index)
    for col in comp.columns:
        out[col] = comp[col]
    return out


def load_popdep_grid(path: str) -> pd.DataFrame | None:
    """Load a TIGER popdepth grid (plain, gzip, or BGZF) for LOO comparison."""
    df = load_popdep_data(path)
    if df is None:
        return None
    if "Position" not in df.columns:
        print(f"[Error] Popdep grid missing Position column: {path}")
        return None
    keep = ["Position", *POPDEP_LOO_METRIC_COLS]
    missing = [c for c in keep if c not in df.columns]
    if missing:
        print(f"[Error] Popdep grid missing columns {missing}: {path}")
        return None
    out = df[keep].copy()
    for col in POPDEP_LOO_METRIC_COLS:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def _resolve_popdep_grid_path(case_dir: str, chr_id: str = "chr001") -> Path | None:
    """Prefer BGZF publish layout ``variant/chrNNN.popdep.txt.bgz`` under *case_dir*."""
    variant = Path(case_dir) / "variant"
    for name in (
        f"{chr_id}.popdep.txt.bgz",
        f"{chr_id}.popdep.txt.gz",
        f"{chr_id}.popdep.txt",
    ):
        candidate = variant / name
        if candidate.exists():
            return candidate
    return None


def compare_popdep_to_baseline(
    baseline_df: pd.DataFrame,
    loo_df: pd.DataFrame,
    *,
    min_baseline_mean: float = 0.0,
) -> dict[str, float | int]:
    """
    Compare one leave-one-out grid to the baseline on overlapping positions.

    Returns per-metric max/mean abs diff, RMSE, and fraction of sites with
    relative change above 1% on any metric (among sites with baseline mean > 0).
    """
    merged = pd.merge(
        baseline_df,
        loo_df,
        on="Position",
        how="inner",
        suffixes=("_base", "_loo"),
    )
    if merged.empty:
        return {"n_positions": 0}

    if min_baseline_mean > 0:
        mask = merged["Depth_Mean_base"] > min_baseline_mean
        merged = merged.loc[mask].copy()
    n_pos = len(merged)
    if n_pos == 0:
        return {"n_positions": 0}

    stats: dict[str, float | int] = {"n_positions": n_pos}
    rel_change_any = np.zeros(n_pos, dtype=bool)

    for col in POPDEP_LOO_METRIC_COLS:
        base = merged[f"{col}_base"].to_numpy(dtype=float)
        loo = merged[f"{col}_loo"].to_numpy(dtype=float)
        diff = loo - base
        abs_diff = np.abs(diff)
        stats[f"max_abs_{col}"] = float(np.nanmax(abs_diff))
        stats[f"mean_abs_{col}"] = float(np.nanmean(abs_diff))
        stats[f"rmse_{col}"] = float(np.sqrt(np.nanmean(diff ** 2)))
        denom = np.maximum(np.abs(base), 1e-12)
        rel = abs_diff / denom
        stats[f"max_rel_{col}"] = float(np.nanmax(rel))
        stats[f"mean_rel_{col}"] = float(np.nanmean(rel))
        stats[f"pct_rel_gt_1pct_{col}"] = float(np.nanmean(rel > 0.01) * 100.0)
        rel_change_any |= rel > 0.01

    stats["pct_any_metric_rel_gt_1pct"] = float(np.mean(rel_change_any) * 100.0)
    append_influence_components(stats)
    return stats


def summarize_popdep_leave_one_out(
    baseline_grid: str,
    loo_root: str,
    output_prefix: str,
    *,
    chr_id: str = "chr001",
    min_baseline_mean: float = 0.0,
) -> pd.DataFrame | None:
    """
    Rank leave-one-out PopDep cases by influence vs a baseline grid.

    Expects ``loo_root/loo_<taxa>/variant/chrNNN.popdep.txt.bgz`` (or .gz/.txt).
    Writes ``{output_prefix}.info.tsv`` sorted by ``influence_score`` descending.
    """
    baseline_df = load_popdep_grid(baseline_grid)
    if baseline_df is None:
        return None

    loo_root_path = Path(loo_root)
    rows: list[dict] = []
    for case_dir in sorted(loo_root_path.glob("loo_*")):
        if not case_dir.is_dir():
            continue
        left_taxa = case_dir.name.replace("loo_", "", 1)
        grid_path = _resolve_popdep_grid_path(str(case_dir), chr_id=chr_id)
        if grid_path is None:
            print(f"[Warning] No popdep grid under {case_dir}")
            continue
        loo_df = load_popdep_grid(str(grid_path))
        if loo_df is None:
            continue
        stats = compare_popdep_to_baseline(
            baseline_df,
            loo_df,
            min_baseline_mean=min_baseline_mean,
        )
        if stats.get("n_positions", 0) == 0:
            print(f"[Warning] No overlapping positions for {left_taxa}")
            continue
        stats["left_out_taxa"] = left_taxa
        stats["loo_grid"] = str(grid_path)
        rows.append(stats)

    if not rows:
        print("[Error] No leave-one-out cases compared.")
        return None

    summary = pd.DataFrame(rows).sort_values(
        "influence_score",
        ascending=False,
        kind="mergesort",
    )
    out_path = f"{output_prefix}.info.tsv"
    save_df_to_tsv(summary, out_path)
    print(f"[Info] Wrote LOO influence summary ({len(summary)} cases) to {out_path}")
    top = summary.head(10)
    for _, row in top.iterrows():
        print(
            f"  {row['left_out_taxa']}: influence_score={row['influence_score']:.6f} "
            f"pct_any_rel>1%={row['pct_any_metric_rel_gt_1pct']:.2f}"
        )
    return summary


def _influence_lookup_maps(chr1: pd.DataFrame) -> dict[str, dict[str, float]]:
    """Build per-sample influence component maps from a chr1 LOO info table."""
    df = ensure_influence_components_df(chr1)
    maps: dict[str, dict[str, float]] = {
        "absolute": df.set_index("left_out_taxa")["influence_absolute"].to_dict(),
        "pattern": df.set_index("left_out_taxa")["influence_pattern"].to_dict(),
        "combined": df.set_index("left_out_taxa")["influence_score"].to_dict(),
    }
    return maps


def _safe_ratio(num: float, denom: float) -> float:
    if pd.notna(num) and pd.notna(denom) and denom > 0:
        return float(num / denom)
    return float("nan")


def compare_aneuploidy_samechr_loo_influence(
    manifest_tsv: str,
    chr1_info_tsv: str,
    samechr_publish_root: str,
    output_prefix: str,
    *,
    bench_max_bp: int = 10_000_000,
) -> pd.DataFrame | None:
    """
    Compare LOO influence on each sample's aneuploidy-matched chromosome vs chr1.

    Reports split metrics: ``influence_absolute`` (depth scale) and
    ``influence_pattern`` (relative coverage + pct-site bonus). Pattern ratio is
    the primary same-chr vs chr1 comparison for anomaly screening.
    """
    manifest = pd.read_csv(manifest_tsv, sep="\t")
    if manifest.empty:
        print("[Error] Empty aneuploidy same-chr LOO manifest.")
        return None

    chr1 = ensure_influence_components_df(pd.read_csv(chr1_info_tsv, sep="\t"))
    chr1_maps = _influence_lookup_maps(chr1)

    publish_root = Path(samechr_publish_root)
    rows: list[dict] = []
    for _, entry in manifest.iterrows():
        sample = str(entry["sample"])
        popdep_chr = int(entry["popdep_chr"])
        prefix = publish_root / f"chr{popdep_chr}" / "stats" / (
            f"n500_aneuploidy_samechr_chr{popdep_chr}_{bench_max_bp}_loo"
        )
        info_path = Path(f"{prefix}.info.tsv")
        if not info_path.exists():
            alt = publish_root / "stats" / (
                f"n500_aneuploidy_samechr_chr{popdep_chr}_{bench_max_bp}_loo.info.tsv"
            )
            info_path = alt if alt.exists() else info_path

        abs_same = pat_same = comb_same = float("nan")
        if info_path.exists():
            same_df = ensure_influence_components_df(pd.read_csv(info_path, sep="\t"))
            hit = same_df.loc[same_df["left_out_taxa"] == sample]
            if not hit.empty:
                row = hit.iloc[0]
                abs_same = float(row["influence_absolute"])
                pat_same = float(row["influence_pattern"])
                comb_same = float(row["influence_score"])

        abs_chr1 = float(chr1_maps["absolute"].get(sample, float("nan")))
        pat_chr1 = float(chr1_maps["pattern"].get(sample, float("nan")))
        comb_chr1 = float(chr1_maps["combined"].get(sample, float("nan")))

        rows.append(
            {
                "sample": sample,
                "popdep_chr": popdep_chr,
                "abnormal_chr": entry.get("abnormal_chr", ""),
                "influence_absolute_chr1": abs_chr1,
                "influence_absolute_samechr": abs_same,
                "influence_absolute_ratio_samechr_over_chr1": _safe_ratio(abs_same, abs_chr1),
                "influence_absolute_delta_samechr_minus_chr1": (
                    abs_same - abs_chr1
                    if pd.notna(abs_same) and pd.notna(abs_chr1)
                    else float("nan")
                ),
                "influence_pattern_chr1": pat_chr1,
                "influence_pattern_samechr": pat_same,
                "influence_pattern_ratio_samechr_over_chr1": _safe_ratio(pat_same, pat_chr1),
                "influence_pattern_delta_samechr_minus_chr1": (
                    pat_same - pat_chr1
                    if pd.notna(pat_same) and pd.notna(pat_chr1)
                    else float("nan")
                ),
                "influence_chr1": comb_chr1,
                "influence_samechr": comb_same,
                "influence_ratio_samechr_over_chr1": _safe_ratio(comb_same, comb_chr1),
                "influence_delta_samechr_minus_chr1": (
                    comb_same - comb_chr1
                    if pd.notna(comb_same) and pd.notna(comb_chr1)
                    else float("nan")
                ),
                "samechr_info_tsv": str(info_path) if info_path.exists() else "",
            }
        )

    out = pd.DataFrame(rows).sort_values(
        "influence_pattern_ratio_samechr_over_chr1",
        ascending=False,
        kind="mergesort",
    )
    save_df_to_tsv(out, f"{output_prefix}.info.tsv")
    print(f"[Info] Wrote aneuploidy same-chr LOO comparison ({len(out)} rows) to {output_prefix}.info.tsv")
    for _, row in out.iterrows():
        print(
            f"  {row['sample']} chr{int(row['popdep_chr'])}: "
            f"pattern same={row['influence_pattern_samechr']:.4f} chr1={row['influence_pattern_chr1']:.4f} "
            f"ratio={row['influence_pattern_ratio_samechr_over_chr1']:.2f} | "
            f"absolute same={row['influence_absolute_samechr']:.4f} chr1={row['influence_absolute_chr1']:.4f}"
        )
    return out


def _plot_aneuploidy_samechr_loo_single_case(
    info_tsv: str,
    output_prefix: str,
    *,
    sample_label: str | None = None,
) -> None:
    """RMSE breakdown and split influence for a one-row same-chr LOO table."""
    df = ensure_influence_components_df(load_df_from_tsv(info_tsv))
    if df is None or df.empty:
        print(f"[Warning] Skipping single-case plot — empty {info_tsv}")
        return
    row = df.iloc[0]
    taxa = sample_label or str(row.get("left_out_taxa", "sample"))
    rmse_cols = [f"rmse_{c}" for c in POPDEP_LOO_METRIC_COLS]
    pct_cols = [f"pct_rel_gt_1pct_{c}" for c in POPDEP_LOO_METRIC_COLS]

    abs_vals = {
        c.replace("rmse_", ""): float(row[c])
        for c in rmse_cols
        if c in row.index and c.replace("rmse_", "") in POPDEP_LOO_ABSOLUTE_METRIC_COLS
    }
    pat_vals = {
        c.replace("rmse_", ""): float(row[c])
        for c in rmse_cols
        if c in row.index and c.replace("rmse_", "") in POPDEP_LOO_PATTERN_METRIC_COLS
    }
    infl_abs = float(row.get("influence_absolute", float("nan")))
    infl_pat = float(row.get("influence_pattern", float("nan")))
    pct_any = float(row.get("pct_any_metric_rel_gt_1pct", 0.0) or 0.0)

    plot_bar_chart(
        list(abs_vals.keys()),
        list(abs_vals.values()),
        title=f"Absolute RMSE — {taxa} (influence_absolute={infl_abs:.4f})",
        ylabel="RMSE",
        filename=f"{output_prefix}.single_case_absolute_rmse.png",
        ylim=None,
        color="coral",
        rotate_xlabels=15,
    )
    plot_bar_chart(
        list(pat_vals.keys()),
        list(pat_vals.values()),
        title=(
            f"Pattern RMSE — {taxa} (influence_pattern={infl_pat:.4f}, "
            f"pct sites rel>1%={pct_any:.1f}%)"
        ),
        ylabel="RMSE",
        filename=f"{output_prefix}.single_case_pattern_rmse.png",
        ylim=None,
        color="mediumpurple",
        rotate_xlabels=15,
    )
    plot_bar_chart(
        ["influence_absolute", "influence_pattern"],
        [infl_abs, infl_pat],
        title=f"Split influence — {taxa}",
        ylabel="Score",
        filename=f"{output_prefix}.single_case_split_influence.png",
        ylim=None,
        color=["coral", "mediumpurple"],
        rotate_xlabels=15,
    )


def _plot_chr1_vs_samechr_influence_metric(
    plot_df: pd.DataFrame,
    output_prefix: str,
    *,
    chr1_col: str,
    samechr_col: str,
    ratio_col: str,
    delta_col: str,
    metric_slug: str,
    metric_label: str,
) -> None:
    """Shared chr1 vs samechr comparison plots for one influence metric."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    from infra.utils.graph import (
        LEGEND_FONT_SIZE,
        TICK_FONT_SIZE,
        TITLE_FONT_SIZE,
        X_LABEL_FONT_SIZE,
        Y_LABEL_FONT_SIZE,
        plot_correlation_with_regression,
        plot_slope_chart,
    )

    df = plot_df.copy()
    df["label"] = df.apply(
        lambda r: f"{r['sample']} (chr{int(r['popdep_chr'])})", axis=1
    )
    df["higher_on_samechr"] = df[ratio_col] >= 1.0
    df["comparison_tier"] = np.where(
        df["higher_on_samechr"], "samechr >= chr1", "samechr < chr1"
    )

    fig_path = f"{output_prefix}.{metric_slug}_chr1_vs_samechr_scatter.png"
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(9, 8))
    for tier, sub in df.groupby("comparison_tier"):
        color = "crimson" if tier.startswith("samechr >=") else "steelblue"
        ax.scatter(
            sub[chr1_col],
            sub[samechr_col],
            c=color,
            alpha=0.75,
            s=80,
            edgecolors="w",
            label=tier,
        )
    lim_hi = max(df[chr1_col].max(), df[samechr_col].max()) * 1.08
    ax.plot([0, lim_hi], [0, lim_hi], "k--", alpha=0.5, label="y = x")
    for _, row in df.iterrows():
        ax.annotate(
            row["sample"],
            (row[chr1_col], row[samechr_col]),
            fontsize=8,
            xytext=(4, 4),
            textcoords="offset points",
        )
    ax.set_xlim(0, lim_hi)
    ax.set_ylim(0, lim_hi)
    ax.set_xlabel(f"{metric_label} — chr1 LOO (10 Mb)", fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(f"{metric_label} — same-chr LOO (10 Mb)", fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(
        f"Aneuploidy LOO {metric_label}: chr1 vs same abnormal chromosome",
        fontsize=TITLE_FONT_SIZE,
    )
    ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.12), fontsize=LEGEND_FONT_SIZE)
    ax.grid(axis="both", linestyle="--", alpha=0.35)
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    plot_correlation_with_regression(
        df,
        x_col=chr1_col,
        y_col=samechr_col,
        title=f"Aneuploidy LOO {metric_label}: chr1 vs same-chr (regression)",
        filename=f"{output_prefix}.{metric_slug}_chr1_vs_samechr_regression.png",
        x_label=f"{metric_label} — chr1 LOO",
        y_label=f"{metric_label} — same-chr LOO",
        group_col="comparison_tier",
        figure_size=(9, 7),
    )

    slope_df = df.rename(columns={chr1_col: "chr1_val", samechr_col: "samechr_val"})
    plot_slope_chart(
        slope_df,
        col1="chr1_val",
        col2="samechr_val",
        xlabel1="chr1 LOO",
        xlabel2="Same abnormal chr LOO",
        title=f"Per-sample {metric_label}: chr1 → same abnormal chromosome",
        filename=f"{output_prefix}.{metric_slug}_chr1_to_samechr_slope.png",
        group_col="comparison_tier",
        ylabel=metric_label,
        line_alpha=0.55,
        figure_size=(7, 9),
    )

    ratio_sorted = df.sort_values(ratio_col, ascending=True, kind="mergesort")
    ratio_colors = [
        "crimson" if r >= 1.0 else "steelblue" for r in ratio_sorted[ratio_col]
    ]
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.barh(ratio_sorted["label"], ratio_sorted[ratio_col], color=ratio_colors, alpha=0.85)
    ax.axvline(1.0, color="black", linestyle="--", linewidth=1.2, label="ratio = 1")
    ax.set_xlabel(f"{metric_label} ratio (samechr / chr1)", fontsize=X_LABEL_FONT_SIZE)
    ax.set_title(
        f"Same-chr vs chr1 LOO {metric_label} ratio (aneuploidy samples)",
        fontsize=TITLE_FONT_SIZE,
    )
    ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax.legend(fontsize=LEGEND_FONT_SIZE)
    fig.savefig(
        f"{output_prefix}.{metric_slug}_ratio_samechr_over_chr1.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    delta_sorted = df.sort_values(delta_col, ascending=True, kind="mergesort")
    delta_colors = [
        "crimson" if d >= 0 else "steelblue" for d in delta_sorted[delta_col]
    ]
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.barh(delta_sorted["label"], delta_sorted[delta_col], color=delta_colors, alpha=0.85)
    ax.axvline(0.0, color="black", linestyle="--", linewidth=1.2, label="delta = 0")
    ax.set_xlabel(f"{metric_label} delta (samechr − chr1)", fontsize=X_LABEL_FONT_SIZE)
    ax.set_title(
        f"Same-chr minus chr1 LOO {metric_label} (aneuploidy samples)",
        fontsize=TITLE_FONT_SIZE,
    )
    ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax.legend(fontsize=LEGEND_FONT_SIZE)
    fig.savefig(
        f"{output_prefix}.{metric_slug}_delta_samechr_minus_chr1.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, max(5, 0.45 * len(df))))
    y_pos = np.arange(len(df))
    bar_h = 0.35
    ordered = df.sort_values(chr1_col, ascending=False, kind="mergesort")
    ax.barh(
        y_pos + bar_h / 2,
        ordered[chr1_col],
        height=bar_h,
        label="chr1 LOO",
        color="steelblue",
        alpha=0.85,
    )
    ax.barh(
        y_pos - bar_h / 2,
        ordered[samechr_col],
        height=bar_h,
        label="same abnormal chr LOO",
        color="coral",
        alpha=0.85,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(ordered["label"], fontsize=TICK_FONT_SIZE)
    ax.set_xlabel(metric_label, fontsize=X_LABEL_FONT_SIZE)
    ax.set_title(
        f"Horizontal comparison: chr1 vs same-chr {metric_label} per aneuploidy sample",
        fontsize=TITLE_FONT_SIZE,
    )
    ax.legend(fontsize=LEGEND_FONT_SIZE)
    fig.savefig(
        f"{output_prefix}.{metric_slug}_chr1_vs_samechr_grouped_bars.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_aneuploidy_samechr_loo_comparison(
    compare_info_tsv: str,
    output_prefix: str,
    *,
    chr1_info_tsv: str | None = None,
) -> pd.DataFrame | None:
    """
    Horizontal comparison: same aneuploidy sample LOO on matched chr vs chr1.

    Primary screening uses ``influence_pattern``; ``influence_absolute`` is plotted
    separately as depth-scale leverage.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    from infra.utils.graph import (
        LEGEND_FONT_SIZE,
        TICK_FONT_SIZE,
        TITLE_FONT_SIZE,
        X_LABEL_FONT_SIZE,
        Y_LABEL_FONT_SIZE,
        plot_correlation_with_regression,
    )

    df = load_df_from_tsv(compare_info_tsv)
    if df is None or df.empty:
        print(f"[Error] Empty aneuploidy same-chr comparison table: {compare_info_tsv}")
        return None

    # Backfill split columns from legacy combined-only tables
    if "influence_pattern_chr1" not in df.columns:
        print("[Warning] Compare table lacks split columns; cannot plot pattern vs absolute.")
        return None

    plot_df = df.copy()
    plot_df["label"] = plot_df.apply(
        lambda r: f"{r['sample']} (chr{int(r['popdep_chr'])})", axis=1
    )

    for slug, chr1_col, same_col, ratio_col, delta_col, label in (
        (
            "pattern",
            "influence_pattern_chr1",
            "influence_pattern_samechr",
            "influence_pattern_ratio_samechr_over_chr1",
            "influence_pattern_delta_samechr_minus_chr1",
            "Pattern influence",
        ),
        (
            "absolute",
            "influence_absolute_chr1",
            "influence_absolute_samechr",
            "influence_absolute_ratio_samechr_over_chr1",
            "influence_absolute_delta_samechr_minus_chr1",
            "Absolute influence",
        ),
    ):
        _plot_chr1_vs_samechr_influence_metric(
            plot_df,
            output_prefix,
            chr1_col=chr1_col,
            samechr_col=same_col,
            ratio_col=ratio_col,
            delta_col=delta_col,
            metric_slug=slug,
            metric_label=label,
        )

    # Absolute vs pattern on chr1 for aneuploidy samples (scale vs profile)
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.scatter(
        plot_df["influence_absolute_chr1"],
        plot_df["influence_pattern_chr1"],
        c="steelblue",
        s=90,
        alpha=0.8,
        edgecolors="w",
    )
    for _, row in plot_df.iterrows():
        ax.annotate(
            row["sample"],
            (row["influence_absolute_chr1"], row["influence_pattern_chr1"]),
            fontsize=8,
            xytext=(4, 4),
            textcoords="offset points",
        )
    ax.set_xlabel("Absolute influence — chr1 LOO", fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel("Pattern influence — chr1 LOO", fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title(
        "Aneuploidy samples: depth scale vs coverage pattern (chr1 LOO)",
        fontsize=TITLE_FONT_SIZE,
    )
    ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax.grid(axis="both", linestyle="--", alpha=0.35)
    fig.savefig(
        f"{output_prefix}.absolute_vs_pattern_chr1_scatter.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    plot_correlation_with_regression(
        plot_df,
        x_col="influence_absolute_chr1",
        y_col="influence_pattern_chr1",
        title="Aneuploidy chr1 LOO: absolute vs pattern influence",
        filename=f"{output_prefix}.absolute_vs_pattern_chr1_regression.png",
        x_label="Absolute influence — chr1 LOO",
        y_label="Pattern influence — chr1 LOO",
        figure_size=(9, 7),
    )

    n_pat_higher = int((plot_df["influence_pattern_ratio_samechr_over_chr1"] >= 1.0).sum())
    n_pat_lower = int((plot_df["influence_pattern_ratio_samechr_over_chr1"] < 1.0).sum())
    n_abs_higher = int((plot_df["influence_absolute_ratio_samechr_over_chr1"] >= 1.0).sum())
    summary = pd.DataFrame(
        [
            {
                "n_cases": len(plot_df),
                "n_pattern_higher_on_samechr": n_pat_higher,
                "n_pattern_lower_on_samechr": n_pat_lower,
                "median_pattern_ratio": float(
                    plot_df["influence_pattern_ratio_samechr_over_chr1"].median()
                ),
                "max_pattern_ratio_sample": str(
                    plot_df.loc[
                        plot_df["influence_pattern_ratio_samechr_over_chr1"].idxmax(),
                        "sample",
                    ]
                ),
                "max_pattern_ratio": float(
                    plot_df["influence_pattern_ratio_samechr_over_chr1"].max()
                ),
                "n_absolute_higher_on_samechr": n_abs_higher,
                "median_absolute_ratio": float(
                    plot_df["influence_absolute_ratio_samechr_over_chr1"].median()
                ),
            }
        ]
    )
    save_df_to_tsv(summary, f"{output_prefix}.comparison_summary.info.tsv")

    if chr1_info_tsv and Path(chr1_info_tsv).exists():
        chr1_df = ensure_influence_components_df(load_df_from_tsv(chr1_info_tsv))
        if chr1_df is not None and not chr1_df.empty:
            for col, slug, title, chr1_key_col, same_key_col in (
                (
                    "influence_pattern",
                    "pattern",
                    "Pattern influence — chr1 LOO distribution",
                    "influence_pattern_chr1",
                    "influence_pattern_samechr",
                ),
                (
                    "influence_absolute",
                    "absolute",
                    "Absolute influence — chr1 LOO distribution",
                    "influence_absolute_chr1",
                    "influence_absolute_samechr",
                ),
            ):
                sns.set_style("white")
                fig, ax = plt.subplots(figsize=(10, 5))
                ax.hist(chr1_df[col], bins=40, color="lightgray", alpha=0.85)
                for _, row in plot_df.iterrows():
                    ax.axvline(
                        row[chr1_key_col],
                        color="steelblue",
                        linestyle=":",
                        alpha=0.5,
                        linewidth=1,
                    )
                    ax.axvline(
                        row[same_key_col],
                        color="crimson",
                        linestyle="-",
                        alpha=0.7,
                        linewidth=1.2,
                    )
                ax.set_xlabel(col.replace("_", " "), fontsize=X_LABEL_FONT_SIZE)
                ax.set_ylabel("Count", fontsize=Y_LABEL_FONT_SIZE)
                ax.set_title(
                    f"{title} with aneuploidy chr1 (blue) and same-chr (red) values",
                    fontsize=TITLE_FONT_SIZE,
                )
                fig.savefig(
                    f"{output_prefix}.chr1_loo_{slug}_distribution_with_aneuploidy.png",
                    dpi=300,
                    bbox_inches="tight",
                )
                plt.close(fig)

    print(
        f"[Info] Split comparison: pattern {n_pat_higher}/{len(plot_df)} higher on same chr; "
        f"absolute {n_abs_higher}/{len(plot_df)} higher on same chr"
    )
    return summary


def plot_aneuploidy_samechr_loo_batch(
    publish_root: str,
    compare_info_tsv: str,
    output_prefix: str,
    *,
    bench_max_bp: int = 10_000_000,
    chr1_info_tsv: str | None = None,
    per_chr_plot_root: str | Path | None = None,
    manifest_tsv: str | None = None,
) -> None:
    """
    Per-chr single-case breakdowns plus cross-sample chr1 vs samechr comparison plots.

    When ``manifest_tsv`` is set, regenerates the compare TSV with split influence
    columns before plotting.
    """
    compare_path = compare_info_tsv
    if manifest_tsv and chr1_info_tsv:
        compare_aneuploidy_samechr_loo_influence(
            manifest_tsv,
            chr1_info_tsv,
            publish_root,
            output_prefix,
            bench_max_bp=bench_max_bp,
        )
        compare_path = f"{output_prefix}.info.tsv"

    read_root = Path(publish_root)
    write_root = Path(per_chr_plot_root) if per_chr_plot_root else Path(".")
    write_root.mkdir(parents=True, exist_ok=True)

    for chr_dir in sorted(read_root.glob("chr[0-9]*")):
        chr_num = chr_dir.name.replace("chr", "")
        info_path = chr_dir / "stats" / (
            f"n500_aneuploidy_samechr_chr{chr_num}_{bench_max_bp}_loo.info.tsv"
        )
        if not info_path.exists():
            print(f"[Warning] Missing per-chr info: {info_path}")
            continue
        df = load_df_from_tsv(str(info_path))
        if df is None or df.empty:
            continue
        chr_plot_dir = write_root / f"chr{chr_num}" / "plots"
        chr_plot_dir.mkdir(parents=True, exist_ok=True)
        out_prefix = chr_plot_dir / (
            f"n500_aneuploidy_samechr_chr{chr_num}_{bench_max_bp}_loo"
        )
        if len(df) == 1:
            _plot_aneuploidy_samechr_loo_single_case(
                str(info_path),
                str(out_prefix),
                sample_label=str(df.iloc[0].get("left_out_taxa", "")),
            )
        else:
            plot_popdep_loo_influence_summary(str(info_path), str(out_prefix))

    plot_aneuploidy_samechr_loo_comparison(
        compare_path,
        output_prefix,
        chr1_info_tsv=chr1_info_tsv,
    )


def plot_popdep_loo_influence_summary(
    info_tsv: str,
    output_prefix: str,
    *,
    top_n: int = 15,
) -> pd.DataFrame | None:
    """
    Plot leave-one-out PopDep influence metrics from an existing ``*.info.tsv``.

    Writes PNG figures and a metric-decomposition TSV under ``{output_prefix}.*``.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    from infra.utils.graph import (
        LEGEND_FONT_SIZE,
        TITLE_FONT_SIZE,
        plot_bar_chart,
        plot_distribution_with_stats,
    )

    df = ensure_influence_components_df(load_df_from_tsv(info_tsv))
    if df is None or df.empty:
        print(f"[Error] Empty or missing LOO influence table: {info_tsv}")
        return None
    if "influence_score" not in df.columns:
        print(f"[Error] Missing influence_score column in {info_tsv}")
        return None

    rmse_cols = [f"rmse_{c}" for c in POPDEP_LOO_METRIC_COLS]
    pct_cols = [f"pct_rel_gt_1pct_{c}" for c in POPDEP_LOO_METRIC_COLS]
    missing = [c for c in rmse_cols if c not in df.columns]
    if missing:
        print(f"[Error] LOO influence table missing columns {missing}")
        return None

    df = df.copy()
    df["rmse_sum"] = df[rmse_cols].sum(axis=1)
    df["pct_bonus"] = PATTERN_INFLUENCE_PCT_WEIGHT * df["pct_any_metric_rel_gt_1pct"]

    for col, title, color in (
        (
            "influence_absolute",
            "PopDep LOO absolute influence (depth scale)",
            "coral",
        ),
        (
            "influence_pattern",
            "PopDep LOO pattern influence (relative coverage + pct bonus)",
            "mediumpurple",
        ),
    ):
        plot_distribution_with_stats(
            df,
            col=col,
            title=f"{title} — distribution",
            filename=f"{output_prefix}.{col}_distribution.png",
            x_label=title,
            y_label="Count",
            bins=40,
            color=color,
        )

    from infra.utils.graph import plot_correlation_with_regression

    plot_correlation_with_regression(
        df,
        x_col="influence_absolute",
        y_col="influence_pattern",
        title="LOO absolute vs pattern influence (all samples)",
        filename=f"{output_prefix}.influence_absolute_vs_pattern_scatter.png",
        x_label="Absolute influence (RMSE Depth_Mean + Depth_SD)",
        y_label="Pattern influence (RMSE Rel + 0.01 × pct sites rel>1%)",
        figure_size=(10, 8),
    )

    plot_distribution_with_stats(
        df,
        col="influence_score",
        title="PopDep LOO influence score distribution",
        filename=f"{output_prefix}.influence_score_distribution.png",
        x_label="Influence score (sum of 4 RMSE + 0.01 × pct sites rel>1%)",
        y_label="Count",
        bins=40,
        color="steelblue",
    )
    plot_distribution_with_stats(
        df,
        col="influence_score",
        title="PopDep LOO influence score distribution (log y)",
        filename=f"{output_prefix}.influence_score_distribution_logy.png",
        x_label="Influence score",
        y_label="Count (log scale)",
        bins=40,
        log_scale=True,
        color="steelblue",
    )

    mean_rmse = {c.replace("rmse_", ""): float(df[c].mean()) for c in rmse_cols}
    plot_bar_chart(
        list(mean_rmse.keys()),
        list(mean_rmse.values()),
        title="Mean RMSE contribution to influence score",
        ylabel="Mean RMSE",
        filename=f"{output_prefix}.mean_rmse_contribution.png",
        ylim=None,
        color="coral",
        rotate_xlabels=15,
    )

    top_rank = min(top_n, len(df))
    top = df.head(top_rank).copy()
    top["dominant_metric"] = top[rmse_cols].idxmax(axis=1).str.replace("rmse_", "")
    dom_counts = top["dominant_metric"].value_counts()
    plot_bar_chart(
        dom_counts.index.tolist(),
        dom_counts.values.tolist(),
        title=f"Dominant RMSE metric among top-{top_rank} influential taxa",
        ylabel="Count",
        filename=f"{output_prefix}.top{top_rank}_dominant_metric.png",
        ylim=(0, float(dom_counts.max()) + 1),
        color="mediumpurple",
    )

    top_bar = df.head(min(15, len(df)))
    plot_bar_chart(
        top_bar["left_out_taxa"].tolist(),
        top_bar["influence_score"].tolist(),
        title="Top influential taxa — combined score (legacy)",
        ylabel="Influence score",
        filename=f"{output_prefix}.top15_influence_score.png",
        ylim=None,
        color="indianred",
        rotate_xlabels=60,
        figure_size=(12, 5),
    )

    top_pat = df.sort_values("influence_pattern", ascending=False, kind="mergesort").head(
        min(15, len(df))
    )
    plot_bar_chart(
        top_pat["left_out_taxa"].tolist(),
        top_pat["influence_pattern"].tolist(),
        title="Top influential taxa — pattern influence (screening)",
        ylabel="Pattern influence",
        filename=f"{output_prefix}.top15_influence_pattern.png",
        ylim=None,
        color="mediumpurple",
        rotate_xlabels=60,
        figure_size=(12, 5),
    )

    top_abs = df.sort_values("influence_absolute", ascending=False, kind="mergesort").head(
        min(15, len(df))
    )
    plot_bar_chart(
        top_abs["left_out_taxa"].tolist(),
        top_abs["influence_absolute"].tolist(),
        title="Top influential taxa — absolute depth influence",
        ylabel="Absolute influence",
        filename=f"{output_prefix}.top15_influence_absolute.png",
        ylim=None,
        color="coral",
        rotate_xlabels=60,
        figure_size=(12, 5),
    )

    pct_mean = {
        c.replace("pct_rel_gt_1pct_", ""): float(df[c].mean())
        for c in pct_cols
        if c in df.columns
    }
    if pct_mean:
        plot_bar_chart(
            list(pct_mean.keys()),
            list(pct_mean.values()),
            title="Mean % sites with >1% relative change per metric",
            ylabel="% positions",
            filename=f"{output_prefix}.mean_pct_sites_rel_gt_1pct.png",
            ylim=(0, max(pct_mean.values()) * 1.15),
            color="seagreen",
        )

    corr_cols = rmse_cols + ["pct_any_metric_rel_gt_1pct", "influence_score"]
    corr = df[corr_cols].corr(method="spearman")
    corr.index = [
        c.replace("rmse_", "RMSE ") if c.startswith("rmse_") else c for c in corr.index
    ]
    corr.columns = [
        c.replace("rmse_", "RMSE ") if c.startswith("rmse_") else c for c in corr.columns
    ]
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(
        corr,
        annot=True,
        fmt=".3f",
        cmap="RdBu_r",
        center=0,
        ax=ax,
        annot_kws={"size": 10},
    )
    ax.set_title(
        "Spearman correlation — LOO metrics vs influence score",
        fontsize=TITLE_FONT_SIZE,
    )
    plt.tight_layout()
    fig.savefig(
        f"{output_prefix}.metric_correlation_heatmap.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()

    summary_rows = []
    for col in rmse_cols:
        metric = col.replace("rmse_", "")
        pct_col = f"pct_rel_gt_1pct_{metric}"
        summary_rows.append(
            {
                "metric": metric,
                "mean_rmse": float(df[col].mean()),
                "median_rmse": float(df[col].median()),
                "max_rmse": float(df[col].max()),
                "max_taxa": str(df.loc[df[col].idxmax(), "left_out_taxa"]),
                "mean_pct_rel_gt_1pct": float(df[pct_col].mean()) if pct_col in df.columns else np.nan,
                "pct_of_total_rmse": float(df[col].sum() / df["rmse_sum"].sum() * 100.0),
                "mean_fraction_of_influence_score": float(
                    (df[col] / df["influence_score"]).mean()
                ),
            }
        )
    summary = pd.DataFrame(summary_rows).sort_values(
        "pct_of_total_rmse", ascending=False, kind="mergesort"
    )
    save_df_to_tsv(summary, f"{output_prefix}.metric_decomposition.info.tsv")

    print(f"[Info] Wrote LOO influence plots with prefix {output_prefix}")
    for _, row in summary.iterrows():
        print(
            f"  {row['metric']}: {row['pct_of_total_rmse']:.1f}% of total RMSE, "
            f"mean RMSE={row['mean_rmse']:.4f}"
        )
    return summary


LOO_INFLUENCE_TIER_BOUNDS = (0.1, 0.3, 0.5)
LOO_INFLUENCE_TIER_LABELS = ("normal", "mild", "moderate", "high")


def classify_loo_influence_tier(score: float) -> str:
    """Map influence_score to review tier (<0.1 normal … >0.5 high)."""
    for bound, label in zip(LOO_INFLUENCE_TIER_BOUNDS, LOO_INFLUENCE_TIER_LABELS[:-1]):
        if score < bound:
            return label
    return LOO_INFLUENCE_TIER_LABELS[-1]


def _classify_loo_influence_driver(row: pd.Series) -> str:
    """Classify why a LOO case has high influence (depth vs relative vs site-wide shift)."""
    pct_any = float(row.get("pct_any_metric_rel_gt_1pct", 0.0) or 0.0)
    if pct_any >= 10.0:
        return "coverage_pattern"
    rmse_cols = [f"rmse_{c}" for c in POPDEP_LOO_METRIC_COLS]
    present = [c for c in rmse_cols if c in row.index]
    if not present:
        return "unknown"
    dominant = row[present].idxmax().replace("rmse_", "")
    if dominant in ("Depth_Mean", "Depth_SD"):
        return "depth_shift"
    return "relative_depth"


def _prepare_loo_qc_review_table(
    influence_df: pd.DataFrame,
    sample_qc_df: pd.DataFrame,
    *,
    top_pct: float = 5.0,
) -> pd.DataFrame:
    """Merge LOO influence with sample QC and assign review tiers / drivers."""
    qc = sample_qc_df.copy()
    if "Sample" not in qc.columns:
        print("[Error] sample QC table missing Sample column")
        return pd.DataFrame()

    fmiss_cols = [c for c in qc.columns if c.endswith("_F_MISS")]
    depth_cols = [c for c in qc.columns if c.endswith("_log10_Depth")]
    if fmiss_cols:
        qc["mean_fmiss"] = qc[fmiss_cols].mean(axis=1)
        qc["max_fmiss"] = qc[fmiss_cols].max(axis=1)
    if depth_cols:
        qc["depth_log10_dispersion"] = qc[depth_cols].std(axis=1)

    merged = influence_df.merge(
        qc,
        left_on="left_out_taxa",
        right_on="Sample",
        how="left",
    )
    merged["influence_tier"] = merged["influence_score"].map(classify_loo_influence_tier)
    merged["influence_driver"] = merged.apply(_classify_loo_influence_driver, axis=1)
    n_top = max(1, int(np.ceil(len(merged) * top_pct / 100.0)))
    top_cut = merged["influence_score"].nlargest(n_top).min()
    merged["review_priority"] = np.where(
        (merged["influence_score"] >= 0.5)
        | (merged["influence_score"] >= top_cut)
        | (merged.get("pct_any_metric_rel_gt_1pct", 0) >= 10.0),
        "required",
        np.where(merged["influence_score"] >= 0.3, "recommended", "optional"),
    )
    return merged.sort_values("influence_score", ascending=False, kind="mergesort")


def plot_popdep_loo_outlier_review(
    info_tsv: str,
    sample_qc_tsv: str,
    output_prefix: str,
    *,
    top_n_highlight: int = 15,
    top_pct: float = 5.0,
) -> pd.DataFrame | None:
    """
    Diagnostic plots for high-influence LOO samples merged with BAM / depth QC.

    Writes tiered distribution, influence-vs-QC scatters, driver summary, and a
    review candidate table under ``{output_prefix}.*``.
    """
    influence_df = load_df_from_tsv(info_tsv)
    sample_qc_df = load_df_from_tsv(sample_qc_tsv)
    if influence_df is None or influence_df.empty:
        print(f"[Error] Empty or missing LOO influence table: {info_tsv}")
        return None
    if sample_qc_df is None or sample_qc_df.empty:
        print(f"[Error] Empty or missing sample QC table: {sample_qc_tsv}")
        return None

    review = _prepare_loo_qc_review_table(
        influence_df,
        sample_qc_df,
        top_pct=top_pct,
    )
    if review.empty:
        return None

    tier_thresholds = [
        {"value": b, "color": c, "linestyle": "--", "label": lbl}
        for b, c, lbl in zip(
            LOO_INFLUENCE_TIER_BOUNDS,
            ("gray", "orange", "crimson"),
            ("0.1 normal", "0.3 moderate", "0.5 high"),
        )
    ]
    plot_distribution_with_stats(
        review,
        col="influence_score",
        title="PopDep LOO influence score — review tiers",
        filename=f"{output_prefix}.influence_score_tiered_distribution.png",
        x_label="Influence score",
        y_label="Count",
        bins=40,
        thresholds=tier_thresholds,
        color="steelblue",
    )

    top_ids = review.head(top_n_highlight)["left_out_taxa"].tolist()
    scatter_specs = [
        ("Coverage", "WGS coverage (×)", "influence_vs_coverage"),
        ("mean_fmiss", "Mean subgenome F_MISS", "influence_vs_mean_fmiss"),
        ("max_fmiss", "Max subgenome F_MISS", "influence_vs_max_fmiss"),
        ("depth_log10_dispersion", "A/B/D log10 depth dispersion", "influence_vs_depth_dispersion"),
    ]
    for x_col, x_label, slug in scatter_specs:
        if x_col not in review.columns:
            continue
        plot_scatter_with_thresholds(
            review,
            x_col=x_col,
            y_col="influence_score",
            title=f"Influence score vs {x_label}",
            filename=f"{output_prefix}.{slug}.png",
            xlabel=x_label,
            ylabel="Influence score",
            thresholds_h=[
                {"value": 0.1, "color": "gray", "linestyle": ":", "label": "0.1"},
                {"value": 0.3, "color": "orange", "linestyle": "--", "label": "0.3"},
                {"value": 0.5, "color": "crimson", "linestyle": "--", "label": "0.5 high"},
            ],
            group_col="Group" if "Group" in review.columns else None,
            highlight_ids=top_ids,
            id_col="left_out_taxa",
            alpha=0.55,
            s=30,
            figure_size=(10, 8),
        )

    driver_counts = review.loc[
        review["influence_score"] >= 0.3, "influence_driver"
    ].value_counts()
    if not driver_counts.empty:
        plot_bar_chart(
            driver_counts.index.tolist(),
            driver_counts.values.tolist(),
            title="Influence driver type (score ≥ 0.3)",
            ylabel="Sample count",
            filename=f"{output_prefix}.influence_driver_counts.png",
            ylim=(0, float(driver_counts.max()) + 1),
            color="mediumpurple",
        )

    top_review = review.head(top_n_highlight).copy()
    plot_bar_chart(
        top_review["left_out_taxa"].tolist(),
        top_review["influence_score"].tolist(),
        title=f"Top {top_n_highlight} samples for manual review",
        ylabel="Influence score",
        filename=f"{output_prefix}.top_review_candidates.png",
        ylim=None,
        color="indianred",
        rotate_xlabels=60,
        figure_size=(12, 5),
    )

    out_cols = [
        "left_out_taxa",
        "influence_score",
        "influence_tier",
        "influence_driver",
        "review_priority",
        "pct_any_metric_rel_gt_1pct",
        "rmse_Depth_Mean",
        "rmse_Depth_SD",
        "rmse_RelativeDepth_Mean",
        "rmse_RelativeDepth_SD",
        "Coverage",
        "Group",
        "mean_fmiss",
        "max_fmiss",
        "depth_log10_dispersion",
    ]
    present_cols = [c for c in out_cols if c in review.columns]
    save_df_to_tsv(review[present_cols], f"{output_prefix}.review_candidates.info.tsv")

    n_high = int((review["influence_score"] >= 0.5).sum())
    n_required = int((review["review_priority"] == "required").sum())
    print(f"[Info] Wrote LOO outlier review plots with prefix {output_prefix}")
    print(f"  review required: {n_required} samples (incl. {n_high} with score ≥ 0.5)")
    for _, row in review.head(10).iterrows():
        print(
            f"  {row['left_out_taxa']}: score={row['influence_score']:.4f} "
            f"tier={row['influence_tier']} driver={row['influence_driver']} "
            f"priority={row['review_priority']}"
        )
    return review


LOO_QUADRANT_ACTIONS = {
    "keep_high_influence": "keep",
    "delete_candidate": "delete_candidate",
    "normal": "keep",
    "check_low_quality": "check",
    "review_coverage_pattern": "review",
}


def _clip_unit_interval(values: pd.Series | float) -> pd.Series | float:
    return np.clip(values, 0.0, 1.0)


def compute_loo_sample_quality_score(df: pd.DataFrame) -> pd.Series:
    """
    Composite sample QC score in [0, 1]; higher means more reliable for PopDep.

    Uses genotype missing (F_MISS), subgenome depth uniformity, and positive GAM
    missing residuals when present. Does not use influence_score.
    """
    parts: list[pd.Series] = []
    if "mean_fmiss" in df.columns:
        parts.append(1.0 - _clip_unit_interval(df["mean_fmiss"] / 0.15))
    if "max_fmiss" in df.columns:
        parts.append(1.0 - _clip_unit_interval(df["max_fmiss"] / 0.20))
    if "depth_log10_dispersion" in df.columns:
        parts.append(1.0 - _clip_unit_interval(df["depth_log10_dispersion"] / 0.10))

    residual_cols = [
        c for c in df.columns if c.endswith("_residual") and not c.startswith("influence")
    ]
    if residual_cols:
        max_pos_residual = df[residual_cols].clip(lower=0).max(axis=1)
        parts.append(1.0 - _clip_unit_interval(max_pos_residual / 0.05))

    if not parts:
        return pd.Series(0.5, index=df.index, dtype=float)
    return sum(parts) / len(parts)


def classify_loo_influence_quality_quadrant(
    influence: float,
    quality: float,
    *,
    influence_threshold: float,
    quality_threshold: float,
    influence_driver: str = "",
    pct_any_rel_gt_1pct: float = 0.0,
    coverage_pattern_pct: float = 10.0,
) -> str:
    """Map influence × quality to quadrant label."""
    if (
        influence_driver == "coverage_pattern"
        and pct_any_rel_gt_1pct >= coverage_pattern_pct
    ):
        return "review_coverage_pattern"
    high_influence = influence >= influence_threshold
    high_quality = quality >= quality_threshold
    if high_influence and high_quality:
        return "keep_high_influence"
    if high_influence and not high_quality:
        return "delete_candidate"
    if not high_influence and high_quality:
        return "normal"
    return "check_low_quality"


LOO_QUADRANT_PALETTE = {
    "keep_high_influence": "#2ca02c",
    "delete_candidate": "#d62728",
    "normal": "#1f77b4",
    "check_low_quality": "#ff7f0e",
    "review_coverage_pattern": "#9467bd",
}


def _load_aneuploidy_manifest_samples(manifest_tsv: str) -> set[str]:
    manifest = pd.read_csv(manifest_tsv, sep="\t")
    if manifest.empty or "sample" not in manifest.columns:
        return set()
    return {str(s) for s in manifest["sample"].astype(str)}


def _assign_loo_quadrant_columns(
    review: pd.DataFrame,
    influence_col: str,
    *,
    influence_threshold: float | None = None,
    quality_threshold: float | None = None,
) -> tuple[pd.DataFrame, float, float]:
    """Add sample_quality_score, thresholds, quadrant, and recommended_action."""
    review = review.copy()
    review["sample_quality_score"] = compute_loo_sample_quality_score(review)
    inf_vals = review[influence_col]
    inf_thr = (
        float(influence_threshold)
        if influence_threshold is not None
        else float(max(0.1, inf_vals.quantile(0.75)))
    )
    qual_thr = (
        float(quality_threshold)
        if quality_threshold is not None
        else float(review["sample_quality_score"].median())
    )
    review["influence_metric"] = influence_col
    review["influence_value"] = inf_vals
    review["influence_threshold"] = inf_thr
    review["quality_threshold"] = qual_thr
    review["quadrant"] = review.apply(
        lambda row: classify_loo_influence_quality_quadrant(
            float(row[influence_col]),
            float(row["sample_quality_score"]),
            influence_threshold=inf_thr,
            quality_threshold=qual_thr,
            influence_driver=str(row.get("influence_driver", "")),
            pct_any_rel_gt_1pct=float(row.get("pct_any_metric_rel_gt_1pct", 0.0) or 0.0),
        ),
        axis=1,
    )
    review["recommended_action"] = review["quadrant"].map(LOO_QUADRANT_ACTIONS)
    return review, inf_thr, qual_thr


def _draw_influence_quality_quadrant_scatter(
    review: pd.DataFrame,
    influence_col: str,
    inf_thr: float,
    qual_thr: float,
    filename: str,
    *,
    title: str,
    x_label: str,
    highlight_ids: set[str] | None = None,
    highlight_label: str = "aneuploidy (manifest)",
    label_col: str = "left_out_taxa",
    figure_size: tuple[float, float] = (10, 8),
) -> None:
    """Draw quadrant scatter; overlay star markers for ``highlight_ids``."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    from infra.utils.graph import LEGEND_FONT_SIZE, TITLE_FONT_SIZE

    plot_df = review.replace([np.inf, -np.inf], np.nan).dropna(
        subset=[influence_col, "sample_quality_score"]
    )
    highlight_ids = highlight_ids or set()
    slug = {
        "influence_score": "influence",
        "influence_pattern": "pattern",
        "influence_absolute": "absolute",
    }.get(influence_col, influence_col.replace("influence_", ""))

    sns.set_style("white")
    fig, ax = plt.subplots(figsize=figure_size)
    bg = plot_df[~plot_df[label_col].astype(str).isin(highlight_ids)]
    hi = plot_df[plot_df[label_col].astype(str).isin(highlight_ids)]

    for quadrant, sub in bg.groupby("quadrant", sort=False):
        ax.scatter(
            sub[influence_col],
            sub["sample_quality_score"],
            label=quadrant,
            c=LOO_QUADRANT_PALETTE.get(quadrant, "gray"),
            alpha=0.55,
            s=32,
            edgecolors="white",
            linewidths=0.3,
            zorder=2,
        )
    if not hi.empty:
        for _, sub in hi.groupby("quadrant", sort=False):
            ax.scatter(
                sub[influence_col],
                sub["sample_quality_score"],
                marker="*",
                s=220,
                c="#FFD700",
                edgecolors="#B22222",
                linewidths=1.2,
                zorder=5,
            )
        if highlight_label:
            ax.scatter(
                [], [], marker="*", s=220, c="#FFD700", edgecolors="#B22222",
                linewidths=1.2, label=highlight_label,
            )
        for _, row in hi.iterrows():
            ax.annotate(
                str(row[label_col]),
                (row[influence_col], row["sample_quality_score"]),
                fontsize=8,
                fontweight="bold",
                color="#8B0000",
                xytext=(5, 5),
                textcoords="offset points",
                zorder=6,
            )

    ax.axvline(inf_thr, color="crimson", linestyle="--", linewidth=1.2, label=f"{slug} ≥ {inf_thr:.3f}")
    ax.axhline(qual_thr, color="steelblue", linestyle="--", linewidth=1.2, label=f"quality ≥ {qual_thr:.3f}")
    ax.set_xlabel(x_label, fontsize=14)
    ax.set_ylabel("Sample quality score (genotype / depth QC)", fontsize=14)
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.set_ylim(0, 1.02)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.14), ncol=2, fontsize=LEGEND_FONT_SIZE, frameon=False)
    plt.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_aneuploidy_samechr_pattern_quality_quadrant(
    manifest_tsv: str,
    samechr_publish_root: str,
    sample_qc_tsv: str,
    output_prefix: str,
    *,
    bench_max_bp: int = 10_000_000,
    influence_threshold: float | None = None,
    quality_threshold: float | None = None,
) -> pd.DataFrame | None:
    """
    Pattern influence × sample quality for aneuploidy same-chr LOO cases.

    Writes a combined plot (all manifest rows) and per-``popdep_chr`` plots under
    ``{output_prefix}.aneuploidy_samechr_*``.
    """
    from infra.utils.graph import TITLE_FONT_SIZE, plot_bar_chart

    manifest = pd.read_csv(manifest_tsv, sep="\t")
    if manifest.empty:
        print("[Error] Empty aneuploidy manifest.")
        return None

    sample_qc_df = load_df_from_tsv(sample_qc_tsv)
    if sample_qc_df is None or sample_qc_df.empty:
        print("[Error] Missing sample QC table for aneuploidy samechr quadrant.")
        return None

    publish_root = Path(samechr_publish_root)
    rows: list[dict] = []
    for _, entry in manifest.iterrows():
        sample = str(entry["sample"])
        popdep_chr = int(entry["popdep_chr"])
        info_path = publish_root / f"chr{popdep_chr}" / "stats" / (
            f"n500_aneuploidy_samechr_chr{popdep_chr}_{bench_max_bp}_loo.info.tsv"
        )
        if not info_path.exists():
            print(f"[Warning] Missing samechr info: {info_path}")
            continue
        inf_df = ensure_influence_components_df(pd.read_csv(info_path, sep="\t"))
        hit = inf_df.loc[inf_df["left_out_taxa"] == sample]
        if hit.empty:
            continue
        row = hit.iloc[0].to_dict()
        row["left_out_taxa"] = sample
        row["plot_label"] = f"{sample} (chr{popdep_chr})"
        row["popdep_chr"] = popdep_chr
        row["abnormal_chr"] = entry.get("abnormal_chr", "")
        rows.append(row)

    if not rows:
        print("[Error] No aneuploidy samechr LOO rows loaded.")
        return None

    influence_df = pd.DataFrame(rows)
    qc = sample_qc_df.copy()
    fmiss_cols = [c for c in qc.columns if c.endswith("_F_MISS")]
    depth_cols = [c for c in qc.columns if c.endswith("_log10_Depth")]
    if fmiss_cols:
        qc["mean_fmiss"] = qc[fmiss_cols].mean(axis=1)
        qc["max_fmiss"] = qc[fmiss_cols].max(axis=1)
    if depth_cols:
        qc["depth_log10_dispersion"] = qc[depth_cols].std(axis=1)

    review = influence_df.merge(qc, left_on="left_out_taxa", right_on="Sample", how="left")
    review["influence_driver"] = review.apply(_classify_loo_influence_driver, axis=1)
    review, inf_thr, qual_thr = _assign_loo_quadrant_columns(
        review,
        "influence_pattern",
        influence_threshold=influence_threshold,
        quality_threshold=quality_threshold,
    )

    _draw_influence_quality_quadrant_scatter(
        review,
        "influence_pattern",
        inf_thr,
        qual_thr,
        f"{output_prefix}.aneuploidy_samechr_pattern_quality_quadrant.png",
        title="Aneuploidy same-chr LOO — pattern influence × sample quality",
        x_label="Pattern influence — same abnormal chromosome LOO (10 Mb)",
        highlight_ids=set(review["plot_label"].astype(str)),
        highlight_label="aneuploidy case",
        label_col="plot_label",
        figure_size=(11, 9),
    )

    counts = review["quadrant"].value_counts()
    plot_bar_chart(
        counts.index.tolist(),
        counts.values.tolist(),
        title="Aneuploidy same-chr pattern × quality quadrant counts",
        ylabel="Case count",
        filename=f"{output_prefix}.aneuploidy_samechr_pattern_quadrant_counts.png",
        ylim=(0, float(counts.max()) + 1),
        color="steelblue",
        rotate_xlabels=20,
    )

    for popdep_chr, sub in review.groupby("popdep_chr", sort=True):
        sub = sub.copy()
        sub_review, chr_inf_thr, chr_qual_thr = _assign_loo_quadrant_columns(
            sub,
            "influence_pattern",
            influence_threshold=inf_thr,
            quality_threshold=qual_thr,
        )
        _draw_influence_quality_quadrant_scatter(
            sub_review,
            "influence_pattern",
            chr_inf_thr,
            chr_qual_thr,
            f"{output_prefix}.aneuploidy_samechr_chr{int(popdep_chr)}_pattern_quality_quadrant.png",
            title=f"Aneuploidy chr{int(popdep_chr)} same-chr LOO — pattern × quality",
            x_label=f"Pattern influence — chr{int(popdep_chr)} LOO (10 Mb)",
            highlight_ids=set(sub_review["plot_label"].astype(str)),
            highlight_label="aneuploidy case",
            label_col="plot_label",
            figure_size=(9, 7),
        )

    out_cols = [
        "plot_label",
        "left_out_taxa",
        "popdep_chr",
        "abnormal_chr",
        "influence_pattern",
        "influence_absolute",
        "influence_score",
        "sample_quality_score",
        "influence_threshold",
        "quality_threshold",
        "quadrant",
        "recommended_action",
        "influence_driver",
        "pct_any_metric_rel_gt_1pct",
        "Coverage",
        "Group",
        "mean_fmiss",
        "max_fmiss",
    ]
    present = [c for c in out_cols if c in review.columns]
    save_df_to_tsv(review[present], f"{output_prefix}.aneuploidy_samechr_pattern_quadrant_review.info.tsv")
    print(f"[Info] Wrote aneuploidy samechr pattern quadrants ({len(review)} cases) prefix={output_prefix}")
    return review


def plot_popdep_loo_influence_quality_quadrant(
    info_tsv: str,
    sample_qc_tsv: str,
    output_prefix: str,
    *,
    influence_col: str = "influence_score",
    influence_threshold: float | None = None,
    quality_threshold: float | None = None,
    top_n_highlight: int = 15,
    highlight_manifest_tsv: str | None = None,
) -> pd.DataFrame | None:
    """
    Influence × sample-quality quadrant review for LOO cases.

    X = ``influence_col`` (default combined ``influence_score``; use
    ``influence_pattern`` for coverage-pattern screening); Y = sample_quality_score.

    When ``highlight_manifest_tsv`` is set, manifest samples are overlaid with
    gold star / crimson edge markers on the scatter plot.
    """
    from infra.utils.graph import plot_bar_chart

    influence_df = ensure_influence_components_df(load_df_from_tsv(info_tsv))
    sample_qc_df = load_df_from_tsv(sample_qc_tsv)
    if influence_df is None or influence_df.empty or sample_qc_df is None:
        print("[Error] Missing LOO influence or sample QC table for quadrant review.")
        return None
    if influence_col not in influence_df.columns:
        print(f"[Error] Missing {influence_col} in LOO influence table.")
        return None

    review = _prepare_loo_qc_review_table(influence_df, sample_qc_df)
    if review.empty:
        return None

    highlight_ids: set[str] = set()
    if highlight_manifest_tsv and Path(highlight_manifest_tsv).exists():
        highlight_ids = _load_aneuploidy_manifest_samples(highlight_manifest_tsv)
        review["is_aneuploidy_manifest"] = review["left_out_taxa"].astype(str).isin(highlight_ids)

    review, inf_thr, qual_thr = _assign_loo_quadrant_columns(
        review,
        influence_col,
        influence_threshold=influence_threshold,
        quality_threshold=quality_threshold,
    )

    slug = {
        "influence_score": "influence",
        "influence_pattern": "pattern",
        "influence_absolute": "absolute",
    }.get(influence_col, influence_col.replace("influence_", ""))
    x_label = (
        "Pattern influence (LOO coverage-profile impact)"
        if influence_col == "influence_pattern"
        else "Influence score (LOO impact on PopDep)"
    )
    title_metric = "pattern influence" if influence_col == "influence_pattern" else "influence"
    title = f"LOO {title_metric} × sample quality — quadrant review"
    if highlight_ids:
        title += " (aneuploidy ★ highlighted)"

    _draw_influence_quality_quadrant_scatter(
        review,
        influence_col,
        inf_thr,
        qual_thr,
        f"{output_prefix}.{slug}_quality_quadrant.png",
        title=title,
        x_label=x_label,
        highlight_ids=highlight_ids,
        highlight_label="aneuploidy (manifest)" if highlight_ids else "",
    )

    counts = review["quadrant"].value_counts()
    plot_bar_chart(
        counts.index.tolist(),
        counts.values.tolist(),
        title=f"LOO quadrant counts ({title_metric} × quality)",
        ylabel="Sample count",
        filename=f"{output_prefix}.{slug}_quadrant_counts.png",
        ylim=(0, float(counts.max()) + 2),
        color="steelblue",
        rotate_xlabels=20,
    )

    out_cols = [
        "left_out_taxa",
        "is_aneuploidy_manifest",
        influence_col,
        "influence_score",
        "influence_absolute",
        "influence_pattern",
        "influence_value",
        "influence_metric",
        "sample_quality_score",
        "influence_threshold",
        "quality_threshold",
        "quadrant",
        "recommended_action",
        "influence_driver",
        "influence_tier",
        "pct_any_metric_rel_gt_1pct",
        "Coverage",
        "Group",
        "mean_fmiss",
        "max_fmiss",
        "depth_log10_dispersion",
        "rmse_Depth_Mean",
        "rmse_Depth_SD",
    ]
    residual_cols = [c for c in review.columns if c.endswith("_residual")]
    present = [c for c in out_cols + residual_cols if c in review.columns]
    save_df_to_tsv(review[present], f"{output_prefix}.{slug}_quadrant_review.info.tsv")

    delete_n = int((review["recommended_action"] == "delete_candidate").sum())
    keep_n = int((review["quadrant"] == "keep_high_influence").sum())
    review_n = int((review["quadrant"] == "review_coverage_pattern").sum())
    print(f"[Info] Wrote LOO quadrant review ({influence_col}) with prefix {output_prefix}")
    print(f"  thresholds: {influence_col}≥{inf_thr:.4f}, quality≥{qual_thr:.4f}")
    print(f"  keep_high_influence={keep_n}, delete_candidate={delete_n}, review_coverage_pattern={review_n}")
    if highlight_ids:
        aneu = review.loc[review["left_out_taxa"].astype(str).isin(highlight_ids)]
        print(f"  aneuploidy manifest highlighted: {len(aneu)} samples in N500 chr1 LOO table")
        for _, row in aneu.sort_values(influence_col, ascending=False).iterrows():
            print(
                f"    {row['left_out_taxa']}: {row['quadrant']} "
                f"{influence_col}={row[influence_col]:.4f} quality={row['sample_quality_score']:.4f}"
            )
    flagged = review.loc[
        review["recommended_action"].isin(["delete_candidate", "review_coverage_pattern"])
    ].sort_values(influence_col, ascending=False, kind="mergesort")
    for _, row in flagged.head(10).iterrows():
        print(
            f"  {row['left_out_taxa']}: action={row['recommended_action']} "
            f"{influence_col}={row[influence_col]:.4f} quality={row['sample_quality_score']:.4f} "
            f"Coverage={float(row.get('Coverage', float('nan'))):.2f} "
            f"mean_fmiss={float(row.get('mean_fmiss', float('nan'))):.4f}"
        )
    return review

