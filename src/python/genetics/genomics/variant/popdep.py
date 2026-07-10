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
from infra.utils.graph import (
    DEFAULT_AXIS_PADDING_FRACTION,
    axis_limits_with_adaptive_upper,
    padded_axis_limits,
    plot_distribution_with_stats,
    plot_mean_variance_fit,
    plot_qq_residuals,
    plot_regression_comparison,
    plot_scatter_with_outliers,
)
from infra.utils.io import load_df_from_space_sep, save_df_to_tsv

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


