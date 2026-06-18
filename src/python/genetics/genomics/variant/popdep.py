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


def _tabix_popdep_lookup(
    chrom: int,
    positions: list[int],
    bgz_path: Path,
) -> dict[int, tuple[float, float]]:
    """Batch tabix fetch for sorted variant positions on one chromosome grid."""
    if not positions:
        return {}

    result: dict[int, tuple[float, float]] = {}
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
            parts = line.rstrip().split("\t")
            if len(parts) >= 4:
                pos_s, mean_s, sd_s = parts[1], parts[2], parts[3]
            else:
                pos_s, mean_s, sd_s = parts[0], parts[1], parts[2]
            result[int(pos_s)] = (float(mean_s), float(sd_s))
    finally:
        os.unlink(regions_file.name)

    return result


def _stream_popdep_lookup(
    popdep_path: str,
    wanted_positions: list[int],
) -> dict[int, tuple[float, float]]:
    """Fallback merge-join when tabix sidecar is unavailable."""
    if not wanted_positions:
        return {}
    wanted_set = set(wanted_positions)
    max_pos = max(wanted_positions)
    result: dict[int, tuple[float, float]] = {}
    if not os.path.isfile(popdep_path):
        print(f"[Warn] Popdep reference missing: {popdep_path}")
        return result

    open_fn = gzip.open if popdep_path.endswith((".gz", ".bgz")) else open
    with open_fn(popdep_path, "rt", encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            return result
        header_cols = header.rstrip("\n").split("\t")
        pos_idx = 1 if header_cols and header_cols[0] == "Chrom" else 0
        mean_idx = pos_idx + 1
        sd_idx = pos_idx + 2
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= sd_idx:
                continue
            pos = int(parts[pos_idx])
            if pos in wanted_set:
                result[pos] = (float(parts[mean_idx]), float(parts[sd_idx]))
            if pos >= max_pos and len(result) == len(wanted_set):
                break
    return result


def save_variant_popdep_info(df: pd.DataFrame, output_path: str) -> None:
    """Write variant popdepth table with PLINK-style ``#CHROM`` header."""
    out = df.copy()
    if "CHROM" in out.columns and "#CHROM" not in out.columns:
        out = out.rename(columns={"CHROM": "#CHROM"})
    cols = ["#CHROM", "ID", "POS", "Depth_Mean", "Depth_SD", "Depth_CV"]
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
        if depth is not None:
            depth_mean, depth_sd = depth
            if depth_mean > 0:
                depth_cv = depth_sd / depth_mean
        rows.append(
            {
                "CHROM": chrom,
                "ID": variant_id,
                "POS": pos,
                "Depth_Mean": depth_mean,
                "Depth_SD": depth_sd,
                "Depth_CV": depth_cv,
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
    ``#CHROM, ID, POS, Depth_Mean, Depth_SD, Depth_CV`` (NaN when reference depth is missing).
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
    if "Depth_CV" not in df.columns:
        df["Depth_CV"] = np.where(
            df["Depth_Mean"] > 0,
            df["Depth_SD"] / df["Depth_Mean"],
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
        
        plot_regression_comparison(
            merged, x_col, y_col, 
            x_label=x_label, y_label=y_label, 
            filename=f"{output_prefix}.reg_{y_col}_vs_{x_col}.png",
            title=f"{y_label} vs {x_label}"
        )

# ==================================================================================
# Analysis 3: Mahalanobis Distance for Variant Depth Distribution
# ==================================================================================

def ana_popdep_mahalanobis(popdep_file, output_prefix, threshold_p_value=0.99):
    """
    Calculates Mahalanobis distance for variant depth distribution (Log_Depth_Mean vs Log_Depth_CV),
    sets an outlier threshold, saves the result, and plots it.
    """
    print(f"[Info] Running Mahalanobis distance analysis on: {popdep_file}")
    df = load_popdep_data(popdep_file)
    if df is None: return

    # Filter out missing/zero variants to perform log transform safely
    df_valid = df[(df['Depth_Mean'] > 0) & (df['Depth_CV'] > 0)].copy()
    if df_valid.empty:
        print("[Error] No valid data for Mahalanobis calculation.")
        return

    df_valid['Log_Depth_Mean'] = np.log(df_valid['Depth_Mean'])
    df_valid['Log_Depth_CV'] = np.log(df_valid['Depth_CV'])
    
    data = df_valid[['Log_Depth_Mean', 'Log_Depth_CV']].values
    
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

    # Save to file
    out_file = f"{output_prefix}.info.tsv"
    out_df = df_valid[['Position', 'Depth_Mean', 'Depth_SD', 'Depth_CV', 'Mahalanobis', 'Mahalanobis_Sq', 'Outlier']]
    save_df_to_tsv(out_df, out_file)

    # Plot
    plot_file = f"{output_prefix}.mahalanobis.png"
    plot_scatter_with_outliers(
        data=df_valid,
        x_col='Log_Depth_Mean',
        y_col='Log_Depth_CV',
        outlier_col='Outlier',
        title=f'Mahalanobis Distance Outlier Detection (p={threshold_p_value})',
        filename=plot_file,
        xlabel='Log Depth Mean',
        ylabel='Log Depth CV',
        color_normal='blue',
        color_outlier='red',
        s=2
    )


