from collections import deque
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from typing import Dict, Optional

import numpy as np
import pandas as pd

from infra.utils.graph import (
    plot_distribution_with_stats,
    plot_multi_line_series,
    plot_scatter_with_thresholds,
)
from infra.utils.io import save_df_to_tsv, save_thresholds


def _normalize_ld_chunk(chunk: pd.DataFrame) -> Optional[pd.DataFrame]:
    """Normalize one PLINK2 LD chunk to BP_A/BP_B/R2 columns."""
    col_map = {
        "#CHROM_A": "CHR_A",
        "CHROM_A": "CHR_A",
        "#CHROM_B": "CHR_B",
        "CHROM_B": "CHR_B",
        "POS_A": "BP_A",
        "POS_B": "BP_B",
        "UNPHASED_R2": "R2",
    }
    chunk = chunk.rename(columns=col_map)
    required_cols = {"BP_A", "BP_B", "R2"}
    if not required_cols.issubset(chunk.columns):
        return None

    local = chunk[["BP_A", "BP_B", "R2"]].copy()
    local["R2"] = pd.to_numeric(local["R2"], errors="coerce")
    local["BP_A"] = pd.to_numeric(local["BP_A"], errors="coerce")
    local["BP_B"] = pd.to_numeric(local["BP_B"], errors="coerce")
    local = local.dropna(subset=["R2", "BP_A", "BP_B"])
    if local.empty:
        return None
    return local


def _reservoir_merge(
    current: np.ndarray, incoming: np.ndarray, max_points: int, rng: np.random.Generator
) -> np.ndarray:
    """Merge two arrays with cap using uniform random downsampling."""
    if len(incoming) == 0:
        return current
    if len(current) == 0:
        if len(incoming) <= max_points:
            return incoming
        idx = rng.choice(len(incoming), size=max_points, replace=False)
        return incoming[idx]

    merged = np.concatenate([current, incoming])
    if len(merged) <= max_points:
        return merged
    idx = rng.choice(len(merged), size=max_points, replace=False)
    return merged[idx]


def _process_chunk_arrays(
    distance_kb: np.ndarray,
    r2: np.ndarray,
    n_bins: int,
    bin_size_kb: float,
    median_sample_per_bin: int,
    scatter_max_points: int,
    seed: int,
) -> Dict[str, object]:
    """Compute partial LD bin stats for one chunk."""
    local_sum = np.zeros(n_bins, dtype=np.float64)
    local_cnt = np.zeros(n_bins, dtype=np.int64)
    local_median: Dict[int, np.ndarray] = {i: np.array([], dtype=np.float32) for i in range(n_bins)}
    local_scatter = np.empty((0, 2), dtype=np.float32)
    rng = np.random.default_rng(seed)

    if len(r2) == 0:
        return {
            "sum_r2": local_sum,
            "cnt_r2": local_cnt,
            "median_samples": local_median,
            "scatter_pairs": local_scatter,
            "total_pairs": 0,
        }

    bin_idx = np.floor(distance_kb / float(bin_size_kb)).astype(np.int64)
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)

    np.add.at(local_sum, bin_idx, r2)
    np.add.at(local_cnt, bin_idx, 1)

    for b in np.unique(bin_idx):
        vals = r2[bin_idx == b].astype(np.float32)
        local_median[int(b)] = _reservoir_merge(
            local_median[int(b)], vals, median_sample_per_bin, rng
        )

    pair_chunk = np.column_stack([distance_kb.astype(np.float32), r2.astype(np.float32)])
    if len(pair_chunk) <= scatter_max_points:
        local_scatter = pair_chunk
    else:
        idx = rng.choice(len(pair_chunk), size=scatter_max_points, replace=False)
        local_scatter = pair_chunk[idx]

    return {
        "sum_r2": local_sum,
        "cnt_r2": local_cnt,
        "median_samples": local_median,
        "scatter_pairs": local_scatter,
        "total_pairs": int(len(r2)),
    }


def ana_ld_decay(
    input_file: str,
    output_prefix: str,
    max_distance_kb: int = 5000,
    bin_size_kb: int = 50,
    scatter_max_points: int = 200000,
    chunksize: int = 2_000_000,
    median_sample_per_bin: int = 20000,
    workers: int = 1,
    random_seed: int = 42,
) -> None:
    """
    Summarize LD versus distance from PLINK2 --r2 output and write TSV plus
    distribution plots (pairwise distance, r2, per-bin mean r2) via infra plotting.
    """
    if max_distance_kb <= 0 or bin_size_kb <= 0:
        print("[Error] Invalid bin settings for LD decay.")
        return
    n_bins = int(np.ceil(float(max_distance_kb) / float(bin_size_kb)))
    if n_bins <= 0:
        print("[Error] Invalid number of bins for LD decay.")
        return

    sum_r2 = np.zeros(n_bins, dtype=np.float64)
    cnt_r2 = np.zeros(n_bins, dtype=np.int64)
    median_samples: Dict[int, np.ndarray] = {i: np.array([], dtype=np.float32) for i in range(n_bins)}
    scatter_pairs = np.empty((0, 2), dtype=np.float32)
    rng = np.random.default_rng(random_seed)
    total_pairs = 0

    try:
        reader = pd.read_csv(input_file, sep=r"\s+", chunksize=chunksize)
    except Exception as e:
        print(f"[Error] Failed to read {input_file}: {e}")
        return

    futures = deque()
    pending_results: Dict[int, Dict[str, object]] = {}
    next_merge_idx = 0
    max_queue = max(1, int(workers) * 2)

    def merge_partial(partial: Dict[str, object]) -> None:
        nonlocal scatter_pairs, total_pairs
        sum_r2[:] = sum_r2 + partial["sum_r2"]
        cnt_r2[:] = cnt_r2 + partial["cnt_r2"]
        total_pairs += int(partial["total_pairs"])
        for b, vals in partial["median_samples"].items():
            if len(vals) == 0:
                continue
            median_samples[int(b)] = _reservoir_merge(
                median_samples[int(b)], vals, median_sample_per_bin, rng
            )
        pairs = partial["scatter_pairs"]
        if len(pairs) > 0:
            if len(scatter_pairs) == 0:
                scatter_pairs = pairs if len(pairs) <= scatter_max_points else pairs[:scatter_max_points]
            else:
                merged_pairs = np.vstack([scatter_pairs, pairs])
                if len(merged_pairs) <= scatter_max_points:
                    scatter_pairs = merged_pairs
                else:
                    idx = rng.choice(len(merged_pairs), size=scatter_max_points, replace=False)
                    scatter_pairs = merged_pairs[idx]

    use_parallel = int(workers) > 1
    chunk_idx = 0
    pool = ProcessPoolExecutor(max_workers=int(workers)) if use_parallel else None

    def flush_pending_in_order() -> None:
        nonlocal next_merge_idx
        while next_merge_idx in pending_results:
            merge_partial(pending_results.pop(next_merge_idx))
            next_merge_idx += 1

    for chunk in reader:
        local = _normalize_ld_chunk(chunk)
        if local is None or local.empty:
            continue

        distance_kb = (local["BP_B"] - local["BP_A"]).abs().to_numpy(dtype=np.float64) / 1000.0
        r2 = local["R2"].to_numpy(dtype=np.float64)

        valid_mask = (distance_kb >= 0.0) & (distance_kb <= float(max_distance_kb)) & np.isfinite(r2)
        if not np.any(valid_mask):
            continue

        distance_kb = distance_kb[valid_mask]
        r2 = r2[valid_mask]
        if use_parallel:
            futures.append(
                (
                    chunk_idx,
                    pool.submit(
                    _process_chunk_arrays,
                    distance_kb,
                    r2,
                    n_bins,
                    float(bin_size_kb),
                    int(median_sample_per_bin),
                    int(scatter_max_points),
                    int(random_seed) + chunk_idx,
                ),
                )
            )
            chunk_idx += 1
            if len(futures) >= max_queue:
                done, _ = wait([fut for _, fut in futures], return_when=FIRST_COMPLETED)
                for pair in list(futures):
                    idx, fut = pair
                    if fut in done:
                        futures.remove(pair)
                        pending_results[idx] = fut.result()
                flush_pending_in_order()
        else:
            partial = _process_chunk_arrays(
                distance_kb,
                r2,
                n_bins,
                float(bin_size_kb),
                int(median_sample_per_bin),
                int(scatter_max_points),
                int(random_seed) + chunk_idx,
            )
            chunk_idx += 1
            merge_partial(partial)

    if use_parallel:
        wait([fut for _, fut in futures])
        for idx, fut in futures:
            pending_results[idx] = fut.result()
        flush_pending_in_order()
        pool.shutdown(wait=True)

    used_bins = cnt_r2 > 0
    if not np.any(used_bins):
        print(f"[Warning] No LD pairs within max distance {max_distance_kb} kb.")
        empty_grouped = pd.DataFrame(
            columns=["bin_start_kb", "bin_end_kb", "bin_mid_kb", "mean_r2", "median_r2", "n_pairs"]
        )
        save_df_to_tsv(empty_grouped, f"{output_prefix}.info.tsv")
        save_thresholds(
            {
                "Total_LD_pairs": int(total_pairs),
                "Pairs_Used_for_Decay": 0,
                "Max_Distance_kb": float(max_distance_kb),
                "Bin_Size_kb": float(bin_size_kb),
                "First_Bin_Mean_R2": np.nan,
                "Half_Decay_R2_Threshold": np.nan,
                "Half_Decay_Distance_kb": np.nan,
            },
            f"{output_prefix}.th.tsv",
        )
        empty_scatter = pd.DataFrame(columns=["distance_kb", "R2"])
        plot_multi_line_series(
            data=empty_grouped,
            x_col="bin_mid_kb",
            y_specs=[
                {"y_col": "mean_r2", "label": "Mean r2", "color": "steelblue", "linestyle": "-", "linewidth": 2},
                {"y_col": "median_r2", "label": "Median r2", "color": "darkorange", "linestyle": "--", "linewidth": 1.5},
            ],
            title="LD decay (binned mean and median r2 vs distance, PLINK2)",
            filename=f"{output_prefix}.png",
            x_label="Distance (kb)",
            y_label="LD (r2)",
        )
        plot_scatter_with_thresholds(
            data=empty_scatter,
            x_col="distance_kb",
            y_col="R2",
            title="LD vs distance (sampled pairs, PLINK2)",
            filename=f"{output_prefix}.scatter.png",
            xlabel="Distance (kb)",
            ylabel="LD (r2)",
            color="steelblue",
            alpha=0.25,
            s=6,
            figure_size=(10, 6),
            xlim=(0, float(max_distance_kb)),
            ylim=(0, 1),
        )
        plot_distribution_with_stats(
            data=empty_scatter,
            col="R2",
            title="LD r2 Distribution (PLINK2)",
            filename=f"{output_prefix}.r2_dist.png",
            mean_val=np.nan,
            median_val=np.nan,
            x_label="LD (r2)",
            y_label="Count",
            bins=100,
            xlim=(0, 1),
        )
        plot_distribution_with_stats(
            data=empty_scatter,
            col="R2",
            title="LD r2 Distribution (PLINK2) - Log Y",
            filename=f"{output_prefix}.r2_dist.log.png",
            mean_val=np.nan,
            median_val=np.nan,
            x_label="LD (r2)",
            y_label="Count",
            bins=100,
            xlim=(0, 1),
            log_scale=True,
        )
        return

    rows = []
    for i in range(n_bins):
        if cnt_r2[i] == 0:
            continue
        bin_start = float(i * bin_size_kb)
        bin_end = float((i + 1) * bin_size_kb)
        bin_mid = (bin_start + bin_end) / 2.0
        mean_r2 = float(sum_r2[i] / cnt_r2[i])
        sampled = median_samples[i]
        median_r2 = float(np.median(sampled)) if len(sampled) > 0 else np.nan
        rows.append([bin_start, bin_end, bin_mid, mean_r2, median_r2, int(cnt_r2[i])])

    grouped = pd.DataFrame(
        rows,
        columns=["bin_start_kb", "bin_end_kb", "bin_mid_kb", "mean_r2", "median_r2", "n_pairs"],
    )
    save_df_to_tsv(grouped, f"{output_prefix}.info.tsv")

    # Keep threshold summary consistent with existing stats outputs.
    first_mean_r2 = float(grouped["mean_r2"].iloc[0]) if len(grouped) > 0 else np.nan
    half_threshold = first_mean_r2 * 0.5 if pd.notna(first_mean_r2) else np.nan
    half_decay_rows = grouped[grouped["mean_r2"] <= half_threshold] if pd.notna(half_threshold) else pd.DataFrame()
    half_decay_kb = (
        float(half_decay_rows["bin_mid_kb"].iloc[0])
        if len(half_decay_rows) > 0
        else np.nan
    )
    th_dict = {
        "Total_LD_pairs": int(total_pairs),
        "Pairs_Used_for_Decay": int(grouped["n_pairs"].sum()),
        "Max_Distance_kb": float(max_distance_kb),
        "Bin_Size_kb": float(bin_size_kb),
        "First_Bin_Mean_R2": first_mean_r2,
        "Half_Decay_R2_Threshold": half_threshold,
        "Half_Decay_Distance_kb": half_decay_kb,
    }
    save_thresholds(th_dict, f"{output_prefix}.th.tsv")

    scatter_df = pd.DataFrame(scatter_pairs, columns=["distance_kb", "R2"])
    r2_mean = float(scatter_df["R2"].mean()) if len(scatter_df) > 0 else np.nan
    r2_median = float(scatter_df["R2"].median()) if len(scatter_df) > 0 else np.nan

    plot_multi_line_series(
        data=grouped,
        x_col="bin_mid_kb",
        y_specs=[
            {"y_col": "mean_r2", "label": "Mean r2", "color": "steelblue", "linestyle": "-", "linewidth": 2},
            {"y_col": "median_r2", "label": "Median r2", "color": "darkorange", "linestyle": "--", "linewidth": 1.5},
        ],
        title="LD decay (binned mean and median r2 vs distance, PLINK2)",
        filename=f"{output_prefix}.png",
        x_label="Distance (kb)",
        y_label="LD (r2)",
    )
    plot_scatter_with_thresholds(
        data=scatter_df,
        x_col="distance_kb",
        y_col="R2",
        title="LD vs distance (sampled pairs, PLINK2)",
        filename=f"{output_prefix}.scatter.png",
        xlabel="Distance (kb)",
        ylabel="LD (r2)",
        color="steelblue",
        alpha=0.25,
        s=6,
        figure_size=(10, 6),
        xlim=(0, float(max_distance_kb)),
        ylim=(0, 1),
    )
    plot_distribution_with_stats(
        data=scatter_df,
        col="R2",
        title="LD r2 Distribution (chunk data, parallel)",
        filename=f"{output_prefix}.r2_dist.png",
        mean_val=r2_mean,
        median_val=r2_median,
        x_label="LD (r2)",
        y_label="Count",
        bins=100,
        xlim=(0, 1),
    )
    plot_distribution_with_stats(
        data=scatter_df,
        col="R2",
        title="LD r2 Distribution (chunk data, parallel) - Log Y",
        filename=f"{output_prefix}.r2_dist.log.png",
        mean_val=r2_mean,
        median_val=r2_median,
        x_label="LD (r2)",
        y_label="Count",
        bins=100,
        xlim=(0, 1),
        log_scale=True,
    )


def ana_ld_crosschr_baseline(
    input_file: str,
    output_prefix: str,
    chunksize: int = 2_000_000,
) -> None:
    """
    Summarize cross-chromosome LD baseline from PLINK2 .vcor output.
    """
    r2_values = []
    total_pairs = 0
    try:
        reader = pd.read_csv(input_file, sep=r"\s+", chunksize=chunksize)
    except Exception as e:
        print(f"[Error] Failed to read cross-chr LD file {input_file}: {e}")
        return

    for chunk in reader:
        local = _normalize_ld_chunk(chunk)
        if local is None or local.empty:
            continue
        r2 = local["R2"].to_numpy(dtype=np.float64)
        r2 = r2[np.isfinite(r2)]
        if len(r2) == 0:
            continue
        r2_values.append(r2)
        total_pairs += int(len(r2))

    if len(r2_values) == 0:
        empty_df = pd.DataFrame(columns=["metric", "value"])
        save_df_to_tsv(empty_df, f"{output_prefix}.info.tsv")
        save_thresholds(
            {
                "Total_CrossChr_Pairs": 0,
                "Mean_R2": np.nan,
                "Median_R2": np.nan,
                "P95_R2": np.nan,
                "P99_R2": np.nan,
            },
            f"{output_prefix}.th.tsv",
        )
        empty_plot_df = pd.DataFrame(columns=["R2"])
        plot_distribution_with_stats(
            data=empty_plot_df,
            col="R2",
            title="Cross-chromosome LD r2 Distribution (PLINK2)",
            filename=f"{output_prefix}.png",
            mean_val=np.nan,
            median_val=np.nan,
            x_label="LD (r2)",
            y_label="Count",
            bins=100,
            xlim=(0, 1),
        )
        plot_distribution_with_stats(
            data=empty_plot_df,
            col="R2",
            title="Cross-chromosome LD r2 Distribution (PLINK2) - Log Y",
            filename=f"{output_prefix}.dist.log.png",
            mean_val=np.nan,
            median_val=np.nan,
            x_label="LD (r2)",
            y_label="Count",
            bins=100,
            xlim=(0, 1),
            log_scale=True,
        )
        return

    r2_all = np.concatenate(r2_values)
    mean_r2 = float(np.mean(r2_all))
    median_r2 = float(np.median(r2_all))
    p95_r2 = float(np.quantile(r2_all, 0.95))
    p99_r2 = float(np.quantile(r2_all, 0.99))
    gt01 = int(np.sum(r2_all >= 0.1))
    gt02 = int(np.sum(r2_all >= 0.2))

    info_df = pd.DataFrame(
        [
            ["Total_CrossChr_Pairs", int(total_pairs)],
            ["Mean_R2", mean_r2],
            ["Median_R2", median_r2],
            ["P95_R2", p95_r2],
            ["P99_R2", p99_r2],
            ["Pairs_R2_ge_0.1", gt01],
            ["Pairs_R2_ge_0.2", gt02],
        ],
        columns=["metric", "value"],
    )
    save_df_to_tsv(info_df, f"{output_prefix}.info.tsv")
    save_thresholds(
        {
            "Total_CrossChr_Pairs": int(total_pairs),
            "Mean_R2": mean_r2,
            "Median_R2": median_r2,
            "P95_R2": p95_r2,
            "P99_R2": p99_r2,
            "Pairs_R2_ge_0.1": gt01,
            "Pairs_R2_ge_0.2": gt02,
        },
        f"{output_prefix}.th.tsv",
    )
    r2_df = pd.DataFrame({"R2": r2_all})
    plot_distribution_with_stats(
        data=r2_df,
        col="R2",
        title="Cross-chromosome LD r2 Distribution (PLINK2)",
        filename=f"{output_prefix}.png",
        mean_val=mean_r2,
        median_val=median_r2,
        x_label="LD (r2)",
        y_label="Count",
        bins=100,
        xlim=(0, 1),
    )
    plot_distribution_with_stats(
        data=r2_df,
        col="R2",
        title="Cross-chromosome LD r2 Distribution (PLINK2) - Log Y",
        filename=f"{output_prefix}.dist.log.png",
        mean_val=mean_r2,
        median_val=median_r2,
        x_label="LD (r2)",
        y_label="Count",
        bins=100,
        xlim=(0, 1),
        log_scale=True,
    )
