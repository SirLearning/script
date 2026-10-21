"""Local relatedness ranking from published SINGER ARG texts (no tskit)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from genetics.genomics.variant.singer_arg import summarize_singer_windows


def load_keep_iids(keep_path: str | Path) -> list[str]:
    iids: list[str] = []
    for line in Path(keep_path).open(encoding="utf-8"):
        parts = line.strip().split()
        if len(parts) >= 2:
            iids.append(parts[1])
        elif parts:
            iids.append(parts[0])
    if not iids:
        raise ValueError(f"no sample IDs in {keep_path}")
    return iids


def load_singer_nodes(path: str | Path) -> np.ndarray:
    return np.loadtxt(path, dtype=np.float64)


def load_singer_branches(path: str | Path) -> np.ndarray:
    """Rows are left, right, parent, child. Parent < 0 (virtual root) is dropped."""
    edges = np.loadtxt(path, dtype=np.float64)
    if edges.ndim == 1:
        edges = edges.reshape(1, -1)
    return edges[edges[:, 2] >= 0.0]


def parent_at_position(edges: np.ndarray, n_nodes: int, position: float) -> np.ndarray:
    parent = np.full(n_nodes, -1, dtype=np.int32)
    left = edges[:, 0]
    right = edges[:, 1]
    cover = (left <= position) & (position < right)
    for row in edges[cover]:
        parent[int(row[3])] = int(row[2])
    return parent


def hap_tmrca_matrix(parent: np.ndarray, times: np.ndarray, n_hap: int) -> np.ndarray:
    """Pairwise TMRCA for haplotypes 0..n_hap-1 on one local tree."""
    children: list[list[int]] = [[] for _ in range(len(times))]
    for child, par in enumerate(parent):
        if par >= 0:
            children[par].append(int(child))
    tmrca = np.zeros((n_hap, n_hap), dtype=np.float64)
    roots = [i for i, par in enumerate(parent) if par < 0 and children[i]]
    if not roots:
        roots = [int(np.argmax(times))]
    for root in roots:
        stack: list[tuple[int, bool]] = [(root, False)]
        leaf_sets: dict[int, list[int]] = {}
        while stack:
            node, seen = stack.pop()
            kids = children[node]
            if not kids:
                leaf_sets[node] = [node] if node < n_hap else []
                continue
            if not seen:
                stack.append((node, True))
                for kid in reversed(kids):
                    stack.append((kid, False))
                continue
            merged: list[int] = []
            for kid in kids:
                leafs = leaf_sets[kid]
                if merged and leafs:
                    left = np.asarray(merged, dtype=np.int32)
                    right = np.asarray(leafs, dtype=np.int32)
                    tmrca[np.ix_(left, right)] = times[node]
                    tmrca[np.ix_(right, left)] = times[node]
                merged.extend(leafs)
            leaf_sets[node] = merged
    np.fill_diagonal(tmrca, 0.0)
    return tmrca


def sample_tmrca_from_hap(hap_tmrca: np.ndarray) -> np.ndarray:
    """Diploid sample TMRCA = mean of the four inter-sample haplotype pairs."""
    n_hap = hap_tmrca.shape[0]
    if n_hap % 2:
        raise ValueError("haplotype count must be even (two haps per sample)")
    n_s = n_hap // 2
    blocks = hap_tmrca.reshape(n_s, 2, n_s, 2)
    sample = blocks.mean(axis=(1, 3))
    np.fill_diagonal(sample, 0.0)
    return sample


def _window_positions(edges: np.ndarray, n_positions: int) -> np.ndarray:
    right = float(edges[:, 1].max()) if len(edges) else 1.0
    left = float(edges[:, 0].min()) if len(edges) else 0.0
    if right <= left:
        return np.array([left], dtype=np.float64)
    width = right - left
    if n_positions <= 1:
        return np.array([left + 0.5 * width], dtype=np.float64)
    # Midpoints of n equal bins in [left, right).
    return left + (np.arange(n_positions, dtype=np.float64) + 0.5) * (width / n_positions)


def sample_tmrca_from_arg(
    singer_dir: str | Path,
    window_id: str,
    *,
    sample_ids: list[str],
    n_positions: int = 11,
    arg_sample: int = 0,
) -> np.ndarray:
    """Mean sample TMRCA across ``n_positions`` local trees of one ARG sample."""
    folder = Path(singer_dir)
    nodes_path = folder / f"singer_{window_id}_nodes_{arg_sample}.txt"
    branches_path = folder / f"singer_{window_id}_branches_{arg_sample}.txt"
    if not nodes_path.is_file() or not branches_path.is_file():
        raise FileNotFoundError(f"missing SINGER texts for {window_id} sample {arg_sample}")
    times = load_singer_nodes(nodes_path)
    edges = load_singer_branches(branches_path)
    n_hap = int(np.sum(times == 0.0))
    if n_hap != 2 * len(sample_ids):
        raise ValueError(
            f"{window_id}: {n_hap} haplotypes vs 2×{len(sample_ids)} keep IDs"
        )
    acc = np.zeros((len(sample_ids), len(sample_ids)), dtype=np.float64)
    for pos in _window_positions(edges, n_positions):
        parent = parent_at_position(edges, len(times), float(pos))
        hap = hap_tmrca_matrix(parent, times, n_hap)
        acc += sample_tmrca_from_hap(hap)
    acc /= max(n_positions, 1)
    np.fill_diagonal(acc, 0.0)
    return acc


def mean_tmrca_across_windows(
    windows: list[dict] | str | Path,
    singer_dir: str | Path,
    sample_ids: list[str],
    *,
    n_positions: int = 11,
    arg_sample: int = 0,
) -> tuple[np.ndarray, list[str]]:
    """Average sample TMRCA over windows that have ``nodes_0``."""
    rows = summarize_singer_windows(windows, singer_dir)
    used: list[str] = []
    acc = np.zeros((len(sample_ids), len(sample_ids)), dtype=np.float64)
    for row in rows:
        if not row["arg_ok"]:
            continue
        window_id = str(row["window_id"])
        acc += sample_tmrca_from_arg(
            singer_dir,
            window_id,
            sample_ids=sample_ids,
            n_positions=n_positions,
            arg_sample=arg_sample,
        )
        used.append(window_id)
    if not used:
        raise ValueError(f"no ARG windows under {singer_dir}")
    acc /= len(used)
    np.fill_diagonal(acc, 0.0)
    return acc, used


def rank_nearest_neighbors(
    matrix: np.ndarray,
    sample_ids: list[str],
    *,
    top_k: int = 5,
) -> list[dict]:
    """Per-sample nearest others by smaller TMRCA."""
    n = len(sample_ids)
    k = min(top_k, n - 1)
    if k <= 0:
        return []
    work = np.array(matrix, dtype=np.float64, copy=True)
    np.fill_diagonal(work, np.inf)
    rows: list[dict] = []
    for i, sample in enumerate(sample_ids):
        order = np.argsort(work[i], kind="stable")[:k]
        for rank, j in enumerate(order, start=1):
            rows.append(
                {
                    "sample": sample,
                    "rank": rank,
                    "neighbor": sample_ids[int(j)],
                    "tmrca": float(work[i, j]),
                }
            )
    return rows


def closest_pairs(matrix: np.ndarray, sample_ids: list[str], *, top_pairs: int = 200) -> list[dict]:
    n = len(sample_ids)
    iu = np.triu_indices(n, k=1)
    vals = matrix[iu]
    take = min(top_pairs, vals.size)
    order = np.argsort(vals, kind="stable")[:take]
    rows: list[dict] = []
    for rank, idx in enumerate(order, start=1):
        i = int(iu[0][idx])
        j = int(iu[1][idx])
        rows.append(
            {
                "rank": rank,
                "sample_a": sample_ids[i],
                "sample_b": sample_ids[j],
                "tmrca": float(matrix[i, j]),
            }
        )
    return rows


def plot_singer_arg_relatedness(
    windows: str | Path,
    singer_dir: str | Path,
    keep_path: str | Path,
    output_prefix: str | Path,
    *,
    n_positions: int = 11,
    top_k: int = 5,
    top_pairs: int = 200,
    arg_sample: int = 0,
) -> dict:
    """Write neighbor/pair ranks and TMRCA histograms from published ARGs."""
    import pandas as pd

    from infra.utils.graph import plot_distribution_with_stats
    from infra.utils.io import save_df_to_tsv, save_thresholds

    sample_ids = load_keep_iids(keep_path)
    matrix, used = mean_tmrca_across_windows(
        windows,
        singer_dir,
        sample_ids,
        n_positions=n_positions,
        arg_sample=arg_sample,
    )
    neighbors = rank_nearest_neighbors(matrix, sample_ids, top_k=top_k)
    pairs = closest_pairs(matrix, sample_ids, top_pairs=top_pairs)
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    save_df_to_tsv(pd.DataFrame(neighbors), f"{prefix}.neighbor.info.tsv", float_format="%.6g")
    save_df_to_tsv(pd.DataFrame(pairs), f"{prefix}.pair_rank.info.tsv", float_format="%.6g")

    iu = np.triu_indices(len(sample_ids), k=1)
    pair_t = matrix[iu]
    nn_t = np.array([row["tmrca"] for row in neighbors if row["rank"] == 1], dtype=np.float64)
    stats = {
        "n_samples": len(sample_ids),
        "n_windows": len(used),
        "n_positions": n_positions,
        "arg_sample": arg_sample,
        "median_pair_tmrca": float(np.median(pair_t)) if pair_t.size else "",
        "median_nearest_tmrca": float(np.median(nn_t)) if nn_t.size else "",
        "min_pair_tmrca": float(np.min(pair_t)) if pair_t.size else "",
    }
    save_thresholds(stats, f"{prefix}.relatedness.th.tsv")

    pair_df = pd.DataFrame({"tmrca": pair_t})
    plot_distribution_with_stats(
        pair_df,
        "tmrca",
        title="Pairwise mean TMRCA across ARG windows",
        filename=f"{prefix}.pair_tmrca.dist.png",
        x_label="Mean TMRCA (SINGER units)",
        y_label="Sample pairs",
        bins=min(80, max(pair_t.size, 1)),
        mean_val=float(np.mean(pair_t)) if pair_t.size else None,
        median_val=float(np.median(pair_t)) if pair_t.size else None,
        log_scale=True,
    )
    if nn_t.size:
        nn_df = pd.DataFrame({"tmrca": nn_t})
        plot_distribution_with_stats(
            nn_df,
            "tmrca",
            title="Nearest-neighbor TMRCA per sample",
            filename=f"{prefix}.nearest_tmrca.dist.png",
            x_label="Nearest-neighbor TMRCA (SINGER units)",
            y_label="Samples",
            bins=min(40, max(nn_t.size, 1)),
            mean_val=float(np.mean(nn_t)),
            median_val=float(np.median(nn_t)),
        )
    return {
        "n_samples": len(sample_ids),
        "n_windows": len(used),
        "used_windows": used,
    }
