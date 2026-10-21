"""Along-chromosome LD: merge PLINK arms, then distance where mean r2 falls to a threshold."""

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.variant.ld_decay import _iter_plink2_vcor_chunks, _normalize_ld_chunk
from infra.utils.graph import plot_scatter_with_thresholds
from infra.utils.io import load_df_from_tsv, save_df_to_tsv, save_thresholds
from infra.wheat.ref_v1 import (
    get_centromere_end,
    get_centromere_start,
    get_ref_chromosome_length,
    is_ref_v1_chr_name,
    plink_chr_ids_for_ref_chrom,
    plink_segment_to_ref_position,
)

R2_DECAY_THRESHOLD = 0.2
DECAY_PLOT_YLIM_PAD = 1.15
DECAY_PLOT_YLIM_MIN_KB = 50.0


def ld_along_chr_output_stem(r2_threshold: float = R2_DECAY_THRESHOLD) -> str:
    """Plot/table stem. Default r2=0.2 keeps the historical filename."""
    value = float(r2_threshold)
    if abs(value - R2_DECAY_THRESHOLD) < 1e-12:
        return "ld_r2_decay_along"
    text = f"{value:g}".replace(".", "p")
    return f"ld_r2_decay_along.r2_{text}"


def _r2_label(r2_threshold: float) -> str:
    return f"{float(r2_threshold):g}"


def _ref_chrom_token(chromosome: str) -> str:
    token = str(chromosome)
    return token[3:] if token.startswith("chr") else token


def subgenome_ref_chroms(subgenome: str) -> list[str]:
    """Return chr1–chr7 names for one wheat subgenome (A, B, or D)."""
    token = _ref_chrom_token(str(subgenome)).strip().upper()
    sg = token[-1] if token else ""
    if sg not in {"A", "B", "D"}:
        raise ValueError(f"Unknown subgenome: {subgenome}")
    return [f"{i}{sg}" for i in range(1, 8)]


def load_pvar_ref_variants(pvar_path: str) -> pd.DataFrame:
    """Load a PLINK2 pvar and attach linear vmap4 reference coordinates.

    Accepts every coding that exists on disk:

    - arm-numbered source pfiles (``#CHROM`` 0-44, POS local to the arm);
    - ref-chr pfiles recoded before the naming change, which carry the first arm's
      id with whole-chromosome POS -- harmless, because a first arm has ``ref_start``
      0 and the segment-to-ref mapping is the identity there;
    - ref-chr pfiles recoded after it, which carry ``chr2A`` with whole-chromosome
      POS and must not be passed through the segment mapping at all.
    """
    rows: list[dict] = []
    with open(pvar_path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            raw_chrom = parts[0].strip()
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            if is_ref_v1_chr_name(raw_chrom):
                plink_chr: int | str = raw_chrom
                ref_chrom = raw_chrom[3:] if raw_chrom.startswith("chr") else raw_chrom
                ref_pos = pos
            else:
                try:
                    plink_chr = int(raw_chrom)
                except ValueError:
                    continue
                try:
                    ref_chrom, ref_pos = plink_segment_to_ref_position(plink_chr, pos)
                except ValueError:
                    continue
            rows.append(
                {
                    "ID": parts[2],
                    "plink_chr": plink_chr,
                    "pos": pos,
                    "ref_chrom": ref_chrom,
                    "ref_pos": ref_pos,
                }
            )
    if not rows:
        return pd.DataFrame(
            columns=["ID", "plink_chr", "pos", "ref_chrom", "ref_pos"]
        )
    return pd.DataFrame(rows)


def ref_merge_chr_label(chromosome: str) -> str | None:
    """Merged ``#CHROM`` label for a vmap4 chromosome, or None if unknown.

    The recode writes whole-chromosome names (``chr2A``) rather than the first arm's
    numeric id, so a merged pfile is self-describing and the recode is idempotent.
    Organelles keep their bare names (``Mit``, ``Chl``); everything else gains ``chr``.
    """
    token = _ref_chrom_token(chromosome)
    for candidate in (token, f"chr{token}"):
        if is_ref_v1_chr_name(candidate):
            return candidate
    return None


def plink_chr_tokens_for_ref_chrom(
    variants: pd.DataFrame,
    chromosome: str,
) -> list[str]:
    """``--chr`` tokens present for one ref chromosome, for either pvar coding.

    A merged pfile contributes the single whole-chromosome label; an arm-numbered
    pfile contributes the arm ids that actually carry variants.
    """
    chrom = _ref_chrom_token(chromosome)
    if variants is None or variants.empty:
        return []
    present = set(variants.loc[variants["ref_chrom"] == chrom, "plink_chr"].astype(str))
    merged = ref_merge_chr_label(chrom)
    if merged is not None and merged in present:
        return [merged]
    return [cid for cid in plink_chr_ids_for_ref_chrom(chrom) if cid in present]


def write_ref_merge_update_maps(
    pvar_path: str,
    chromosome: str,
    chr_map_path: str,
    pos_map_path: str,
) -> pd.DataFrame:
    """Write plink2 ``--update-chr`` / ``--update-map`` files for one ref chromosome.

    Both PLINK arms (e.g. 13+14 for chr3A) are recoded onto the whole-chromosome
    label with linear reference POS so ``--r2-unphased`` treats them as one chromosome.
    """
    chrom = _ref_chrom_token(chromosome)
    variants = load_pvar_ref_variants(pvar_path)
    sub = variants[variants["ref_chrom"] == chrom].copy()
    if sub.empty:
        raise ValueError(f"No variants for chromosome {chrom} in {pvar_path}")
    merge_chr = ref_merge_chr_label(chrom)
    if merge_chr is None:
        raise ValueError(f"No merge chromosome label for {chrom}")
    Path(chr_map_path).parent.mkdir(parents=True, exist_ok=True)
    Path(pos_map_path).parent.mkdir(parents=True, exist_ok=True)
    with open(chr_map_path, "w", encoding="utf-8") as chr_out, open(
        pos_map_path, "w", encoding="utf-8"
    ) as pos_out:
        for row in sub.itertuples(index=False):
            chr_out.write(f"{row.ID}\t{merge_chr}\n")
            pos_out.write(f"{row.ID}\t{int(row.ref_pos)}\n")
    sub["merge_chr"] = merge_chr
    print(
        f"[Info] Ref-merge maps for chr{chrom}: {len(sub)} variants, "
        f"PLINK arms {','.join(plink_chr_ids_for_ref_chrom(chrom))} -> chr {merge_chr}"
    )
    return sub


def write_all_ref_merge_update_maps(
    pvar_path: str,
    chr_map_path: str,
    pos_map_path: str,
) -> pd.DataFrame:
    """Write plink2 update maps for every wheat chromosome that still has split arms.

    Both arms are recoded onto the whole-chromosome label (``chr3A``) with linear
    reference POS. Identity rows (already merged, Un, Mit, Chl) are omitted so
    plink2 leaves those sites unchanged, which makes a second pass a no-op.
    """
    variants = load_pvar_ref_variants(pvar_path)
    Path(chr_map_path).parent.mkdir(parents=True, exist_ok=True)
    Path(pos_map_path).parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    with open(chr_map_path, "w", encoding="utf-8") as chr_out, open(
        pos_map_path, "w", encoding="utf-8"
    ) as pos_out:
        for row in variants.itertuples(index=False):
            merge_chr = ref_merge_chr_label(str(row.ref_chrom))
            if merge_chr is None:
                continue
            new_pos = int(row.ref_pos)
            if str(row.plink_chr) == merge_chr and int(row.pos) == new_pos:
                continue
            chr_out.write(f"{row.ID}\t{merge_chr}\n")
            pos_out.write(f"{row.ID}\t{new_pos}\n")
            rows.append(
                {
                    "ID": row.ID,
                    "plink_chr": str(row.plink_chr),
                    "pos": int(row.pos),
                    "ref_chrom": str(row.ref_chrom),
                    "ref_pos": new_pos,
                    "merge_chr": merge_chr,
                }
            )
    print(
        f"[Info] Genome-wide ref-merge maps: {len(rows)} variants recoded "
        f"from {pvar_path}"
    )
    return pd.DataFrame(rows)


def evenly_sample_variants_along_chromosome(
    variants: pd.DataFrame,
    n: int,
    *,
    pos_col: str = "ref_pos",
    id_col: str = "ID",
) -> pd.DataFrame:
    """Pick ``n`` variants nearest to evenly spaced target positions (no replacement)."""
    if variants is None or variants.empty or int(n) <= 0:
        return pd.DataFrame(columns=variants.columns if variants is not None else [])
    df = variants.dropna(subset=[pos_col, id_col]).drop_duplicates(subset=[id_col])
    df = df.sort_values(pos_col).reset_index(drop=True)
    if len(df) <= int(n):
        return df
    pos = df[pos_col].to_numpy(dtype=np.float64)
    targets = np.linspace(pos[0], pos[-1], int(n))
    used: set[int] = set()
    picked: list = []
    for target in targets:
        for idx in np.argsort(np.abs(pos - target), kind="stable"):
            i = int(idx)
            if i in used:
                continue
            used.add(i)
            picked.append(df.iloc[i])
            break
    return pd.DataFrame(picked).sort_values(pos_col).reset_index(drop=True)


def _spatial_bin_bounds(bin_idx: int, bin_bp: int, chrom_length_bp: int) -> tuple[int, int, float]:
    start = int(bin_idx) * int(bin_bp) + 1
    end = min((int(bin_idx) + 1) * int(bin_bp), int(chrom_length_bp))
    mid_mb = (start + end) / 2.0 / 1e6
    return start, end, mid_mb


def sample_focals_per_spatial_bin(
    variants: pd.DataFrame,
    *,
    n_per_bin: int = 100,
    bin_mb: float = 5.0,
    chrom_length_bp: int | None = None,
    pos_col: str = "ref_pos",
) -> pd.DataFrame:
    """Evenly sample up to ``n_per_bin`` focals in each ``bin_mb`` window along the chromosome.

    LD partners are *not* restricted to the window; this only chooses focals. Bins with
    fewer than ``n_per_bin`` sites keep every variant.
    """
    if variants is None or variants.empty or int(n_per_bin) <= 0 or float(bin_mb) <= 0:
        return pd.DataFrame()
    df = variants.dropna(subset=[pos_col, "ID"]).drop_duplicates(subset=["ID"]).copy()
    df = df.sort_values(pos_col).reset_index(drop=True)
    bin_bp = int(round(float(bin_mb) * 1_000_000))
    pos = df[pos_col].to_numpy(dtype=np.int64)
    length = int(chrom_length_bp) if chrom_length_bp is not None else int(pos.max())
    length = max(length, int(pos.max()) if len(pos) else 1)
    df["bin_idx"] = np.floor((np.maximum(pos, 1) - 1) / bin_bp).astype(int)
    frames: list[pd.DataFrame] = []
    for bin_idx, sub in df.groupby("bin_idx", sort=True):
        sampled = evenly_sample_variants_along_chromosome(
            sub.drop(columns=["bin_idx"]), int(n_per_bin), pos_col=pos_col
        )
        if sampled.empty:
            continue
        start, end, mid_mb = _spatial_bin_bounds(int(bin_idx), bin_bp, length)
        sampled = sampled.copy()
        sampled["bin_idx"] = int(bin_idx)
        sampled["bin_start_bp"] = start
        sampled["bin_end_bp"] = end
        sampled["bin_mid_mb"] = mid_mb
        sampled["n_in_bin"] = int(len(sub))
        frames.append(sampled)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def sample_focals_per_spatial_bin_from_pvar(
    pvar_path: str,
    chromosomes: list[str],
    *,
    n_per_bin: int = 100,
    bin_mb: float = 5.0,
    output_prefix: str | None = None,
) -> pd.DataFrame:
    """Sample focals in ``bin_mb`` windows on each vmap4 chromosome (both PLINK arms)."""
    variants = load_pvar_ref_variants(pvar_path)
    frames: list[pd.DataFrame] = []
    for raw in chromosomes:
        chrom = _ref_chrom_token(raw)
        sub = variants[variants["ref_chrom"] == chrom]
        try:
            length = get_ref_chromosome_length(chrom)
        except ValueError:
            length = int(sub["ref_pos"].max()) if not sub.empty else 1
        sampled = sample_focals_per_spatial_bin(
            sub,
            n_per_bin=int(n_per_bin),
            bin_mb=float(bin_mb),
            chrom_length_bp=length,
        )
        if sampled.empty:
            print(f"[Warning] No variants to sample on chr{chrom}")
            continue
        frames.append(sampled)
        n_bins = sampled["bin_idx"].nunique()
        print(
            f"[Info] chr{chrom}: {len(sampled)} focals in {n_bins} x {bin_mb:g} Mb bins "
            f"(target {n_per_bin}/bin; LD window is not limited to the bin)"
        )
    if not frames:
        out = pd.DataFrame(
            columns=[
                "ID",
                "plink_chr",
                "pos",
                "ref_chrom",
                "ref_pos",
                "bin_idx",
                "bin_start_bp",
                "bin_end_bp",
                "bin_mid_mb",
                "n_in_bin",
            ]
        )
    else:
        out = pd.concat(frames, ignore_index=True)
    if output_prefix:
        save_df_to_tsv(out, f"{output_prefix}.focals.info.tsv")
        with open(f"{output_prefix}.focals.id", "w", encoding="utf-8") as handle:
            for vid in out["ID"].astype(str):
                handle.write(f"{vid}\n")
        ids: list[str] = []
        for raw in chromosomes:
            ids.extend(plink_chr_tokens_for_ref_chrom(variants, raw))
        Path(f"{output_prefix}.plink_chr_ids.txt").write_text(
            ",".join(dict.fromkeys(ids)),
            encoding="utf-8",
        )
    return out


def sample_even_focals_from_pvar(
    pvar_path: str,
    chromosomes: list[str],
    n_per_chr: int,
    output_prefix: str | None = None,
) -> pd.DataFrame:
    """Evenly sample variants along each vmap4 chromosome (both PLINK arms)."""
    variants = load_pvar_ref_variants(pvar_path)
    frames: list[pd.DataFrame] = []
    for raw in chromosomes:
        chrom = _ref_chrom_token(raw)
        sub = variants[variants["ref_chrom"] == chrom]
        sampled = evenly_sample_variants_along_chromosome(sub, int(n_per_chr))
        if sampled.empty:
            print(f"[Warning] No variants to sample on chr{chrom}")
            continue
        frames.append(sampled)
    if not frames:
        out = pd.DataFrame(columns=["ID", "plink_chr", "pos", "ref_chrom", "ref_pos"])
    else:
        out = pd.concat(frames, ignore_index=True)
    if output_prefix:
        save_df_to_tsv(out, f"{output_prefix}.focals.info.tsv")
        with open(f"{output_prefix}.focals.id", "w", encoding="utf-8") as handle:
            for vid in out["ID"].astype(str):
                handle.write(f"{vid}\n")
        ids: list[str] = []
        for raw in chromosomes:
            ids.extend(plink_chr_tokens_for_ref_chrom(variants, raw))
        Path(f"{output_prefix}.plink_chr_ids.txt").write_text(
            ",".join(dict.fromkeys(ids)),
            encoding="utf-8",
        )
    return out


def mean_r2_decay_distance_kb(
    distance_kb,
    r2,
    *,
    threshold: float = R2_DECAY_THRESHOLD,
    bin_size_kb: float = 50.0,
    min_pairs_per_bin: int = 5,
    max_distance_kb: float | None = None,
) -> tuple[float, bool, int]:
    """Distance (kb) where binned mean r2 first falls to ``threshold`` (interpolated)."""
    dist = np.asarray(distance_kb, dtype=np.float64)
    r2_arr = np.asarray(r2, dtype=np.float64)
    ok = np.isfinite(dist) & np.isfinite(r2_arr) & (dist >= 0.0)
    if max_distance_kb is not None:
        ok &= dist <= float(max_distance_kb)
    dist = dist[ok]
    r2_arr = r2_arr[ok]
    n_pairs = int(dist.size)
    if n_pairs == 0 or float(bin_size_kb) <= 0:
        return float("nan"), False, n_pairs
    max_d = float(max_distance_kb) if max_distance_kb is not None else float(np.max(dist))
    n_bins = int(np.ceil(max_d / float(bin_size_kb)))
    if n_bins <= 0:
        return float("nan"), False, n_pairs
    bin_idx = np.floor(dist / float(bin_size_kb)).astype(int)
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)
    sums = np.zeros(n_bins, dtype=np.float64)
    cnts = np.zeros(n_bins, dtype=np.int64)
    np.add.at(sums, bin_idx, r2_arr)
    np.add.at(cnts, bin_idx, 1)
    mids = (np.arange(n_bins) + 0.5) * float(bin_size_kb)
    means = np.full(n_bins, np.nan, dtype=np.float64)
    valid = cnts >= int(min_pairs_per_bin)
    means[valid] = sums[valid] / cnts[valid]
    prev_mid = None
    prev_mean = None
    thr = float(threshold)
    for i in range(n_bins):
        if not np.isfinite(means[i]):
            continue
        if means[i] <= thr:
            if prev_mid is None:
                return float(mids[i]), True, n_pairs
            y0 = float(prev_mean)
            y1 = float(means[i])
            x0 = float(prev_mid)
            x1 = float(mids[i])
            if y0 == y1:
                return x1, True, n_pairs
            t = (y0 - thr) / (y0 - y1)
            t = min(max(t, 0.0), 1.0)
            return x0 + t * (x1 - x0), True, n_pairs
        prev_mid = float(mids[i])
        prev_mean = float(means[i])
    return float("nan"), False, n_pairs


def _load_focal_pairs(vcor_file: str, focal_ids: set[str]) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for chunk in _iter_plink2_vcor_chunks(vcor_file, chunksize=500_000):
        local = _normalize_ld_chunk(chunk, keep_ids=True)
        if local is None or local.empty:
            continue
        id_a = local["ID_A"].astype(str)
        id_b = local["ID_B"].astype(str)
        mask = id_a.isin(focal_ids) | id_b.isin(focal_ids)
        if not mask.any():
            continue
        rows.append(local.loc[mask])
    if not rows:
        return pd.DataFrame(columns=["ID_A", "ID_B", "BP_A", "BP_B", "R2"])
    return pd.concat(rows, ignore_index=True)


def _collect_focal_pair_arrays(
    vcor_file: str, focal_ids: set[str]
) -> tuple[dict[str, list], dict[str, list]]:
    """Stream vcor chunks into per-focal distance / r2 arrays (full LD window)."""
    dist_map: dict[str, list] = defaultdict(list)
    r2_map: dict[str, list] = defaultdict(list)
    if not focal_ids:
        return dist_map, r2_map
    for chunk in _iter_plink2_vcor_chunks(vcor_file, chunksize=500_000):
        local = _normalize_ld_chunk(chunk, keep_ids=True)
        if local is None or local.empty:
            continue
        local = local.copy()
        local["ID_A"] = local["ID_A"].astype(str)
        local["ID_B"] = local["ID_B"].astype(str)
        dist = (
            np.abs(
                local["BP_B"].to_numpy(dtype=np.float64)
                - local["BP_A"].to_numpy(dtype=np.float64)
            )
            / 1000.0
        )
        local["_dist"] = dist
        local["_r2"] = local["R2"].to_numpy(dtype=np.float64)
        for id_col in ("ID_A", "ID_B"):
            mask = local[id_col].isin(focal_ids)
            if not mask.any():
                continue
            sub = local.loc[mask]
            for vid, grp in sub.groupby(id_col, sort=False):
                dist_map[str(vid)].append(grp["_dist"].to_numpy(dtype=np.float64))
                r2_map[str(vid)].append(grp["_r2"].to_numpy(dtype=np.float64))
    return dist_map, r2_map


def focal_mean_r2_decay_table(
    pairs: pd.DataFrame,
    focals: pd.DataFrame,
    *,
    r2_threshold: float = R2_DECAY_THRESHOLD,
    bin_size_kb: float = 50.0,
    min_pairs_per_bin: int = 5,
    max_distance_kb: float = 5000.0,
) -> pd.DataFrame:
    """Per-focal distance where mean r2 falls to ``r2_threshold``."""
    if pairs is None or pairs.empty or focals is None or focals.empty:
        return pd.DataFrame(
            columns=[
                "ID",
                "ref_chrom",
                "ref_pos",
                "ref_pos_mb",
                "decay_kb",
                "plot_decay_kb",
                "reached",
                "n_pairs",
            ]
        )
    id_a = pairs["ID_A"].astype(str).to_numpy()
    id_b = pairs["ID_B"].astype(str).to_numpy()
    dist_kb = (
        np.abs(
            pairs["BP_B"].to_numpy(dtype=np.float64)
            - pairs["BP_A"].to_numpy(dtype=np.float64)
        )
        / 1000.0
    )
    r2 = pairs["R2"].to_numpy(dtype=np.float64)
    dist_map: dict[str, np.ndarray] = {}
    r2_map: dict[str, np.ndarray] = {}
    for vid in focals["ID"].astype(str):
        mask = (id_a == vid) | (id_b == vid)
        dist_map[vid] = dist_kb[mask]
        r2_map[vid] = r2[mask]
    return _focal_decay_rows_from_maps(
        focals,
        dist_map,
        r2_map,
        r2_threshold=float(r2_threshold),
        bin_size_kb=float(bin_size_kb),
        min_pairs_per_bin=int(min_pairs_per_bin),
        max_distance_kb=float(max_distance_kb),
    )


def _focal_decay_rows_from_maps(
    focals: pd.DataFrame,
    dist_map,
    r2_map,
    *,
    r2_threshold: float,
    bin_size_kb: float,
    min_pairs_per_bin: int,
    max_distance_kb: float,
) -> pd.DataFrame:
    out_rows: list[dict] = []
    for row in focals.itertuples(index=False):
        vid = str(row.ID)
        dist_parts = dist_map.get(vid, [])
        r2_parts = r2_map.get(vid, [])
        if isinstance(dist_parts, np.ndarray):
            dist_arr = dist_parts
            r2_arr = r2_parts
        elif dist_parts:
            dist_arr = np.concatenate(dist_parts)
            r2_arr = np.concatenate(r2_parts)
        else:
            dist_arr = np.asarray([], dtype=np.float64)
            r2_arr = np.asarray([], dtype=np.float64)
        decay_kb, reached, n_pairs = mean_r2_decay_distance_kb(
            dist_arr,
            r2_arr,
            threshold=float(r2_threshold),
            bin_size_kb=float(bin_size_kb),
            min_pairs_per_bin=int(min_pairs_per_bin),
            max_distance_kb=float(max_distance_kb),
        )
        plot_decay = decay_kb if reached else float(max_distance_kb)
        row_out = {
            "ID": vid,
            "ref_chrom": str(row.ref_chrom),
            "ref_pos": int(row.ref_pos),
            "ref_pos_mb": float(row.ref_pos) / 1e6,
            "decay_kb": decay_kb,
            "plot_decay_kb": plot_decay,
            "reached": bool(reached),
            "n_pairs": int(n_pairs),
        }
        for col in (
            "bin_idx",
            "bin_start_bp",
            "bin_end_bp",
            "bin_mid_mb",
            "n_in_bin",
        ):
            if hasattr(row, col):
                row_out[col] = getattr(row, col)
        out_rows.append(row_out)
    return pd.DataFrame(out_rows)


def aggregate_decay_by_spatial_bin(
    focal_table: pd.DataFrame,
    *,
    max_distance_kb: float = 20000.0,
) -> pd.DataFrame:
    """Mean per-focal decay-to-r2 distance in each spatial bin (partners may lie outside the bin).

    Bins with no focal that reaches the r2 threshold keep ``mean_decay_kb`` as
    NaN and set ``plot_decay_kb`` to ``max_distance_kb`` (window cap).
    """
    needed = {"ref_chrom", "bin_idx", "bin_mid_mb", "decay_kb", "reached"}
    if focal_table is None or focal_table.empty or not needed.issubset(focal_table.columns):
        return pd.DataFrame()
    cap = float(max_distance_kb)
    group_cols = ["ref_chrom", "bin_idx"]
    rows: list[dict] = []
    for key, grp in focal_table.groupby(group_cols, sort=True):
        chrom, bin_idx = key
        decay = pd.to_numeric(grp["decay_kb"], errors="coerce").to_numpy(dtype=np.float64)
        reached = grp["reached"].to_numpy(dtype=bool)
        ok = reached & np.isfinite(decay)
        n_ok = int(ok.sum())
        mean_d = float(np.mean(decay[ok])) if n_ok else float("nan")
        median_d = float(np.median(decay[ok])) if n_ok else float("nan")
        first = grp.iloc[0]
        rows.append(
            {
                "ref_chrom": str(chrom),
                "bin_idx": int(bin_idx),
                "bin_start_bp": int(first["bin_start_bp"]) if "bin_start_bp" in grp.columns else int(bin_idx),
                "bin_end_bp": int(first["bin_end_bp"]) if "bin_end_bp" in grp.columns else int(bin_idx),
                "bin_mid_mb": float(first["bin_mid_mb"]),
                "ref_pos_mb": float(first["bin_mid_mb"]),
                "n_focals": int(len(grp)),
                "n_reached": n_ok,
                "n_in_bin": int(first["n_in_bin"]) if "n_in_bin" in grp.columns else int(len(grp)),
                "mean_decay_kb": mean_d,
                "median_decay_kb": median_d,
                "decay_kb": mean_d,
                "plot_decay_kb": mean_d if n_ok else cap,
                "reached": n_ok > 0,
            }
        )
    return pd.DataFrame(rows)


def summarize_decay_distance_by_chromosome(
    focal_table: pd.DataFrame,
    r2_threshold: float = R2_DECAY_THRESHOLD,
) -> pd.DataFrame:
    """Median/mean decay_kb among focals that reached r2_threshold, one row per chromosome."""
    needed = {"ref_chrom", "reached", "decay_kb"}
    if focal_table is None or focal_table.empty or not needed.issubset(focal_table.columns):
        return pd.DataFrame()
    rows: list[dict] = []
    work = focal_table.copy()
    work["ref_chrom"] = work["ref_chrom"].map(_ref_chrom_token)
    for chrom, grp in work.groupby("ref_chrom", sort=True):
        decay = pd.to_numeric(grp["decay_kb"], errors="coerce").to_numpy(dtype=np.float64)
        reached = grp["reached"].to_numpy(dtype=bool)
        ok = reached & np.isfinite(decay)
        n_ok = int(ok.sum())
        rows.append(
            {
                "chromosome": f"chr{chrom}" if not str(chrom).startswith("chr") else str(chrom),
                "r2_threshold": float(r2_threshold),
                "n_focals": int(len(grp)),
                "n_reached": n_ok,
                "median_decay_kb": float(np.median(decay[ok])) if n_ok else float("nan"),
                "mean_decay_kb": float(np.mean(decay[ok])) if n_ok else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def _centromere_x_spans(chrom: str) -> list[dict]:
    token = _ref_chrom_token(chrom)
    start = get_centromere_start(token)
    end = get_centromere_end(token)
    if start is None or end is None:
        return []
    return [
        {
            "start": float(start) / 1e6,
            "end": float(end) / 1e6,
            "color": "#2ca02c",
            "alpha": 0.18,
            "label": "Centromere",
        }
    ]


def decay_plot_ylim(
    decay_kb,
    *,
    pad: float = DECAY_PLOT_YLIM_PAD,
    min_span_kb: float = DECAY_PLOT_YLIM_MIN_KB,
) -> tuple[float, float]:
    """Y-limit from observed decay distances, not the LD window cap."""
    vals = np.asarray(decay_kb, dtype=np.float64)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return (0.0, float(min_span_kb))
    y_hi = float(np.max(vals)) * float(pad)
    return (0.0, max(y_hi, float(min_span_kb)))


def _prepare_decay_plot_frame(
    table: pd.DataFrame,
    chrom: str,
    max_distance_kb: float,
    *,
    r2_threshold: float = R2_DECAY_THRESHOLD,
) -> tuple[pd.DataFrame, tuple[float, float], str, str]:
    """Return plot frame, ylim, x column, and y-axis label for one chromosome."""
    sub = table[table["ref_chrom"].astype(str).map(_ref_chrom_token) == _ref_chrom_token(chrom)].copy()
    x_col = "bin_mid_mb" if "mean_decay_kb" in table.columns and "bin_mid_mb" in table.columns else "ref_pos_mb"
    r2_txt = _r2_label(r2_threshold)
    ylabel = (
        f"Mean distance to mean r2 = {r2_txt} (kb)"
        if "mean_decay_kb" in table.columns
        else f"Distance to mean r2 = {r2_txt} (kb)"
    )
    if sub.empty:
        plot_df = pd.DataFrame({x_col: [], "plot_decay_kb": [], "reached": []})
        return plot_df, decay_plot_ylim([]), x_col, ylabel
    plot_df = sub.copy()
    y_src = plot_df["plot_decay_kb"] if "plot_decay_kb" in plot_df.columns else plot_df["decay_kb"]
    y_num = pd.to_numeric(y_src, errors="coerce").to_numpy(dtype=np.float64)
    if "reached" in plot_df.columns:
        reached = plot_df["reached"].to_numpy(dtype=bool)
        y_num = np.where(reached & np.isfinite(y_num), y_num, float(max_distance_kb))
    else:
        y_num = np.where(np.isfinite(y_num), y_num, float(max_distance_kb))
    plot_df["plot_decay_kb"] = y_num
    ylim = decay_plot_ylim(plot_df["plot_decay_kb"])
    plot_df["plot_decay_kb"] = np.clip(
        plot_df["plot_decay_kb"].to_numpy(dtype=np.float64), ylim[0], ylim[1]
    )
    return plot_df, ylim, x_col, ylabel


def _plot_decay_along_chromosome(
    table: pd.DataFrame,
    chrom: str,
    output_prefix: str,
    max_distance_kb: float,
    *,
    ax=None,
    save: bool = True,
    title: str | None = None,
    xlabel: str | None = "Position (Mb)",
    ylabel: str | None = None,
    ylim: tuple[float, float] | None = None,
    show_legend: bool = True,
    legend_outside: bool = True,
    x_span_label: bool = True,
    r2_threshold: float = R2_DECAY_THRESHOLD,
):
    plot_df, data_ylim, x_col, default_ylabel = _prepare_decay_plot_frame(
        table, chrom, max_distance_kb, r2_threshold=r2_threshold
    )
    token = _ref_chrom_token(chrom)
    try:
        chr_len_mb = get_ref_chromosome_length(token) / 1e6
    except ValueError:
        chr_len_mb = float(plot_df[x_col].max()) if len(plot_df) else 1.0
    ylim = data_ylim if ylim is None else ylim
    if len(plot_df):
        plot_df = plot_df.copy()
        plot_df["plot_decay_kb"] = np.clip(
            plot_df["plot_decay_kb"].to_numpy(dtype=np.float64), ylim[0], ylim[1]
        )
    thresholds_h = None
    reached_col = "reached" in plot_df.columns
    if reached_col and len(plot_df) and not bool(plot_df["reached"].all()):
        thresholds_h = [
            {
                "value": float(max_distance_kb),
                "color": "0.4",
                "linestyle": ":",
                "label": f"Unreached ({max_distance_kb:g} kb window)",
            }
        ]
    x_spans = _centromere_x_spans(token)
    if x_spans and not x_span_label:
        x_spans = [{**span, "label": None} for span in x_spans]
    is_spatial_bin = "mean_decay_kb" in table.columns
    r2_txt = _r2_label(r2_threshold)
    if title is None:
        title = f"Distance to mean r2 = {r2_txt} along chr{token}"
        if is_spatial_bin and "bin_start_bp" in table.columns and "bin_end_bp" in table.columns:
            sub = table[table["ref_chrom"].astype(str).map(_ref_chrom_token) == token]
            widths = (
                pd.to_numeric(sub["bin_end_bp"], errors="coerce")
                - pd.to_numeric(sub["bin_start_bp"], errors="coerce")
                + 1
            ) / 1e6
            bin_mb = float(np.nanmedian(widths.to_numpy(dtype=np.float64))) if len(widths) else 5.0
            title = f"Mean distance to mean r2 = {r2_txt} along chr{token} ({bin_mb:g} Mb bins)"
    return plot_scatter_with_thresholds(
        data=plot_df,
        x_col=x_col,
        y_col="plot_decay_kb",
        title=title,
        filename=None if not save else f"{output_prefix}.png",
        xlabel=xlabel,
        ylabel=default_ylabel if ylabel is None else ylabel,
        color="steelblue",
        alpha=0.9,
        s=28,
        figure_size=(12, 4.5),
        xlim=(0.0, float(chr_len_mb)),
        ylim=ylim,
        x_spans=x_spans,
        connect=True,
        connect_flag_col="reached" if reached_col else None,
        connect_alt_color="crimson",
        connect_alt_linestyle="--",
        connect_flag_true_label=f"Reached r2={r2_txt}",
        connect_flag_false_label=f"Unreached ({max_distance_kb:g} kb)",
        thresholds_h=thresholds_h,
        ax=ax,
        save=save,
        show_legend=show_legend,
        legend_outside=legend_outside,
    )


def plot_decay_along_subgenome(
    table: pd.DataFrame,
    subgenome: str,
    output_prefix: str,
    max_distance_kb: float,
    *,
    r2_threshold: float = R2_DECAY_THRESHOLD,
):
    """Stack chr1–chr7 along-chromosome decay panels for one subgenome."""
    import matplotlib.pyplot as plt

    chroms = subgenome_ref_chroms(subgenome)
    y_vals = []
    frames = {}
    ylabel = f"Mean distance to mean r2 = {_r2_label(r2_threshold)} (kb)"
    for chrom in chroms:
        plot_df, _, _, ylabel = _prepare_decay_plot_frame(
            table, chrom, max_distance_kb, r2_threshold=r2_threshold
        )
        frames[chrom] = plot_df
        if len(plot_df) and "plot_decay_kb" in plot_df.columns:
            y_vals.append(plot_df["plot_decay_kb"].to_numpy(dtype=np.float64))
    ylim = decay_plot_ylim(np.concatenate(y_vals) if y_vals else [])
    n = len(chroms)
    fig, axes = plt.subplots(n, 1, figsize=(12, 2.15 * n), sharex=False, sharey=True)
    if n == 1:
        axes = [axes]
    for i, (chrom, ax) in enumerate(zip(chroms, axes)):
        last = i == n - 1
        _plot_decay_along_chromosome(
            table,
            chrom,
            output_prefix,
            max_distance_kb,
            ax=ax,
            save=False,
            title=f"chr{chrom}",
            xlabel="Position (Mb)" if last else "",
            ylabel=ylabel,
            ylim=ylim,
            show_legend=last,
            legend_outside=False,
            x_span_label=last,
            r2_threshold=r2_threshold,
        )
    fig.suptitle(
        f"Mean distance to mean r2 = {_r2_label(r2_threshold)} along "
        f"{subgenome_ref_chroms(subgenome)[0][-1]} genome "
        "(5 Mb bins)",
        fontsize=18,
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    filename = f"{output_prefix}.png"
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Saved scatter plot to {filename}")
    return fig


def collect_chrom_ld_files(work_dir: str, sample_id: str) -> list[tuple[str, str, str]]:
    """Pair ``{id}.{chrom}.focal.ld.vcor`` with matching focals TSV files."""
    root = Path(work_dir)
    suffix = ".focal.ld.vcor"
    prefix = f"{sample_id}."
    pairs: list[tuple[str, str, str]] = []
    for vcor in sorted(root.glob(f"{prefix}*{suffix}")):
        chrom = vcor.name[len(prefix) : -len(suffix)]
        focals = root / f"{sample_id}.{chrom}.focals.info.tsv"
        if focals.is_file():
            pairs.append((chrom, str(vcor), str(focals)))
    return pairs


def ana_ld_r2_decay_along_subgenome(
    chrom_ld_files: list[tuple[str, str, str]],
    output_prefix: str,
    *,
    subgenome: str,
    r2_threshold: float = R2_DECAY_THRESHOLD,
    max_distance_kb: float = 20000.0,
    bin_size_kb: float = 100.0,
    min_pairs_per_bin: int = 5,
) -> pd.DataFrame:
    """Compute per-chromosome along-chr decay tables and one stacked subgenome PNG."""
    plot_frames: list[pd.DataFrame] = []
    focal_frames: list[pd.DataFrame] = []
    for chrom, vcor, focals in chrom_ld_files:
        chrom_prefix = f"{output_prefix}.{chrom}"
        table = ana_ld_r2_decay_along_chr(
            vcor,
            focals,
            chrom_prefix,
            r2_threshold=r2_threshold,
            max_distance_kb=max_distance_kb,
            bin_size_kb=bin_size_kb,
            min_pairs_per_bin=min_pairs_per_bin,
            plot=False,
        )
        if table is not None and not table.empty:
            focal_frames.append(table)
        plot_table = table
        info_path = Path(f"{chrom_prefix}.info.tsv")
        if info_path.is_file():
            loaded = load_df_from_tsv(str(info_path))
            if loaded is not None and not loaded.empty:
                plot_table = loaded
        if plot_table is not None and not plot_table.empty:
            plot_frames.append(plot_table)
    combined = pd.concat(plot_frames, ignore_index=True) if plot_frames else pd.DataFrame()
    if focal_frames:
        summary = summarize_decay_distance_by_chromosome(
            pd.concat(focal_frames, ignore_index=True),
            r2_threshold=float(r2_threshold),
        )
        if not summary.empty:
            save_df_to_tsv(summary, f"{output_prefix}.chr_summary.info.tsv")
            print(f"[Info] Wrote chromosome decay summary ({len(summary)} chromosomes)")
    fig = plot_decay_along_subgenome(
        combined,
        subgenome,
        output_prefix,
        float(max_distance_kb),
        r2_threshold=float(r2_threshold),
    )
    plt_close = getattr(fig, "clf", None)
    if plt_close is not None:
        import matplotlib.pyplot as plt

        plt.close(fig)
    return combined


def ana_ld_r2_decay_along_chr(
    input_file: str,
    focals_file: str,
    output_prefix: str,
    *,
    r2_threshold: float = R2_DECAY_THRESHOLD,
    max_distance_kb: float = 20000.0,
    bin_size_kb: float = 100.0,
    min_pairs_per_bin: int = 5,
    plot: bool = True,
) -> pd.DataFrame:
    """Compute per-focal mean-r2 decay distance and plot it along the chromosome.

    When focals carry spatial-bin columns, the published plot and ``*.info.tsv``
    are the per-bin *mean* of those decay distances. Pairwise r2 still uses the
    full LD window (not the spatial bin).
    """
    focals = load_df_from_tsv(focals_file)
    if focals is None or focals.empty:
        print(f"[Error] Empty focals table: {focals_file}")
        return pd.DataFrame()
    focal_ids = set(focals["ID"].astype(str))
    dist_map, r2_map = _collect_focal_pair_arrays(input_file, focal_ids)
    table = _focal_decay_rows_from_maps(
        focals,
        dist_map,
        r2_map,
        r2_threshold=float(r2_threshold),
        bin_size_kb=float(bin_size_kb),
        min_pairs_per_bin=int(min_pairs_per_bin),
        max_distance_kb=float(max_distance_kb),
    )
    plot_table = table
    if not table.empty and "bin_mid_mb" in table.columns:
        plot_table = aggregate_decay_by_spatial_bin(
            table, max_distance_kb=float(max_distance_kb)
        )
        save_df_to_tsv(table, f"{output_prefix}.focal_decay.info.tsv")
        save_df_to_tsv(plot_table, f"{output_prefix}.info.tsv")
        n_reached = int(table["reached"].sum())
        n_bins = int(len(plot_table))
        n_bins_ok = int(plot_table["reached"].sum()) if n_bins else 0
        print(
            f"[Info] Aggregated {len(table)} focals into {n_bins} spatial bins "
            f"({n_bins_ok} with at least one reached decay distance)"
        )
    else:
        save_df_to_tsv(table, f"{output_prefix}.info.tsv")
        n_reached = int(table["reached"].sum()) if not table.empty else 0
        n_bins = 0
        n_bins_ok = 0
    save_thresholds(
        {
            "N_Focals": int(len(table)),
            "N_Reached_R2_Threshold": n_reached,
            "N_Spatial_Bins": n_bins,
            "N_Spatial_Bins_Reached": n_bins_ok,
            "R2_Threshold": float(r2_threshold),
            "Max_Distance_kb": float(max_distance_kb),
            "Bin_Size_kb": float(bin_size_kb),
            "Min_Pairs_Per_Bin": int(min_pairs_per_bin),
            "Median_Decay_kb": float(table["decay_kb"].median()) if n_reached else np.nan,
            "Median_Bin_Mean_Decay_kb": (
                float(plot_table["mean_decay_kb"].median())
                if n_bins_ok and "mean_decay_kb" in plot_table.columns
                else np.nan
            ),
        },
        f"{output_prefix}.th.tsv",
    )
    if plot:
        chroms = (
            list(plot_table["ref_chrom"].astype(str).unique()) if not plot_table.empty else []
        )
        if len(chroms) == 1:
            _plot_decay_along_chromosome(
                plot_table,
                chroms[0],
                output_prefix,
                float(max_distance_kb),
                r2_threshold=float(r2_threshold),
            )
        else:
            for chrom in chroms:
                _plot_decay_along_chromosome(
                    plot_table,
                    chrom,
                    f"{output_prefix}.{chrom}",
                    float(max_distance_kb),
                    r2_threshold=float(r2_threshold),
                )
    print(
        f"[Info] Along-chr LD r2 decay: {len(table)} focals, "
        f"{n_reached} reached r2={r2_threshold} within {max_distance_kb} kb"
    )
    return table
