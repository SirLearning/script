"""Three-ring Circos map: rare density, common density, and r2/D-prime LD blocks.

LD blocks are contiguous SNP runs whose adjacent unphased r2 meets ``r2_threshold``.
Adjacent blocks are then merged when the lower quartile of |D'| among junction
pairs exceeds ``dprime_q1_threshold`` (Wu-style Q1(D') > 0.98 merge).
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.variant.ld_along_chr import load_pvar_ref_variants
from infra.utils.io import load_df_from_tsv, save_df_to_tsv, save_thresholds
from infra.wheat.ref_v1 import get_chromosomes, get_ref_chromosome_length

logger = logging.getLogger(__name__)

MAF_RARE_THRESHOLD = 0.05
R2_BLOCK_THRESHOLD = 0.8
DPRIME_Q1_MERGE = 0.98
DEFAULT_BIN_BP = 1_000_000
DEFAULT_BLOCK_COVERAGE_FRAC = 0.50
DEFAULT_EDGE_SNPS = 10
DEFAULT_LD_CHUNKSIZE = 2_000_000
WHEAT_CIRCOS_CHR_ORDER = tuple(f"{i}{sg}" for i in range(1, 8) for sg in "ABD")
_LD_COL_MAP = {
    "SNP_A": "id_a",
    "ID_A": "id_a",
    "SNP_B": "id_b",
    "ID_B": "id_b",
    "R2": "r2",
    "UNPHASED_R2": "r2",
    "DP": "dprime",
    "DPRIME": "dprime",
    "D_PRIME": "dprime",
}
RARE_COLOR = "#4c72b0"
COMMON_COLOR = "#8172b3"
LD_BLOCK_COLOR = "#c44e52"
RARE_TRACK_LABEL = "Rare variants"
COMMON_TRACK_LABEL = "Common variants"
LD_BLOCK_TRACK_LABEL = "LD block"


def _norm_chrom(value) -> str:
    token = str(value).strip()
    if token.startswith("chr"):
        token = token[3:]
    return token


def wheat_circos_chrom_lengths() -> dict[str, int]:
    """IWGSC v1 lengths for the 21 A/B/D chromosomes, Circos order (1A, 1B, 1D, …)."""
    known = set(get_chromosomes())
    lengths = {}
    for chrom in WHEAT_CIRCOS_CHR_ORDER:
        if chrom not in known:
            continue
        lengths[chrom] = int(get_ref_chromosome_length(chrom))
    return lengths


def load_sites_from_pvar_afreq(pvar_path: str, afreq_path: str) -> pd.DataFrame:
    """Join pvar coordinates to MAF from a PLINK2 ``.afreq`` table."""
    variants = load_pvar_ref_variants(pvar_path)
    afreq = pd.read_csv(afreq_path, sep="\t")
    afreq.columns = [str(c).lstrip("#") for c in afreq.columns]
    if "ID" not in afreq.columns or "ALT_FREQS" not in afreq.columns:
        raise ValueError(f"{afreq_path} must contain ID and ALT_FREQS")
    af = pd.to_numeric(afreq["ALT_FREQS"], errors="coerce")
    maf = pd.Series(
        np.minimum(af.to_numpy(dtype=np.float64), 1.0 - af.to_numpy(dtype=np.float64)),
        index=afreq["ID"].astype(str),
    )
    out = variants.copy()
    out["id"] = out["ID"].astype(str)
    out["chrom"] = out["ref_chrom"].map(_norm_chrom)
    out["pos"] = pd.to_numeric(out["ref_pos"], errors="coerce")
    out["maf"] = out["id"].map(maf)
    return out.dropna(subset=["chrom", "pos", "maf"]).reset_index(drop=True)


def load_block_snps_from_pvar(pvar_path: str) -> pd.DataFrame:
    """Reference-chromosome coordinates from a pvar (no MAF join)."""
    variants = load_pvar_ref_variants(pvar_path)
    out = variants.copy()
    out["id"] = out["ID"].astype(str)
    out["chrom"] = out["ref_chrom"].map(_norm_chrom)
    out["pos"] = pd.to_numeric(out["ref_pos"], errors="coerce")
    return out.dropna(subset=["chrom", "pos", "id"]).reset_index(drop=True)


def _normalize_ld_frame(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()
    work.columns = [str(c).lstrip("#") for c in work.columns]
    rename = {col: _LD_COL_MAP[col.upper()] for col in work.columns if col.upper() in _LD_COL_MAP}
    work = work.rename(columns=rename)
    if "id_a" not in work.columns or "id_b" not in work.columns:
        raise ValueError(f"LD table missing SNP ids; got {list(work.columns)}")
    out = pd.DataFrame(
        {
            "id_a": work["id_a"].astype(str),
            "id_b": work["id_b"].astype(str),
            "r2": pd.to_numeric(work["r2"], errors="coerce") if "r2" in work.columns else np.nan,
            "dprime": pd.to_numeric(work["dprime"], errors="coerce")
            if "dprime" in work.columns
            else np.nan,
        }
    )
    if "r2" in work.columns:
        return out.dropna(subset=["r2"]).reset_index(drop=True)
    return out.reset_index(drop=True)


def _plink_ld_usecols(ld_path: str) -> list[str]:
    header = pd.read_csv(ld_path, sep=r"\s+", nrows=0)
    names = [str(c).lstrip("#") for c in header.columns]
    keep = set(_LD_COL_MAP)
    cols = [name for name in names if name.upper() in keep]
    if not cols:
        raise ValueError(f"{ld_path} missing SNP_A/SNP_B/R2 columns; got {names}")
    return cols


def iter_plink_ld_chunks(ld_path: str, chunksize: int = DEFAULT_LD_CHUNKSIZE):
    """Yield normalized PLINK ``--r2 dprime`` chunks (id_a/id_b/r2/dprime)."""
    usecols = _plink_ld_usecols(ld_path)
    reader = pd.read_csv(
        ld_path,
        sep=r"\s+",
        usecols=usecols,
        chunksize=max(int(chunksize), 1),
        low_memory=True,
    )
    for chunk in reader:
        yield _normalize_ld_frame(chunk)


def load_plink_ld_pairs(ld_path: str) -> pd.DataFrame:
    """Load PLINK 1.9 ``--r2 dprime`` pairs (SNP_A/SNP_B/R2/DP)."""
    frames = list(iter_plink_ld_chunks(ld_path, chunksize=DEFAULT_LD_CHUNKSIZE))
    if not frames:
        return pd.DataFrame(columns=["id_a", "id_b", "r2", "dprime"])
    return pd.concat(frames, ignore_index=True)


def attach_ld_values(
    query: pd.DataFrame,
    ld_source: pd.DataFrame | str,
    columns: Sequence[str] = ("r2", "dprime"),
    chunksize: int = DEFAULT_LD_CHUNKSIZE,
) -> pd.DataFrame:
    """Left-join LD values onto ``query`` rows keyed by ``id_a``/``id_b``.

    ``ld_source`` may be an in-memory pair table or a PLINK ``.ld`` path. The
    file path is scanned in chunks so a multi-GB windowed LD table is never
    fully materialized. Either pair orientation is accepted.
    """
    cols = tuple(columns)
    out = query.copy().reset_index(drop=True)
    n = len(out)
    filled = {col: np.full(n, np.nan, dtype=np.float64) for col in cols}
    found = np.zeros(n, dtype=bool)
    if n == 0:
        for col in cols:
            out[col] = filled[col]
        return out

    def ingest(pairs: pd.DataFrame) -> None:
        if found.all() or pairs is None or pairs.empty:
            return
        src_cols = ["id_a", "id_b"] + [c for c in cols if c in pairs.columns]
        if "id_a" not in pairs.columns or "id_b" not in pairs.columns:
            return
        src = pairs[src_cols].drop_duplicates(subset=["id_a", "id_b"], keep="first")
        pending = np.flatnonzero(~found)
        q = pd.DataFrame(
            {
                "id_a": out["id_a"].to_numpy()[pending],
                "id_b": out["id_b"].to_numpy()[pending],
                "_i": pending,
            }
        )
        for orient in (src, src.rename(columns={"id_a": "id_b", "id_b": "id_a"})):
            if q.empty:
                break
            hit = q.merge(orient, on=["id_a", "id_b"], how="inner")
            if hit.empty:
                continue
            ii = hit["_i"].to_numpy()
            for col in cols:
                if col not in hit.columns:
                    continue
                vals = pd.to_numeric(hit[col], errors="coerce").to_numpy(dtype=np.float64)
                dest = filled[col]
                unset = ~np.isfinite(dest[ii])
                dest[ii[unset]] = vals[unset]
            found[ii] = True
            hit_set = set(ii.tolist())
            q = q[~q["_i"].isin(hit_set)]

    if isinstance(ld_source, str):
        for chunk in iter_plink_ld_chunks(ld_source, chunksize=chunksize):
            ingest(chunk)
            if found.all():
                break
    else:
        ingest(_normalize_ld_frame(ld_source) if "SNP_A" in ld_source.columns else ld_source)

    for col in cols:
        out[col] = filled[col]
    return out


def bin_variant_density(
    sites: pd.DataFrame,
    bin_bp: int = DEFAULT_BIN_BP,
    maf_threshold: float = MAF_RARE_THRESHOLD,
) -> pd.DataFrame:
    """Count rare (MAF < threshold) and common sites in fixed bp bins."""
    work = sites.dropna(subset=["chrom", "pos", "maf"]).copy()
    if work.empty:
        return pd.DataFrame(
            columns=["chrom", "bin_start", "bin_end", "bin_mid", "n_rare", "n_common", "n_sites"]
        )
    bin_bp = int(bin_bp)
    work["bin_start"] = (work["pos"].astype(int) // bin_bp) * bin_bp
    work["is_rare"] = pd.to_numeric(work["maf"], errors="coerce") < float(maf_threshold)
    grouped = work.groupby(["chrom", "bin_start"], as_index=False).agg(
        n_rare=("is_rare", "sum"),
        n_sites=("maf", "size"),
    )
    grouped["n_rare"] = grouped["n_rare"].astype(int)
    grouped["n_common"] = grouped["n_sites"].astype(int) - grouped["n_rare"]
    grouped["bin_end"] = grouped["bin_start"] + bin_bp
    grouped["bin_mid"] = grouped["bin_start"] + bin_bp / 2.0
    return grouped.sort_values(["chrom", "bin_start"]).reset_index(drop=True)


def _ordered_snps(snps: pd.DataFrame) -> pd.DataFrame:
    ordered = snps.dropna(subset=["chrom", "id", "pos"]).copy()
    ordered["id"] = ordered["id"].astype(str)
    ordered["pos"] = pd.to_numeric(ordered["pos"], errors="coerce")
    return ordered.dropna(subset=["pos"]).sort_values(["chrom", "pos"]).reset_index(drop=True)


def _empty_blocks() -> pd.DataFrame:
    return pd.DataFrame(columns=["chrom", "start_id", "end_id", "start_pos", "end_pos", "n_snps"])


def _block_row(chrom: str, ids: Sequence[str], poss: Sequence[int], i0: int, i1: int) -> dict:
    return {
        "chrom": chrom,
        "start_id": ids[i0],
        "end_id": ids[i1],
        "start_pos": int(poss[i0]),
        "end_pos": int(poss[i1]),
        "n_snps": int(i1 - i0 + 1),
    }


def ld_blocks_from_pairs(
    snps: pd.DataFrame,
    pairs: pd.DataFrame | str,
    r2_threshold: float = R2_BLOCK_THRESHOLD,
    chunksize: int = DEFAULT_LD_CHUNKSIZE,
) -> pd.DataFrame:
    """Maximal runs of consecutive SNPs whose adjacent r2 meets ``r2_threshold``.

    Only consecutive-neighbour r2 is used. A PLINK ``.ld`` path is scanned in
    chunks and joined onto those edges; the full windowed pair table is not
    loaded into a Python dict.
    """
    ordered = _ordered_snps(snps)
    if ordered.empty:
        return _empty_blocks()
    nxt = ordered.groupby("chrom", sort=False)["id"].shift(-1)
    keep = nxt.notna()
    edges = pd.DataFrame(
        {
            "id_a": ordered.loc[keep, "id"].to_numpy(),
            "id_b": nxt.loc[keep].to_numpy(),
        }
    )
    if not edges.empty:
        edges = attach_ld_values(edges, pairs, columns=("r2",), chunksize=chunksize)
        r2_vals = pd.to_numeric(edges["r2"], errors="coerce").to_numpy(dtype=np.float64)
        r2_vals = np.where(np.isfinite(r2_vals), r2_vals, 0.0)
    else:
        r2_vals = np.array([], dtype=np.float64)
    rows = []
    edge_i = 0
    thr = float(r2_threshold)
    for chrom, grp in ordered.groupby("chrom", sort=False):
        ids = grp["id"].tolist()
        poss = grp["pos"].astype(int).tolist()
        if not ids:
            continue
        start_i = 0
        for i in range(len(ids) - 1):
            r2 = float(r2_vals[edge_i]) if edge_i < len(r2_vals) else 0.0
            edge_i += 1
            if r2 < thr:
                rows.append(_block_row(str(chrom), ids, poss, start_i, i))
                start_i = i + 1
        rows.append(_block_row(str(chrom), ids, poss, start_i, len(ids) - 1))
    if not rows:
        return _empty_blocks()
    return pd.DataFrame(rows)


def _snps_in_block(block: Mapping, snp_df: pd.DataFrame) -> list[str]:
    mask = (snp_df["pos"] >= int(block["start_pos"])) & (snp_df["pos"] <= int(block["end_pos"]))
    return snp_df.loc[mask, "id"].astype(str).tolist()


def _merge_block_row(left: Mapping, right: Mapping) -> dict:
    return {
        "chrom": left["chrom"],
        "start_id": left["start_id"],
        "end_id": right["end_id"],
        "start_pos": int(left["start_pos"]),
        "end_pos": int(right["end_pos"]),
        "n_snps": int(left["n_snps"]) + int(right["n_snps"]),
    }


def _junction_pair_query(
    blocks: pd.DataFrame,
    snp_df: pd.DataFrame,
    n_edge: int,
) -> tuple[pd.DataFrame, list[tuple[str, list[dict]]]]:
    """SNP pairs at each neighbouring-block junction, plus per-chromosome block lists."""
    edge = max(int(n_edge), 1)
    by_chrom = {str(chrom): grp for chrom, grp in snp_df.groupby("chrom", sort=False)}
    chrom_blocks: list[tuple[str, list[dict]]] = []
    rows: list[dict] = []
    jid = 0
    for chrom, grp in blocks.groupby("chrom", sort=False):
        ordered = grp.sort_values("start_pos").to_dict("records")
        chrom_blocks.append((str(chrom), ordered))
        snp_chr = by_chrom.get(str(chrom), pd.DataFrame(columns=snp_df.columns))
        for left, right in zip(ordered, ordered[1:]):
            left_ids = _snps_in_block(left, snp_chr)
            right_ids = _snps_in_block(right, snp_chr)
            if left_ids and right_ids:
                for a in left_ids[-edge:]:
                    for b in right_ids[:edge]:
                        if a != b:
                            rows.append({"id_a": a, "id_b": b, "jid": jid})
            jid += 1
    query = (
        pd.DataFrame(rows, columns=["id_a", "id_b", "jid"])
        if rows
        else pd.DataFrame(columns=["id_a", "id_b", "jid"])
    )
    return query, chrom_blocks


def merge_adjacent_blocks_by_dprime_q1(
    blocks: pd.DataFrame,
    snps: pd.DataFrame,
    pairs: pd.DataFrame | str,
    q1_threshold: float = DPRIME_Q1_MERGE,
    n_edge: int = DEFAULT_EDGE_SNPS,
    chunksize: int = DEFAULT_LD_CHUNKSIZE,
) -> pd.DataFrame:
    """Merge neighbouring r2 blocks when Q1(|D'|) at the junction exceeds the cutoff."""
    if blocks is None or blocks.empty:
        return blocks.copy() if blocks is not None else pd.DataFrame()
    snp_df = _ordered_snps(snps)
    query, chrom_blocks = _junction_pair_query(blocks, snp_df, n_edge)
    q1_by_jid: dict[int, float] = {}
    if not query.empty:
        attached = attach_ld_values(query, pairs, columns=("dprime",), chunksize=chunksize)
        attached["dprime"] = pd.to_numeric(attached["dprime"], errors="coerce")
        for jid, sub in attached.dropna(subset=["dprime"]).groupby("jid"):
            q1_by_jid[int(jid)] = float(
                np.nanpercentile(np.abs(sub["dprime"].to_numpy(dtype=float)), 25)
            )
    out_rows: list[dict] = []
    jid = 0
    for _chrom, ordered in chrom_blocks:
        if not ordered:
            continue
        current = dict(ordered[0])
        for nxt in ordered[1:]:
            q1 = q1_by_jid.get(jid)
            jid += 1
            if q1 is not None and q1 > float(q1_threshold):
                current = _merge_block_row(current, nxt)
            else:
                out_rows.append(current)
                current = dict(nxt)
        out_rows.append(current)
    return pd.DataFrame(out_rows)


def _merge_closed_intervals(starts: np.ndarray, ends: np.ndarray) -> list[tuple[int, int]]:
    """Union overlapping or abutting closed intervals ``[start, end]``."""
    if len(starts) == 0:
        return []
    order = np.argsort(starts, kind="mergesort")
    merged: list[list[int]] = []
    for s, e in zip(starts[order].tolist(), ends[order].tolist()):
        if not merged or s > merged[-1][1] + 1:
            merged.append([int(s), int(e)])
        else:
            merged[-1][1] = max(merged[-1][1], int(e))
    return [(a, b) for a, b in merged]


def ld_block_coverage_intervals(
    blocks: pd.DataFrame,
    chrom_lengths: Mapping[str, int],
    bin_bp: int = DEFAULT_BIN_BP,
    coverage_frac: float = DEFAULT_BLOCK_COVERAGE_FRAC,
) -> pd.DataFrame:
    """Merge ``bin_bp`` windows whose LD-block base coverage meets ``coverage_frac``.

    Zero-width (single-SNP) blocks are ignored. Overlapping blocks are unioned
    so covered bases in a window are not double-counted. Adjacent passing
    windows are concatenated into one interval.
    """
    empty = pd.DataFrame(columns=["chrom", "start_pos", "end_pos", "n_bins"])
    if blocks is None or blocks.empty:
        return empty
    bin_bp = int(bin_bp)
    if bin_bp <= 0:
        raise ValueError("bin_bp must be positive")
    thr = float(coverage_frac)
    work = blocks.dropna(subset=["chrom", "start_pos", "end_pos"]).copy()
    work["start_pos"] = pd.to_numeric(work["start_pos"], errors="coerce")
    work["end_pos"] = pd.to_numeric(work["end_pos"], errors="coerce")
    work = work.dropna(subset=["start_pos", "end_pos"])
    work = work.loc[work["end_pos"] > work["start_pos"]]
    if work.empty:
        return empty
    rows: list[dict] = []
    for chrom, length in chrom_lengths.items():
        chrom_len = int(length)
        if chrom_len <= 0:
            continue
        grp = work.loc[work["chrom"].astype(str) == str(chrom)]
        nbin = (chrom_len + bin_bp - 1) // bin_bp
        covered = np.zeros(nbin, dtype=np.int64)
        if not grp.empty:
            intervals = _merge_closed_intervals(
                grp["start_pos"].to_numpy(dtype=np.int64),
                grp["end_pos"].to_numpy(dtype=np.int64),
            )
            for start, end in intervals:
                hi = min(int(end) + 1, chrom_len)
                lo = max(int(start), 0)
                if hi <= lo:
                    continue
                b0 = lo // bin_bp
                b1 = (hi - 1) // bin_bp
                for b in range(max(b0, 0), min(b1 + 1, nbin)):
                    bin_lo = b * bin_bp
                    bin_hi = min(bin_lo + bin_bp, chrom_len)
                    ov_lo = max(lo, bin_lo)
                    ov_hi = min(hi, bin_hi)
                    if ov_hi > ov_lo:
                        covered[b] += ov_hi - ov_lo
        run_start = None
        for b in range(nbin):
            bin_lo = b * bin_bp
            bin_hi = min(bin_lo + bin_bp, chrom_len)
            bin_len = bin_hi - bin_lo
            frac = (covered[b] / bin_len) if bin_len > 0 else 0.0
            if frac >= thr:
                if run_start is None:
                    run_start = b
            elif run_start is not None:
                rows.append(
                    {
                        "chrom": str(chrom),
                        "start_pos": run_start * bin_bp,
                        "end_pos": min(b * bin_bp, chrom_len),
                        "n_bins": b - run_start,
                    }
                )
                run_start = None
        if run_start is not None:
            rows.append(
                {
                    "chrom": str(chrom),
                    "start_pos": run_start * bin_bp,
                    "end_pos": chrom_len,
                    "n_bins": nbin - run_start,
                }
            )
    if not rows:
        return empty
    return pd.DataFrame(rows)


def _histogram_track(density: pd.DataFrame, value_col: str) -> pd.DataFrame:
    if density.empty:
        return pd.DataFrame(columns=["chrom", "pos", "value", "bin_start", "bin_end"])
    out = density.copy()
    out["pos"] = out["bin_mid"]
    out["value"] = pd.to_numeric(out[value_col], errors="coerce").fillna(0.0)
    return out[["chrom", "pos", "value", "bin_start", "bin_end"]]


def plot_circos_ld_map_from_tables(
    density_table: str | Path | pd.DataFrame,
    blocks_table: str | Path | pd.DataFrame,
    output_prefix: str,
    chrom_lengths: Mapping[str, int] | None = None,
    bin_bp: int = DEFAULT_BIN_BP,
    coverage_frac: float = DEFAULT_BLOCK_COVERAGE_FRAC,
) -> None:
    """Draw the three-ring Circos PNG from published density and LD-block tables.

    The outer ring uses ``bin_bp`` windows whose unioned LD-block coverage is
    at least ``coverage_frac``; adjacent passing windows are merged.
    """
    from infra.utils.graph import plot_circos_chromosome_tracks

    if isinstance(density_table, (str, Path)):
        density = load_df_from_tsv(str(density_table))
    else:
        density = density_table.copy()
    if isinstance(blocks_table, (str, Path)):
        blocks = load_df_from_tsv(str(blocks_table))
    else:
        blocks = blocks_table.copy()

    lengths = dict(chrom_lengths) if chrom_lengths is not None else wheat_circos_chrom_lengths()
    chrom_order = [c for c in WHEAT_CIRCOS_CHR_ORDER if c in lengths] or list(lengths.keys())
    painted = ld_block_coverage_intervals(
        blocks,
        chrom_lengths=lengths,
        bin_bp=bin_bp,
        coverage_frac=coverage_frac,
    )
    save_df_to_tsv(painted, f"{output_prefix}.blocks_1mb.info.tsv")
    print(
        f"[circos] {len(painted)} LD intervals after {int(bin_bp)} bp "
        f"coverage>={coverage_frac:g}",
        flush=True,
    )
    plot_circos_chromosome_tracks(
        chrom_lengths=lengths,
        chrom_order=chrom_order,
        filename=f"{output_prefix}.png",
        title="Rare / common variant density and LD blocks",
        tracks=[
            {
                "kind": "histogram",
                "label": RARE_TRACK_LABEL,
                "data": _histogram_track(density, "n_rare"),
                "color": RARE_COLOR,
                "r_inner": 0.30,
                "r_outer": 0.50,
            },
            {
                "kind": "histogram",
                "label": COMMON_TRACK_LABEL,
                "data": _histogram_track(density, "n_common"),
                "color": COMMON_COLOR,
                "r_inner": 0.52,
                "r_outer": 0.70,
            },
            {
                "kind": "karyotype",
                "label": LD_BLOCK_TRACK_LABEL,
                "intervals": painted,
                "color": LD_BLOCK_COLOR,
                "interval_color": LD_BLOCK_COLOR,
                "base_color": "#1a1a1a",
                "r_inner": 0.74,
                "r_outer": 0.90,
            },
        ],
    )


def _process_subgenome_blocks(payload: tuple) -> pd.DataFrame:
    """Call LD blocks from one subgenome pvar plus its windowed ``.ld`` file."""
    pvar_path, ld_path, r2_threshold, dprime_q1_threshold, n_edge, chunksize = payload
    print(f"[circos] block SNPs {pvar_path}", flush=True)
    snps = load_block_snps_from_pvar(pvar_path)[["chrom", "id", "pos"]].drop_duplicates(
        subset=["chrom", "id"]
    )
    print(f"[circos] adjacent r2 scan {ld_path} ({len(snps)} sites)", flush=True)
    blocks = ld_blocks_from_pairs(
        snps, ld_path, r2_threshold=r2_threshold, chunksize=chunksize
    )
    print(f"[circos] {len(blocks)} r2 blocks; D' merge scan {ld_path}", flush=True)
    blocks = merge_adjacent_blocks_by_dprime_q1(
        blocks,
        snps,
        ld_path,
        q1_threshold=dprime_q1_threshold,
        n_edge=n_edge,
        chunksize=chunksize,
    )
    print(f"[circos] {len(blocks)} blocks after D' merge ({ld_path})", flush=True)
    return blocks


def ana_circos_ld_map(
    pvar_files: Sequence[str],
    afreq_files: Sequence[str],
    ld_files: Sequence[str],
    output_prefix: str,
    maf_threshold: float = MAF_RARE_THRESHOLD,
    r2_threshold: float = R2_BLOCK_THRESHOLD,
    dprime_q1_threshold: float = DPRIME_Q1_MERGE,
    bin_bp: int = DEFAULT_BIN_BP,
    chrom_lengths: Mapping[str, int] | None = None,
    n_jobs: int = 1,
    ld_chunksize: int = DEFAULT_LD_CHUNKSIZE,
    ld_pvar_files: Sequence[str] | None = None,
    coverage_frac: float = DEFAULT_BLOCK_COVERAGE_FRAC,
) -> None:
    """Write rare/common density, LD-block tables, and a three-ring Circos PNG.

    Density uses ``pvar_files`` + ``afreq_files``. LD blocks use ``ld_files``
    with SNP coordinates from ``ld_pvar_files`` when given, otherwise from
    ``pvar_files``. ``n_jobs`` runs A/B/D block scans in separate processes.
    After writing the density and block tables, the PNG is drawn by
    ``plot_circos_ld_map_from_tables``.
    """
    if len(pvar_files) != len(afreq_files):
        raise ValueError("pvar_files and afreq_files must be paired")
    block_pvars = list(ld_pvar_files) if ld_pvar_files is not None else list(pvar_files)
    if len(block_pvars) != len(ld_files):
        raise ValueError("ld_pvar_files (or pvar_files) and ld_files must be paired")
    print(f"[circos] density from {len(pvar_files)} pvar/afreq pairs", flush=True)
    site_frames = [
        load_sites_from_pvar_afreq(pvar, afreq)
        for pvar, afreq in zip(pvar_files, afreq_files)
    ]
    sites = pd.concat(site_frames, ignore_index=True) if site_frames else pd.DataFrame()
    payloads = [
        (
            pvar,
            ld_path,
            float(r2_threshold),
            float(dprime_q1_threshold),
            DEFAULT_EDGE_SNPS,
            int(ld_chunksize),
        )
        for pvar, ld_path in zip(block_pvars, ld_files)
    ]
    workers = max(1, min(int(n_jobs), len(payloads))) if payloads else 1
    print(f"[circos] {len(payloads)} LD subgenomes, n_jobs={workers}", flush=True)
    if workers == 1:
        block_frames = [_process_subgenome_blocks(item) for item in payloads]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            block_frames = list(pool.map(_process_subgenome_blocks, payloads))
    blocks = (
        pd.concat(block_frames, ignore_index=True)
        if block_frames
        else pd.DataFrame(columns=["chrom", "start_id", "end_id", "start_pos", "end_pos", "n_snps"])
    )

    density = bin_variant_density(sites, bin_bp=bin_bp, maf_threshold=maf_threshold)
    save_df_to_tsv(density, f"{output_prefix}.density.info.tsv")
    save_df_to_tsv(blocks, f"{output_prefix}.blocks.info.tsv")

    lengths = dict(chrom_lengths) if chrom_lengths is not None else wheat_circos_chrom_lengths()
    n_rare = int((pd.to_numeric(sites["maf"], errors="coerce") < float(maf_threshold)).sum())
    n_common = int(len(sites) - n_rare)
    stats = {
        "N_sites": int(len(sites)),
        "N_rare": n_rare,
        "N_common": n_common,
        "MAF_threshold": float(maf_threshold),
        "R2_block_threshold": float(r2_threshold),
        "Dprime_Q1_merge": float(dprime_q1_threshold),
        "N_ld_blocks": int(len(blocks)),
        "Bin_bp": int(bin_bp),
    }
    save_thresholds(stats, f"{output_prefix}.th.tsv")
    print(stats, flush=True)
    plot_circos_ld_map_from_tables(
        density,
        blocks,
        output_prefix,
        chrom_lengths=lengths,
        bin_bp=bin_bp,
        coverage_frac=coverage_frac,
    )
    logger.info("Wrote Circos LD map under %s", output_prefix)
