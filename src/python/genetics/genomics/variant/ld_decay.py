from collections import deque
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path
from typing import Dict, Optional, Sequence
import re
import subprocess

import numpy as np
import pandas as pd

from infra.utils.graph import (
    SEQUENTIAL_WHITE_RED,
    TWO_CATEGORY_BLUE,
    TWO_CATEGORY_RED,
    plot_distribution_with_stats,
    plot_faceted_multi_line_series,
    plot_faceted_stacked_bars,
    plot_grouped_density_lines,
    plot_multi_line_series,
    plot_parallel_coordinates,
    plot_scatter_continuous_color,
    plot_scatter_two_categories,
    plot_scatter_with_thresholds,
    plot_scatter_facet_grid,
    plot_stacked_histogram,
)
from infra.utils.io import load_df_generic, save_df_to_tsv, save_thresholds
from infra.wheat.ref_v1 import get_chromosome

PAIR_RARE_RARE = "rare-rare"
PAIR_RARE_COMMON = "rare-common"
PAIR_COMMON_COMMON = "common-common"
PAIR_CLASSES = (PAIR_RARE_RARE, PAIR_RARE_COMMON, PAIR_COMMON_COMMON)
PAIR_TYPES = PAIR_CLASSES
PAIR_CLASS_COLORS = {
    PAIR_RARE_RARE: TWO_CATEGORY_BLUE,
    PAIR_RARE_COMMON: TWO_CATEGORY_RED,
    PAIR_COMMON_COMMON: "#8c8c8c",
}
# Decay lines only: brighter scarlet so translucent rare-common still contrasts with blue.
PAIR_CLASS_DECAY_COLORS = {
    PAIR_RARE_RARE: TWO_CATEGORY_BLUE,
    PAIR_RARE_COMMON: "#ff4d4d",
    PAIR_COMMON_COMMON: "#8c8c8c",
}
PAIR_CLASS_LABELS = {
    PAIR_RARE_RARE: "Rare-rare",
    PAIR_RARE_COMMON: "Rare-common",
    PAIR_COMMON_COMMON: "Common-common",
}
# Only rare-common is translucent so it does not hide the opaque rare-rare line.
PAIR_CLASS_LINE_ALPHA = {
    PAIR_RARE_RARE: 1.0,
    PAIR_RARE_COMMON: 0.6,
    PAIR_COMMON_COMMON: 0.9,
}
PAIR_CLASS_LINE_ZORDER = {
    PAIR_COMMON_COMMON: 1,
    PAIR_RARE_RARE: 2,
    PAIR_RARE_COMMON: 3,
}

# Within-bin MAF classes for LD decay (both sites in the same interval).
# Edges are left-open, right-closed: (0, 0.01], (0.01, 0.05], …
MAF_LD_BIN_EDGES = (0.0, 0.01, 0.05, 0.10, 0.20, 0.50)
MAF_LD_BIN_LABELS = (
    "(0, 0.01]",
    "(0.01, 0.05]",
    "(0.05, 0.1]",
    "(0.1, 0.2]",
    "(0.2, 0.5]",
)
MAF_LD_BIN_COLORS = {
    "(0, 0.01]": "#08306b",
    "(0.01, 0.05]": TWO_CATEGORY_BLUE,
    "(0.05, 0.1]": "#8172b3",
    "(0.1, 0.2]": TWO_CATEGORY_RED,
    "(0.2, 0.5]": "#dd8452",
}
DPRIME_CLASS_POS1 = "Dprime_pos1"
DPRIME_CLASS_NEG1 = "Dprime_neg1"
DPRIME_CLASS_POSITIVE = "Dprime_positive"
DPRIME_CLASS_NEGATIVE = "Dprime_negative"
DPRIME_CLASS_ZERO = "Dprime_zero"
DPRIME_CLASS_EPS = 1e-6
SIGNED_DPRIME_PROP_COLORS = {
    "p_plus1": TWO_CATEGORY_BLUE,
    "p_minus1": TWO_CATEGORY_RED,
    "p_other": "#7f7f7f",
}
SIGNED_DPRIME_BOUNDARY_COLORS = {
    DPRIME_CLASS_POS1: TWO_CATEGORY_BLUE,
    DPRIME_CLASS_NEG1: TWO_CATEGORY_RED,
}
MAF_SCATTER_HIGH_R2 = 0.2
MAF_SCATTER_ULTRARARE_MAX = 0.01
R2_MAX_CONTOUR_LEVELS = (0.8, 0.6, 0.2)
R2_MAX_CONTOUR_STYLES = (
    {"color": "#111111", "linestyle": "-", "linewidth": 1.5},
    {"color": "#111111", "linestyle": "--", "linewidth": 1.5},
    {"color": "#111111", "linestyle": ":", "linewidth": 1.8},
)
ALL_BIN_LABEL = "all"
# Homoeolog 1–7 (Okabe–Ito minus black; black is reserved for the overall mean).
LD_DECAY_HOMOEOLOG_COLORS = (
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#CC79A7",
    "#56B4E9",
    "#D55E00",
    "#882255",
)
LD_DECAY_OVERALL_COLOR = "#111111"
_CHR_HOMOEOLOG_RE = re.compile(r"^chr([1-7])([ABD])$")


def ld_decay_chrom_label(chr_id) -> str:
    """Map a PLINK/vmap4 chrom token to a plot label (``chr1A``, ``Mit``, …)."""
    if chr_id is None:
        return "unknown"
    if isinstance(chr_id, float) and not np.isfinite(chr_id):
        return "unknown"
    if isinstance(chr_id, str):
        text = chr_id.strip()
        if not text:
            return "unknown"
        try:
            cid = int(text)
        except ValueError:
            if text.startswith("chr"):
                return text
            if len(text) >= 2 and text[0].isdigit():
                return f"chr{text}"
            return text
    else:
        try:
            cid = int(chr_id)
        except (TypeError, ValueError):
            return str(chr_id)
    name = get_chromosome(cid)
    if not name:
        return f"chr{cid}"
    if name in ("Mit", "Chl"):
        return name
    if name.startswith("chr"):
        return name
    return f"chr{name}"


def ld_decay_chrom_labels(chr_values) -> np.ndarray:
    """Vectorized ``ld_decay_chrom_label`` via unique-value lookup."""
    ser = pd.Series(chr_values, dtype=object)
    mapping = {val: ld_decay_chrom_label(val) for val in ser.unique()}
    return ser.map(mapping).to_numpy()


def _ld_decay_chrom_sort_key(label: str):
    match = _CHR_HOMOEOLOG_RE.match(str(label))
    if match:
        return (0, int(match.group(1)), str(label))
    return (1, 99, str(label))


def _ld_decay_chrom_color(label: str, index: int) -> str:
    match = _CHR_HOMOEOLOG_RE.match(str(label))
    if match:
        return LD_DECAY_HOMOEOLOG_COLORS[int(match.group(1)) - 1]
    return LD_DECAY_HOMOEOLOG_COLORS[int(index) % len(LD_DECAY_HOMOEOLOG_COLORS)]


def parse_wheat_homoeolog_label(label: str) -> Optional[tuple]:
    """Return ``(homoeolog_number, subgenome)`` for labels like ``chr4A``."""
    match = _CHR_HOMOEOLOG_RE.match(str(label).strip())
    if not match:
        return None
    return int(match.group(1)), match.group(2)


def summarize_ld_decay_chromosomes(by_chr: pd.DataFrame) -> pd.DataFrame:
    """One row per wheat homoeolog chromosome with pair-weighted mean r2."""
    columns = ["chromosome", "homoeolog", "subgenome", "mean_r2", "half_decay_kb", "n_pairs"]
    if by_chr is None or by_chr.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    for chrom, sub in by_chr.groupby("chromosome", sort=False):
        parsed = parse_wheat_homoeolog_label(str(chrom))
        if parsed is None:
            continue
        homoeolog, subgenome = parsed
        sub = sub.sort_values("bin_mid_kb")
        weights = sub["n_pairs"].to_numpy(dtype=np.float64)
        r2 = sub["mean_r2"].to_numpy(dtype=np.float64)
        finite = np.isfinite(r2) & np.isfinite(weights) & (weights > 0)
        wmean = (
            float(np.average(r2[finite], weights=weights[finite])) if np.any(finite) else np.nan
        )
        _, _, half_kb = _half_decay_kb(sub)
        rows.append(
            {
                "chromosome": str(chrom),
                "homoeolog": homoeolog,
                "subgenome": subgenome,
                "mean_r2": wmean,
                "half_decay_kb": half_kb,
                "n_pairs": int(np.nansum(weights)),
            }
        )
    out = pd.DataFrame(rows, columns=columns)
    if out.empty:
        return out
    return out.sort_values(["homoeolog", "subgenome"]).reset_index(drop=True)


def homoeolog_subgenome_ld_table(summaries: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Wide A/B/D table of pair-weighted mean r2 and within-subgenome rank."""
    frames = []
    for sg, summary in summaries.items():
        if summary is None or summary.empty:
            continue
        work = summary.copy()
        work["subgenome"] = str(sg)
        frames.append(work)
    empty_cols = ["homoeolog"] + [
        f"{prefix}_{sg}"
        for sg in ("A", "B", "D")
        for prefix in ("mean_r2", "rank", "n_pairs")
    ]
    if not frames:
        return pd.DataFrame(columns=empty_cols)
    long = pd.concat(frames, ignore_index=True)
    long["homoeolog"] = long["homoeolog"].astype(int)
    mean_wide = long.pivot_table(
        index="homoeolog", columns="subgenome", values="mean_r2", aggfunc="first"
    )
    n_wide = long.pivot_table(
        index="homoeolog", columns="subgenome", values="n_pairs", aggfunc="first"
    )
    for sg in ("A", "B", "D"):
        if sg not in mean_wide.columns:
            mean_wide[sg] = np.nan
        if sg not in n_wide.columns:
            n_wide[sg] = np.nan
    mean_wide = mean_wide[["A", "B", "D"]]
    n_wide = n_wide[["A", "B", "D"]]
    rank_wide = mean_wide.rank(axis=0, ascending=False, method="min")
    out = pd.DataFrame({"homoeolog": mean_wide.index.astype(int)})
    for sg in ("A", "B", "D"):
        out[f"mean_r2_{sg}"] = mean_wide[sg].to_numpy()
        out[f"rank_{sg}"] = rank_wide[sg].to_numpy()
        out[f"n_pairs_{sg}"] = n_wide[sg].to_numpy()
    return out.sort_values("homoeolog").reset_index(drop=True)


def plot_ld_homoeolog_subgenome_compare(
    a_by_chr: str,
    b_by_chr: str,
    d_by_chr: str,
    output_prefix: str,
) -> None:
    """Compare homoeolog LD size across A/B/D as parallel coordinates plus ranks."""
    summaries = {
        "A": summarize_ld_decay_chromosomes(pd.read_csv(a_by_chr, sep="\t")),
        "B": summarize_ld_decay_chromosomes(pd.read_csv(b_by_chr, sep="\t")),
        "D": summarize_ld_decay_chromosomes(pd.read_csv(d_by_chr, sep="\t")),
    }
    table = homoeolog_subgenome_ld_table(summaries)
    save_df_to_tsv(table, f"{output_prefix}.info.tsv")
    palette = {i + 1: LD_DECAY_HOMOEOLOG_COLORS[i] for i in range(len(LD_DECAY_HOMOEOLOG_COLORS))}
    plot_parallel_coordinates(
        data=table,
        value_cols=["mean_r2_A", "mean_r2_B", "mean_r2_D"],
        series_col="homoeolog",
        palette_map=palette,
        x_labels=["A", "B", "D"],
        title="Homoeolog LD across subgenomes (pair-weighted mean r2)",
        filename=f"{output_prefix}.png",
        ylabel="Pair-weighted mean r2 (0-5 Mb)",
        annotate=True,
        annotate_template="{series}{axis}",
    )
    rank_vals = table[["rank_A", "rank_B", "rank_D"]].to_numpy(dtype=np.float64)
    n_rank = int(np.nanmax(rank_vals)) if rank_vals.size and np.isfinite(rank_vals).any() else 7
    plot_parallel_coordinates(
        data=table,
        value_cols=["rank_A", "rank_B", "rank_D"],
        series_col="homoeolog",
        palette_map=palette,
        x_labels=["A", "B", "D"],
        title="Homoeolog LD rank across subgenomes (1 = highest mean r2)",
        filename=f"{output_prefix}.rank.png",
        ylabel="Rank (1 = highest LD)",
        invert_yaxis=True,
        yticks=list(range(1, n_rank + 1)),
        annotate=True,
        annotate_template="{series}{axis}",
    )


COMMON_Q_BIN_EDGES = (0.05, 0.10, 0.20, 0.30, 0.40, 0.5000001)
RARE_Q_BIN_EDGES = (0.0, 0.005, 0.01, 0.02, 0.03, 0.0500001)
RARE_MAC_BIN_EDGES = (0.0, 2.0, 5.0, 20.0, float("inf"))
RARE_MAC_BIN_LABELS = ("MAC 1-2", "MAC 3-5", "MAC 6-20", "MAC >20")
DISTANCE_KB_BIN_EDGES = (0.0, 100.0, 1000.0, 5000.0)
DISTANCE_KB_BIN_LABELS = ("0-100 kb", "100-1000 kb", "1000-5000 kb")


def _normalize_ld_chunk(
    chunk: pd.DataFrame,
    keep_ids: bool = False,
    keep_chrom: bool = False,
) -> Optional[pd.DataFrame]:
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
    if keep_ids:
        required_cols = required_cols | {"ID_A", "ID_B"}
    if keep_chrom:
        required_cols = required_cols | {"CHR_A"}
    if not required_cols.issubset(chunk.columns):
        return None

    cols = ["BP_A", "BP_B", "R2"]
    if keep_ids:
        cols = ["ID_A", "ID_B"] + cols
    if keep_chrom:
        cols = ["CHR_A"] + cols
    local = chunk[cols]
    already_numeric = (
        pd.api.types.is_numeric_dtype(local["R2"])
        and pd.api.types.is_numeric_dtype(local["BP_A"])
        and pd.api.types.is_numeric_dtype(local["BP_B"])
    )
    if not already_numeric:
        local = local.copy()
        local["R2"] = pd.to_numeric(local["R2"], errors="coerce")
        local["BP_A"] = pd.to_numeric(local["BP_A"], errors="coerce")
        local["BP_B"] = pd.to_numeric(local["BP_B"], errors="coerce")
    drop_cols = ["R2", "BP_A", "BP_B"]
    local = local.dropna(subset=drop_cols)
    if local.empty:
        return None
    return local


_VCOR_NUMERIC_DTYPES = {
    "POS_A": "int32",
    "POS_B": "int32",
    "UNPHASED_R2": "float32",
}


def _ids_to_maf(ids: pd.Series, maf_by_id: pd.Series) -> np.ndarray:
    """Map variant IDs to MAF via unique-ID lookup (avoids per-row Series.map)."""
    return _ids_to_values(ids, maf_by_id)


def _ids_to_values(ids: pd.Series, values: pd.Series) -> np.ndarray:
    """Map variant IDs to a numeric series via unique-ID lookup."""
    codes, uniques = pd.factorize(ids, sort=False)
    lookup = values.reindex(uniques).to_numpy(dtype=np.float64)
    return lookup[codes]


def r2_max_from_maf(maf_a: np.ndarray, maf_b: np.ndarray) -> np.ndarray:
    """Theoretical max haplotype r2 given two MAFs (p <= q <= 0.5): p(1-q)/(q(1-p))."""
    p = np.minimum(maf_a, maf_b).astype(np.float64, copy=False)
    q = np.maximum(maf_a, maf_b).astype(np.float64, copy=False)
    out = np.full(p.shape, np.nan, dtype=np.float64)
    denom = q * (1.0 - p)
    ok = np.isfinite(p) & np.isfinite(q) & (p > 0.0) & (q > 0.0) & (denom > 0.0)
    out[ok] = (p[ok] * (1.0 - q[ok])) / denom[ok]
    equal = np.isfinite(p) & np.isfinite(q) & (p == q)
    out[equal] = 1.0
    return out


def abs_dprime_from_maf_r2(
    maf_a: np.ndarray, maf_b: np.ndarray, r2: np.ndarray
) -> np.ndarray:
    """|D'| from MAFs and r2 **assuming D >= 0** (minor-allele coupling).

    Equals sqrt(r2 / r_max^2) with r_max^2 = p(1-q)/(q(1-p)). This is **not**
    general |D'|: r2 drops the sign of D, and D < 0 uses denominator pq.
    Prefer ``pairwise_complete_ld_stats`` when dosages are available.
    """
    r2_arr = np.asarray(r2, dtype=np.float64)
    rmax = r2_max_from_maf(
        np.asarray(maf_a, dtype=np.float64),
        np.asarray(maf_b, dtype=np.float64),
    )
    out = np.full(r2_arr.shape, np.nan, dtype=np.float64)
    ok = np.isfinite(r2_arr) & np.isfinite(rmax) & (rmax > 0.0) & (r2_arr >= 0.0)
    out[ok] = np.sqrt(np.clip(r2_arr[ok] / rmax[ok], 0.0, 1.0))
    return out


def pair_distance_kb_from_ids(id_a, id_b) -> np.ndarray:
    """Physical distance in kb from ``{chrom}-{pos}`` IDs; NaN if chroms differ."""
    a = pd.Series(np.asarray(id_a), dtype=str)
    b = pd.Series(np.asarray(id_b), dtype=str)
    n = len(a)
    out = np.full(n, np.nan, dtype=np.float64)
    pa = a.str.rsplit("-", n=1, expand=True)
    pb = b.str.rsplit("-", n=1, expand=True)
    if pa.shape[1] < 2 or pb.shape[1] < 2:
        return out
    chrom_a = pa[0].to_numpy()
    chrom_b = pb[0].to_numpy()
    pos_a = pd.to_numeric(pa[1], errors="coerce").to_numpy(dtype=np.float64)
    pos_b = pd.to_numeric(pb[1], errors="coerce").to_numpy(dtype=np.float64)
    same = (chrom_a == chrom_b) & (chrom_a != "") & np.isfinite(pos_a) & np.isfinite(pos_b)
    out[same] = np.abs(pos_a[same] - pos_b[same]) / 1000.0
    return out


def annotate_pairs_abs_dprime(pairs: pd.DataFrame) -> pd.DataFrame:
    """Add ``DISTANCE_KB`` from ``{chrom}-{pos}`` IDs. Does not infer |D'| from r2."""
    out = pairs.copy()
    if {"ID_A", "ID_B"} <= set(out.columns):
        out["DISTANCE_KB"] = pair_distance_kb_from_ids(out["ID_A"], out["ID_B"])
    return out


def r2_max_contour_branches(r2: float, n: int = 400) -> tuple[np.ndarray, np.ndarray]:
    """MAF pairs (p, q) with p <= q <= 0.5 and r_max^2(p, q) = r2.

    Invert p(1-q)/(q(1-p)) = r2 to q = p / (p + r2*(1-p)).
    p runs from near 0 to r2/(1+r2), where q reaches 0.5.
    """
    r2 = float(r2)
    if r2 <= 0.0 or r2 >= 1.0:
        raise ValueError(f"r2 contour level must be in (0, 1), got {r2}")
    n = max(2, int(n))
    p_hi = r2 / (1.0 + r2)
    p = np.linspace(p_hi / float(n), p_hi, n, dtype=np.float64)
    q = p / (p + r2 * (1.0 - p))
    return p, q


def r2_max_scatter_line_overlays(
    levels: Sequence[float] = R2_MAX_CONTOUR_LEVELS,
) -> list[dict]:
    """Six MAF-vs-MAF curves: each r_max^2 level mirrored across the diagonal."""
    overlays: list[dict] = []
    for i, r2 in enumerate(levels):
        style = R2_MAX_CONTOUR_STYLES[i % len(R2_MAX_CONTOUR_STYLES)]
        p, q = r2_max_contour_branches(float(r2))
        overlays.append(
            {
                **style,
                "x": p,
                "y": q,
                "label": f"$r_{{max}}^2$ = {float(r2):g}",
            }
        )
        overlays.append({**style, "x": q, "y": p, "label": None})
    return overlays


def pairwise_complete_ld_stats(
    dosage_a: np.ndarray, dosage_b: np.ndarray
) -> tuple[float, float, int, float, float]:
    """Pairwise-complete MAF, signed D, and signed D' (dosage/2 = allele fraction).

    Minor-allele coding so p, q <= 0.5, then D = x - pq with x = mean(A*B).

        D' = D / min[p(1-q), (1-p)q]     if D > 0
        D' = D / min[pq, (1-p)(1-q)]     if D < 0
        D' = 0                            if D = 0

    Clipped to [-1, 1]. Fifth return value is signed D', not |D'|.
    """
    a = np.asarray(dosage_a, dtype=np.float64)
    b = np.asarray(dosage_b, dtype=np.float64)
    both = np.isfinite(a) & np.isfinite(b)
    n = int(np.sum(both))
    nan5 = (float("nan"), float("nan"), n, float("nan"), float("nan"))
    if n == 0:
        return nan5
    fa = a[both] / 2.0
    fb = b[both] / 2.0
    af_a = float(np.mean(fa))
    af_b = float(np.mean(fb))
    if af_a > 0.5:
        fa = 1.0 - fa
        af_a = 1.0 - af_a
    if af_b > 0.5:
        fb = 1.0 - fb
        af_b = 1.0 - af_b
    p = af_a
    q = af_b
    maf_a = min(p, 1.0 - p)
    maf_b = min(q, 1.0 - q)
    x = float(np.mean(fa * fb))
    d = x - p * q
    if d > 0.0:
        d_den = min(p * (1.0 - q), (1.0 - p) * q)
    elif d < 0.0:
        d_den = min(p * q, (1.0 - p) * (1.0 - q))
    else:
        return maf_a, maf_b, n, d, 0.0
    if not np.isfinite(d_den) or d_den <= 0.0:
        return maf_a, maf_b, n, d, float("nan")
    signed = float(np.clip(d / d_den, -1.0, 1.0))
    return maf_a, maf_b, n, d, signed


def _minor_allele_count(maf, *, n_pairwise=None, obs_ct=None) -> np.ndarray:
    """Minor-allele count: pairwise ``MAF * 2 * N`` when sample N is known."""
    maf = np.asarray(maf, dtype=np.float64)
    out = np.full(np.shape(maf), np.nan, dtype=np.float64)
    if n_pairwise is not None:
        n = np.asarray(n_pairwise, dtype=np.float64)
        if n.shape != maf.shape:
            n = np.broadcast_to(n, maf.shape)
        ok = np.isfinite(maf) & np.isfinite(n) & (n > 0)
        out[ok] = maf[ok] * 2.0 * n[ok]
        return out
    if obs_ct is not None:
        obs = np.asarray(obs_ct, dtype=np.float64)
        ok = np.isfinite(maf) & np.isfinite(obs) & (obs > 0)
        out[ok] = maf[ok] * obs[ok]
    return out


def pair_haplotype_background_roles(
    maf_a,
    maf_b,
    maf_threshold: float = 0.05,
    *,
    n_pairwise=None,
    obs_ct_a=None,
    obs_ct_b=None,
) -> pd.DataFrame:
    """Background MAF ``q`` and rare MAC for the haplotype-background D' model.

    Rare-common: ``q`` is the common-site MAF; rare MAC is from the rare site.
    Rare-rare: ``q`` is the higher MAF (the compared / background rare);
    rare MAC is from the lower-MAF rare.
    """
    maf_a = np.asarray(maf_a, dtype=np.float64)
    maf_b = np.asarray(maf_b, dtype=np.float64)
    n = len(maf_a)
    masks = pair_type_masks(maf_a, maf_b, float(maf_threshold))
    pair_class = np.full(n, "", dtype=object)
    background_maf = np.full(n, np.nan, dtype=np.float64)
    rare_maf = np.full(n, np.nan, dtype=np.float64)
    rare_is_a = np.zeros(n, dtype=bool)
    rr = masks[PAIR_RARE_RARE]
    rc = masks[PAIR_RARE_COMMON]
    pair_class[rr] = PAIR_RARE_RARE
    pair_class[rc] = PAIR_RARE_COMMON
    a_ge_b = maf_a >= maf_b
    background_maf[rr] = np.where(a_ge_b[rr], maf_a[rr], maf_b[rr])
    rare_maf[rr] = np.where(a_ge_b[rr], maf_b[rr], maf_a[rr])
    rare_is_a[rr] = ~a_ge_b[rr]
    rare_a = maf_a < float(maf_threshold)
    background_maf[rc] = np.where(rare_a[rc], maf_b[rc], maf_a[rc])
    rare_maf[rc] = np.where(rare_a[rc], maf_a[rc], maf_b[rc])
    rare_is_a[rc] = rare_a[rc]
    mac_a = _minor_allele_count(maf_a, n_pairwise=n_pairwise, obs_ct=obs_ct_a)
    mac_b = _minor_allele_count(maf_b, n_pairwise=n_pairwise, obs_ct=obs_ct_b)
    rare_mac = np.where(rare_is_a, mac_a, mac_b)
    return pd.DataFrame(
        {
            "pair_class": pair_class,
            "background_maf": background_maf,
            "rare_maf": rare_maf,
            "rare_mac": rare_mac,
        }
    )


def _cut_with_labels(values, edges, labels, *, right: bool) -> pd.Series:
    return pd.cut(
        np.asarray(values, dtype=np.float64),
        bins=list(edges),
        labels=list(labels),
        right=right,
        include_lowest=True,
    )


def _q_bin_series(values, edges) -> pd.Series:
    return pd.cut(
        np.asarray(values, dtype=np.float64),
        bins=list(edges),
        right=False,
        include_lowest=True,
    )


def _aggregate_pplus_bins(work: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "pair_class",
        "mac_bin",
        "dist_bin",
        "q_bin",
        "q_mean",
        "q_mid",
        "n_pos1",
        "n_neg1",
        "n_bound",
        "p_plus",
        "expected_p_plus",
    ]
    if work.empty:
        return pd.DataFrame(columns=cols)
    grouped = work.groupby(
        ["pair_class", "mac_bin", "dist_bin", "q_bin"],
        observed=True,
        dropna=True,
    )
    rows: list[dict] = []
    for (pair_class, mac_bin, dist_bin, q_bin), g in grouped:
        n_pos1 = int((g["dprime_class"] == DPRIME_CLASS_POS1).sum())
        n_neg1 = int((g["dprime_class"] == DPRIME_CLASS_NEG1).sum())
        n_bound = n_pos1 + n_neg1
        q_mean = float(g["background_maf"].mean()) if n_bound else float("nan")
        if isinstance(q_bin, pd.Interval):
            q_mid = 0.5 * (float(q_bin.left) + float(q_bin.right))
            q_label = f"{float(q_bin.left):.3g}-{float(q_bin.right):.3g}"
        else:
            q_mid = float("nan")
            q_label = str(q_bin)
        p_plus = (n_pos1 / n_bound) if n_bound else float("nan")
        rows.append(
            {
                "pair_class": pair_class,
                "mac_bin": str(mac_bin),
                "dist_bin": str(dist_bin),
                "q_bin": q_label,
                "q_mean": q_mean,
                "q_mid": q_mid,
                "n_pos1": n_pos1,
                "n_neg1": n_neg1,
                "n_bound": n_bound,
                "p_plus": p_plus,
                "expected_p_plus": q_mean,
            }
        )
    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(rows)


def bin_signed_dprime_pplus_vs_q(
    pairs: pd.DataFrame,
    maf_threshold: float = 0.05,
    max_distance_kb: float = 5000,
    eps: float = DPRIME_CLASS_EPS,
) -> pd.DataFrame:
    """Bin |D'|≈1 pairs: observed P(D'=+1 | |D'|=1) vs background MAF ``q``.

    Rare-common uses the common MAF as ``q``. Rare-rare uses the higher
    (compared) rare MAF as ``q``. Also splits by rare MAC and distance.
    """
    pairs = _attach_signed_dprime(pairs)
    need = {"MAF_A", "MAF_B", "D_PRIME_SIGNED"}
    empty_cols = [
        "pair_class",
        "mac_bin",
        "dist_bin",
        "q_bin",
        "q_mean",
        "q_mid",
        "n_pos1",
        "n_neg1",
        "n_bound",
        "p_plus",
        "expected_p_plus",
    ]
    empty = pd.DataFrame(columns=empty_cols)
    if pairs.empty or not need <= set(pairs.columns):
        return empty
    maf_a = pairs["MAF_A"].to_numpy(dtype=np.float64)
    maf_b = pairs["MAF_B"].to_numpy(dtype=np.float64)
    signed = pairs["D_PRIME_SIGNED"].to_numpy(dtype=np.float64)
    if "DISTANCE_KB" in pairs.columns:
        distance_kb = pairs["DISTANCE_KB"].to_numpy(dtype=np.float64)
    else:
        distance_kb = np.full(len(pairs), np.nan, dtype=np.float64)
    n_pairwise = (
        pairs["N_PAIRWISE"].to_numpy(dtype=np.float64)
        if "N_PAIRWISE" in pairs.columns
        else None
    )
    obs_a = (
        pairs["OBS_CT_A"].to_numpy(dtype=np.float64)
        if "OBS_CT_A" in pairs.columns
        else None
    )
    obs_b = (
        pairs["OBS_CT_B"].to_numpy(dtype=np.float64)
        if "OBS_CT_B" in pairs.columns
        else None
    )
    roles = pair_haplotype_background_roles(
        maf_a,
        maf_b,
        maf_threshold=float(maf_threshold),
        n_pairwise=n_pairwise,
        obs_ct_a=obs_a,
        obs_ct_b=obs_b,
    )
    labels = classify_signed_dprime(signed, eps=float(eps))
    boundary = (labels == DPRIME_CLASS_POS1) | (labels == DPRIME_CLASS_NEG1)
    in_range = np.isfinite(distance_kb) & (distance_kb >= 0.0) & (
        distance_kb <= float(max_distance_kb)
    )
    keep = (
        boundary
        & in_range
        & roles["pair_class"].isin([PAIR_RARE_RARE, PAIR_RARE_COMMON]).to_numpy()
        & np.isfinite(roles["background_maf"].to_numpy(dtype=np.float64))
    )
    if not np.any(keep):
        return empty
    work = pd.DataFrame(
        {
            "pair_class": roles.loc[keep, "pair_class"].to_numpy(),
            "background_maf": roles.loc[keep, "background_maf"].to_numpy(dtype=np.float64),
            "rare_mac": np.rint(roles.loc[keep, "rare_mac"].to_numpy(dtype=np.float64)),
            "distance_kb": distance_kb[keep],
            "dprime_class": labels[keep],
        }
    )
    mac_cut = _cut_with_labels(
        work["rare_mac"], RARE_MAC_BIN_EDGES, RARE_MAC_BIN_LABELS, right=True
    )
    dist_cut = _cut_with_labels(
        work["distance_kb"],
        DISTANCE_KB_BIN_EDGES,
        DISTANCE_KB_BIN_LABELS,
        right=True,
    )
    work["mac_bin"] = mac_cut.astype(str)
    work["dist_bin"] = dist_cut.astype(str)
    q_all = np.empty(len(work), dtype=object)
    rr = work["pair_class"] == PAIR_RARE_RARE
    rc = work["pair_class"] == PAIR_RARE_COMMON
    if bool(rr.any()):
        q_all[rr.to_numpy()] = _q_bin_series(
            work.loc[rr, "background_maf"], RARE_Q_BIN_EDGES
        ).to_numpy()
    if bool(rc.any()):
        q_all[rc.to_numpy()] = _q_bin_series(
            work.loc[rc, "background_maf"], COMMON_Q_BIN_EDGES
        ).to_numpy()
    work["q_bin"] = q_all
    faceted = _aggregate_pplus_bins(work)
    overall = work.copy()
    overall["mac_bin"] = ALL_BIN_LABEL
    overall["dist_bin"] = ALL_BIN_LABEL
    overall_bins = _aggregate_pplus_bins(overall)
    return pd.concat([overall_bins, faceted], ignore_index=True)


def classify_signed_dprime(signed, eps: float = DPRIME_CLASS_EPS) -> np.ndarray:
    """Five-class labels for signed D' with tolerance ``eps`` (default 1e-6)."""
    s = np.asarray(signed, dtype=np.float64)
    out = np.full(s.shape, "", dtype=object)
    finite = np.isfinite(s)
    eps = float(eps)
    out[finite & (s >= 1.0 - eps)] = DPRIME_CLASS_POS1
    out[finite & (s <= -1.0 + eps)] = DPRIME_CLASS_NEG1
    out[finite & (s > eps) & (s < 1.0 - eps)] = DPRIME_CLASS_POSITIVE
    out[finite & (s < -eps) & (s > -1.0 + eps)] = DPRIME_CLASS_NEGATIVE
    out[finite & (np.abs(s) <= eps)] = DPRIME_CLASS_ZERO
    return out


def pairwise_complete_allele_stats(
    dosage_a: np.ndarray, dosage_b: np.ndarray
) -> tuple[float, float, int]:
    """MAF of two dosage vectors on samples called at both sites (NaN = missing)."""
    maf_a, maf_b, n, _d, _dprime = pairwise_complete_ld_stats(dosage_a, dosage_b)
    return maf_a, maf_b, n


def annotate_pairs_pairwise_maf(
    pairs: pd.DataFrame, dosage_by_id: dict[str, np.ndarray]
) -> pd.DataFrame:
    """Replace plotted MAF/R2_MAX with pairwise-complete values; keep global copies.

    Also writes signed haplotype ``D``, signed ``D_PRIME_SIGNED`` (D/D_den,
    minor-allele coding), and ``D_PRIME`` = |D'|.
    """
    out = pairs.copy()
    if "MAF_A_GLOBAL" not in out.columns:
        out["MAF_A_GLOBAL"] = out["MAF_A"]
        out["MAF_B_GLOBAL"] = out["MAF_B"]
    if "R2_MAX_GLOBAL" not in out.columns:
        out["R2_MAX_GLOBAL"] = out["R2_MAX"]
    maf_a_pair = np.full(len(out), np.nan, dtype=np.float64)
    maf_b_pair = np.full(len(out), np.nan, dtype=np.float64)
    n_pair = np.zeros(len(out), dtype=np.int32)
    d_hap = np.full(len(out), np.nan, dtype=np.float64)
    d_prime = np.full(len(out), np.nan, dtype=np.float64)
    for i, (id_a, id_b) in enumerate(zip(out["ID_A"].astype(str), out["ID_B"].astype(str))):
        da = dosage_by_id.get(id_a)
        db = dosage_by_id.get(id_b)
        if da is None or db is None:
            continue
        maf_a_pair[i], maf_b_pair[i], n_pair[i], d_hap[i], d_prime[i] = (
            pairwise_complete_ld_stats(da, db)
        )
    out["MAF_A"] = maf_a_pair
    out["MAF_B"] = maf_b_pair
    out["N_PAIRWISE"] = n_pair
    out["D"] = d_hap
    out["D_PRIME_SIGNED"] = d_prime
    out["D_PRIME"] = np.abs(d_prime)
    out["R2_MAX"] = r2_max_from_maf(maf_a_pair, maf_b_pair)
    out["R2_OVER"] = out["R2"].to_numpy(dtype=np.float64) - out["R2_MAX"].to_numpy(dtype=np.float64)
    return out


def load_dosages_from_pfile(
    pfile_prefix: str,
    variant_ids: Sequence[str],
    work_dir: str,
    threads: int = 1,
) -> dict[str, np.ndarray]:
    """Extract unique variant IDs from a PLINK2 pfile and return 0/1/2 dosage vectors."""
    ids = sorted({str(v) for v in variant_ids if str(v)})
    if not ids:
        return {}
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    extract_path = work / "pairwise_maf.extract.ids"
    extract_path.write_text("\n".join(ids) + "\n", encoding="utf-8")
    bed_prefix = work / "pairwise_maf_bed"
    cmd = [
        "plink2",
        "--pfile",
        pfile_prefix,
        "--allow-extra-chr",
        "--extract",
        str(extract_path),
        "--make-bed",
        "--threads",
        str(max(1, int(threads))),
        "--out",
        str(bed_prefix),
    ]
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"plink2 --make-bed failed for pairwise MAF extract: {proc.stderr or proc.stdout}"
        )
    return _load_bed_dosages(str(bed_prefix))


def _load_bed_dosages(bed_prefix: str) -> dict[str, np.ndarray]:
    """Read PLINK1 SNP-major .bed dosages (NaN = missing) keyed by BIM ID."""
    prefix = Path(bed_prefix)
    bim = pd.read_csv(
        prefix.with_suffix(".bim"),
        sep=r"\s+",
        engine="python",
        header=None,
        names=["chrom", "id", "cm", "pos", "a1", "a2"],
        dtype={"id": str},
    )
    fam_n = sum(1 for _ in prefix.with_suffix(".fam").open(encoding="utf-8"))
    n_var = len(bim)
    n_bytes = (fam_n + 3) // 4
    raw = np.fromfile(prefix.with_suffix(".bed"), dtype=np.uint8)
    if raw.size < 3 + n_bytes * n_var:
        raise ValueError(f"{prefix}.bed is shorter than expected for {n_var} SNPs")
    if raw[0] != 0x6C or raw[1] != 0x1B:
        raise ValueError(f"{prefix}.bed is not a PLINK bed file")
    if int(raw[2]) != 1:
        raise ValueError(f"{prefix}.bed is not SNP-major")
    body = raw[3:].reshape(n_var, n_bytes)
    dosages: dict[str, np.ndarray] = {}
    for i, vid in enumerate(bim["id"].tolist()):
        dosages[str(vid)] = _decode_bed_snp_bytes(body[i], fam_n)
    return dosages


def _decode_bed_snp_bytes(snp_bytes: np.ndarray, n_samples: int) -> np.ndarray:
    """Decode one SNP-major PLINK bed row to 0/1/2 dosages with NaN missing."""
    n_bytes = (n_samples + 3) // 4
    packed = np.asarray(snp_bytes[:n_bytes], dtype=np.uint8)
    geno = np.empty((n_bytes, 4), dtype=np.uint8)
    geno[:, 0] = packed & np.uint8(3)
    geno[:, 1] = (packed >> 2) & np.uint8(3)
    geno[:, 2] = (packed >> 4) & np.uint8(3)
    geno[:, 3] = (packed >> 6) & np.uint8(3)
    codes = geno.reshape(-1)[:n_samples]
    out = np.full(n_samples, np.nan, dtype=np.float32)
    out[codes == 0] = 0.0
    out[codes == 2] = 1.0
    out[codes == 3] = 2.0
    return out


def _iter_plink2_vcor_chunks(
    input_file: str,
    chunksize: int,
    usecols: list[str] | None = None,
):
    """Stream PLINK2 ``--r2-unphased`` ``.vcor`` in row batches.

    Prefers Polars byte-window chunks (avoids the 32-bit full-file row cap),
    then PyArrow, then pandas C-engine. Files are TAB-delimited; do not treat
    ``#`` as a comment (header starts with ``#CHROM_A``).
    """
    chunksize = int(chunksize)
    try:
        yield from _iter_plink2_vcor_chunks_polars(input_file, chunksize, usecols)
        return
    except ImportError:
        pass
    try:
        yield from _iter_plink2_vcor_chunks_pyarrow(input_file, chunksize, usecols)
        return
    except ImportError:
        pass
    yield from _iter_plink2_vcor_chunks_pandas(input_file, chunksize, usecols)


def _iter_plink2_vcor_chunks_pandas(
    input_file: str,
    chunksize: int,
    usecols: list[str] | None = None,
):
    kwargs: dict = {
        "sep": "\t",
        "chunksize": int(chunksize),
        "engine": "c",
        "dtype": _VCOR_NUMERIC_DTYPES,
    }
    if usecols is not None:
        kwargs["usecols"] = usecols
    print("[Info] Reading vcor with pandas C-engine")
    yield from pd.read_csv(input_file, **kwargs)


def _iter_plink2_vcor_chunks_polars(
    input_file: str,
    chunksize: int,
    usecols: list[str] | None = None,
):
    """Stream a ``.vcor`` with Polars without counting the whole file.

    ``scan_csv().collect_batches()`` still materialises a file-wide row
    index and fails above ``2**32-1`` pairs. File-handle ``n_rows`` also
    fails because Polars does not advance the stream. Read sequential
    byte windows and parse each window with ``read_csv``.
    """
    from io import BytesIO

    import polars as pl

    schema_overrides = {
        "POS_A": pl.Int32,
        "POS_B": pl.Int32,
        "UNPHASED_R2": pl.Float32,
        "ID_A": pl.String,
        "ID_B": pl.String,
    }
    if usecols is not None:
        schema_overrides = {k: v for k, v in schema_overrides.items() if k in usecols}
    chunksize = max(1, int(chunksize))
    avg_row_bytes = 64
    print("[Info] Reading vcor with Polars byte windows")
    with open(input_file, "rb") as fh:
        header_line = fh.readline()
        if not header_line:
            return
        while True:
            buf = fh.read(chunksize * avg_row_bytes)
            if not buf:
                return
            if not buf.endswith(b"\n"):
                buf += fh.readline()
            df = pl.read_csv(
                BytesIO(header_line + buf),
                separator="\t",
                quote_char=None,
                schema_overrides=schema_overrides,
                infer_schema_length=64,
            )
            if usecols is not None:
                df = df.select(usecols)
            n_rows = int(df.height)
            if n_rows == 0:
                return
            avg_row_bytes = max(32, len(buf) // n_rows)
            yield df.to_pandas(use_pyarrow_extension_array=True)


def _iter_plink2_vcor_chunks_pyarrow(
    input_file: str,
    chunksize: int,
    usecols: list[str] | None = None,
):
    import pyarrow as pa
    import pyarrow.csv as pacsv

    column_types = {
        "POS_A": pa.int32(),
        "POS_B": pa.int32(),
        "UNPHASED_R2": pa.float32(),
    }
    convert_opts = pacsv.ConvertOptions(
        include_columns=usecols,
        column_types=column_types if usecols is None else {k: v for k, v in column_types.items() if k in usecols},
    )
    avg_row_bytes = 64
    block_size = max(1 << 20, int(chunksize) * avg_row_bytes)
    reader = pacsv.open_csv(
        input_file,
        read_options=pacsv.ReadOptions(use_threads=True, block_size=block_size),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
        convert_options=convert_opts,
    )
    print("[Info] Reading vcor with PyArrow")
    pending: list = []
    n_pending = 0
    target = max(1, int(chunksize))
    for batch in reader:
        pending.append(batch)
        n_pending += batch.num_rows
        if n_pending >= target:
            yield pa.Table.from_batches(pending).to_pandas()
            pending = []
            n_pending = 0
    if pending:
        yield pa.Table.from_batches(pending).to_pandas()


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


def _reservoir_add(
    reservoir: np.ndarray,
    incoming: np.ndarray,
    n_seen: int,
    cap: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, int]:
    """Unbiased Algorithm R sample; never concatenates the full incoming chunk."""
    if cap <= 0 or len(incoming) == 0:
        return reservoir, n_seen
    if n_seen < cap:
        space = cap - int(n_seen)
        take = min(space, len(incoming))
        fill = incoming[:take]
        reservoir = fill.copy() if len(reservoir) == 0 else np.vstack([reservoir, fill])
        n_seen += take
        incoming = incoming[take:]
        if len(incoming) == 0:
            return reservoir, n_seen
    m = len(incoming)
    counts = n_seen + np.arange(1, m + 1, dtype=np.float64)
    keep = rng.random(m) < (float(cap) / counts)
    kept = np.flatnonzero(keep)
    if len(kept):
        slots = rng.integers(0, cap, size=len(kept))
        if not reservoir.flags.writeable:
            reservoir = reservoir.copy()
        reservoir[slots] = incoming[kept]
    return reservoir, n_seen + m


def _reservoir_add_with_ids(
    reservoir: np.ndarray,
    reservoir_ids: np.ndarray,
    incoming: np.ndarray,
    incoming_ids: np.ndarray,
    n_seen: int,
    cap: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Algorithm R on a numeric reservoir plus an aligned ID table."""
    if cap <= 0 or len(incoming) == 0:
        return reservoir, reservoir_ids, n_seen
    if n_seen < cap:
        space = cap - int(n_seen)
        take = min(space, len(incoming))
        fill = incoming[:take]
        fill_ids = incoming_ids[:take]
        if len(reservoir) == 0:
            reservoir = fill.copy()
            reservoir_ids = fill_ids.copy()
        else:
            reservoir = np.vstack([reservoir, fill])
            reservoir_ids = np.vstack([reservoir_ids, fill_ids])
        n_seen += take
        incoming = incoming[take:]
        incoming_ids = incoming_ids[take:]
        if len(incoming) == 0:
            return reservoir, reservoir_ids, n_seen
    m = len(incoming)
    counts = n_seen + np.arange(1, m + 1, dtype=np.float64)
    keep = rng.random(m) < (float(cap) / counts)
    kept = np.flatnonzero(keep)
    if len(kept):
        slots = rng.integers(0, cap, size=len(kept))
        if not reservoir.flags.writeable:
            reservoir = reservoir.copy()
        if not reservoir_ids.flags.writeable:
            reservoir_ids = reservoir_ids.copy()
        reservoir[slots] = incoming[kept]
        reservoir_ids[slots] = incoming_ids[kept]
    return reservoir, reservoir_ids, n_seen + m


def _process_chunk_arrays(
    distance_kb: np.ndarray,
    r2: np.ndarray,
    n_bins: int,
    bin_size_kb: float,
    median_sample_per_bin: int,
    scatter_max_points: int,
    seed: int,
    chrom_labels: Optional[np.ndarray] = None,
) -> Dict[str, object]:
    """Compute partial LD bin stats for one chunk."""
    local_sum = np.zeros(n_bins, dtype=np.float64)
    local_cnt = np.zeros(n_bins, dtype=np.int64)
    local_median: Dict[int, np.ndarray] = {i: np.array([], dtype=np.float32) for i in range(n_bins)}
    local_scatter = np.empty((0, 2), dtype=np.float32)
    rng = np.random.default_rng(seed)
    chrom_sum: Dict[str, np.ndarray] = {}
    chrom_cnt: Dict[str, np.ndarray] = {}

    if len(r2) == 0:
        return {
            "sum_r2": local_sum,
            "cnt_r2": local_cnt,
            "median_samples": local_median,
            "scatter_pairs": local_scatter,
            "chrom_sum": chrom_sum,
            "chrom_cnt": chrom_cnt,
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

    if chrom_labels is not None and len(chrom_labels) == len(r2):
        labels = np.asarray(chrom_labels, dtype=object)
        for lab in np.unique(labels):
            mask = labels == lab
            csum = np.zeros(n_bins, dtype=np.float64)
            ccnt = np.zeros(n_bins, dtype=np.int64)
            np.add.at(csum, bin_idx[mask], r2[mask])
            np.add.at(ccnt, bin_idx[mask], 1)
            chrom_sum[str(lab)] = csum
            chrom_cnt[str(lab)] = ccnt

    return {
        "sum_r2": local_sum,
        "cnt_r2": local_cnt,
        "median_samples": local_median,
        "scatter_pairs": local_scatter,
        "chrom_sum": chrom_sum,
        "chrom_cnt": chrom_cnt,
            "total_pairs": int(len(r2)),
    }


def _ld_decay_by_chr_table(
    chrom_sum_r2: Dict[str, np.ndarray],
    chrom_cnt_r2: Dict[str, np.ndarray],
    n_bins: int,
    bin_size_kb: float,
) -> pd.DataFrame:
    """One row per chromosome × distance bin with pair-weighted mean r2."""
    rows = []
    for chrom in sorted(chrom_cnt_r2.keys(), key=_ld_decay_chrom_sort_key):
        cnt = chrom_cnt_r2[chrom]
        sm = chrom_sum_r2[chrom]
        for i in range(n_bins):
            if cnt[i] == 0:
                continue
            bin_start = float(i * bin_size_kb)
            bin_end = float((i + 1) * bin_size_kb)
            rows.append(
                [
                    chrom,
                    bin_start,
                    bin_end,
                    (bin_start + bin_end) / 2.0,
                    float(sm[i] / cnt[i]),
                    int(cnt[i]),
                ]
            )
    return pd.DataFrame(
        rows,
        columns=["chromosome", "bin_start_kb", "bin_end_kb", "bin_mid_kb", "mean_r2", "n_pairs"],
    )


def _ld_decay_plot_frame(grouped: pd.DataFrame, by_chr: pd.DataFrame):
    """Wide table + y_specs: per-chromosome mean curves, overall mean last."""
    if by_chr is None or by_chr.empty:
        return grouped, [
            {"y_col": "mean_r2", "label": "Mean r2", "color": "steelblue", "linestyle": "-", "linewidth": 2},
            {"y_col": "median_r2", "label": "Median r2", "color": "darkorange", "linestyle": "--", "linewidth": 1.5},
        ]
    plot_df = grouped[["bin_mid_kb", "mean_r2"]].copy()
    y_specs = []
    chroms = sorted(by_chr["chromosome"].unique(), key=_ld_decay_chrom_sort_key)
    for i, chrom in enumerate(chroms):
        col = f"mean_r2_{chrom}"
        sub = by_chr.loc[by_chr["chromosome"] == chrom, ["bin_mid_kb", "mean_r2"]].rename(
            columns={"mean_r2": col}
        )
        plot_df = plot_df.merge(sub, on="bin_mid_kb", how="outer")
        y_specs.append(
            {
                "y_col": col,
                "label": str(chrom),
                "color": _ld_decay_chrom_color(str(chrom), i),
                "linestyle": "-",
                "linewidth": 1.8,
            }
        )
    plot_df = plot_df.sort_values("bin_mid_kb")
    y_specs.append(
        {
            "y_col": "mean_r2",
            "label": "All chromosomes (mean)",
            "color": LD_DECAY_OVERALL_COLOR,
            "linestyle": "--",
            "linewidth": 2.8,
            "zorder": 5,
        }
    )
    return plot_df, y_specs


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

    The decay PNG overlays one mean-r2 curve per reference chromosome and the
    overall mean as a thick dashed line.
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
    chrom_sum_r2: Dict[str, np.ndarray] = {}
    chrom_cnt_r2: Dict[str, np.ndarray] = {}

    try:
        reader = _iter_plink2_vcor_chunks(input_file, chunksize)
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
        for lab, arr in partial.get("chrom_sum", {}).items():
            if lab not in chrom_sum_r2:
                chrom_sum_r2[lab] = np.array(arr, dtype=np.float64, copy=True)
                chrom_cnt_r2[lab] = np.array(partial["chrom_cnt"][lab], dtype=np.int64, copy=True)
            else:
                chrom_sum_r2[lab] = chrom_sum_r2[lab] + arr
                chrom_cnt_r2[lab] = chrom_cnt_r2[lab] + partial["chrom_cnt"][lab]

    use_parallel = int(workers) > 1
    chunk_idx = 0
    pool = ProcessPoolExecutor(max_workers=int(workers)) if use_parallel else None

    def flush_pending_in_order() -> None:
        nonlocal next_merge_idx
        while next_merge_idx in pending_results:
            merge_partial(pending_results.pop(next_merge_idx))
            next_merge_idx += 1

    for chunk in reader:
        local = _normalize_ld_chunk(chunk, keep_chrom=True)
        if local is None:
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
        chrom_labels = None
        if "CHR_A" in local.columns:
            chrom_labels = ld_decay_chrom_labels(local["CHR_A"].to_numpy()[valid_mask])
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
                    chrom_labels,
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
                chrom_labels,
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
        save_df_to_tsv(
            pd.DataFrame(
                columns=["chromosome", "bin_start_kb", "bin_end_kb", "bin_mid_kb", "mean_r2", "n_pairs"]
            ),
            f"{output_prefix}.by_chr.info.tsv",
        )
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
    by_chr = _ld_decay_by_chr_table(chrom_sum_r2, chrom_cnt_r2, n_bins, float(bin_size_kb))
    save_df_to_tsv(by_chr, f"{output_prefix}.by_chr.info.tsv")

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

    plot_df, y_specs = _ld_decay_plot_frame(grouped, by_chr)
    plot_multi_line_series(
        data=plot_df,
        x_col="bin_mid_kb",
        y_specs=y_specs,
        title="LD decay (mean r2 vs distance by chromosome, PLINK2)",
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


def ana_ld_decay_mean_r2_threshold(
    input_file: str,
    output_prefix: str,
    *,
    r2_threshold: float = 0.5,
    max_distance_kb: float = 100.0,
    bin_size_kb: float = 1.0,
    scatter_max_points: int = 200000,
    chunksize: int = 2_000_000,
    median_sample_per_bin: int = 20000,
    workers: int = 1,
    random_seed: int = 42,
) -> pd.DataFrame:
    """Fine-bin all-pairs mean r2 and the per-chromosome distance to ``r2_threshold``.

    Writes the usual ``ana_ld_decay`` tables/plots under ``output_prefix``, plus
    ``{prefix}.chr_summary.info.tsv``. Use a tagged prefix (for example
    ``*.ld_decay.r2_0p5.1kb``) so default 50 kb decay artefacts stay.
    """
    ana_ld_decay(
        input_file,
        output_prefix,
        max_distance_kb=int(max_distance_kb) if float(max_distance_kb) >= 1 else 1,
        bin_size_kb=int(bin_size_kb) if float(bin_size_kb) >= 1 else 1,
        scatter_max_points=int(scatter_max_points),
        chunksize=int(chunksize),
        median_sample_per_bin=int(median_sample_per_bin),
        workers=int(workers),
        random_seed=int(random_seed),
    )
    by_chr = load_df_generic(f"{output_prefix}.by_chr.info.tsv")
    summary = distance_mean_r2_to_threshold(by_chr, threshold=float(r2_threshold))
    save_df_to_tsv(summary, f"{output_prefix}.chr_summary.info.tsv")
    print(
        f"[Info] Mean-r2 threshold {r2_threshold:g}: "
        f"{int(summary['reached'].sum()) if not summary.empty else 0}/"
        f"{len(summary)} chromosomes reached within {max_distance_kb:g} kb"
    )
    return summary


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
        reader = _iter_plink2_vcor_chunks(input_file, chunksize)
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


def load_afreq_ld_metrics(afreq_file: str) -> pd.DataFrame:
    """Return ID-indexed MAF and OBS_CT from a PLINK2 ``.afreq`` table."""
    df = pd.read_csv(afreq_file, sep="\t")
    df.columns = [str(c).lstrip("#") for c in df.columns]
    if "ID" not in df.columns or "ALT_FREQS" not in df.columns:
        raise ValueError(f"{afreq_file} must contain ID and ALT_FREQS columns")
    af = pd.to_numeric(df["ALT_FREQS"], errors="coerce")
    maf = np.minimum(af.to_numpy(dtype=np.float64), 1.0 - af.to_numpy(dtype=np.float64))
    obs = (
        pd.to_numeric(df["OBS_CT"], errors="coerce")
        if "OBS_CT" in df.columns
        else pd.Series(np.nan, index=df.index)
    )
    out = pd.DataFrame(
        {"MAF": maf, "OBS_CT": obs.to_numpy(dtype=np.float64)},
        index=df["ID"].astype(str),
    )
    if out.index.duplicated().any():
        out = out[~out.index.duplicated(keep="first")]
    return out


def load_maf_by_variant_id(afreq_file: str) -> pd.Series:
    """Return a unique-ID MAF series from a PLINK2 ``.afreq`` table."""
    return load_afreq_ld_metrics(afreq_file)["MAF"]


def pair_class_masks(
    maf_a: np.ndarray,
    maf_b: np.ndarray,
    maf_threshold: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Boolean masks for rare-rare and rare-common pairs (MAF < threshold is rare)."""
    rare_a = maf_a < maf_threshold
    rare_b = maf_b < maf_threshold
    valid = np.isfinite(maf_a) & np.isfinite(maf_b)
    rare_rare = valid & rare_a & rare_b
    rare_common = valid & (rare_a != rare_b)
    return rare_rare, rare_common


def pair_type_masks(
    maf_a: np.ndarray,
    maf_b: np.ndarray,
    maf_threshold: float,
) -> dict[str, np.ndarray]:
    """Boolean masks for rare-rare, rare-common, and common-common pairs."""
    rare_a = maf_a < maf_threshold
    rare_b = maf_b < maf_threshold
    valid = np.isfinite(maf_a) & np.isfinite(maf_b)
    return {
        PAIR_RARE_RARE: valid & rare_a & rare_b,
        PAIR_RARE_COMMON: valid & (rare_a != rare_b),
        PAIR_COMMON_COMMON: valid & ~rare_a & ~rare_b,
    }


def maf_ld_bin_index(maf: np.ndarray) -> np.ndarray:
    """Map MAF to within-bin class index, or -1 if outside ``(0, 0.5]``."""
    values = np.asarray(maf, dtype=np.float64)
    cuts = np.asarray(MAF_LD_BIN_EDGES[1:], dtype=np.float64)
    idx = np.digitize(values, cuts, right=True)
    out = np.full(values.shape, -1, dtype=np.int64)
    ok = (
        np.isfinite(values)
        & (values > float(MAF_LD_BIN_EDGES[0]))
        & (values <= float(MAF_LD_BIN_EDGES[-1]))
    )
    out[ok] = idx[ok]
    return out


def within_maf_bin_masks(
    maf_a: np.ndarray,
    maf_b: np.ndarray,
) -> dict[str, np.ndarray]:
    """Boolean masks for pairs whose two MAFs fall in the same LD MAF bin."""
    ia = maf_ld_bin_index(maf_a)
    ib = maf_ld_bin_index(maf_b)
    same = (ia >= 0) & (ia == ib)
    return {
        label: same & (ia == i) for i, label in enumerate(MAF_LD_BIN_LABELS)
    }


def _process_two_class_chunk(
    distance_kb: np.ndarray,
    r2: np.ndarray,
    mask_rare_rare: np.ndarray,
    mask_rare_common: np.ndarray,
    mask_common_common: np.ndarray,
    n_bins: int,
    bin_size_kb: float,
    median_sample_per_bin: int,
    scatter_max_points: int,
    seed: int,
) -> Dict[str, object]:
    """Bin stats, per-class scatter, and proportion-preserving r2 sample for one chunk."""
    class_masks = (
        (PAIR_RARE_RARE, mask_rare_rare),
        (PAIR_RARE_COMMON, mask_rare_common),
        (PAIR_COMMON_COMMON, mask_common_common),
    )
    by_class: Dict[str, object] = {}
    r2_parts: list[np.ndarray] = []
    cls_parts: list[np.ndarray] = []
    for code, (cls, mask) in enumerate(class_masks):
        by_class[cls] = _process_chunk_arrays(
            distance_kb[mask],
            r2[mask],
            n_bins,
            bin_size_kb,
            median_sample_per_bin,
            scatter_max_points,
            seed + code,
        )
        r2_c = r2[mask]
        if len(r2_c) == 0:
            continue
        r2_parts.append(r2_c.astype(np.float32))
        cls_parts.append(np.full(len(r2_c), float(code), dtype=np.float32))
    if not r2_parts:
        dist_pairs = np.empty((0, 2), dtype=np.float32)
    else:
        dist_pairs = np.column_stack(
            [np.concatenate(r2_parts), np.concatenate(cls_parts)]
        )
        if len(dist_pairs) > scatter_max_points:
            rng = np.random.default_rng(seed)
            idx = rng.choice(len(dist_pairs), size=scatter_max_points, replace=False)
            dist_pairs = dist_pairs[idx]
    return {"by_class": by_class, "dist_pairs": dist_pairs}


def _binned_stats_rows(
    sum_r2: np.ndarray,
    cnt_r2: np.ndarray,
    median_samples: Dict[int, np.ndarray],
    n_bins: int,
    bin_size_kb: float,
    pair_class: str,
) -> list[list]:
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
        rows.append(
            [pair_class, bin_start, bin_end, bin_mid, mean_r2, median_r2, int(cnt_r2[i])]
        )
    return rows


def _half_decay_kb(grouped: pd.DataFrame) -> tuple[float, float, float]:
    if grouped.empty:
        return np.nan, np.nan, np.nan
    first_mean = float(grouped["mean_r2"].iloc[0])
    half_thr = first_mean * 0.5 if pd.notna(first_mean) else np.nan
    if pd.isna(half_thr):
        return first_mean, half_thr, np.nan
    below = grouped[grouped["mean_r2"] <= half_thr]
    half_kb = float(below["bin_mid_kb"].iloc[0]) if len(below) > 0 else np.nan
    return first_mean, half_thr, half_kb


def distance_mean_r2_to_threshold(
    by_chr: pd.DataFrame,
    threshold: float = 0.5,
) -> pd.DataFrame:
    """Per-chromosome distance where pair-weighted mean r2 first falls to ``threshold``.

    Uses bin midpoints. If the first populated bin is already at or below the
    threshold, ``decay_kb`` is that bin's midpoint (the crossing is inside the
    first bin, not later along the chromosome).
    """
    columns = [
        "chromosome",
        "r2_threshold",
        "bin_size_kb",
        "first_bin_mean_r2",
        "first_bin_n_pairs",
        "reached",
        "decay_kb",
    ]
    if by_chr is None or by_chr.empty:
        return pd.DataFrame(columns=columns)
    needed = {"chromosome", "bin_mid_kb", "mean_r2"}
    if not needed.issubset(by_chr.columns):
        return pd.DataFrame(columns=columns)
    thr = float(threshold)
    rows: list[dict] = []
    for chrom, sub in by_chr.groupby("chromosome", sort=False):
        g = sub.sort_values("bin_mid_kb")
        first = g.iloc[0]
        first_mean = float(first["mean_r2"])
        first_n = int(first["n_pairs"]) if "n_pairs" in g.columns else 0
        bin_size = (
            float(first["bin_end_kb"] - first["bin_start_kb"])
            if {"bin_start_kb", "bin_end_kb"}.issubset(g.columns)
            else float("nan")
        )
        prev_mid = None
        prev_mean = None
        decay_kb = float("nan")
        reached = False
        for _, rec in g.iterrows():
            y = float(rec["mean_r2"])
            x = float(rec["bin_mid_kb"])
            if not np.isfinite(y):
                continue
            if y <= thr:
                reached = True
                if prev_mid is None:
                    decay_kb = x
                else:
                    y0 = float(prev_mean)
                    y1 = y
                    x0 = float(prev_mid)
                    x1 = x
                    if y0 == y1:
                        decay_kb = x1
                    else:
                        t = min(max((y0 - thr) / (y0 - y1), 0.0), 1.0)
                        decay_kb = x0 + t * (x1 - x0)
                break
            prev_mid, prev_mean = x, y
        rows.append(
            {
                "chromosome": str(chrom),
                "r2_threshold": thr,
                "bin_size_kb": bin_size,
                "first_bin_mean_r2": first_mean,
                "first_bin_n_pairs": first_n,
                "reached": reached,
                "decay_kb": decay_kb,
            }
        )
    out = pd.DataFrame(rows, columns=columns)
    if out.empty:
        return out
    return out.sort_values("chromosome", key=lambda s: s.map(_ld_decay_chrom_sort_key)).reset_index(
        drop=True
    )


def plot_ld_rare_common_decay_from_table(
    info_table: str | Path | pd.DataFrame,
    output_prefix: str,
    max_distance_kb: float | None = None,
) -> None:
    """Redraw class mean-r2 decay from ``*.variant.ld_decay.rare_common.info.tsv``.

    Rare-common is translucent so overlap with opaque rare-rare stays visible.
    Writes ``{prefix}.png`` only.
    """
    if isinstance(info_table, (str, Path)):
        grouped = load_df_generic(info_table)
    else:
        grouped = info_table.copy()
    if grouped is None:
        grouped = pd.DataFrame()
    mean_cols = {
        PAIR_RARE_RARE: "mean_r2_rare_rare",
        PAIR_RARE_COMMON: "mean_r2_rare_common",
        PAIR_COMMON_COMMON: "mean_r2_common_common",
    }
    mean_wide = pd.DataFrame(columns=["bin_mid_kb", *mean_cols.values()])
    if not grouped.empty:
        pivoted = grouped.pivot_table(
            index="bin_mid_kb",
            columns="pair_class",
            values="mean_r2",
            aggfunc="first",
        )
        pivoted = pivoted.reindex(columns=list(PAIR_CLASSES))
        mean_wide = pivoted.rename(columns=mean_cols).reset_index()
        for col in mean_cols.values():
            if col not in mean_wide.columns:
                mean_wide[col] = np.nan
    xlim_max = max_distance_kb
    if xlim_max is None:
        if not grouped.empty and "bin_end_kb" in grouped.columns:
            xlim_max = float(pd.to_numeric(grouped["bin_end_kb"], errors="coerce").max())
        elif not mean_wide.empty:
            xlim_max = float(pd.to_numeric(mean_wide["bin_mid_kb"], errors="coerce").max())
        else:
            xlim_max = 5000.0
    plot_multi_line_series(
        data=mean_wide,
        x_col="bin_mid_kb",
        y_specs=[
            {
                "y_col": mean_cols[cls],
                "label": f"{PAIR_CLASS_LABELS[cls]} mean r2",
                "color": PAIR_CLASS_DECAY_COLORS[cls],
                "linestyle": "-",
                "linewidth": 2,
                "marker": "o",
                "alpha": PAIR_CLASS_LINE_ALPHA[cls],
                "zorder": PAIR_CLASS_LINE_ZORDER[cls],
            }
            for cls in PAIR_CLASSES
        ],
        title="LD decay by pair class (rare-rare, rare-common, common-common)",
        filename=f"{output_prefix}.png",
        x_label="Distance (kb)",
        y_label="LD (r2)",
        xlim=(0, float(xlim_max)),
    )


def plot_ld_maf_bin_decay_from_table(
    info_table: str | Path | pd.DataFrame,
    output_prefix: str,
    max_distance_kb: float | None = None,
    title: str = "LD decay by within-bin MAF class",
) -> None:
    """Redraw within-bin MAF mean-r2 decay from ``*.variant.ld_decay.maf_bin.info.tsv``."""
    if isinstance(info_table, (str, Path)):
        grouped = load_df_generic(info_table)
    else:
        grouped = info_table.copy()
    if grouped is None:
        grouped = pd.DataFrame()
    mean_cols = {label: f"mean_r2_{i}" for i, label in enumerate(MAF_LD_BIN_LABELS)}
    mean_wide = pd.DataFrame(columns=["bin_mid_kb", *mean_cols.values()])
    if not grouped.empty:
        pivoted = grouped.pivot_table(
            index="bin_mid_kb",
            columns="pair_class",
            values="mean_r2",
            aggfunc="first",
        )
        pivoted = pivoted.reindex(columns=list(MAF_LD_BIN_LABELS))
        mean_wide = pivoted.rename(columns=mean_cols).reset_index()
        for col in mean_cols.values():
            if col not in mean_wide.columns:
                mean_wide[col] = np.nan
    xlim_max = max_distance_kb
    if xlim_max is None:
        if not grouped.empty and "bin_end_kb" in grouped.columns:
            xlim_max = float(pd.to_numeric(grouped["bin_end_kb"], errors="coerce").max())
        elif not mean_wide.empty:
            xlim_max = float(pd.to_numeric(mean_wide["bin_mid_kb"], errors="coerce").max())
        else:
            xlim_max = 5000.0
    n_class = len(MAF_LD_BIN_LABELS)
    plot_multi_line_series(
        data=mean_wide,
        x_col="bin_mid_kb",
        y_specs=[
            {
                "y_col": mean_cols[label],
                "label": f"MAF {label} mean $r^2$",
                "color": MAF_LD_BIN_COLORS[label],
                "linestyle": "-",
                "linewidth": 2,
                "marker": "o",
                "alpha": 1.0,
                "zorder": n_class - i,
            }
            for i, label in enumerate(MAF_LD_BIN_LABELS)
        ],
        title=title,
        filename=f"{output_prefix}.png",
        x_label="Distance (kb)",
        y_label=r"LD ($r^2$)",
        xlim=(0, float(xlim_max)),
    )


def combine_ld_maf_bin_decay_tables(
    info_tables: Sequence[str | Path | pd.DataFrame],
    output_prefix: str,
    max_distance_kb: float | None = None,
    title: str = "LD decay by within-bin MAF class (A+B+D)",
) -> pd.DataFrame:
    """Pool per-subgenome MAF-bin decay tables by pair-count-weighted mean r2.

    Equivalent to concatenating all A/B/D within-bin LD pairs into the same
    distance bins. Writes ``{prefix}.info.tsv``, ``{prefix}.th.tsv``, and
    ``{prefix}.png``. Median r2 is left empty (cannot be pooled from means).
    """
    frames: list[pd.DataFrame] = []
    for table in info_tables:
        if isinstance(table, pd.DataFrame):
            df = table.copy()
        else:
            df = load_df_generic(str(table))
        if df is None or df.empty:
            continue
        need = {"pair_class", "bin_mid_kb", "mean_r2", "n_pairs"}
        if not need.issubset(set(df.columns)):
            continue
        frames.append(df)
    if not frames:
        grouped = pd.DataFrame(
            columns=[
                "pair_class",
                "bin_start_kb",
                "bin_end_kb",
                "bin_mid_kb",
                "mean_r2",
                "median_r2",
                "n_pairs",
            ]
        )
    else:
        stacked = pd.concat(frames, ignore_index=True)
        stacked["mean_r2"] = pd.to_numeric(stacked["mean_r2"], errors="coerce")
        stacked["n_pairs"] = pd.to_numeric(stacked["n_pairs"], errors="coerce").fillna(0)
        stacked["bin_mid_kb"] = pd.to_numeric(stacked["bin_mid_kb"], errors="coerce")
        stacked["weighted_r2"] = stacked["mean_r2"] * stacked["n_pairs"]
        keys = ["pair_class", "bin_mid_kb"]
        extra = [c for c in ("bin_start_kb", "bin_end_kb") if c in stacked.columns]
        grouped = stacked.groupby(keys, dropna=False, as_index=False).agg(
            weighted_r2=("weighted_r2", "sum"),
            n_pairs=("n_pairs", "sum"),
            **{c: (c, "first") for c in extra},
        )
        grouped["mean_r2"] = np.where(
            grouped["n_pairs"] > 0,
            grouped["weighted_r2"] / grouped["n_pairs"],
            np.nan,
        )
        grouped["median_r2"] = np.nan
        if "bin_start_kb" not in grouped.columns:
            grouped["bin_start_kb"] = np.nan
        if "bin_end_kb" not in grouped.columns:
            grouped["bin_end_kb"] = np.nan
        grouped["bin_start_kb"] = pd.to_numeric(grouped["bin_start_kb"], errors="coerce")
        grouped["bin_end_kb"] = pd.to_numeric(grouped["bin_end_kb"], errors="coerce")
        order = {label: i for i, label in enumerate(MAF_LD_BIN_LABELS)}
        grouped["_ord"] = grouped["pair_class"].map(order).fillna(len(order))
        grouped = grouped.sort_values(["_ord", "bin_mid_kb"]).drop(
            columns=["weighted_r2", "_ord"]
        )
        grouped = grouped[
            [
                "pair_class",
                "bin_start_kb",
                "bin_end_kb",
                "bin_mid_kb",
                "mean_r2",
                "median_r2",
                "n_pairs",
            ]
        ]
        grouped = grouped.loc[grouped["n_pairs"] > 0].reset_index(drop=True)

    save_df_to_tsv(grouped, f"{output_prefix}.info.tsv")
    bin_size_kb = np.nan
    xlim_max = max_distance_kb
    if not grouped.empty:
        widths = grouped["bin_end_kb"] - grouped["bin_start_kb"]
        if widths.notna().any():
            bin_size_kb = float(widths.dropna().iloc[0])
        if xlim_max is None and grouped["bin_end_kb"].notna().any():
            xlim_max = float(grouped["bin_end_kb"].max())
    if xlim_max is None:
        xlim_max = 5000.0
    th_dict: Dict[str, object] = {
        "N_Subgenome_Tables": int(len(frames)),
        "Total_Within_Bin_Pairs": int(grouped["n_pairs"].sum()) if not grouped.empty else 0,
        "Max_Distance_kb": float(xlim_max),
        "Bin_Size_kb": float(bin_size_kb) if pd.notna(bin_size_kb) else np.nan,
    }
    for label in MAF_LD_BIN_LABELS:
        sub = grouped[grouped["pair_class"] == label].sort_values("bin_mid_kb")
        th_dict[f"Pairs_{label}"] = int(sub["n_pairs"].sum()) if not sub.empty else 0
        first_mean, half_thr, half_kb = _half_decay_kb(sub)
        th_dict[f"{label}_First_Bin_Mean_R2"] = first_mean
        th_dict[f"{label}_Half_Decay_R2_Threshold"] = half_thr
        th_dict[f"{label}_Half_Decay_kb"] = half_kb
    save_thresholds(th_dict, f"{output_prefix}.th.tsv")
    plot_ld_maf_bin_decay_from_table(
        grouped,
        output_prefix,
        max_distance_kb=float(xlim_max),
        title=title,
    )
    return grouped


def _plot_rare_common_figures(
    grouped: pd.DataFrame,
    scatter_df: pd.DataFrame,
    dist_df: pd.DataFrame,
    output_prefix: str,
    max_distance_kb: float,
    maf_scatter_df: pd.DataFrame | None = None,
    pairwise_maf: bool = False,
    maf_threshold: float = 0.05,
    bin_size_kb: float = 50,
) -> None:
    plot_ld_rare_common_decay_from_table(
        grouped,
        output_prefix,
        max_distance_kb=float(max_distance_kb),
    )
    plot_scatter_two_categories(
        data=scatter_df,
        x_col="distance_kb",
        y_col="R2",
        category_col="pair_class",
        category_labels=PAIR_CLASS_LABELS,
        category_colors=PAIR_CLASS_COLORS,
        title="LD vs distance by pair class (sampled pairs)",
        filename=f"{output_prefix}.scatter.png",
        xlabel="Distance (kb)",
        ylabel="LD (r2)",
        alpha=0.25,
        s=6,
        figure_size=(10, 6),
        xlim=(0, float(max_distance_kb)),
        ylim=(0, 1),
        legend_below=True,
    )
    dist_plot = dist_df if not dist_df.empty else scatter_df
    if dist_plot.empty:
        dist_plot = pd.DataFrame(
            {"R2": pd.Series(dtype=float), "pair_class": pd.Series(dtype=str)}
        )
    plot_stacked_histogram(
        data=dist_plot,
        x_col="R2",
        hue_col="pair_class",
        hue_order=list(PAIR_CLASSES),
        palette_map=PAIR_CLASS_COLORS,
        title="LD r2 distribution by pair class",
        filename=f"{output_prefix}.r2_dist.png",
        xlabel="LD (r2)",
        ylabel="Count",
        xlim=(0, 1),
        legend_title="Pair class",
    )
    plot_stacked_histogram(
        data=dist_plot,
        x_col="R2",
        hue_col="pair_class",
        hue_order=list(PAIR_CLASSES),
        palette_map=PAIR_CLASS_COLORS,
        title="LD r2 distribution by pair class - Log Y",
        filename=f"{output_prefix}.r2_dist.log.png",
        xlabel="LD (r2)",
        ylabel="Count",
        xlim=(0, 1),
        log_scale=True,
        legend_title="Pair class",
    )
    if maf_scatter_df is not None:
        plot_ld_maf_r2_scatter_from_table(
            maf_scatter_df,
            output_prefix,
            pairwise_maf=pairwise_maf,
            maf_max=float(maf_threshold),
            maf_threshold=float(maf_threshold),
            max_distance_kb=float(max_distance_kb),
            bin_size_kb=float(bin_size_kb),
        )


def _plot_maf_r2_scatter(
    data: pd.DataFrame,
    filename: str,
    title: str,
    s: float,
    alpha: float,
    xlabel: str = "MAF (site A)",
    ylabel: str = "MAF (site B)",
    xlim: tuple[float, float] = (0.0, 0.5),
    ylim: tuple[float, float] = (0.0, 0.5),
) -> None:
    plot_scatter_continuous_color(
        data=data,
        x_col="MAF_A",
        y_col="MAF_B",
        color_col="R2",
        title=title,
        filename=filename,
        xlabel=xlabel,
        ylabel=ylabel,
        cbar_label="LD (r2)",
        cmap=SEQUENTIAL_WHITE_RED,
        center=None,
        vmin=0.0,
        vmax=1.0,
        xlim=xlim,
        ylim=ylim,
        alpha=alpha,
        s=s,
        figure_size=(8, 8),
        line_overlays=r2_max_scatter_line_overlays(),
    )


def _plot_maf_r2_scatter_window(
    data: pd.DataFrame,
    output_prefix: str,
    maf_label: str,
    maf_max: float,
    suffix: str,
    n_total: int,
) -> None:
    """MAF-window scatter plus r2>0.2 subset, both axes limited to ``maf_max``."""
    maf_max = float(maf_max)
    window = data[(data["MAF_A"] < maf_max) & (data["MAF_B"] < maf_max)]
    if window.empty:
        print(f"[Info] No sampled pairs with both MAFs < {maf_max:g} ({suffix})")
        return
    print(
        f"[Info] {suffix} MAF scatter: {len(window)}/{n_total} pairs "
        f"with both MAFs < {maf_max:g}"
    )
    _plot_maf_r2_scatter(
        window,
        filename=f"{output_prefix}.maf_scatter_{suffix}.png",
        title=(
            f"LD r2 by pair MAF ({maf_label}, sampled pairs, "
            f"both MAF < {maf_max:g})"
        ),
        s=8,
        alpha=0.45,
        xlabel=f"MAF (site A, {maf_label})",
        ylabel=f"MAF (site B, {maf_label})",
        xlim=(0.0, maf_max),
        ylim=(0.0, maf_max),
    )
    high_r2 = window[window["R2"] > MAF_SCATTER_HIGH_R2]
    if high_r2.empty:
        print(
            f"[Info] No sampled pairs with both MAFs < {maf_max:g} "
            f"and r2 > {MAF_SCATTER_HIGH_R2} ({suffix})"
        )
        return
    print(
        f"[Info] {suffix} r2>{MAF_SCATTER_HIGH_R2} MAF scatter: "
        f"{len(high_r2)}/{len(window)} pairs"
    )
    _plot_maf_r2_scatter(
        high_r2,
        filename=f"{output_prefix}.maf_scatter_0.2_{suffix}.png",
        title=(
            f"LD r2 by pair MAF ({maf_label}, sampled pairs, "
            f"both MAF < {maf_max:g}, r2 > {MAF_SCATTER_HIGH_R2})"
        ),
        s=36,
        alpha=0.75,
        xlabel=f"MAF (site A, {maf_label})",
        ylabel=f"MAF (site B, {maf_label})",
        xlim=(0.0, maf_max),
        ylim=(0.0, maf_max),
    )


def _bin_abs_dprime_decay(
    pairs: pd.DataFrame,
    maf_threshold: float,
    max_distance_kb: float,
    bin_size_kb: float,
) -> pd.DataFrame:
    """Bin |D'| vs distance for rare-rare, rare-common, and common-common pairs."""
    cols = [
        "pair_class",
        "bin_start_kb",
        "bin_end_kb",
        "bin_mid_kb",
        "mean_dprime",
        "median_dprime",
        "n_pairs",
    ]
    empty = pd.DataFrame(columns=cols)
    if pairs.empty or not {"MAF_A", "MAF_B", "D_PRIME", "DISTANCE_KB"} <= set(pairs.columns):
        return empty
    if max_distance_kb <= 0 or bin_size_kb <= 0:
        return empty
    n_bins = int(np.ceil(float(max_distance_kb) / float(bin_size_kb)))
    if n_bins <= 0:
        return empty
    maf_a = pairs["MAF_A"].to_numpy(dtype=np.float64)
    maf_b = pairs["MAF_B"].to_numpy(dtype=np.float64)
    dprime = pairs["D_PRIME"].to_numpy(dtype=np.float64)
    distance_kb = pairs["DISTANCE_KB"].to_numpy(dtype=np.float64)
    valid = (
        np.isfinite(dprime)
        & np.isfinite(distance_kb)
        & (distance_kb >= 0.0)
        & (distance_kb <= float(max_distance_kb))
    )
    if not np.any(valid):
        return empty
    maf_a = maf_a[valid]
    maf_b = maf_b[valid]
    dprime = dprime[valid]
    distance_kb = distance_kb[valid]
    masks = pair_type_masks(maf_a, maf_b, float(maf_threshold))
    rows: list[list] = []
    for cls, mask in masks.items():
        if not np.any(mask):
            continue
        partial = _process_chunk_arrays(
            distance_kb[mask],
            dprime[mask],
            n_bins,
            float(bin_size_kb),
            median_sample_per_bin=max(len(dprime[mask]), 1),
            scatter_max_points=max(len(dprime[mask]), 1),
            seed=0,
        )
        rows.extend(
            _binned_stats_rows(
                partial["sum_r2"],
                partial["cnt_r2"],
                partial["median_samples"],
                n_bins,
                float(bin_size_kb),
                cls,
            )
        )
    if not rows:
        return empty
    grouped = pd.DataFrame(
        rows,
        columns=[
            "pair_class",
            "bin_start_kb",
            "bin_end_kb",
            "bin_mid_kb",
            "mean_r2",
            "median_r2",
            "n_pairs",
        ],
    )
    return grouped.rename(columns={"mean_r2": "mean_dprime", "median_r2": "median_dprime"})


def _attach_signed_dprime(pairs: pd.DataFrame) -> pd.DataFrame:
    """Use ``D_PRIME_SIGNED`` if present; otherwise warn (sign is not recoverable from |D'|)."""
    out = pairs.copy()
    if "D_PRIME_SIGNED" in out.columns and not bool(out["D_PRIME_SIGNED"].isna().all()):
        return out
    print(
        "[Warning] No signed D' on the pair table; skip signed-D' plots. "
        "Do not recover +1/-1 from |D'|."
    )
    return out


def bin_signed_dprime_windows(
    pairs: pd.DataFrame,
    maf_threshold: float = 0.05,
    max_distance_kb: float = 5000,
    bin_size_kb: float = 50,
    eps: float = DPRIME_CLASS_EPS,
) -> pd.DataFrame:
    """Per-distance-window signed-D' class counts and proportions.

    ``p_plus1`` / ``p_minus1`` / ``p_other`` are fractions of all pairs in the
    window. ``p_plus1_given_boundary`` and ``asymmetry`` use only |D'|≈1 pairs.
    """
    cols = [
        "pair_class",
        "bin_start_kb",
        "bin_end_kb",
        "bin_mid_kb",
        "n_all",
        "n_pos1",
        "n_neg1",
        "n_positive",
        "n_negative",
        "n_zero",
        "n_other",
        "p_plus1",
        "p_minus1",
        "p_other",
        "p_plus1_given_boundary",
        "asymmetry",
    ]
    empty = pd.DataFrame(columns=cols)
    pairs = _attach_signed_dprime(pairs)
    need = {"MAF_A", "MAF_B", "D_PRIME_SIGNED", "DISTANCE_KB"}
    if pairs.empty or not need <= set(pairs.columns):
        return empty
    if max_distance_kb <= 0 or bin_size_kb <= 0:
        return empty
    n_bins = int(np.ceil(float(max_distance_kb) / float(bin_size_kb)))
    if n_bins <= 0:
        return empty
    maf_a = pairs["MAF_A"].to_numpy(dtype=np.float64)
    maf_b = pairs["MAF_B"].to_numpy(dtype=np.float64)
    signed = pairs["D_PRIME_SIGNED"].to_numpy(dtype=np.float64)
    distance_kb = pairs["DISTANCE_KB"].to_numpy(dtype=np.float64)
    valid = (
        np.isfinite(signed)
        & np.isfinite(distance_kb)
        & (distance_kb >= 0.0)
        & (distance_kb <= float(max_distance_kb))
    )
    if not np.any(valid):
        return empty
    maf_a = maf_a[valid]
    maf_b = maf_b[valid]
    signed = signed[valid]
    distance_kb = distance_kb[valid]
    labels = classify_signed_dprime(signed, eps=eps)
    masks = pair_type_masks(maf_a, maf_b, float(maf_threshold))
    bin_idx = np.floor(distance_kb / float(bin_size_kb)).astype(np.int64)
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)
    rows: list[list] = []
    for cls in PAIR_TYPES:
        mask = masks[cls]
        if not np.any(mask):
            continue
        idx = bin_idx[mask]
        lab = labels[mask]
        for i in range(n_bins):
            in_bin = idx == i
            n_all = int(np.sum(in_bin))
            if n_all == 0:
                continue
            bin_lab = lab[in_bin]
            n_pos1 = int(np.sum(bin_lab == DPRIME_CLASS_POS1))
            n_neg1 = int(np.sum(bin_lab == DPRIME_CLASS_NEG1))
            n_positive = int(np.sum(bin_lab == DPRIME_CLASS_POSITIVE))
            n_negative = int(np.sum(bin_lab == DPRIME_CLASS_NEGATIVE))
            n_zero = int(np.sum(bin_lab == DPRIME_CLASS_ZERO))
            n_other = n_all - n_pos1 - n_neg1
            n_bound = n_pos1 + n_neg1
            start = i * float(bin_size_kb)
            end = start + float(bin_size_kb)
            mid = start + 0.5 * float(bin_size_kb)
            p_given = (n_pos1 / n_bound) if n_bound > 0 else np.nan
            asym = ((n_pos1 - n_neg1) / n_bound) if n_bound > 0 else np.nan
            rows.append(
                [
                    cls,
                    start,
                    end,
                    mid,
                    n_all,
                    n_pos1,
                    n_neg1,
                    n_positive,
                    n_negative,
                    n_zero,
                    n_other,
                    n_pos1 / n_all,
                    n_neg1 / n_all,
                    n_other / n_all,
                    p_given,
                    asym,
                ]
            )
    if not rows:
        return empty
    return pd.DataFrame(rows, columns=cols)


def _dprime_scatter_pairs(
    pairs: pd.DataFrame,
    maf_threshold: float,
    max_distance_kb: float,
) -> pd.DataFrame:
    """Rare-rare / rare-common / common-common |D'| vs distance rows."""
    empty = pd.DataFrame(columns=["distance_kb", "D_PRIME", "pair_class"])
    need = {"MAF_A", "MAF_B", "D_PRIME", "DISTANCE_KB"}
    if pairs.empty or not need <= set(pairs.columns):
        return empty
    maf_a = pairs["MAF_A"].to_numpy(dtype=np.float64)
    maf_b = pairs["MAF_B"].to_numpy(dtype=np.float64)
    dprime = pairs["D_PRIME"].to_numpy(dtype=np.float64)
    distance_kb = pairs["DISTANCE_KB"].to_numpy(dtype=np.float64)
    masks = pair_type_masks(maf_a, maf_b, float(maf_threshold))
    in_range = (
        np.isfinite(dprime)
        & np.isfinite(distance_kb)
        & (distance_kb >= 0.0)
        & (distance_kb <= float(max_distance_kb))
    )
    cls = np.full(len(pairs), "", dtype=object)
    for name, mask in masks.items():
        cls[mask] = name
    keep = in_range & (cls != "")
    if not np.any(keep):
        return empty
    return pd.DataFrame(
        {
            "distance_kb": distance_kb[keep],
            "D_PRIME": dprime[keep],
            "pair_class": cls[keep],
        }
    )


def plot_ld_dprime_decay_from_table(
    pair_table: str | Path | pd.DataFrame,
    output_prefix: str,
    maf_threshold: float = 0.05,
    max_distance_kb: float = 5000,
    bin_size_kb: float = 50,
) -> None:
    """Bin sampled-pair |D'| vs distance and plot decay by pair class.

    Expects dosage-based ``D_PRIME`` (sign-aware denominator). Does not infer
    |D'| from r2. Writes ``{prefix}.dprime.png``, ``{prefix}.dprime.scatter.png``,
    ``{prefix}.dprime.dist.png`` (fine-bin density lines),
    ``{prefix}.dprime.dist.log.png`` (``log10(density + 1)``),
    signed-D' figures A/B/C (``*.dprime.signed_prop.png``,
    ``*.dprime.signed_boundary_stack.png``, ``*.dprime.signed_interior_dist.png``),
    and ``{prefix}.dprime.info.tsv``.
    """
    if isinstance(pair_table, (str, Path)):
        df = pd.read_csv(pair_table, sep="\t")
    else:
        df = pair_table.copy()
    df = annotate_pairs_abs_dprime(df)
    if "D_PRIME" not in df.columns:
        print(
            "[Warning] No dosage-based D_PRIME; skip |D'| decay "
            "(r2 cannot recover D sign / D<0 denominator)."
        )
    grouped = _bin_abs_dprime_decay(
        df,
        maf_threshold=float(maf_threshold),
        max_distance_kb=float(max_distance_kb),
        bin_size_kb=float(bin_size_kb),
    )
    save_df_to_tsv(grouped, f"{output_prefix}.dprime.info.tsv")
    dprime_cols = {
        PAIR_RARE_RARE: "mean_dprime_rare_rare",
        PAIR_RARE_COMMON: "mean_dprime_rare_common",
        PAIR_COMMON_COMMON: "mean_dprime_common_common",
    }
    mean_wide = pd.DataFrame(columns=["bin_mid_kb", *dprime_cols.values()])
    if not grouped.empty:
        pivoted = grouped.pivot_table(
            index="bin_mid_kb",
            columns="pair_class",
            values="mean_dprime",
            aggfunc="first",
        )
        pivoted = pivoted.reindex(columns=list(PAIR_CLASSES))
        mean_wide = pivoted.rename(columns=dprime_cols).reset_index()
        for col in dprime_cols.values():
            if col not in mean_wide.columns:
                mean_wide[col] = np.nan
    plot_multi_line_series(
        data=mean_wide,
        x_col="bin_mid_kb",
        y_specs=[
            {
                "y_col": dprime_cols[cls],
                "label": f"{PAIR_CLASS_LABELS[cls]} mean |D'|",
                "color": PAIR_CLASS_COLORS[cls],
                "linestyle": "-",
                "linewidth": 2,
                "marker": "o",
            }
            for cls in PAIR_CLASSES
        ],
        title="LD |D'| decay by pair class (rare-rare, rare-common, common-common)",
        filename=f"{output_prefix}.dprime.png",
        x_label="Distance (kb)",
        y_label="LD (|D'|)",
        xlim=(0, float(max_distance_kb)),
        ylim=(0, 1),
    )
    scatter_df = _dprime_scatter_pairs(
        df, maf_threshold=float(maf_threshold), max_distance_kb=float(max_distance_kb)
    )
    plot_scatter_two_categories(
        data=scatter_df,
        x_col="distance_kb",
        y_col="D_PRIME",
        category_col="pair_class",
        category_labels=PAIR_CLASS_LABELS,
        category_colors=PAIR_CLASS_COLORS,
        title="LD |D'| vs distance by pair class (sampled pairs)",
        filename=f"{output_prefix}.dprime.scatter.png",
        xlabel="Distance (kb)",
        ylabel="LD (|D'|)",
        alpha=0.25,
        s=6,
        figure_size=(10, 6),
        xlim=(0, float(max_distance_kb)),
        ylim=(0, 1),
        legend_below=True,
    )
    plot_grouped_density_lines(
        data=scatter_df,
        x_col="D_PRIME",
        hue_col="pair_class",
        hue_order=list(PAIR_CLASSES),
        palette_map=PAIR_CLASS_COLORS,
        hue_labels=PAIR_CLASS_LABELS,
        title="LD |D'| distribution by pair class",
        filename=f"{output_prefix}.dprime.dist.png",
        xlabel="LD (|D'|)",
        ylabel="Density",
        xlim=(0, 1),
        n_bins=400,
        density=True,
        log_scale=False,
    )
    plot_grouped_density_lines(
        data=scatter_df,
        x_col="D_PRIME",
        hue_col="pair_class",
        hue_order=list(PAIR_CLASSES),
        palette_map=PAIR_CLASS_COLORS,
        hue_labels=PAIR_CLASS_LABELS,
        title="LD |D'| distribution by pair class",
        filename=f"{output_prefix}.dprime.dist.log.png",
        xlabel="LD (|D'|)",
        ylabel="Density",
        xlim=(0, 1),
        n_bins=400,
        density=True,
        log_scale=True,
    )
    n_used = int(grouped["n_pairs"].sum()) if not grouped.empty else 0
    print(
        f"[Info] |D'| decay: {n_used} classified pairs in "
        f"{len(grouped)} bins -> {output_prefix}.dprime.png; "
        f"scatter n={len(scatter_df)} -> {output_prefix}.dprime.scatter.png; "
        f"density -> {output_prefix}.dprime.dist.png"
    )
    _plot_signed_dprime_figures(
        df,
        output_prefix,
        maf_threshold=float(maf_threshold),
        max_distance_kb=float(max_distance_kb),
        bin_size_kb=float(bin_size_kb),
    )


def _plot_signed_dprime_figures(
    pairs: pd.DataFrame,
    output_prefix: str,
    maf_threshold: float,
    max_distance_kb: float,
    bin_size_kb: float,
) -> None:
    """Figures A/B/C for signed D' vs distance (minor-allele coding)."""
    pairs = _attach_signed_dprime(pairs)
    windows = bin_signed_dprime_windows(
        pairs,
        maf_threshold=float(maf_threshold),
        max_distance_kb=float(max_distance_kb),
        bin_size_kb=float(bin_size_kb),
    )
    save_df_to_tsv(windows, f"{output_prefix}.dprime.signed.info.tsv")
    facet_order = list(PAIR_CLASSES)
    plot_faceted_multi_line_series(
        data=windows,
        x_col="bin_mid_kb",
        facet_col="pair_class",
        facet_order=facet_order,
        y_specs=[
            {
                "y_col": "p_plus1",
                "label": r"$P(D'\approx +1)$",
                "color": SIGNED_DPRIME_PROP_COLORS["p_plus1"],
                "linestyle": "-",
                "linewidth": 1.8,
            },
            {
                "y_col": "p_minus1",
                "label": r"$P(D'\approx -1)$",
                "color": SIGNED_DPRIME_PROP_COLORS["p_minus1"],
                "linestyle": "-",
                "linewidth": 1.8,
            },
            {
                "y_col": "p_other",
                "label": r"$P(|D'|<1)$",
                "color": SIGNED_DPRIME_PROP_COLORS["p_other"],
                "linestyle": "--",
                "linewidth": 1.5,
            },
        ],
        title="Signed D' class proportions vs distance",
        filename=f"{output_prefix}.dprime.signed_prop.png",
        xlabel="Distance (kb)",
        ylabel="Pair proportion",
        xlim=(0, float(max_distance_kb)),
        ylim=(0, 1),
        facet_labels=PAIR_CLASS_LABELS,
    )
    stack_df = windows.copy()
    if not stack_df.empty:
        n_bound = stack_df["n_pos1"] + stack_df["n_neg1"]
        stack_df[DPRIME_CLASS_POS1] = np.where(n_bound > 0, stack_df["n_pos1"] / n_bound, np.nan)
        stack_df[DPRIME_CLASS_NEG1] = np.where(n_bound > 0, stack_df["n_neg1"] / n_bound, np.nan)
    plot_faceted_stacked_bars(
        data=stack_df,
        x_col="bin_mid_kb",
        facet_col="pair_class",
        facet_order=facet_order,
        stack_cols=[DPRIME_CLASS_POS1, DPRIME_CLASS_NEG1],
        stack_colors=SIGNED_DPRIME_BOUNDARY_COLORS,
        stack_labels={
            DPRIME_CLASS_POS1: r"$D'\approx +1$",
            DPRIME_CLASS_NEG1: r"$D'\approx -1$",
        },
        title=r"Composition of $|D'|=1$ pairs ($+1$ vs $-1$)",
        filename=f"{output_prefix}.dprime.signed_boundary_stack.png",
        xlabel="Distance (kb)",
        ylabel="Proportion of |D'|=1 pairs",
        facet_labels=PAIR_CLASS_LABELS,
        bar_width=float(bin_size_kb),
        xlim=(0, float(max_distance_kb)),
        ylim=(0.0, 1.0),
    )
    interior = pd.DataFrame(columns=["D_PRIME_SIGNED", "pair_class"])
    if "D_PRIME_SIGNED" in pairs.columns and not pairs.empty:
        signed = pairs["D_PRIME_SIGNED"].to_numpy(dtype=np.float64)
        maf_a = pairs["MAF_A"].to_numpy(dtype=np.float64)
        maf_b = pairs["MAF_B"].to_numpy(dtype=np.float64)
        labels = classify_signed_dprime(signed)
        masks = pair_type_masks(maf_a, maf_b, float(maf_threshold))
        keep = (labels == DPRIME_CLASS_POSITIVE) | (labels == DPRIME_CLASS_NEGATIVE) | (
            labels == DPRIME_CLASS_ZERO
        )
        cls = np.full(len(pairs), "", dtype=object)
        for name, mask in masks.items():
            cls[mask] = name
        interior = pd.DataFrame(
            {
                "D_PRIME_SIGNED": signed[keep & (cls != "")],
                "pair_class": cls[keep & (cls != "")],
            }
        )
        interior = interior[interior["pair_class"].isin(facet_order)]
    plot_grouped_density_lines(
        data=interior,
        x_col="D_PRIME_SIGNED",
        hue_col="pair_class",
        hue_order=facet_order,
        palette_map=PAIR_CLASS_COLORS,
        hue_labels=PAIR_CLASS_LABELS,
        title=r"Signed $D'$ density excluding $|D'|\approx 1$",
        filename=f"{output_prefix}.dprime.signed_interior_dist.png",
        xlabel="Signed D'",
        ylabel="Density",
        xlim=(-1.0, 1.0),
        n_bins=400,
        density=True,
        log_scale=False,
    )
    _plot_signed_dprime_pplus_vs_q(
        pairs,
        output_prefix,
        maf_threshold=float(maf_threshold),
        max_distance_kb=float(max_distance_kb),
    )
    print(
        f"[Info] signed D' windows n_rows={len(windows)} -> "
        f"{output_prefix}.dprime.signed_prop.png / "
        f"{output_prefix}.dprime.signed_boundary_stack.png / "
        f"{output_prefix}.dprime.signed_interior_dist.png / "
        f"{output_prefix}.dprime.signed_pplus_vs_q.png"
    )


def _draw_pplus_vs_q_overall(
    table: pd.DataFrame,
    filename: str,
    maf_threshold: float,
    title: str,
) -> None:
    """Overall P+ vs q (no MAC×distance facets)."""
    overall = pd.DataFrame()
    if not table.empty:
        overall = table[
            (table["mac_bin"] == ALL_BIN_LABEL) & (table["dist_bin"] == ALL_BIN_LABEL)
        ].copy()
        overall["panel"] = "all"
        overall = overall.sort_values(["pair_class", "q_mean"])
    plot_scatter_facet_grid(
        data=overall,
        x_col="q_mean",
        y_col="p_plus",
        row_col="pair_class",
        col_col="panel",
        row_order=list(PAIR_CLASSES),
        col_order=["all"],
        title=title,
        filename=filename,
        xlabel=r"Background MAF $q$",
        ylabel=r"$P_{+}$",
        ylim=(0.0, 1.0),
        identity_line=True,
        identity_label=r"$P_{+}=q$",
        size_col="n_bound",
        row_labels=PAIR_CLASS_LABELS,
        col_labels={"all": ""},
        sharex=False,
        sharey=False,
        row_xlim={
            PAIR_RARE_RARE: (0.0, float(maf_threshold)),
            PAIR_RARE_COMMON: (0.0, 0.5),
            PAIR_COMMON_COMMON: (0.0, 0.5),
        },
        row_ylim={
            PAIR_RARE_RARE: (0.0, float(maf_threshold)),
            PAIR_RARE_COMMON: (0.0, 0.5),
            PAIR_COMMON_COMMON: (0.0, 0.5),
        },
    )


def plot_signed_dprime_pplus_vs_q_overall(
    pair_table: str | Path | pd.DataFrame,
    output_prefix: str,
    maf_threshold: float = 0.05,
    max_distance_kb: float = 1000,
    filename_suffix: str = "1mb",
) -> None:
    """Overall P+ vs q using only pairs with distance ≤ ``max_distance_kb``.

    Writes ``{prefix}.dprime.signed_pplus_vs_q.{suffix}.png`` and
    ``{prefix}.dprime.signed_pplus_vs_q.{suffix}.info.tsv``. Does not redraw
    the default (wider-window) overall figure or MAC×distance grids.
    """
    if isinstance(pair_table, (str, Path)):
        df = pd.read_csv(pair_table, sep="\t")
    else:
        df = pair_table.copy()
    df = annotate_pairs_abs_dprime(df)
    suffix = str(filename_suffix).strip() or "1mb"
    table = bin_signed_dprime_pplus_vs_q(
        df,
        maf_threshold=float(maf_threshold),
        max_distance_kb=float(max_distance_kb),
    )
    save_df_to_tsv(table, f"{output_prefix}.dprime.signed_pplus_vs_q.{suffix}.info.tsv")
    mb = float(max_distance_kb) / 1000.0
    dist_txt = f"{int(round(mb))} Mb" if abs(mb - round(mb)) < 1e-9 else f"{mb:g} Mb"
    _draw_pplus_vs_q_overall(
        table,
        filename=f"{output_prefix}.dprime.signed_pplus_vs_q.{suffix}.png",
        maf_threshold=float(maf_threshold),
        title=(
            r"Observed $P(D'\approx +1 \mid |D'|=1)$ vs background MAF $q$"
            rf" (distance $\leq$ {dist_txt})"
        ),
    )


def _plot_signed_dprime_pplus_vs_q(
    pairs: pd.DataFrame,
    output_prefix: str,
    maf_threshold: float,
    max_distance_kb: float,
) -> None:
    """Observed P(D'=+1 | |D'|=1) vs background MAF q, overall and MAC×distance."""
    table = bin_signed_dprime_pplus_vs_q(
        pairs,
        maf_threshold=float(maf_threshold),
        max_distance_kb=float(max_distance_kb),
    )
    save_df_to_tsv(table, f"{output_prefix}.dprime.signed_pplus_vs_q.info.tsv")
    _draw_pplus_vs_q_overall(
        table,
        filename=f"{output_prefix}.dprime.signed_pplus_vs_q.png",
        maf_threshold=float(maf_threshold),
        title=r"Observed $P(D'\approx +1 \mid |D'|=1)$ vs background MAF $q$",
    )
    faceted = pd.DataFrame()
    if not table.empty:
        faceted = table[
            (table["mac_bin"] != ALL_BIN_LABEL)
            & (table["dist_bin"] != ALL_BIN_LABEL)
            & (~table["mac_bin"].isin(["nan", "NaN"]))
            & (~table["dist_bin"].isin(["nan", "NaN"]))
        ].copy()
    for cls, slug, lim in (
        (PAIR_RARE_RARE, "rare_rare", (0.0, float(maf_threshold))),
        (PAIR_RARE_COMMON, "rare_common", (0.0, 0.5)),
        (PAIR_COMMON_COMMON, "common_common", (0.0, 0.5)),
    ):
        sub = (
            faceted[faceted["pair_class"] == cls].sort_values(
                ["mac_bin", "dist_bin", "q_mean"]
            )
            if not faceted.empty
            else faceted
        )
        plot_scatter_facet_grid(
            data=sub,
            x_col="q_mean",
            y_col="p_plus",
            row_col="mac_bin",
            col_col="dist_bin",
            row_order=list(RARE_MAC_BIN_LABELS),
            col_order=list(DISTANCE_KB_BIN_LABELS),
            title=(
                rf"{PAIR_CLASS_LABELS[cls]}: $P_{{+}}$ vs $q$ by rare MAC and distance"
            ),
            filename=f"{output_prefix}.dprime.signed_pplus_vs_q.{slug}.mac_dist.png",
            xlabel=r"Background MAF $q$",
            ylabel=r"$P_{+}$",
            xlim=lim,
            ylim=lim,
            identity_line=True,
            identity_label=r"$P_{+}=q$",
            size_col="n_bound",
        )


def plot_ld_maf_r2_scatter_from_table(
    pair_table: str | Path | pd.DataFrame,
    output_prefix: str,
    pairwise_maf: bool | None = None,
    maf_max: float = 0.05,
    ultra_maf_max: float = MAF_SCATTER_ULTRARARE_MAX,
    maf_threshold: float | None = None,
    max_distance_kb: float = 5000,
    bin_size_kb: float = 50,
    pfile_prefix: str | None = None,
    workers: int = 1,
) -> None:
    """Redraw MAF-vs-MAF r2 scatters (with r_max contours) from a pair info TSV.

    When ``ID_A``/``ID_B`` are present, adds ``DISTANCE_KB``. |D'| decay uses
    dosage-based ``D_PRIME`` (from the table or recomputed via ``pfile_prefix``).
    """
    if isinstance(pair_table, (str, Path)):
        df = pd.read_csv(pair_table, sep="\t")
    else:
        df = pair_table.copy()
    has_ids = {"ID_A", "ID_B"} <= set(df.columns)
    if maf_threshold is None:
        maf_threshold = float(maf_max)
    if not {"MAF_A", "MAF_B", "R2"} <= set(df.columns):
        print(f"[Warning] No MAF scatter rows for {output_prefix}")
        if has_ids:
            plot_ld_dprime_decay_from_table(
                df,
                output_prefix,
                maf_threshold=float(maf_threshold),
                max_distance_kb=float(max_distance_kb),
                bin_size_kb=float(bin_size_kb),
            )
        return
    if df.empty:
        print(f"[Warning] No MAF scatter rows for {output_prefix}")
        if has_ids:
            plot_ld_dprime_decay_from_table(
                df,
                output_prefix,
                maf_threshold=float(maf_threshold),
                max_distance_kb=float(max_distance_kb),
                bin_size_kb=float(bin_size_kb),
            )
        return
    if pairwise_maf is None:
        pairwise_maf = "N_PAIRWISE" in df.columns
    maf_label = "pairwise-complete MAF" if pairwise_maf else "global MAF"
    _plot_maf_r2_scatter(
        df,
        filename=f"{output_prefix}.maf_scatter.png",
        title=f"LD r2 by pair MAF ({maf_label}, sampled pairs)",
        s=8,
        alpha=0.45,
        xlabel=f"MAF (site A, {maf_label})",
        ylabel=f"MAF (site B, {maf_label})",
    )
    high_r2 = df[df["R2"] > MAF_SCATTER_HIGH_R2]
    if not high_r2.empty:
        _plot_maf_r2_scatter(
            high_r2,
            filename=f"{output_prefix}.maf_scatter_0.2.png",
            title=f"LD r2 by pair MAF ({maf_label}, sampled pairs, r2 > 0.2)",
            s=36,
            alpha=0.75,
            xlabel=f"MAF (site A, {maf_label})",
            ylabel=f"MAF (site B, {maf_label})",
        )
    _plot_maf_r2_scatter_window(
        df, output_prefix, maf_label, float(maf_max), "rare", len(df)
    )
    _plot_maf_r2_scatter_window(
        df, output_prefix, maf_label, float(ultra_maf_max), "ultrarare", len(df)
    )
    if pfile_prefix and has_ids and not df.empty:
        need_signed = (
            "D_PRIME_SIGNED" not in df.columns
            or bool(df["D_PRIME_SIGNED"].isna().all())
        )
        if need_signed:
            variant_ids = (
                pd.concat(
                    [df["ID_A"].astype(str), df["ID_B"].astype(str)],
                    ignore_index=True,
                )
                .unique()
                .tolist()
            )
            print(
                f"[Info] Computing signed D' from pairwise pfile dosages "
                f"(minor-allele coding) {pfile_prefix}"
            )
            dosage_by_id = load_dosages_from_pfile(
                pfile_prefix,
                variant_ids,
                work_dir=str(Path.cwd()),
                threads=int(workers),
            )
            df = annotate_pairs_pairwise_maf(df, dosage_by_id)
            pairwise_maf = True
        else:
            print("[Info] Using dosage-based signed D' already on the pair table")
    if has_ids:
        annotated = annotate_pairs_abs_dprime(df)
        save_df_to_tsv(
            annotated,
            f"{output_prefix}.maf_scatter.info.tsv",
            float_format="%.8g",
        )
        plot_ld_dprime_decay_from_table(
            annotated,
            output_prefix,
            maf_threshold=float(maf_threshold),
            max_distance_kb=float(max_distance_kb),
            bin_size_kb=float(bin_size_kb),
        )


def ana_ld_decay_rare_common(
    input_file: str,
    afreq_file: str,
    output_prefix: str,
    maf_threshold: float = 0.05,
    max_distance_kb: int = 5000,
    bin_size_kb: int = 50,
    scatter_max_points: int = 200000,
    maf_scatter_points: int = 20000,
    chunksize: int = 2_000_000,
    median_sample_per_bin: int = 20000,
    workers: int = 1,
    random_seed: int = 42,
    pfile_prefix: str | None = None,
) -> None:
    """
    Classify mixed-MAF PLINK2 LD pairs as rare-rare, rare-common, or
    common-common and write three-colour decay, stacked r2 distribution, and scatter.
    Rare is MAF < ``maf_threshold``.
    Also writes a sampled MAF-vs-MAF scatter coloured by r2 (includes all pair classes),
    a scatter of those points with r2 > 0.2, rare-rare zooms (both MAF <
    ``maf_threshold``), and ultra-rare zooms (both MAF < 0.01) of both panels.
    When ``pfile_prefix`` is set, sampled-pair MAF axes use pairwise-complete
    frequencies from the pfile; the pair table is ``{prefix}.maf_scatter.info.tsv``
    (adds signed ``D``, dosage-based ``D_PRIME``, and ``DISTANCE_KB``). |D'| uses
    the sign-aware Lewontin denominator (not r2 / r_max^2). Also writes
    ``{prefix}.dprime.png`` and ``{prefix}.dprime.scatter.png``.
    """
    if max_distance_kb <= 0 or bin_size_kb <= 0 or maf_threshold <= 0:
        print("[Error] Invalid LD rare/common settings.")
        return
    n_bins = int(np.ceil(float(max_distance_kb) / float(bin_size_kb)))
    if n_bins <= 0:
        print("[Error] Invalid number of bins for LD decay.")
        return

    try:
        afreq_metrics = load_afreq_ld_metrics(afreq_file)
    except Exception as e:
        print(f"[Error] Failed to read afreq {afreq_file}: {e}")
        return
    maf_by_id = afreq_metrics["MAF"]
    obs_by_id = afreq_metrics["OBS_CT"]
    print(
        f"[Info] Loaded MAF for {len(maf_by_id)} variants; "
        f"rare = MAF < {maf_threshold}"
    )

    sum_r2 = {cls: np.zeros(n_bins, dtype=np.float64) for cls in PAIR_CLASSES}
    cnt_r2 = {cls: np.zeros(n_bins, dtype=np.int64) for cls in PAIR_CLASSES}
    median_samples: Dict[str, Dict[int, np.ndarray]] = {
        cls: {i: np.array([], dtype=np.float32) for i in range(n_bins)}
        for cls in PAIR_CLASSES
    }
    scatter_pairs = {cls: np.empty((0, 2), dtype=np.float32) for cls in PAIR_CLASSES}
    dist_pairs = np.empty((0, 2), dtype=np.float32)
    maf_scatter = np.empty((0, 6), dtype=np.float32)
    maf_scatter_ids = np.empty((0, 2), dtype=object)
    maf_scatter_seen = 0
    rng = np.random.default_rng(random_seed)
    total_pairs = 0
    used_counts = {cls: 0 for cls in PAIR_CLASSES}
    dropped_unmapped = 0
    n_r2_gt_r2max = 0

    try:
        reader = _iter_plink2_vcor_chunks(
            input_file,
            chunksize,
            usecols=["POS_A", "ID_A", "POS_B", "ID_B", "UNPHASED_R2"],
        )
    except Exception as e:
        print(f"[Error] Failed to read {input_file}: {e}")
        return

    futures = deque()
    pending_results: Dict[int, Dict[str, object]] = {}
    next_merge_idx = 0
    max_queue = max(1, int(workers) * 2)

    def merge_scatter(current: np.ndarray, incoming: np.ndarray, max_points: int | None = None) -> np.ndarray:
        cap = int(scatter_max_points if max_points is None else max_points)
        if cap <= 0 or len(incoming) == 0:
            return current
        if len(current) == 0:
            if len(incoming) <= cap:
                return incoming
            idx = rng.choice(len(incoming), size=cap, replace=False)
            return incoming[idx]
        merged = np.vstack([current, incoming])
        if len(merged) <= cap:
            return merged
        idx = rng.choice(len(merged), size=cap, replace=False)
        return merged[idx]

    def merge_partial(partial: Dict[str, object]) -> None:
        nonlocal dist_pairs
        by_class = partial["by_class"]
        for cls, stats in by_class.items():
            sum_r2[cls][:] = sum_r2[cls] + stats["sum_r2"]
            cnt_r2[cls][:] = cnt_r2[cls] + stats["cnt_r2"]
            used_counts[cls] += int(stats["total_pairs"])
            for b, vals in stats["median_samples"].items():
                if len(vals) == 0:
                    continue
                median_samples[cls][int(b)] = _reservoir_merge(
                    median_samples[cls][int(b)], vals, median_sample_per_bin, rng
                )
            scatter_pairs[cls] = merge_scatter(scatter_pairs[cls], stats["scatter_pairs"])
        dist_pairs = merge_scatter(dist_pairs, partial["dist_pairs"])

    use_parallel = int(workers) > 1
    chunk_idx = 0
    pool = ProcessPoolExecutor(max_workers=int(workers)) if use_parallel else None

    def flush_pending_in_order() -> None:
        nonlocal next_merge_idx
        while next_merge_idx in pending_results:
            merge_partial(pending_results.pop(next_merge_idx))
            next_merge_idx += 1

    for chunk in reader:
        local = _normalize_ld_chunk(chunk, keep_ids=True)
        if local is None or local.empty:
            continue

        distance_kb = (local["BP_B"] - local["BP_A"]).abs().to_numpy(dtype=np.float64) / 1000.0
        r2 = local["R2"].to_numpy(dtype=np.float64)
        valid_mask = (distance_kb >= 0.0) & (distance_kb <= float(max_distance_kb)) & np.isfinite(r2)
        if not np.any(valid_mask):
            continue

        local = local.loc[valid_mask]
        distance_kb = distance_kb[valid_mask]
        r2 = r2[valid_mask]
        total_pairs += int(len(r2))

        maf_a = _ids_to_values(local["ID_A"], maf_by_id)
        maf_b = _ids_to_values(local["ID_B"], maf_by_id)
        obs_a = _ids_to_values(local["ID_A"], obs_by_id)
        obs_b = _ids_to_values(local["ID_B"], obs_by_id)
        unmapped = ~np.isfinite(maf_a) | ~np.isfinite(maf_b)
        dropped_unmapped += int(np.sum(unmapped))

        mapped = np.isfinite(maf_a) & np.isfinite(maf_b)
        r2_max = r2_max_from_maf(maf_a, maf_b)
        n_r2_gt_r2max += int(np.sum(mapped & np.isfinite(r2_max) & (r2 > r2_max)))
        if np.any(mapped) and int(maf_scatter_points) > 0:
            incoming_maf = np.column_stack(
                [
                    maf_a[mapped].astype(np.float32),
                    maf_b[mapped].astype(np.float32),
                    r2[mapped].astype(np.float32),
                    r2_max[mapped].astype(np.float32),
                    obs_a[mapped].astype(np.float32),
                    obs_b[mapped].astype(np.float32),
                ]
            )
            incoming_ids = np.column_stack(
                [
                    local["ID_A"].to_numpy()[mapped],
                    local["ID_B"].to_numpy()[mapped],
                ]
            )
            maf_scatter, maf_scatter_ids, maf_scatter_seen = _reservoir_add_with_ids(
                maf_scatter,
                maf_scatter_ids,
                incoming_maf,
                incoming_ids,
                maf_scatter_seen,
                int(maf_scatter_points),
                rng,
            )

        masks = pair_type_masks(maf_a, maf_b, maf_threshold)
        mask_rr = masks[PAIR_RARE_RARE]
        mask_rc = masks[PAIR_RARE_COMMON]
        mask_cc = masks[PAIR_COMMON_COMMON]
        if not (np.any(mask_rr) or np.any(mask_rc) or np.any(mask_cc)):
            continue

        if use_parallel:
            futures.append(
                (
                    chunk_idx,
                    pool.submit(
                        _process_two_class_chunk,
                        distance_kb,
                        r2,
                        mask_rr,
                        mask_rc,
                        mask_cc,
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
            partial = _process_two_class_chunk(
                distance_kb,
                r2,
                mask_rr,
                mask_rc,
                mask_cc,
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

    rows = []
    for cls in PAIR_CLASSES:
        rows.extend(
            _binned_stats_rows(
                sum_r2[cls],
                cnt_r2[cls],
                median_samples[cls],
                n_bins,
                float(bin_size_kb),
                cls,
            )
        )
    grouped = pd.DataFrame(
        rows,
        columns=[
            "pair_class",
            "bin_start_kb",
            "bin_end_kb",
            "bin_mid_kb",
            "mean_r2",
            "median_r2",
            "n_pairs",
        ],
    )
    save_df_to_tsv(grouped, f"{output_prefix}.info.tsv")

    th_dict: Dict[str, object] = {
        "Total_LD_pairs": int(total_pairs),
        "Pairs_Rare_Rare": int(used_counts[PAIR_RARE_RARE]),
        "Pairs_Rare_Common": int(used_counts[PAIR_RARE_COMMON]),
        "Pairs_Common_Common": int(used_counts[PAIR_COMMON_COMMON]),
        "Pairs_Dropped_Unmapped": int(dropped_unmapped),
        "MAF_Threshold": float(maf_threshold),
        "Max_Distance_kb": float(max_distance_kb),
        "Bin_Size_kb": float(bin_size_kb),
        "N_Variants_With_MAF": int(len(maf_by_id)),
        "Maf_Scatter_Pairs": int(len(maf_scatter)),
        "Maf_Scatter_Pairs_R2_gt_0.2": (
            int(np.sum(maf_scatter[:, 2] > MAF_SCATTER_HIGH_R2)) if len(maf_scatter) else 0
        ),
        "Pairs_R2_gt_R2max_GlobalMAF": int(n_r2_gt_r2max),
        "Maf_Scatter_Pairs_R2_gt_R2max": (
            int(np.sum(maf_scatter[:, 2] > maf_scatter[:, 3])) if len(maf_scatter) else 0
        ),
    }
    for cls, key_prefix in (
        (PAIR_RARE_RARE, "Rare_Rare"),
        (PAIR_RARE_COMMON, "Rare_Common"),
        (PAIR_COMMON_COMMON, "Common_Common"),
    ):
        sub = grouped[grouped["pair_class"] == cls].sort_values("bin_mid_kb")
        first_mean, half_thr, half_kb = _half_decay_kb(sub)
        th_dict[f"{key_prefix}_First_Bin_Mean_R2"] = first_mean
        th_dict[f"{key_prefix}_Half_Decay_R2_Threshold"] = half_thr
        th_dict[f"{key_prefix}_Half_Decay_Distance_kb"] = half_kb
    save_thresholds(th_dict, f"{output_prefix}.th.tsv")
    print(
        f"[Info] rare-rare={used_counts[PAIR_RARE_RARE]:,} "
        f"rare-common={used_counts[PAIR_RARE_COMMON]:,} "
        f"common-common={used_counts[PAIR_COMMON_COMMON]:,} "
        f"unmapped={dropped_unmapped:,} "
        f"r2>r2_max(global MAF)={n_r2_gt_r2max:,}"
    )

    scatter_frames = []
    for cls in PAIR_CLASSES:
        arr = scatter_pairs[cls]
        if len(arr) == 0:
            continue
        scatter_frames.append(
            pd.DataFrame({"distance_kb": arr[:, 0], "R2": arr[:, 1], "pair_class": cls})
        )
    scatter_df = (
        pd.concat(scatter_frames, ignore_index=True)
        if scatter_frames
        else pd.DataFrame(columns=["distance_kb", "R2", "pair_class"])
    )
    if len(dist_pairs) == 0:
        dist_df = scatter_df[["R2", "pair_class"]].copy() if not scatter_df.empty else scatter_df
    else:
        class_codes = np.clip(np.rint(dist_pairs[:, 1]).astype(np.int64), 0, 2)
        cls_labels = np.take(np.array(PAIR_CLASSES, dtype=object), class_codes)
        dist_df = pd.DataFrame({"R2": dist_pairs[:, 0], "pair_class": cls_labels})

    maf_scatter_cols = [
        "ID_A",
        "ID_B",
        "MAF_A",
        "MAF_B",
        "R2",
        "R2_MAX",
        "OBS_CT_A",
        "OBS_CT_B",
        "R2_OVER",
    ]
    if len(maf_scatter) > 0:
        maf_scatter_df = pd.DataFrame(
            {
                "ID_A": maf_scatter_ids[:, 0],
                "ID_B": maf_scatter_ids[:, 1],
                "MAF_A": maf_scatter[:, 0],
                "MAF_B": maf_scatter[:, 1],
                "R2": maf_scatter[:, 2],
                "R2_MAX": maf_scatter[:, 3],
                "OBS_CT_A": maf_scatter[:, 4],
                "OBS_CT_B": maf_scatter[:, 5],
            }
        )
        maf_scatter_df["R2_OVER"] = maf_scatter_df["R2"] - maf_scatter_df["R2_MAX"]
    else:
        maf_scatter_df = pd.DataFrame(columns=maf_scatter_cols)

    used_pairwise = False
    if pfile_prefix and len(maf_scatter_df) > 0:
        variant_ids = pd.concat(
            [maf_scatter_df["ID_A"].astype(str), maf_scatter_df["ID_B"].astype(str)],
            ignore_index=True,
        ).unique().tolist()
        try:
            dosage_by_id = load_dosages_from_pfile(
                pfile_prefix,
                variant_ids,
                work_dir=str(Path.cwd()),
                threads=int(workers),
            )
            maf_scatter_df = annotate_pairs_pairwise_maf(maf_scatter_df, dosage_by_id)
            used_pairwise = True
            n_pair_ok = int((maf_scatter_df["N_PAIRWISE"] > 0).sum())
            print(
                f"[Info] Pairwise-complete MAF for {n_pair_ok}/{len(maf_scatter_df)} "
                f"sampled pairs from {pfile_prefix}"
            )
        except Exception as e:
            print(f"[Error] Pairwise-complete MAF from pfile failed: {e}")
            raise

    save_df_to_tsv(maf_scatter_df, f"{output_prefix}.maf_scatter.info.tsv")
    if used_pairwise and len(maf_scatter_df) > 0:
        th_dict["Maf_Scatter_Pairs_R2_gt_R2max"] = int(
            np.sum(
                np.isfinite(maf_scatter_df["R2_MAX"].to_numpy(dtype=np.float64))
                & (
                    maf_scatter_df["R2"].to_numpy(dtype=np.float64)
                    > maf_scatter_df["R2_MAX"].to_numpy(dtype=np.float64)
                )
            )
        )
        th_dict["Maf_Scatter_Median_N_Pairwise"] = float(maf_scatter_df["N_PAIRWISE"].median())
        save_thresholds(th_dict, f"{output_prefix}.th.tsv")
    _plot_rare_common_figures(
        grouped,
        scatter_df,
        dist_df,
        output_prefix,
        float(max_distance_kb),
        maf_scatter_df=maf_scatter_df,
        pairwise_maf=used_pairwise,
        maf_threshold=float(maf_threshold),
        bin_size_kb=float(bin_size_kb),
    )


def ana_ld_decay_maf_bins(
    input_file: str,
    afreq_file: str,
    output_prefix: str,
    max_distance_kb: int = 5000,
    bin_size_kb: int = 50,
    chunksize: int = 2_000_000,
) -> None:
    """
    LD decay using only pairs whose two sites fall in the same MAF interval.

    Intervals are left-open, right-closed: (0, 0.01], (0.01, 0.05], (0.05, 0.1],
    (0.1, 0.2], (0.2, 0.5]. Writes ``{prefix}.info.tsv``, ``{prefix}.th.tsv``,
    and ``{prefix}.png`` (mean r2 vs distance, one curve per interval).
    """
    if max_distance_kb <= 0 or bin_size_kb <= 0:
        print("[Error] Invalid LD MAF-bin settings.")
        return
    n_bins = int(np.ceil(float(max_distance_kb) / float(bin_size_kb)))
    if n_bins <= 0:
        print("[Error] Invalid number of bins for LD decay.")
        return

    try:
        maf_by_id = load_maf_by_variant_id(afreq_file)
    except Exception as e:
        print(f"[Error] Failed to read afreq {afreq_file}: {e}")
        return
    print(
        f"[Info] Loaded MAF for {len(maf_by_id)} variants; "
        f"within-bin classes {list(MAF_LD_BIN_LABELS)}"
    )

    sum_r2 = {label: np.zeros(n_bins, dtype=np.float64) for label in MAF_LD_BIN_LABELS}
    cnt_r2 = {label: np.zeros(n_bins, dtype=np.int64) for label in MAF_LD_BIN_LABELS}
    used_counts = {label: 0 for label in MAF_LD_BIN_LABELS}
    empty_median = {
        i: np.array([], dtype=np.float32) for i in range(n_bins)
    }
    total_pairs = 0
    dropped_unmapped = 0
    dropped_cross_bin = 0

    try:
        reader = _iter_plink2_vcor_chunks(
            input_file,
            chunksize,
            usecols=["POS_A", "ID_A", "POS_B", "ID_B", "UNPHASED_R2"],
        )
    except Exception as e:
        print(f"[Error] Failed to read {input_file}: {e}")
        return

    for chunk in reader:
        local = _normalize_ld_chunk(chunk, keep_ids=True)
        if local is None or local.empty:
            continue

        distance_kb = (local["BP_B"] - local["BP_A"]).abs().to_numpy(dtype=np.float64) / 1000.0
        r2 = local["R2"].to_numpy(dtype=np.float64)
        valid_mask = (
            (distance_kb >= 0.0)
            & (distance_kb <= float(max_distance_kb))
            & np.isfinite(r2)
        )
        if not np.any(valid_mask):
            continue

        local = local.loc[valid_mask]
        distance_kb = distance_kb[valid_mask]
        r2 = r2[valid_mask]
        total_pairs += int(len(r2))

        maf_a = _ids_to_values(local["ID_A"], maf_by_id)
        maf_b = _ids_to_values(local["ID_B"], maf_by_id)
        unmapped = ~np.isfinite(maf_a) | ~np.isfinite(maf_b)
        dropped_unmapped += int(np.sum(unmapped))
        masks = within_maf_bin_masks(maf_a, maf_b)
        in_any = np.zeros(len(r2), dtype=bool)
        bin_idx = np.floor(distance_kb / float(bin_size_kb)).astype(np.int64)
        bin_idx = np.clip(bin_idx, 0, n_bins - 1)
        for label, mask in masks.items():
            if not np.any(mask):
                continue
            in_any |= mask
            np.add.at(sum_r2[label], bin_idx[mask], r2[mask])
            np.add.at(cnt_r2[label], bin_idx[mask], 1)
            used_counts[label] += int(np.sum(mask))
        dropped_cross_bin += int(np.sum(~unmapped & ~in_any))

    rows = []
    for label in MAF_LD_BIN_LABELS:
        rows.extend(
            _binned_stats_rows(
                sum_r2[label],
                cnt_r2[label],
                empty_median,
                n_bins,
                float(bin_size_kb),
                label,
            )
        )
    grouped = pd.DataFrame(
        rows,
        columns=[
            "pair_class",
            "bin_start_kb",
            "bin_end_kb",
            "bin_mid_kb",
            "mean_r2",
            "median_r2",
            "n_pairs",
        ],
    )
    save_df_to_tsv(grouped, f"{output_prefix}.info.tsv")

    th_dict: Dict[str, object] = {
        "Total_LD_pairs": int(total_pairs),
        "Pairs_Dropped_Unmapped": int(dropped_unmapped),
        "Pairs_Dropped_Cross_Bin": int(dropped_cross_bin),
        "Max_Distance_kb": float(max_distance_kb),
        "Bin_Size_kb": float(bin_size_kb),
        "N_Variants_With_MAF": int(len(maf_by_id)),
    }
    for label in MAF_LD_BIN_LABELS:
        th_dict[f"Pairs_{label}"] = int(used_counts[label])
        sub = grouped[grouped["pair_class"] == label].sort_values("bin_mid_kb")
        first_mean, half_thr, half_kb = _half_decay_kb(sub)
        th_dict[f"{label}_First_Bin_Mean_R2"] = first_mean
        th_dict[f"{label}_Half_Decay_R2_Threshold"] = half_thr
        th_dict[f"{label}_Half_Decay_kb"] = half_kb
    save_thresholds(th_dict, f"{output_prefix}.th.tsv")
    plot_ld_maf_bin_decay_from_table(
        grouped,
        output_prefix,
        max_distance_kb=float(max_distance_kb),
    )
