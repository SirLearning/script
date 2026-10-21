"""Sample-identity metrics that are not raw PLINK IBS.

Rare-site IBS is uninformative (almost everyone is hom-ref). These helpers
score rare-allele sharing, genotype discordance, a somalier-style
relatedness coefficient, and an extra IBS that keeps only biallelic SNPs
and counts a site when both samples are homozygous.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.sample.sample_utils import load_df_from_king, load_ids_file
from genetics.germplasm.sample.cohort_manifest import load_sample_group_table
from genetics.germplasm.sample.mosdepth_aneuploidy import normalize_sample_id
from genetics.germplasm.sample.ploidy import (
    PLOIDY_ORDER,
    PLOIDY_PALETTE,
    PLOIDY_SUBGENOMES,
    assign_corrected_ploidy,
    load_breachr_segments,
    subgenome_ploidies,
)
from infra.utils.graph import (
    plot_distribution_with_stats,
    plot_faceted_histograms,
    plot_grouped_density_lines,
    plot_scatter_density,
    plot_slope_chart,
)
from infra.utils.io import load_df_from_tsv, save_df_to_tsv, save_thresholds

IDENTITY_PLOT_SPECS = {
    "king": {
        "col": "KINSHIP",
        "title": "KING kinship (pairs above threshold)",
        "x_label": "KINSHIP",
        "xlim": (0.15, 0.55),
        "thresholds": [
            {"value": 0.177, "label": "1st degree 0.177", "color": "#4c72b0", "linestyle": "--"},
            {"value": 0.354, "label": "duplicate 0.354", "color": "#c44e52", "linestyle": "--"},
        ],
    },
    "jaccard": {
        "col": "jaccard",
        "title": "Rare-allele Jaccard (pairs above threshold)",
        "x_label": "Jaccard",
        "xlim": (0.48, 1.02),
        "thresholds": [
            {"value": 0.5, "label": "keep 0.5", "color": "#4c72b0", "linestyle": "--"},
        ],
    },
    "discord": {
        "col": "discord_rate",
        "title": "Genotype discordance (pairs below threshold)",
        "x_label": "Discordance rate",
        "xlim": (0.0, 0.055),
        "thresholds": [
            {"value": 0.05, "label": "keep 0.05", "color": "#4c72b0", "linestyle": "--"},
        ],
    },
    "ibd": {
        "col": "PI_HAT",
        "title": "PLINK IBD PI_HAT (pairs above threshold)",
        "x_label": "PI_HAT",
        "xlim": (0.2, 1.05),
        "thresholds": [
            {"value": 0.25, "label": "keep 0.25", "color": "#4c72b0", "linestyle": "--"},
            {"value": 0.9, "label": "duplicate 0.9", "color": "#c44e52", "linestyle": "--"},
        ],
    },
    "fingerprint": {
        "col": "relatedness",
        "title": "Somalier-style relatedness (pairs above threshold)",
        "x_label": "Relatedness",
        "xlim": (0.35, 1.15),
        "thresholds": [
            {"value": 0.4, "label": "keep 0.4", "color": "#4c72b0", "linestyle": "--"},
            {"value": 0.9, "label": "duplicate 0.9", "color": "#c44e52", "linestyle": "--"},
        ],
    },
    "hom_snp": {
        "col": "hom_ibs",
        "title": "Hom-SNP IBS (biallelic SNPs; both samples homozygous)",
        "x_label": "Hom-SNP IBS",
        "xlim": (0.9, 1.02),
        "thresholds": [
            {"value": 0.99, "label": "keep 0.99", "color": "#4c72b0", "linestyle": "--"},
            {"value": 0.995, "label": "merge 0.995", "color": "#c44e52", "linestyle": "--"},
        ],
    },
}

PAIR_COLUMNS = ["Sample1", "Sample2"]


def _order_pair(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def _empty_pair_frame(extra_cols: list[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=PAIR_COLUMNS + extra_cols)


def jaccard_from_carriers(
    carriers: dict[str, set[str]],
    *,
    min_shared: int = 5,
    min_jaccard: float = 0.5,
) -> pd.DataFrame:
    """Pairwise Jaccard on rare-ALT carrier sets.

    ``carriers`` maps sample ID -> set of variant IDs where the sample carries
    at least one ALT allele. Pairs below ``min_shared`` or ``min_jaccard`` are
    omitted. Shared counts are accumulated from an inverted index so the work
    scales with co-occurrence, not n².
    """
    n_carry = {sample: len(variants) for sample, variants in carriers.items()}
    variant_to_samples: dict[str, list[str]] = defaultdict(list)
    for sample, variants in carriers.items():
        for variant in variants:
            variant_to_samples[variant].append(sample)

    shared: dict[tuple[str, str], int] = defaultdict(int)
    for samples in variant_to_samples.values():
        uniq = sorted(set(samples))
        if len(uniq) < 2:
            continue
        for a, b in combinations(uniq, 2):
            shared[(a, b)] += 1

    rows: list[dict] = []
    for (a, b), n_shared in shared.items():
        if n_shared < min_shared:
            continue
        n_union = n_carry[a] + n_carry[b] - n_shared
        jaccard = n_shared / n_union if n_union else 0.0
        if jaccard < min_jaccard:
            continue
        s1, s2 = _order_pair(a, b)
        rows.append(
            {
                "Sample1": s1,
                "Sample2": s2,
                "n_shared": n_shared,
                "n_union": n_union,
                "jaccard": jaccard,
            }
        )
    if not rows:
        return _empty_pair_frame(["n_shared", "n_union", "jaccard"])
    return pd.DataFrame(rows).sort_values(["jaccard", "n_shared"], ascending=False).reset_index(drop=True)


def discord_rate_pairs(
    ids: list[str],
    gt: np.ndarray,
    *,
    max_discord: float = 0.05,
    min_nobs: int = 100,
) -> pd.DataFrame:
    """Pairwise hardcall discordance. Missing calls are skipped, not counted."""
    if gt.ndim != 2 or gt.shape[0] != len(ids):
        raise ValueError("gt must be (n_samples, n_sites) aligned with ids")
    n_samples = len(ids)
    valid = ~np.isnan(gt)
    rows: list[dict] = []
    for i in range(n_samples - 1):
        both = valid[i] & valid[i + 1 :]
        n_obs = both.sum(axis=1)
        n_diff = np.sum((gt[i] != gt[i + 1 :]) & both, axis=1)
        rate = np.divide(n_diff, n_obs, out=np.ones_like(n_diff, dtype=float), where=n_obs > 0)
        keep = (n_obs >= min_nobs) & (rate <= max_discord)
        for offset in np.flatnonzero(keep):
            j = i + 1 + int(offset)
            s1, s2 = _order_pair(str(ids[i]), str(ids[j]))
            rows.append(
                {
                    "Sample1": s1,
                    "Sample2": s2,
                    "n_obs": int(n_obs[offset]),
                    "n_diff": int(n_diff[offset]),
                    "discord_rate": float(rate[offset]),
                }
            )
    if not rows:
        return _empty_pair_frame(["n_obs", "n_diff", "discord_rate"])
    return pd.DataFrame(rows).sort_values(["discord_rate", "n_diff"]).reset_index(drop=True)


def somalier_relatedness_pairs(
    ids: list[str],
    gt: np.ndarray,
    *,
    min_relatedness: float = 0.4,
    min_nobs: int = 100,
) -> pd.DataFrame:
    """Somalier-style relatedness: ``(shared_hets - 2 * IBS0) / min(hets)``.

    This is the VCF-panel form of somalier, not BAM pileup at a sites VCF.
    """
    if gt.ndim != 2 or gt.shape[0] != len(ids):
        raise ValueError("gt must be (n_samples, n_sites) aligned with ids")
    n_samples = len(ids)
    valid = ~np.isnan(gt)
    hets = np.sum((gt == 1) & valid, axis=1)
    rows: list[dict] = []
    for i in range(n_samples - 1):
        both = valid[i] & valid[i + 1 :]
        n_obs = both.sum(axis=1)
        shared_hets = np.sum((gt[i] == 1) & (gt[i + 1 :] == 1) & both, axis=1)
        ibs0 = np.sum(
            (((gt[i] == 0) & (gt[i + 1 :] == 2)) | ((gt[i] == 2) & (gt[i + 1 :] == 0))) & both,
            axis=1,
        )
        denom = np.minimum(hets[i], hets[i + 1 :])
        rel = np.divide(
            shared_hets - 2.0 * ibs0,
            denom,
            out=np.full(denom.shape, np.nan, dtype=float),
            where=denom > 0,
        )
        keep = (n_obs >= min_nobs) & np.isfinite(rel) & (rel >= min_relatedness)
        for offset in np.flatnonzero(keep):
            j = i + 1 + int(offset)
            s1, s2 = _order_pair(str(ids[i]), str(ids[j]))
            rows.append(
                {
                    "Sample1": s1,
                    "Sample2": s2,
                    "n_obs": int(n_obs[offset]),
                    "shared_hets": int(shared_hets[offset]),
                    "ibs0": int(ibs0[offset]),
                    "relatedness": float(rel[offset]),
                }
            )
    if not rows:
        return _empty_pair_frame(["n_obs", "shared_hets", "ibs0", "relatedness"])
    return (
        pd.DataFrame(rows)
        .sort_values(["relatedness", "shared_hets"], ascending=False)
        .reset_index(drop=True)
    )


def cluster_undirected_pairs(pairs: pd.DataFrame) -> list[list[str]]:
    """Connected components from undirected Sample1/Sample2 edges."""
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    if pairs is None or pairs.empty:
        return []
    for a, b in zip(pairs["Sample1"].astype(str), pairs["Sample2"].astype(str)):
        union(a, b)
    groups: dict[str, list[str]] = defaultdict(list)
    for sample in parent:
        groups[find(sample)].append(sample)
    return [sorted(members) for members in groups.values() if len(members) >= 2]


def consensus_pairs(
    method_tables: dict[str, pd.DataFrame],
    *,
    min_methods: int = 2,
) -> pd.DataFrame:
    """Pairs that appear in at least ``min_methods`` identity tables."""
    hits: dict[tuple[str, str], list[str]] = defaultdict(list)
    for method, table in method_tables.items():
        if table is None or table.empty:
            continue
        for a, b in zip(table["Sample1"].astype(str), table["Sample2"].astype(str)):
            hits[_order_pair(a, b)].append(method)
    rows = []
    for (s1, s2), methods in hits.items():
        uniq = sorted(set(methods))
        if len(uniq) < min_methods:
            continue
        rows.append(
            {
                "Sample1": s1,
                "Sample2": s2,
                "n_methods": len(uniq),
                "methods": ",".join(uniq),
            }
        )
    if not rows:
        return _empty_pair_frame(["n_methods", "methods"])
    return pd.DataFrame(rows).sort_values(["n_methods", "Sample1"], ascending=[False, True]).reset_index(drop=True)


def clusters_to_frame(clusters: list[list[str]]) -> pd.DataFrame:
    rows = [
        {
            "cluster_id": idx,
            "n_samples": len(members),
            "samples": ",".join(members),
        }
        for idx, members in enumerate(clusters, start=1)
    ]
    if not rows:
        return pd.DataFrame(columns=["cluster_id", "n_samples", "samples"])
    return pd.DataFrame(rows)


def carriers_from_vcf(vcf_path: str, *, max_carriers: int = 400) -> dict[str, set[str]]:
    """Collect ALT carriers from a VCF. Sites with too many carriers are skipped."""
    carriers: dict[str, set[str]] = defaultdict(set)
    path = Path(vcf_path)
    if not path.exists() or path.stat().st_size == 0:
        return {}
    samples: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                for sample in samples:
                    carriers.setdefault(sample, set())
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or not samples:
                continue
            chrom, pos, vid = fields[0], fields[1], fields[2]
            variant = vid if vid not in {".", ""} else f"{chrom}:{pos}"
            fmt = fields[8].split(":")
            try:
                gt_idx = fmt.index("GT")
            except ValueError:
                continue
            hit: list[str] = []
            for sample, cell in zip(samples, fields[9:]):
                gt = cell.split(":")[gt_idx] if cell else "./."
                alleles = gt.replace("|", "/").split("/")
                if any(a not in {".", "0", ""} for a in alleles):
                    hit.append(sample)
            if not hit or len(hit) > max_carriers:
                continue
            for sample in hit:
                carriers[sample].add(variant)
    return dict(carriers)


def load_plink_export_a(raw_path: str) -> tuple[list[str], np.ndarray]:
    """Load PLINK ``--export A`` additive dosages as (ids, gt). Missing → NaN."""
    df = pd.read_csv(raw_path, sep=r"\s+", engine="python")
    id_col = "IID" if "IID" in df.columns else df.columns[1]
    skip = {c for c in ("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE") if c in df.columns}
    dose = df.drop(columns=list(skip))
    gt = dose.to_numpy(dtype=float)
    gt[gt < 0] = np.nan
    return df[id_col].astype(str).tolist(), gt


_SNP_NT = frozenset("ACGT")


def is_snp_alleles(a1: str, a2: str) -> bool:
    """True when both alleles are single A/C/G/T bases (drop indels)."""
    s1 = str(a1).strip().upper()
    s2 = str(a2).strip().upper()
    return len(s1) == 1 and len(s2) == 1 and s1 in _SNP_NT and s2 in _SNP_NT


def hom_snp_ibs_from_gt(
    ids: list[str],
    gt: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Pairwise IBS on sites where both samples are homozygous (0 or 2).

    Hets and missing calls are skipped, not scored. Returns square
    ``ibs``, ``n_compared``, and ``n_match`` matrices aligned with ``ids``.
    """
    if gt.ndim != 2 or gt.shape[0] != len(ids):
        raise ValueError("gt must be (n_samples, n_sites) aligned with ids")
    n_samples = len(ids)
    h0 = gt == 0
    h2 = gt == 2
    hom = h0 | h2
    n_match = (h0.astype(np.int64) @ h0.T) + (h2.astype(np.int64) @ h2.T)
    n_comp = hom.astype(np.int64) @ hom.T
    ibs = np.divide(
        n_match,
        n_comp,
        out=np.full((n_samples, n_samples), np.nan, dtype=float),
        where=n_comp > 0,
    )
    np.fill_diagonal(ibs, 1.0)
    return ibs.astype(float), n_comp.astype(float), n_match.astype(float)


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


def _plink1_path(bed_prefix: str, suffix: str) -> Path:
    """Join a PLINK1 prefix with ``.bed`` / ``.bim`` / ``.fam``.

    ``Path.with_suffix`` is wrong for prefixes such as ``A_test.plink``
    (it would replace ``.plink`` instead of appending).
    """
    return Path(f"{bed_prefix}{suffix}")


def _fam_iids(bed_prefix: str) -> list[str]:
    fam = pd.read_csv(
        _plink1_path(bed_prefix, ".fam"),
        sep=r"\s+",
        engine="python",
        header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"],
        dtype={"IID": str},
    )
    return fam["IID"].astype(str).tolist()


def _iter_bed_snp_dosages(bed_prefix: str, *, snp_only: bool = True):
    """Yield per-variant dosage vectors from a SNP-major PLINK1 bed."""
    bim = pd.read_csv(
        _plink1_path(bed_prefix, ".bim"),
        sep=r"\s+",
        engine="python",
        header=None,
        names=["chrom", "id", "cm", "pos", "a1", "a2"],
        dtype={"a1": str, "a2": str},
    )
    ids = _fam_iids(bed_prefix)
    n_samples = len(ids)
    n_var = len(bim)
    n_bytes = (n_samples + 3) // 4
    bed_path = _plink1_path(bed_prefix, ".bed")
    raw = np.memmap(bed_path, dtype=np.uint8, mode="r")
    if raw.size < 3 + n_bytes * n_var:
        raise ValueError(f"{bed_path} is shorter than expected for {n_var} SNPs")
    if raw[0] != 0x6C or raw[1] != 0x1B:
        raise ValueError(f"{bed_path} is not a PLINK bed file")
    if int(raw[2]) != 1:
        raise ValueError(f"{bed_path} is not SNP-major")
    body = raw[3:]
    a1 = bim["a1"].astype(str)
    a2 = bim["a2"].astype(str)
    keep = np.ones(n_var, dtype=bool)
    if snp_only:
        keep = np.array([is_snp_alleles(x, y) for x, y in zip(a1, a2)], dtype=bool)
    for i in np.flatnonzero(keep):
        offset = int(i) * n_bytes
        yield _decode_bed_snp_bytes(body[offset : offset + n_bytes], n_samples)


def hom_snp_ibs_from_bed(
    bed_prefix: str,
    *,
    chunk_sites: int = 4096,
    snp_only: bool = True,
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    """Stream a PLINK1 bed and accumulate hom-SNP IBS via chunked GEMM."""
    ids = _fam_iids(bed_prefix)
    n_samples = len(ids)
    n_match = np.zeros((n_samples, n_samples), dtype=np.float64)
    n_comp = np.zeros((n_samples, n_samples), dtype=np.float64)
    buf0: list[np.ndarray] = []
    buf2: list[np.ndarray] = []
    chunk_sites = max(1, int(chunk_sites))

    def flush() -> None:
        nonlocal buf0, buf2, n_match, n_comp
        if not buf0:
            return
        h0 = np.ascontiguousarray(np.stack(buf0, axis=1), dtype=np.float32)
        h2 = np.ascontiguousarray(np.stack(buf2, axis=1), dtype=np.float32)
        n_match += np.asarray(h0 @ h0.T, dtype=np.float64)
        n_match += np.asarray(h2 @ h2.T, dtype=np.float64)
        hom = np.ascontiguousarray(h0 + h2, dtype=np.float32)
        n_comp += np.asarray(hom @ hom.T, dtype=np.float64)
        buf0 = []
        buf2 = []

    for dosage in _iter_bed_snp_dosages(bed_prefix, snp_only=snp_only):
        buf0.append(dosage == 0)
        buf2.append(dosage == 2)
        if len(buf0) >= chunk_sites:
            flush()
    flush()
    ibs = np.divide(
        n_match,
        n_comp,
        out=np.full((n_samples, n_samples), np.nan, dtype=float),
        where=n_comp > 0,
    )
    np.fill_diagonal(ibs, 1.0)
    return ids, ibs, n_comp, n_match


def pairs_from_square_ibs(
    ids: list[str],
    ibs: np.ndarray,
    n_comp: np.ndarray,
    n_match: np.ndarray,
    *,
    min_ibs: float = 0.99,
    min_nobs: int = 100,
) -> pd.DataFrame:
    """Upper-triangle pairs with ``hom_ibs >= min_ibs``."""
    n_samples = len(ids)
    iu = np.triu_indices(n_samples, k=1)
    scores = ibs[iu]
    comps = n_comp[iu]
    matches = n_match[iu]
    keep = np.isfinite(scores) & (scores >= min_ibs) & (comps >= min_nobs)
    rows: list[dict] = []
    for i, j, score, comp, match in zip(
        iu[0][keep], iu[1][keep], scores[keep], comps[keep], matches[keep]
    ):
        s1, s2 = _order_pair(str(ids[i]), str(ids[j]))
        rows.append(
            {
                "Sample1": s1,
                "Sample2": s2,
                "n_hom": int(comp),
                "n_match": int(match),
                "hom_ibs": float(score),
            }
        )
    if not rows:
        return _empty_pair_frame(["n_hom", "n_match", "hom_ibs"])
    return (
        pd.DataFrame(rows)
        .sort_values(["hom_ibs", "n_match"], ascending=False)
        .reset_index(drop=True)
    )


def compare_square_ibs(
    ids: list[str],
    plink_ibs: np.ndarray,
    hom_ibs: np.ndarray,
    n_comp: np.ndarray,
    *,
    min_either: float = 0.95,
    min_nobs: int = 1,
    max_pairs: int = 200_000,
) -> pd.DataFrame:
    """Pairs high on PLINK IBS or hom-SNP IBS, with both scores side by side."""
    n_samples = len(ids)
    iu = np.triu_indices(n_samples, k=1)
    plink = plink_ibs[iu]
    hom = hom_ibs[iu]
    comps = n_comp[iu]
    keep = (
        ((plink >= min_either) | (hom >= min_either))
        & np.isfinite(plink)
        & np.isfinite(hom)
        & (comps >= min_nobs)
    )
    idx = np.flatnonzero(keep)
    if idx.size > max_pairs:
        peak = np.maximum(plink[keep], hom[keep])
        idx = idx[np.argpartition(peak, -max_pairs)[-max_pairs:]]
    rows: list[dict] = []
    for pos in idx:
        i = int(iu[0][pos])
        j = int(iu[1][pos])
        s1, s2 = _order_pair(str(ids[i]), str(ids[j]))
        rows.append(
            {
                "Sample1": s1,
                "Sample2": s2,
                "plink_ibs": float(plink[pos]),
                "hom_ibs": float(hom[pos]),
                "n_hom": int(comps[pos]),
                "ibs_delta": float(hom[pos] - plink[pos]),
            }
        )
    if not rows:
        return _empty_pair_frame(["plink_ibs", "hom_ibs", "n_hom", "ibs_delta"])
    return (
        pd.DataFrame(rows)
        .sort_values(["hom_ibs", "plink_ibs"], ascending=False)
        .reset_index(drop=True)
    )


def pairs_hom_ibs_one(
    ids: list[str],
    hom_ibs: np.ndarray,
    plink_ibs: np.ndarray,
    *,
    atol: float = 5e-7,
) -> pd.DataFrame:
    """Upper-triangle pairs whose hom-SNP IBS is 1, with matching PLINK IBS."""
    n_samples = len(ids)
    extra = ["hom_ibs", "plink_ibs", "ibs_delta"]
    if n_samples < 2:
        return _empty_pair_frame(extra)
    iu = np.triu_indices(n_samples, k=1)
    hom = np.asarray(hom_ibs, dtype=float)[iu]
    plink = np.asarray(plink_ibs, dtype=float)[iu]
    keep = np.isfinite(hom) & np.isfinite(plink) & np.isclose(hom, 1.0, rtol=0.0, atol=atol)
    idx = np.flatnonzero(keep)
    if idx.size == 0:
        return _empty_pair_frame(extra)
    names = np.asarray(ids, dtype=object)
    left = names[iu[0][idx]]
    right = names[iu[1][idx]]
    swap = left > right
    left[swap], right[swap] = right[swap], left[swap]
    return (
        pd.DataFrame(
            {
                "Sample1": left,
                "Sample2": right,
                "hom_ibs": hom[idx],
                "plink_ibs": plink[idx],
                "ibs_delta": hom[idx] - plink[idx],
            }
        )
        .sort_values(["plink_ibs", "Sample1", "Sample2"], ascending=[True, True, True])
        .reset_index(drop=True)
    )


def sample_square_ibs_pairs(
    ids: list[str],
    plink_ibs: np.ndarray,
    hom_ibs: np.ndarray,
    n_comp: np.ndarray | None = None,
    *,
    n_pairs: int = 200_000,
    seed: int = 1,
) -> pd.DataFrame:
    """Random upper-triangle pairs from two square IBS matrices. No score cutoff."""
    n_samples = len(ids)
    if n_samples < 2:
        return _empty_pair_frame(["plink_ibs", "hom_ibs", "n_hom", "ibs_delta"])
    iu = np.triu_indices(n_samples, k=1)
    plink = np.asarray(plink_ibs, dtype=float)[iu]
    hom = np.asarray(hom_ibs, dtype=float)[iu]
    finite = np.isfinite(plink) & np.isfinite(hom)
    idx = np.flatnonzero(finite)
    n_keep = min(int(n_pairs), int(idx.size))
    if n_keep <= 0:
        return _empty_pair_frame(["plink_ibs", "hom_ibs", "n_hom", "ibs_delta"])
    if idx.size > n_keep:
        rng = np.random.default_rng(int(seed))
        idx = rng.choice(idx, size=n_keep, replace=False)
    names = np.asarray(ids, dtype=object)
    left = names[iu[0][idx]]
    right = names[iu[1][idx]]
    swap = left > right
    left[swap], right[swap] = right[swap], left[swap]
    out = pd.DataFrame(
        {
            "Sample1": left,
            "Sample2": right,
            "plink_ibs": plink[idx],
            "hom_ibs": hom[idx],
            "ibs_delta": hom[idx] - plink[idx],
        }
    )
    if n_comp is not None:
        out["n_hom"] = np.asarray(n_comp, dtype=float)[iu][idx]
    return out.reset_index(drop=True)


def histogram_peak_centers(
    values: np.ndarray,
    *,
    n_bins: int = 80,
    min_rel_height: float = 0.05,
    xlim: tuple[float, float] | None = None,
) -> list[float]:
    """Bin centers of histogram local maxima at least ``min_rel_height`` of the tallest bin."""
    scores = np.asarray(values, dtype=float)
    scores = scores[np.isfinite(scores)]
    if scores.size < 3:
        return []
    lo = float(xlim[0]) if xlim is not None else float(scores.min())
    hi = float(xlim[1]) if xlim is not None else float(scores.max())
    if hi <= lo:
        return [lo]
    n_bins = max(8, int(n_bins))
    counts, edges = np.histogram(scores, bins=n_bins, range=(lo, hi))
    peak_h = float(counts.max())
    if peak_h <= 0:
        return []
    floor = float(min_rel_height) * peak_h
    peaks: list[float] = []
    for i in range(1, len(counts) - 1):
        if counts[i] < floor:
            continue
        if counts[i] >= counts[i - 1] and counts[i] >= counts[i + 1]:
            if counts[i] > counts[i - 1] or counts[i] > counts[i + 1]:
                peaks.append(float(0.5 * (edges[i] + edges[i + 1])))
    return peaks


def _xlim_for_ibs_scores(scores: pd.Series) -> tuple[float, float]:
    finite = pd.to_numeric(scores, errors="coerce").dropna()
    if finite.empty:
        return (0.5, 1.02)
    lo = float(finite.min())
    hi = float(finite.max())
    pad = max(0.01, (hi - lo) * 0.05) if hi > lo else 0.02
    return (max(0.0, lo - pad), min(1.05, hi + pad))


def ana_ibs_sampled_pair_dist(
    pairs_path: str,
    output_prefix: str,
    *,
    n_pairs: int | None = None,
    seed: int = 1,
) -> pd.DataFrame:
    """Full-range PLINK and hom-SNP IBS histograms from a proportional pair sample.

    The input table is already a uniform draw from the upper triangle (no 0.99
    cutoff). This does not zoom to the high-IBS tail.
    """
    df = load_df_from_tsv(pairs_path)
    if df is None or df.empty or "plink_ibs" not in df.columns or "hom_ibs" not in df.columns:
        empty = _empty_pair_frame(["plink_ibs", "hom_ibs"])
        save_thresholds({"n_pairs": 0}, f"{output_prefix}.th.tsv")
        return empty
    work = df.dropna(subset=["plink_ibs", "hom_ibs"]).copy()
    if n_pairs is not None and len(work) > int(n_pairs):
        work = work.sample(n=int(n_pairs), random_state=int(seed)).reset_index(drop=True)
    plink_xlim = _xlim_for_ibs_scores(work["plink_ibs"])
    hom_xlim = _xlim_for_ibs_scores(work["hom_ibs"])
    plink_peaks = histogram_peak_centers(work["plink_ibs"].to_numpy(), xlim=plink_xlim)
    hom_peaks = histogram_peak_centers(work["hom_ibs"].to_numpy(), xlim=hom_xlim)
    save_thresholds(
        {
            "n_pairs": int(len(work)),
            "min_plink_ibs": float(work["plink_ibs"].min()),
            "median_plink_ibs": float(work["plink_ibs"].median()),
            "max_plink_ibs": float(work["plink_ibs"].max()),
            "min_hom_ibs": float(work["hom_ibs"].min()),
            "median_hom_ibs": float(work["hom_ibs"].median()),
            "max_hom_ibs": float(work["hom_ibs"].max()),
            "frac_plink_ge_0.99": float((work["plink_ibs"] >= 0.99).mean()),
            "frac_hom_ge_0.99": float((work["hom_ibs"] >= 0.99).mean()),
            "n_plink_peaks": int(len(plink_peaks)),
            "plink_peak_centers": ",".join(f"{p:.4f}" for p in plink_peaks),
            "n_hom_peaks": int(len(hom_peaks)),
            "hom_peak_centers": ",".join(f"{p:.4f}" for p in hom_peaks),
        },
        f"{output_prefix}.th.tsv",
    )
    cut = [{"value": 0.99, "label": "0.99", "color": "#c44e52", "linestyle": "--"}]
    specs = (
        ("plink_ibs", "PLINK IBS (sampled pairs)", "PLINK IBS", "plink", plink_xlim),
        ("hom_ibs", "Hom-SNP IBS (sampled pairs)", "Hom-SNP IBS", "hom", hom_xlim),
    )
    for col, title, x_label, tag, xlim in specs:
        plot_distribution_with_stats(
            data=work,
            col=col,
            title=title,
            filename=f"{output_prefix}.{tag}.dist.png",
            x_label=x_label,
            y_label="Count of pairs",
            mean_val=float(work[col].mean()),
            median_val=float(work[col].median()),
            xlim=xlim,
            thresholds=cut,
        )
        plot_distribution_with_stats(
            data=work,
            col=col,
            title=f"{title} (log count)",
            filename=f"{output_prefix}.{tag}.dist_log.png",
            x_label=x_label,
            y_label="Count of pairs",
            mean_val=float(work[col].mean()),
            median_val=float(work[col].median()),
            xlim=xlim,
            log_scale=True,
            thresholds=cut,
        )
    return work


def sample_square_ibs_pairs_by_label(
    ids: list[str],
    plink_ibs: np.ndarray,
    hom_ibs: np.ndarray,
    labels: dict[str, str],
    *,
    allowed_labels: tuple[str, ...] | None = None,
    n_pairs_per_class: int = 200_000,
    seed: int = 1,
) -> pd.DataFrame:
    """Upper-triangle IBS pairs whose two samples share a label (ploidy class)."""
    extra = ["plink_ibs", "hom_ibs", "ibs_delta", "ploidy"]
    if len(ids) < 2:
        return _empty_pair_frame(extra)
    pos = {sid: i for i, sid in enumerate(ids)}
    order = allowed_labels if allowed_labels is not None else tuple(
        dict.fromkeys(labels[sid] for sid in ids if sid in labels)
    )
    rng = np.random.default_rng(int(seed))
    frames: list[pd.DataFrame] = []
    names = np.asarray(ids, dtype=object)
    plink = np.asarray(plink_ibs, dtype=float)
    hom = np.asarray(hom_ibs, dtype=float)
    for label in order:
        members = [sid for sid in ids if labels.get(sid) == label]
        idx = np.array([pos[sid] for sid in members if sid in pos], dtype=int)
        if idx.size < 2:
            continue
        iu = np.triu_indices(int(idx.size), k=1)
        n_total = int(iu[0].size)
        n_keep = min(int(n_pairs_per_class), n_total)
        if n_keep <= 0:
            continue
        pick = np.arange(n_total)
        if n_total > n_keep:
            pick = rng.choice(n_total, size=n_keep, replace=False)
        left_i = idx[iu[0][pick]]
        right_i = idx[iu[1][pick]]
        left = names[left_i].copy()
        right = names[right_i].copy()
        swap = left > right
        left[swap], right[swap] = right[swap], left[swap]
        plink_vals = plink[left_i, right_i]
        hom_vals = hom[left_i, right_i]
        finite = np.isfinite(plink_vals) & np.isfinite(hom_vals)
        if not np.any(finite):
            continue
        frames.append(
            pd.DataFrame(
                {
                    "Sample1": left[finite],
                    "Sample2": right[finite],
                    "plink_ibs": plink_vals[finite],
                    "hom_ibs": hom_vals[finite],
                    "ibs_delta": hom_vals[finite] - plink_vals[finite],
                    "ploidy": label,
                }
            )
        )
    if not frames:
        return _empty_pair_frame(extra)
    return pd.concat(frames, ignore_index=True)


def _ploidy_labels_for_ids(ids: list[str], ploidy_table: pd.DataFrame) -> dict[str, str]:
    by_exact: dict[str, str] = {}
    by_norm: dict[str, str] = {}
    for row in ploidy_table.itertuples(index=False):
        ploidy = str(row.ploidy)
        sample = str(row.Sample)
        by_exact[sample] = ploidy
        by_norm[normalize_sample_id(sample)] = ploidy
        depth = getattr(row, "sample_breachr", None)
        if depth is not None and not (isinstance(depth, float) and pd.isna(depth)):
            depth_s = str(depth)
            if depth_s and depth_s.lower() != "nan":
                by_exact[depth_s] = ploidy
                by_norm[normalize_sample_id(depth_s)] = ploidy
    out: dict[str, str] = {}
    for sid in ids:
        if sid in by_exact:
            out[sid] = by_exact[sid]
        else:
            mapped = by_norm.get(normalize_sample_id(sid))
            if mapped is not None:
                out[sid] = mapped
    return out


def ana_ibs_ploidy_pair_dist(
    plink_mibs: str,
    plink_mibs_id: str,
    hom_mibs: str,
    hom_mibs_id: str,
    group_file: str,
    breachr_segments_path: str,
    output_prefix: str,
    *,
    subgenome: str,
    n_pairs_per_class: int = 200_000,
    seed: int = 1,
) -> pd.DataFrame:
    """Within-ploidy PLINK and hom-SNP IBS histograms for one subgenome matrix."""
    ids, plink = load_square_mibs(plink_mibs, plink_mibs_id)
    hom_ids, hom = load_square_mibs(hom_mibs, hom_mibs_id)
    ids, plink, hom = _align_square_by_ids(ids, plink, hom_ids, hom)
    groups = load_sample_group_table(group_file)
    segments = load_breachr_segments(breachr_segments_path)
    ploidy_table = assign_corrected_ploidy(groups, segments)
    labels = _ploidy_labels_for_ids(ids, ploidy_table)
    allowed = subgenome_ploidies(subgenome)
    pairs = sample_square_ibs_pairs_by_label(
        ids,
        plink,
        hom,
        labels,
        allowed_labels=allowed,
        n_pairs_per_class=n_pairs_per_class,
        seed=seed,
    )
    sample_rows = []
    ploidy_lookup = ploidy_table.drop_duplicates(subset=["Sample"], keep="first").set_index("Sample")
    for sid in ids:
        ploidy = labels.get(sid)
        if ploidy is None or ploidy not in allowed:
            continue
        meta = ploidy_lookup.loc[sid] if sid in ploidy_lookup.index else None
        if meta is None:
            norm = normalize_sample_id(sid)
            hits = ploidy_table[ploidy_table["Sample"].map(normalize_sample_id) == norm]
            meta = hits.iloc[0] if not hits.empty else None
        sample_rows.append(
            {
                "Sample": sid,
                "ploidy": ploidy,
                "Group": "" if meta is None else str(meta["Group"]),
                "nominal_ploidy": "" if meta is None else str(meta["nominal_ploidy"]),
                "reason": "" if meta is None else str(meta["reason"]),
            }
        )
    samples = pd.DataFrame(sample_rows)
    save_df_to_tsv(samples, f"{output_prefix}.samples.tsv")
    save_df_to_tsv(pairs, f"{output_prefix}.pairs.tsv")

    class_rows: list[dict[str, object]] = []
    present = [p for p in allowed if p in set(pairs["ploidy"])] if not pairs.empty else []
    for ploidy in allowed:
        sub = pairs[pairs["ploidy"] == ploidy] if not pairs.empty else pairs
        n_samples = int((samples["ploidy"] == ploidy).sum()) if not samples.empty else 0
        class_rows.append(
            {
                "ploidy": ploidy,
                "n_samples": n_samples,
                "n_pairs": int(len(sub)),
                "median_plink_ibs": float(sub["plink_ibs"].median()) if len(sub) else float("nan"),
                "median_hom_ibs": float(sub["hom_ibs"].median()) if len(sub) else float("nan"),
                "min_plink_ibs": float(sub["plink_ibs"].min()) if len(sub) else float("nan"),
                "max_plink_ibs": float(sub["plink_ibs"].max()) if len(sub) else float("nan"),
                "plink_peak_centers": ",".join(
                    f"{p:.4f}"
                    for p in histogram_peak_centers(
                        sub["plink_ibs"].to_numpy() if len(sub) else np.array([])
                    )
                ),
                "hom_peak_centers": ",".join(
                    f"{p:.4f}"
                    for p in histogram_peak_centers(
                        sub["hom_ibs"].to_numpy() if len(sub) else np.array([])
                    )
                ),
            }
        )
    class_df = pd.DataFrame(class_rows)
    save_df_to_tsv(class_df, f"{output_prefix}.info.tsv")
    save_thresholds(
        {
            "subgenome": subgenome,
            "n_samples": int(len(samples)),
            "n_pairs": int(len(pairs)),
            "n_ploidy_classes": int(len(present)),
            "ploidy_classes": ",".join(present),
            "n_reclassified": int((samples["reason"] == "breachr_ploidy_mismatch").sum())
            if not samples.empty and "reason" in samples.columns
            else 0,
        },
        f"{output_prefix}.th.tsv",
    )
    if pairs.empty or not present:
        return pairs

    facet_labels = {}
    n_by_ploidy = samples.groupby("ploidy").size() if not samples.empty else pd.Series(dtype=int)
    for ploidy in present:
        n_s = int(n_by_ploidy.get(ploidy, 0))
        n_p = int((pairs["ploidy"] == ploidy).sum())
        facet_labels[ploidy] = f"{ploidy} (n={n_s}, pairs={n_p})"
    cut = [{"value": 0.99, "label": "0.99", "color": "#c44e52", "linestyle": "--"}]
    specs = (
        ("plink_ibs", "PLINK IBS", "plink"),
        ("hom_ibs", "Hom-SNP IBS", "hom"),
    )
    for col, x_label, tag in specs:
        xlim = _xlim_for_ibs_scores(pairs[col])
        plot_faceted_histograms(
            data=pairs,
            x_col=col,
            facet_col="ploidy",
            facet_order=present,
            title=f"{x_label} within ploidy ({subgenome} subgenome)",
            filename=f"{output_prefix}.{tag}.dist.png",
            xlabel=x_label,
            ylabel="Count of pairs",
            xlim=xlim,
            facet_colors=PLOIDY_PALETTE,
            facet_labels=facet_labels,
            sharey=False,
            thresholds=cut,
        )
        plot_faceted_histograms(
            data=pairs,
            x_col=col,
            facet_col="ploidy",
            facet_order=present,
            title=f"{x_label} within ploidy ({subgenome} subgenome, log count)",
            filename=f"{output_prefix}.{tag}.dist_log.png",
            xlabel=x_label,
            ylabel="Count of pairs",
            log_scale=True,
            xlim=xlim,
            facet_colors=PLOIDY_PALETTE,
            facet_labels=facet_labels,
            sharey=False,
            thresholds=cut,
        )
        plot_grouped_density_lines(
            data=pairs,
            x_col=col,
            hue_col="ploidy",
            hue_order=present,
            palette_map=PLOIDY_PALETTE,
            title=f"{x_label} within ploidy ({subgenome} subgenome)",
            filename=f"{output_prefix}.{tag}.density.png",
            xlabel=x_label,
            ylabel="Density",
            xlim=xlim,
            hue_labels=facet_labels,
        )
    return pairs


def _norm_index(ids: list[str]) -> dict[str, int]:
    pos: dict[str, int] = {}
    for i, sid in enumerate(ids):
        key = normalize_sample_id(str(sid))
        if key not in pos:
            pos[key] = i
    return pos


def sample_combined_subgenome_ibs_pairs(
    matrices: dict[str, tuple[list[str], np.ndarray, np.ndarray]],
    labels: dict[str, str],
    *,
    n_pairs_per_class: int = 200_000,
    seed: int = 1,
) -> pd.DataFrame:
    """Within-ploidy pairs scored on every subgenome that class carries.

    AA uses A, SS uses B, DD uses D, AABB means A and B, AABBDD means A+B+D.
    Combined PLINK / hom-SNP IBS is the unweighted mean of those matrices.
    """
    extra = [
        "plink_ibs",
        "hom_ibs",
        "ibs_delta",
        "ploidy",
        "subgenomes",
        "plink_ibs_A",
        "plink_ibs_B",
        "plink_ibs_D",
        "hom_ibs_A",
        "hom_ibs_B",
        "hom_ibs_D",
    ]
    loaded: dict[str, tuple[list[str], np.ndarray, np.ndarray, dict[str, int]]] = {}
    for sg, payload in matrices.items():
        ids, plink, hom = payload
        loaded[str(sg).upper()] = (
            list(ids),
            np.asarray(plink, dtype=float),
            np.asarray(hom, dtype=float),
            _norm_index(list(ids)),
        )
    labels_norm: dict[str, str] = {}
    display: dict[str, str] = {}
    for sid, ploidy in labels.items():
        if not ploidy:
            continue
        key = normalize_sample_id(str(sid))
        if key not in labels_norm:
            labels_norm[key] = str(ploidy)
            display[key] = str(sid)
    rng = np.random.default_rng(int(seed))
    frames: list[pd.DataFrame] = []
    for ploidy in PLOIDY_ORDER:
        required = tuple(sg for sg in ("A", "B", "D") if sg in PLOIDY_SUBGENOMES[ploidy])
        if any(sg not in loaded for sg in required):
            continue
        members = [
            key
            for key, label in labels_norm.items()
            if label == ploidy and all(key in loaded[sg][3] for sg in required)
        ]
        members.sort()
        n_mem = len(members)
        if n_mem < 2:
            continue
        iu = np.triu_indices(n_mem, k=1)
        n_total = int(iu[0].size)
        n_keep = min(int(n_pairs_per_class), n_total)
        if n_keep <= 0:
            continue
        pick = np.arange(n_total)
        if n_total > n_keep:
            pick = rng.choice(n_total, size=n_keep, replace=False)
        left_keys = [members[i] for i in iu[0][pick]]
        right_keys = [members[j] for j in iu[1][pick]]
        plink_stack = []
        hom_stack = []
        per_plink: dict[str, np.ndarray] = {}
        per_hom: dict[str, np.ndarray] = {}
        finite = np.ones(len(left_keys), dtype=bool)
        for sg in required:
            _ids, plink, hom, pos = loaded[sg]
            i = np.array([pos[k] for k in left_keys], dtype=int)
            j = np.array([pos[k] for k in right_keys], dtype=int)
            pv = plink[i, j]
            hv = hom[i, j]
            per_plink[sg] = pv
            per_hom[sg] = hv
            plink_stack.append(pv)
            hom_stack.append(hv)
            finite &= np.isfinite(pv) & np.isfinite(hv)
        if not np.any(finite):
            continue
        plink_mean = np.mean(np.vstack(plink_stack), axis=0)
        hom_mean = np.mean(np.vstack(hom_stack), axis=0)
        names1 = np.array([display[k] for k in left_keys], dtype=object)
        names2 = np.array([display[k] for k in right_keys], dtype=object)
        swap = names1 > names2
        names1[swap], names2[swap] = names2[swap], names1[swap]
        frame = {
            "Sample1": names1[finite],
            "Sample2": names2[finite],
            "plink_ibs": plink_mean[finite],
            "hom_ibs": hom_mean[finite],
            "ibs_delta": (hom_mean - plink_mean)[finite],
            "ploidy": ploidy,
            "subgenomes": "".join(required),
            "plink_ibs_A": per_plink["A"][finite] if "A" in per_plink else np.full(int(finite.sum()), np.nan),
            "plink_ibs_B": per_plink["B"][finite] if "B" in per_plink else np.full(int(finite.sum()), np.nan),
            "plink_ibs_D": per_plink["D"][finite] if "D" in per_plink else np.full(int(finite.sum()), np.nan),
            "hom_ibs_A": per_hom["A"][finite] if "A" in per_hom else np.full(int(finite.sum()), np.nan),
            "hom_ibs_B": per_hom["B"][finite] if "B" in per_hom else np.full(int(finite.sum()), np.nan),
            "hom_ibs_D": per_hom["D"][finite] if "D" in per_hom else np.full(int(finite.sum()), np.nan),
        }
        frames.append(pd.DataFrame(frame))
    if not frames:
        return _empty_pair_frame(extra)
    return pd.concat(frames, ignore_index=True)


def ana_ibs_ploidy_group_pair_dist(
    a_plink_mibs: str,
    a_plink_mibs_id: str,
    a_hom_mibs: str,
    a_hom_mibs_id: str,
    b_plink_mibs: str,
    b_plink_mibs_id: str,
    b_hom_mibs: str,
    b_hom_mibs_id: str,
    d_plink_mibs: str,
    d_plink_mibs_id: str,
    d_hom_mibs: str,
    d_hom_mibs_id: str,
    group_file: str,
    breachr_segments_path: str,
    output_prefix: str,
    *,
    n_pairs_per_class: int = 200_000,
    seed: int = 1,
) -> pd.DataFrame:
    """Within-ploidy PLINK and hom-SNP IBS using each class's own subgenomes."""
    matrices: dict[str, tuple[list[str], np.ndarray, np.ndarray]] = {}
    path_sets = {
        "A": (a_plink_mibs, a_plink_mibs_id, a_hom_mibs, a_hom_mibs_id),
        "B": (b_plink_mibs, b_plink_mibs_id, b_hom_mibs, b_hom_mibs_id),
        "D": (d_plink_mibs, d_plink_mibs_id, d_hom_mibs, d_hom_mibs_id),
    }
    for sg, (plink_path, plink_id, hom_path, hom_id) in path_sets.items():
        ids, plink = load_square_mibs(plink_path, plink_id)
        hom_ids, hom = load_square_mibs(hom_path, hom_id)
        ids, plink, hom = _align_square_by_ids(ids, plink, hom_ids, hom)
        matrices[sg] = (ids, plink, hom)
    groups = load_sample_group_table(group_file)
    segments = load_breachr_segments(breachr_segments_path)
    ploidy_table = assign_corrected_ploidy(groups, segments)
    all_ids: list[str] = []
    seen: set[str] = set()
    for sg in ("A", "B", "D"):
        for sid in matrices[sg][0]:
            if sid not in seen:
                seen.add(sid)
                all_ids.append(sid)
    labels = _ploidy_labels_for_ids(all_ids, ploidy_table)
    pairs = sample_combined_subgenome_ibs_pairs(
        matrices,
        labels,
        n_pairs_per_class=n_pairs_per_class,
        seed=seed,
    )
    sg_pos = {sg: _norm_index(matrices[sg][0]) for sg in ("A", "B", "D")}
    sample_rows = []
    ploidy_lookup = ploidy_table.drop_duplicates(subset=["Sample"], keep="first").set_index("Sample")
    for sid in all_ids:
        ploidy = labels.get(sid)
        if ploidy is None or ploidy not in PLOIDY_SUBGENOMES:
            continue
        required = PLOIDY_SUBGENOMES[ploidy]
        key = normalize_sample_id(sid)
        if not all(key in sg_pos[sg] for sg in required):
            continue
        meta = ploidy_lookup.loc[sid] if sid in ploidy_lookup.index else None
        if meta is None:
            hits = ploidy_table[ploidy_table["Sample"].map(normalize_sample_id) == key]
            meta = hits.iloc[0] if not hits.empty else None
        sample_rows.append(
            {
                "Sample": sid,
                "ploidy": ploidy,
                "Group": "" if meta is None else str(meta["Group"]),
                "nominal_ploidy": "" if meta is None else str(meta["nominal_ploidy"]),
                "reason": "" if meta is None else str(meta["reason"]),
                "subgenomes": "".join(sg for sg in ("A", "B", "D") if sg in required),
            }
        )
    samples = pd.DataFrame(sample_rows)
    if not samples.empty:
        samples = samples.drop_duplicates(subset=["Sample"], keep="first").reset_index(drop=True)
    save_df_to_tsv(samples, f"{output_prefix}.samples.tsv")
    save_df_to_tsv(pairs, f"{output_prefix}.pairs.tsv")

    class_rows: list[dict[str, object]] = []
    present = [p for p in PLOIDY_ORDER if not pairs.empty and p in set(pairs["ploidy"])]
    for ploidy in PLOIDY_ORDER:
        sub = pairs[pairs["ploidy"] == ploidy] if not pairs.empty else pairs
        n_samples = int((samples["ploidy"] == ploidy).sum()) if not samples.empty else 0
        required = "".join(sg for sg in ("A", "B", "D") if sg in PLOIDY_SUBGENOMES[ploidy])
        class_rows.append(
            {
                "ploidy": ploidy,
                "subgenomes": required,
                "n_samples": n_samples,
                "n_pairs": int(len(sub)),
                "median_plink_ibs": float(sub["plink_ibs"].median()) if len(sub) else float("nan"),
                "median_hom_ibs": float(sub["hom_ibs"].median()) if len(sub) else float("nan"),
                "min_plink_ibs": float(sub["plink_ibs"].min()) if len(sub) else float("nan"),
                "max_plink_ibs": float(sub["plink_ibs"].max()) if len(sub) else float("nan"),
            }
        )
    class_df = pd.DataFrame(class_rows)
    save_df_to_tsv(class_df, f"{output_prefix}.info.tsv")
    save_thresholds(
        {
            "n_samples": int(len(samples)),
            "n_pairs": int(len(pairs)),
            "n_ploidy_classes": int(len(present)),
            "ploidy_classes": ",".join(present),
            "n_reclassified": int((samples["reason"] == "breachr_ploidy_mismatch").sum())
            if not samples.empty and "reason" in samples.columns
            else 0,
        },
        f"{output_prefix}.th.tsv",
    )
    if pairs.empty or not present:
        return pairs

    facet_labels = {}
    n_by_ploidy = samples.groupby("ploidy").size() if not samples.empty else pd.Series(dtype=int)
    for ploidy in present:
        n_s = int(n_by_ploidy.get(ploidy, 0))
        n_p = int((pairs["ploidy"] == ploidy).sum())
        sgs = "".join(sg for sg in ("A", "B", "D") if sg in PLOIDY_SUBGENOMES[ploidy])
        facet_labels[ploidy] = f"{ploidy} {sgs} (n={n_s}, pairs={n_p})"
    cut = [{"value": 0.99, "label": "0.99", "color": "#c44e52", "linestyle": "--"}]
    specs = (
        ("plink_ibs", "PLINK IBS", "plink"),
        ("hom_ibs", "Hom-SNP IBS", "hom"),
    )
    for col, x_label, tag in specs:
        xlim = _xlim_for_ibs_scores(pairs[col])
        plot_faceted_histograms(
            data=pairs,
            x_col=col,
            facet_col="ploidy",
            facet_order=present,
            title=f"{x_label} within ploidy (combined subgenomes)",
            filename=f"{output_prefix}.{tag}.dist.png",
            xlabel=x_label,
            ylabel="Count of pairs",
            xlim=xlim,
            facet_colors=PLOIDY_PALETTE,
            facet_labels=facet_labels,
            sharey=False,
            thresholds=cut,
        )
        plot_faceted_histograms(
            data=pairs,
            x_col=col,
            facet_col="ploidy",
            facet_order=present,
            title=f"{x_label} within ploidy (combined subgenomes, log count)",
            filename=f"{output_prefix}.{tag}.dist_log.png",
            xlabel=x_label,
            ylabel="Count of pairs",
            log_scale=True,
            xlim=xlim,
            facet_colors=PLOIDY_PALETTE,
            facet_labels=facet_labels,
            sharey=False,
            thresholds=cut,
        )
        plot_grouped_density_lines(
            data=pairs,
            x_col=col,
            hue_col="ploidy",
            hue_order=present,
            palette_map=PLOIDY_PALETTE,
            title=f"{x_label} within ploidy (combined subgenomes)",
            filename=f"{output_prefix}.{tag}.density.png",
            xlabel=x_label,
            ylabel="Density",
            xlim=xlim,
            hue_labels=facet_labels,
        )
    return pairs


def nearest_neighbors_from_square(
    ids: list[str],
    matrix: np.ndarray,
    *,
    k: int = 5,
) -> pd.DataFrame:
    """Top-k off-diagonal neighbors per sample from a square similarity matrix."""
    n_samples = len(ids)
    if n_samples < 2:
        return pd.DataFrame(columns=["Sample", "rank", "Neighbor", "score"])
    k = min(int(k), n_samples - 1)
    rows: list[dict] = []
    for i, sid in enumerate(ids):
        row = np.asarray(matrix[i], dtype=float).copy()
        row[i] = -np.inf
        row[~np.isfinite(row)] = -np.inf
        top = np.argpartition(row, -k)[-k:]
        top = top[np.argsort(row[top])[::-1]]
        for rank, j in enumerate(top, start=1):
            rows.append(
                {
                    "Sample": str(sid),
                    "rank": int(rank),
                    "Neighbor": str(ids[j]),
                    "score": float(matrix[i, j]),
                }
            )
    return pd.DataFrame(rows)


def write_square_mibs(
    matrix_path: str,
    id_path: str,
    ids: list[str],
    matrix: np.ndarray,
) -> None:
    """Write a PLINK-style square IBS matrix and FID/IID id file."""
    np.savetxt(matrix_path, matrix, fmt="%.6f")
    Path(id_path).write_text("".join(f"{sid}\t{sid}\n" for sid in ids), encoding="utf-8")


def load_square_mibs(matrix_path: str, id_path: str) -> tuple[list[str], np.ndarray]:
    """Load a space- or tab-separated square IBS matrix plus its id file."""
    ids = load_ids_file(id_path)
    mat = np.loadtxt(matrix_path, dtype=float)
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        raise ValueError(f"{matrix_path} is not a square matrix: {getattr(mat, 'shape', None)}")
    if ids and mat.shape[0] != len(ids):
        raise ValueError(f"{matrix_path} shape {mat.shape} != {len(ids)} IDs")
    return ids, mat


def _align_square_by_ids(
    ids_a: list[str],
    mat_a: np.ndarray,
    ids_b: list[str],
    mat_b: np.ndarray,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    pos_b = {sid: i for i, sid in enumerate(ids_b)}
    common = [sid for sid in ids_a if sid in pos_b]
    ia = {sid: i for i, sid in enumerate(ids_a)}
    idx_a = [ia[sid] for sid in common]
    idx_b = [pos_b[sid] for sid in common]
    return common, mat_a[np.ix_(idx_a, idx_a)], mat_b[np.ix_(idx_b, idx_b)]


def _write_score_outputs(
    df: pd.DataFrame,
    output_prefix: str,
    *,
    score_col: str,
    title: str,
    x_label: str,
) -> None:
    save_df_to_tsv(df, f"{output_prefix}.pairs.tsv")
    stats = {
        "n_pairs": int(len(df)),
        f"min_{score_col}": float(df[score_col].min()) if len(df) else np.nan,
        f"max_{score_col}": float(df[score_col].max()) if len(df) else np.nan,
        f"median_{score_col}": float(df[score_col].median()) if len(df) else np.nan,
    }
    save_thresholds(stats, f"{output_prefix}.th.tsv")
    if df.empty:
        return
    plot_distribution_with_stats(
        data=df,
        col=score_col,
        title=title,
        filename=f"{output_prefix}.dist.png",
        x_label=x_label,
        y_label="Count of pairs",
        mean_val=float(df[score_col].mean()),
        median_val=float(df[score_col].median()),
    )


def ana_rare_jaccard_vcf(
    vcf_path: str,
    output_prefix: str,
    *,
    min_shared: int = 5,
    min_jaccard: float = 0.5,
    max_carriers: int = 400,
) -> pd.DataFrame:
    carriers = carriers_from_vcf(vcf_path, max_carriers=max_carriers)
    df = jaccard_from_carriers(carriers, min_shared=min_shared, min_jaccard=min_jaccard)
    _write_score_outputs(
        df,
        output_prefix,
        score_col="jaccard",
        title="Rare-allele Jaccard (pairs above threshold)",
        x_label="Jaccard",
    )
    return df


def ana_fingerprint_additive(
    raw_path: str,
    output_prefix: str,
    *,
    max_discord: float = 0.05,
    min_relatedness: float = 0.4,
    min_nobs: int = 100,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    ids, gt = load_plink_export_a(raw_path)
    discord = discord_rate_pairs(ids, gt, max_discord=max_discord, min_nobs=min_nobs)
    related = somalier_relatedness_pairs(
        ids, gt, min_relatedness=min_relatedness, min_nobs=min_nobs
    )
    _write_score_outputs(
        discord,
        f"{output_prefix}.discord",
        score_col="discord_rate",
        title="Genotype discordance (pairs below threshold)",
        x_label="Discordance rate",
    )
    _write_score_outputs(
        related,
        f"{output_prefix}.fingerprint",
        score_col="relatedness",
        title="Somalier-style relatedness (VCF panel)",
        x_label="Relatedness",
    )
    return discord, related


def ana_king_identity(
    king_path: str,
    output_prefix: str,
    *,
    min_kinship: float = 0.177,
) -> pd.DataFrame:
    df = load_df_from_king(king_path)
    if df is None or df.empty or "KINSHIP" not in df.columns:
        out = _empty_pair_frame(["KINSHIP", "IBS0", "NSNP"])
        _write_score_outputs(out, output_prefix, score_col="KINSHIP", title="KING kinship", x_label="KINSHIP")
        return out
    rename = {}
    if "Sample1" not in df.columns:
        for src in ("#IID1", "IID1", "ID1"):
            if src in df.columns:
                rename[src] = "Sample1"
                break
    if "Sample2" not in df.columns:
        for src in ("IID2", "ID2"):
            if src in df.columns:
                rename[src] = "Sample2"
                break
    df = df.rename(columns=rename)
    keep = df[df["KINSHIP"] >= min_kinship].copy()
    ordered = [_order_pair(str(a), str(b)) for a, b in zip(keep["Sample1"], keep["Sample2"])]
    keep["Sample1"] = [p[0] for p in ordered]
    keep["Sample2"] = [p[1] for p in ordered]
    extra = [c for c in ("KINSHIP", "IBS0", "NSNP") if c in keep.columns]
    keep = keep[PAIR_COLUMNS + extra]
    _write_score_outputs(
        keep,
        output_prefix,
        score_col="KINSHIP",
        title="KING kinship (pairs above threshold)",
        x_label="KINSHIP",
    )
    return keep


def ana_genome_ibd(
    genome_path: str,
    output_prefix: str,
    *,
    min_pi_hat: float = 0.25,
) -> pd.DataFrame:
    df = pd.read_csv(genome_path, sep=r"\s+", engine="python")
    if df.empty or "PI_HAT" not in df.columns:
        out = _empty_pair_frame(["PI_HAT", "Z0", "Z1", "Z2", "IBS0"])
        _write_score_outputs(out, output_prefix, score_col="PI_HAT", title="PLINK IBD", x_label="PI_HAT")
        return out
    keep = df[df["PI_HAT"] >= min_pi_hat].copy()
    keep["Sample1"] = keep["IID1"].astype(str)
    keep["Sample2"] = keep["IID2"].astype(str)
    ordered = [_order_pair(a, b) for a, b in zip(keep["Sample1"], keep["Sample2"])]
    keep["Sample1"] = [p[0] for p in ordered]
    keep["Sample2"] = [p[1] for p in ordered]
    extra = [c for c in ("PI_HAT", "Z0", "Z1", "Z2", "IBS0") if c in keep.columns]
    keep = keep[PAIR_COLUMNS + extra]
    _write_score_outputs(
        keep,
        output_prefix,
        score_col="PI_HAT",
        title="PLINK IBD PI_HAT (pairs above threshold)",
        x_label="PI_HAT",
    )
    return keep


def ana_identity_consensus(
    method_tables: dict[str, pd.DataFrame],
    output_prefix: str,
    *,
    min_methods: int = 2,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pairs = consensus_pairs(method_tables, min_methods=min_methods)
    clusters = clusters_to_frame(cluster_undirected_pairs(pairs))
    save_df_to_tsv(pairs, f"{output_prefix}.consensus.pairs.tsv")
    save_df_to_tsv(clusters, f"{output_prefix}.consensus.clusters.tsv")
    save_thresholds(
        {
            "n_consensus_pairs": int(len(pairs)),
            "n_consensus_clusters": int(len(clusters)),
            "n_samples_in_clusters": int(clusters["n_samples"].sum()) if len(clusters) else 0,
            "min_methods": min_methods,
        },
        f"{output_prefix}.consensus.th.tsv",
    )
    return pairs, clusters


def write_hom_snp_ibs_matrix(
    bed_prefix: str,
    matrix_prefix: str,
    *,
    chunk_sites: int = 2048,
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    """Compute hom-SNP IBS from a bed and write a PLINK-style square matrix."""
    ids, ibs, n_comp, n_match = hom_snp_ibs_from_bed(
        bed_prefix, chunk_sites=chunk_sites, snp_only=True
    )
    write_square_mibs(f"{matrix_prefix}.mibs", f"{matrix_prefix}.mibs.id", ids, ibs)
    return ids, ibs, n_comp, n_match


def ana_hom_snp_compare(
    ids: list[str],
    ibs: np.ndarray,
    n_comp: np.ndarray,
    n_match: np.ndarray,
    output_prefix: str,
    *,
    plink_ibs: np.ndarray | None = None,
    plink_ids: list[str] | None = None,
    min_pair_ibs: float = 0.99,
    min_compare_ibs: float = 0.97,
    min_nobs: int = 100,
    neighbor_k: int = 5,
    scatter_n: int = 200_000,
    scatter_seed: int = 1,
) -> pd.DataFrame:
    """Pair tables, neighbor shifts, and PLINK-vs-hom plots from square matrices."""
    pairs = pairs_from_square_ibs(
        ids, ibs, n_comp, n_match, min_ibs=min_pair_ibs, min_nobs=min_nobs
    )
    _write_score_outputs(
        pairs,
        output_prefix,
        score_col="hom_ibs",
        title="Hom-SNP IBS (biallelic SNPs; both samples homozygous)",
        x_label="Hom-SNP IBS",
    )
    clusters = clusters_to_frame(cluster_undirected_pairs(pairs))
    save_df_to_tsv(clusters, f"{output_prefix}.clusters.tsv")
    save_thresholds(
        {
            "n_samples": int(len(ids)),
            "n_hom_pairs": int(len(pairs)),
            "n_hom_clusters": int(len(clusters)),
            "min_pair_ibs": float(min_pair_ibs),
            "median_n_hom": float(np.median(n_comp[np.triu_indices(len(ids), k=1)]))
            if len(ids) > 1
            else 0.0,
        },
        f"{output_prefix}.th.tsv",
    )
    nn_hom = nearest_neighbors_from_square(ids, ibs, k=neighbor_k)
    save_df_to_tsv(nn_hom, f"{output_prefix}.neighbors.tsv")

    if plink_ibs is None or plink_ids is None:
        return pairs
    common, plink_al, hom_al = _align_square_by_ids(plink_ids, plink_ibs, ids, ibs)
    if len(common) < 2:
        return pairs
    pos = {sid: i for i, sid in enumerate(ids)}
    idx = [pos[sid] for sid in common]
    n_al = n_comp[np.ix_(idx, idx)]
    nn_plink = nearest_neighbors_from_square(common, plink_al, k=neighbor_k)
    save_df_to_tsv(nn_plink, f"{output_prefix}.neighbors.plink.tsv")
    rank1_h = nn_hom[nn_hom["rank"] == 1][["Sample", "Neighbor", "score"]].rename(
        columns={"Neighbor": "neighbor_hom", "score": "hom_ibs"}
    )
    rank1_p = nn_plink[nn_plink["rank"] == 1][["Sample", "Neighbor", "score"]].rename(
        columns={"Neighbor": "neighbor_plink", "score": "plink_ibs"}
    )
    shifts = rank1_p.merge(rank1_h, on="Sample")
    shifts = shifts[shifts["neighbor_plink"] != shifts["neighbor_hom"]].copy()
    save_df_to_tsv(shifts, f"{output_prefix}.neighbors.shift.tsv")
    scatter = sample_square_ibs_pairs(
        common, plink_al, hom_al, n_al, n_pairs=scatter_n, seed=scatter_seed
    )
    save_df_to_tsv(scatter, f"{output_prefix}.vs_plink.pairs.tsv")
    if not scatter.empty:
        lo = float(min(scatter["plink_ibs"].min(), scatter["hom_ibs"].min()))
        hi = float(max(scatter["plink_ibs"].max(), scatter["hom_ibs"].max()))
        pad = max(0.002, (hi - lo) * 0.05) if hi > lo else 0.01
        lim = (lo - pad, min(1.02, hi + pad))
        plot_scatter_density(
            scatter,
            "plink_ibs",
            "hom_ibs",
            title="PLINK IBS vs hom-SNP IBS",
            filename=f"{output_prefix}.vs_plink.scatter.png",
            xlabel="PLINK IBS (all genotypes)",
            ylabel="Hom-SNP IBS (both homozygous)",
            cbar_label="Local point density",
            cmap="inferno",
            s=12,
            figure_size=(8, 8),
            xlim=lim,
            ylim=lim,
            aspect="equal",
            line_overlays=[
                {
                    "x": [lim[0], lim[1]],
                    "y": [lim[0], lim[1]],
                    "color": "#7f7f7f",
                    "linestyle": "--",
                    "linewidth": 1.2,
                    "label": r"$y=x$",
                    "zorder": 5,
                }
            ],
        )
        plot_distribution_with_stats(
            data=scatter,
            col="ibs_delta",
            title="Hom-SNP IBS minus PLINK IBS",
            filename=f"{output_prefix}.vs_plink.delta.png",
            x_label="hom_ibs - plink_ibs",
            y_label="Count of pairs",
            mean_val=float(scatter["ibs_delta"].mean()),
            median_val=float(scatter["ibs_delta"].median()),
        )
    return pairs


def ana_hom_snp_ibs_from_mibs(
    hom_mibs: str,
    hom_mibs_id: str,
    output_prefix: str,
    *,
    plink_mibs: str | None = None,
    plink_mibs_id: str | None = None,
    min_pair_ibs: float = 0.99,
    min_compare_ibs: float = 0.97,
    min_nobs: int = 0,
    neighbor_k: int = 5,
    scatter_n: int = 200_000,
    scatter_seed: int = 1,
) -> pd.DataFrame:
    """Compare published hom-SNP and PLINK square IBS matrices (no bed recompute)."""
    ids, ibs = load_square_mibs(hom_mibs, hom_mibs_id)
    if not ids:
        return _empty_pair_frame(["n_hom", "n_match", "hom_ibs"])
    n_comp = np.zeros_like(ibs)
    n_match = np.zeros_like(ibs)
    plink_mat = None
    plink_ids = None
    if plink_mibs and plink_mibs_id:
        plink_ids, plink_mat = load_square_mibs(plink_mibs, plink_mibs_id)
    return ana_hom_snp_compare(
        ids,
        ibs,
        n_comp,
        n_match,
        output_prefix,
        plink_ibs=plink_mat,
        plink_ids=plink_ids,
        min_pair_ibs=min_pair_ibs,
        min_compare_ibs=min_compare_ibs,
        min_nobs=min_nobs,
        neighbor_k=neighbor_k,
        scatter_n=scatter_n,
        scatter_seed=scatter_seed,
    )


def ana_hom_snp_ibs(
    bed_prefix: str,
    output_prefix: str,
    *,
    matrix_prefix: str | None = None,
    plink_mibs: str | None = None,
    plink_mibs_id: str | None = None,
    min_pair_ibs: float = 0.99,
    min_compare_ibs: float = 0.97,
    min_nobs: int = 100,
    chunk_sites: int = 2048,
    neighbor_k: int = 5,
) -> pd.DataFrame:
    """Hom-SNP IBS from a bed, plus optional side-by-side PLINK IBS compare."""
    mx = matrix_prefix or output_prefix
    ids, ibs, n_comp, n_match = write_hom_snp_ibs_matrix(
        bed_prefix, mx, chunk_sites=chunk_sites
    )
    plink_mat = None
    plink_ids = None
    if plink_mibs and plink_mibs_id:
        plink_ids, plink_mat = load_square_mibs(plink_mibs, plink_mibs_id)
    return ana_hom_snp_compare(
        ids,
        ibs,
        n_comp,
        n_match,
        output_prefix,
        plink_ibs=plink_mat,
        plink_ids=plink_ids,
        min_pair_ibs=min_pair_ibs,
        min_compare_ibs=min_compare_ibs,
        min_nobs=min_nobs,
        neighbor_k=neighbor_k,
    )


def plot_identity_pair_table(pairs_path: str, output_prefix: str, method: str) -> None:
    """Histogram a published identity pair table via ``plot_distribution_with_stats``."""
    spec = IDENTITY_PLOT_SPECS.get(method)
    if spec is None:
        raise ValueError(f"Unknown identity method {method}")
    df = load_df_from_tsv(pairs_path)
    col = spec["col"]
    if df is None or df.empty or col not in df.columns:
        return
    scores = pd.to_numeric(df[col], errors="coerce").dropna()
    if scores.empty:
        return
    work = df.loc[scores.index, [col]].copy()
    plot_distribution_with_stats(
        data=work,
        col=col,
        title=spec["title"],
        filename=f"{output_prefix}.dist.png",
        x_label=spec["x_label"],
        y_label="Count of pairs",
        mean_val=float(scores.mean()),
        median_val=float(scores.median()),
        xlim=spec["xlim"],
        thresholds=spec["thresholds"],
    )


def ana_hom_ibs_one_vs_plink(
    hom_mibs: str,
    hom_mibs_id: str,
    plink_mibs: str,
    plink_mibs_id: str,
    output_prefix: str,
    *,
    atol: float = 5e-7,
) -> pd.DataFrame:
    """Pairs with hom-SNP IBS = 1, plus the PLINK IBS distribution of those pairs."""
    hom_ids, hom = load_square_mibs(hom_mibs, hom_mibs_id)
    plink_ids, plink = load_square_mibs(plink_mibs, plink_mibs_id)
    common, hom_al, plink_al = _align_square_by_ids(hom_ids, hom, plink_ids, plink)
    pairs = pairs_hom_ibs_one(common, hom_al, plink_al, atol=atol)
    save_df_to_tsv(pairs, f"{output_prefix}.eq1.pairs.tsv")
    stats = {
        "n_samples": int(len(common)),
        "n_pairs": int(len(pairs)),
        "min_plink_ibs": float(pairs["plink_ibs"].min()) if len(pairs) else np.nan,
        "median_plink_ibs": float(pairs["plink_ibs"].median()) if len(pairs) else np.nan,
        "max_plink_ibs": float(pairs["plink_ibs"].max()) if len(pairs) else np.nan,
        "n_plink_lt_0.99": int((pairs["plink_ibs"] < 0.99).sum()) if len(pairs) else 0,
        "n_plink_eq1": int(np.isclose(pairs["plink_ibs"], 1.0, rtol=0.0, atol=atol).sum())
        if len(pairs)
        else 0,
    }
    save_thresholds(stats, f"{output_prefix}.eq1.th.tsv")
    if pairs.empty:
        return pairs
    lo = float(pairs["plink_ibs"].min())
    pad = max(0.002, (1.0 - lo) * 0.08)
    plot_distribution_with_stats(
        data=pairs,
        col="plink_ibs",
        title="PLINK IBS of pairs with hom-SNP IBS = 1",
        filename=f"{output_prefix}.eq1.plink.dist.png",
        x_label="PLINK IBS (all genotypes)",
        y_label="Count of pairs",
        mean_val=float(pairs["plink_ibs"].mean()),
        median_val=float(pairs["plink_ibs"].median()),
        xlim=(max(0.0, lo - pad), 1.02),
        thresholds=[
            {
                "value": 1.0,
                "label": "IBS = 1",
                "color": "#c44e52",
                "linestyle": "--",
            }
        ],
    )
    plot_slope_chart(
        pairs,
        "plink_ibs",
        "hom_ibs",
        "PLINK IBS",
        "Hom-SNP IBS",
        "Pairs with hom-SNP IBS = 1",
        f"{output_prefix}.eq1.slope.png",
        ylabel="IBS",
        ylim=(max(0.0, lo - pad), 1.02),
        line_alpha=0.35,
    )
    return pairs

