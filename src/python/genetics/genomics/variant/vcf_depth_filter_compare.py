"""Compare two FastCall3 VCFs for MiDR/MaDR genotype masking (same sites, same AD)."""

from __future__ import annotations

import gzip
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

from infra.utils.errors import fail

logger = logging.getLogger(__name__)

MISSING_GT = {".", "./.", ".|."}


@dataclass
class CellClassCounts:
    both_missing: int = 0
    unfiltered_called_filtered_missing: int = 0
    both_called: int = 0
    unfiltered_missing_filtered_called: int = 0
    both_called_gt_mismatch: int = 0
    both_called_ad_mismatch: int = 0
    het_unfiltered: int = 0
    het_filtered: int = 0
    het_to_missing: int = 0
    het_to_hom: int = 0
    midr_explained: int = 0
    madr_explained: int = 0
    called_to_missing_unexplained: int = 0
    leak_still_called: int = 0


@dataclass
class VcfFilterCompareResult:
    n_sites_unfiltered: int = 0
    n_sites_filtered: int = 0
    n_sites_matched: int = 0
    n_sites_only_unfiltered: int = 0
    n_sites_only_filtered: int = 0
    n_sites_allele_mismatch: int = 0
    samples_equal: bool = False
    n_samples: int = 0
    cells: CellClassCounts = field(default_factory=CellClassCounts)
    info_dp_unfiltered: int = 0
    info_dp_filtered: int = 0
    info_nz_unfiltered: int = 0
    info_nz_filtered: int = 0
    missing_unfiltered: int = 0
    missing_filtered: int = 0
    n_cells: int = 0


def parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if not info or info == ".":
        return out
    for item in info.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            out[key] = val
        else:
            out[item] = "true"
    return out


def is_missing_gt(gt: str) -> bool:
    return gt.split(":")[0] in MISSING_GT


def is_het_gt(gt: str) -> bool:
    token = gt.split(":")[0]
    if token in MISSING_GT:
        return False
    sep = "/" if "/" in token else "|" if "|" in token else None
    if sep is None:
        return False
    left, right = token.split(sep, 1)
    return left != right


def parse_ad(sample: str, format_keys: list[str]) -> list[int] | None:
    parts = sample.split(":")
    if "AD" not in format_keys:
        return None
    idx = format_keys.index("AD")
    if idx >= len(parts) or parts[idx] in {".", ""}:
        return None
    try:
        return [int(x) for x in parts[idx].split(",")]
    except ValueError:
        return None


def site_depth(ad: list[int] | None) -> int:
    if not ad:
        return 0
    return int(sum(ad))


def load_taxa_coverage(taxa_bam_map: Path) -> dict[str, float]:
    coverage: dict[str, float] = {}
    with taxa_bam_map.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            fail(f"Empty taxaBamMap: {taxa_bam_map}")
        for line in handle:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            coverage[cols[0]] = float(cols[1])
    return coverage


def _open_vcf(path: Path):
    name = path.name
    if name.endswith(".gz") or name.endswith(".tmp"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def iter_vcf_records(path: Path, max_sites: int | None = None) -> tuple[list[str], Iterator[tuple]]:
    """Yield (chrom, pos, ref, alt, info, format_keys, samples) after reading header."""
    handle = _open_vcf(path)
    samples: list[str] = []

    def _records():
        nonlocal samples
        n = 0
        try:
            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    cols = line.rstrip("\n").split("\t")
                    samples = cols[9:]
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 10:
                    continue
                info = parse_info(cols[7])
                fmt = cols[8].split(":")
                yield cols[0], int(cols[1]), cols[3], cols[4], info, fmt, cols[9:]
                n += 1
                if max_sites is not None and n >= max_sites:
                    break
        finally:
            handle.close()

    # Prime header so samples is populated before returning the iterator.
    rec_iter = _records()
    first = None
    try:
        first = next(rec_iter)
    except StopIteration:
        return samples, iter(())
    def _with_first():
        yield first
        yield from rec_iter
    return samples, _with_first()


def compare_filter_vcfs(
    unfiltered_vcf: Path,
    filtered_vcf: Path,
    taxa_bam_map: Path,
    *,
    midr: float = 0.2,
    madr: float = 3.0,
    max_sites: int | None = 2000,
    max_pos: int | None = None,
) -> VcfFilterCompareResult:
    """Walk two VCFs in POS order and classify per-sample GT masking."""
    u_samples, u_iter = iter_vcf_records(unfiltered_vcf, max_sites=None)
    f_samples, f_iter = iter_vcf_records(filtered_vcf, max_sites=None)
    result = VcfFilterCompareResult()
    result.samples_equal = u_samples == f_samples
    result.n_samples = len(f_samples)
    if not result.samples_equal:
        logger.warning("Sample columns differ; comparing intersection order of filtered VCF")
    coverage = load_taxa_coverage(taxa_bam_map)
    cov_list = [coverage.get(name, float("nan")) for name in f_samples]

    u_rec = next(u_iter, None)
    f_rec = next(f_iter, None)
    cells = result.cells

    while u_rec is not None and f_rec is not None:
        if max_pos is not None and u_rec[1] > max_pos and f_rec[1] > max_pos:
            break
        if max_sites is not None and result.n_sites_matched >= max_sites:
            break
        u_key = (u_rec[0], u_rec[1])
        f_key = (f_rec[0], f_rec[1])
        if u_key < f_key:
            result.n_sites_only_unfiltered += 1
            result.n_sites_unfiltered += 1
            u_rec = next(u_iter, None)
            continue
        if f_key < u_key:
            result.n_sites_only_filtered += 1
            result.n_sites_filtered += 1
            f_rec = next(f_iter, None)
            continue

        result.n_sites_unfiltered += 1
        result.n_sites_filtered += 1
        if u_rec[2] != f_rec[2] or u_rec[3] != f_rec[3]:
            result.n_sites_allele_mismatch += 1
            u_rec = next(u_iter, None)
            f_rec = next(f_iter, None)
            continue

        result.n_sites_matched += 1
        u_info, u_fmt, u_samp = u_rec[4], u_rec[5], u_rec[6]
        f_info, f_fmt, f_samp = f_rec[4], f_rec[5], f_rec[6]
        try:
            result.info_dp_unfiltered += int(float(u_info.get("DP", "0")))
            result.info_dp_filtered += int(float(f_info.get("DP", "0")))
            result.info_nz_unfiltered += int(float(u_info.get("NZ", "0")))
            result.info_nz_filtered += int(float(f_info.get("NZ", "0")))
        except ValueError:
            pass

        n = min(len(u_samp), len(f_samp), len(f_samples))
        for i in range(n):
            u_gt = u_samp[i]
            f_gt = f_samp[i]
            u_miss = is_missing_gt(u_gt)
            f_miss = is_missing_gt(f_gt)
            result.n_cells += 1
            if u_miss:
                result.missing_unfiltered += 1
            if f_miss:
                result.missing_filtered += 1
            u_het = is_het_gt(u_gt)
            f_het = is_het_gt(f_gt)
            if u_het:
                cells.het_unfiltered += 1
            if f_het:
                cells.het_filtered += 1

            u_ad = parse_ad(u_gt, u_fmt)
            f_ad = parse_ad(f_gt, f_fmt)
            dp = site_depth(u_ad)
            cov = cov_list[i]
            ratio = (dp / cov) if cov and cov > 0 else float("inf")
            should_mask = (ratio < midr) or (ratio > madr)

            if u_miss and f_miss:
                cells.both_missing += 1
            elif (not u_miss) and f_miss:
                cells.unfiltered_called_filtered_missing += 1
                if u_het:
                    cells.het_to_missing += 1
                if ratio < midr:
                    cells.midr_explained += 1
                elif ratio > madr:
                    cells.madr_explained += 1
                else:
                    cells.called_to_missing_unexplained += 1
            elif u_miss and (not f_miss):
                cells.unfiltered_missing_filtered_called += 1
            else:
                cells.both_called += 1
                if u_gt.split(":")[0] != f_gt.split(":")[0]:
                    cells.both_called_gt_mismatch += 1
                    if u_het and not f_het:
                        cells.het_to_hom += 1
                if u_ad != f_ad:
                    cells.both_called_ad_mismatch += 1
                if should_mask:
                    cells.leak_still_called += 1

        u_rec = next(u_iter, None)
        f_rec = next(f_iter, None)

    return result


def result_to_rows(result: VcfFilterCompareResult) -> list[dict]:
    n_cells = result.n_cells or 1
    miss_u = result.missing_unfiltered / n_cells
    miss_f = result.missing_filtered / n_cells
    return [
        {"metric": "n_sites_matched", "value": result.n_sites_matched},
        {"metric": "n_sites_only_unfiltered", "value": result.n_sites_only_unfiltered},
        {"metric": "n_sites_only_filtered", "value": result.n_sites_only_filtered},
        {"metric": "n_sites_allele_mismatch", "value": result.n_sites_allele_mismatch},
        {"metric": "samples_equal", "value": int(result.samples_equal)},
        {"metric": "n_samples", "value": result.n_samples},
        {"metric": "n_cells", "value": result.n_cells},
        {"metric": "both_missing", "value": result.cells.both_missing},
        {"metric": "unfiltered_called_filtered_missing", "value": result.cells.unfiltered_called_filtered_missing},
        {"metric": "both_called", "value": result.cells.both_called},
        {"metric": "unfiltered_missing_filtered_called", "value": result.cells.unfiltered_missing_filtered_called},
        {"metric": "both_called_gt_mismatch", "value": result.cells.both_called_gt_mismatch},
        {"metric": "both_called_ad_mismatch", "value": result.cells.both_called_ad_mismatch},
        {"metric": "het_unfiltered", "value": result.cells.het_unfiltered},
        {"metric": "het_filtered", "value": result.cells.het_filtered},
        {"metric": "het_to_missing", "value": result.cells.het_to_missing},
        {"metric": "het_to_hom", "value": result.cells.het_to_hom},
        {"metric": "midr_explained", "value": result.cells.midr_explained},
        {"metric": "madr_explained", "value": result.cells.madr_explained},
        {"metric": "called_to_missing_unexplained", "value": result.cells.called_to_missing_unexplained},
        {"metric": "leak_still_called", "value": result.cells.leak_still_called},
        {"metric": "missing_rate_unfiltered", "value": miss_u},
        {"metric": "missing_rate_filtered", "value": miss_f},
        {"metric": "info_dp_unfiltered", "value": result.info_dp_unfiltered},
        {"metric": "info_dp_filtered", "value": result.info_dp_filtered},
        {"metric": "info_nz_unfiltered", "value": result.info_nz_unfiltered},
        {"metric": "info_nz_filtered", "value": result.info_nz_filtered},
    ]
