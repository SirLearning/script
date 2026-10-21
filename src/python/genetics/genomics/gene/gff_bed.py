"""Merged BED intervals from a GFF3 feature type (plink2 ``--extract bed0``)."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pandas as pd

from genetics.genomics.variant.coding_mask import subgenome_from_gene_id


DEFAULT_FEATURE_TYPES = ("exon",)


def parse_gff_feature_types(raw: str | None) -> tuple[str, ...]:
    """Split a comma-separated GFF type list; empty input means exons only."""

    if raw is None:
        return DEFAULT_FEATURE_TYPES
    types = tuple(part.strip() for part in str(raw).split(",") if part.strip())
    return types or DEFAULT_FEATURE_TYPES


def _merge_half_open(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals.sort()
    merged: list[tuple[int, int]] = []
    start, end = intervals[0]
    for next_start, next_end in intervals[1:]:
        if next_start <= end:
            if next_end > end:
                end = next_end
        else:
            merged.append((start, end))
            start, end = next_start, next_end
    merged.append((start, end))
    return merged


def collect_merged_gff_intervals(
    gff_path: str | Path,
    *,
    feature_types: tuple[str, ...] = DEFAULT_FEATURE_TYPES,
) -> dict[str, list[tuple[int, int]]]:
    """0-based half-open intervals per chromosome, isoforms merged."""

    wanted = frozenset(feature_types)
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    path = Path(gff_path)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[2] not in wanted:
                continue
            chrom = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            if end <= start:
                continue
            by_chrom[chrom].append((start, end))
    return {chrom: _merge_half_open(iv) for chrom, iv in by_chrom.items()}


def _circos_chrom_name(chrom: str) -> str:
    token = str(chrom).strip()
    if token.startswith("chr"):
        return token[3:]
    return token


def bin_gff_feature_counts(
    gff_path: str | Path,
    bin_bp: int,
    *,
    feature_type: str = "gene",
) -> pd.DataFrame:
    """Count GFF features per ``bin_bp`` window; chromosome names drop a ``chr`` prefix."""

    bin_bp = int(bin_bp)
    if bin_bp <= 0:
        raise ValueError("bin_bp must be positive")
    counts: dict[tuple[str, int], int] = defaultdict(int)
    path = Path(gff_path)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[2] != feature_type:
                continue
            start = int(fields[3]) - 1
            if start < 0:
                continue
            chrom = _circos_chrom_name(fields[0])
            bin_start = (start // bin_bp) * bin_bp
            counts[(chrom, bin_start)] += 1
    rows = [
        {
            "chrom": chrom,
            "bin_start": start,
            "bin_end": start + bin_bp,
            "value": n,
        }
        for (chrom, start), n in sorted(counts.items())
    ]
    if not rows:
        return pd.DataFrame(columns=["chrom", "bin_start", "bin_end", "value"])
    return pd.DataFrame(rows)


def write_merged_gff_bed(
    gff_path: str | Path,
    bed_path: str | Path,
    *,
    feature_types: tuple[str, ...] = DEFAULT_FEATURE_TYPES,
) -> int:
    """Write a 3-column BED for ``plink2 --extract bed0``. Returns interval count."""

    intervals = collect_merged_gff_intervals(gff_path, feature_types=feature_types)
    n_written = 0
    out = Path(bed_path)
    with out.open("w", encoding="utf-8") as handle:
        for chrom in sorted(intervals, key=lambda name: (len(name), name)):
            for start, end in intervals[chrom]:
                handle.write(f"{chrom}\t{start}\t{end}\n")
                n_written += 1
    if n_written == 0:
        raise ValueError(
            f"No {','.join(feature_types)} intervals in {gff_path}"
        )
    return n_written


def _gff_attributes(raw: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in str(raw).split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        out[key] = value
    return out


def gene_lengths_from_gff(gff_path: str | Path) -> pd.DataFrame:
    """One row per gene: genomic span, longest-transcript CDS bp, exon count."""

    path = Path(gff_path)
    gene_span: dict[str, int] = {}
    transcript_gene: dict[str, str] = {}
    cds_bp: dict[str, int] = defaultdict(int)
    n_exons: dict[str, int] = defaultdict(int)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            kind = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            if end < start:
                continue
            attrs = _gff_attributes(fields[8])
            if kind == "gene":
                gene_id = attrs.get("ID", "")
                if gene_id:
                    gene_span[gene_id] = end - start + 1
            elif kind == "mRNA":
                tx_id = attrs.get("ID", "")
                parent = attrs.get("Parent", "")
                if tx_id and parent:
                    transcript_gene[tx_id] = parent
            elif kind == "CDS":
                parent = attrs.get("Parent", "")
                if parent:
                    cds_bp[parent] += end - start + 1
            elif kind == "exon":
                parent = attrs.get("Parent", "")
                if parent:
                    n_exons[parent] += 1

    best_tx: dict[str, str] = {}
    best_cds: dict[str, int] = {}
    for tx_id, length in cds_bp.items():
        gene_id = transcript_gene.get(tx_id, tx_id.rsplit(".", 1)[0] if "." in tx_id else tx_id)
        if gene_id not in best_cds or length > best_cds[gene_id]:
            best_cds[gene_id] = length
            best_tx[gene_id] = tx_id
        if gene_id not in gene_span:
            gene_span[gene_id] = 0

    rows = []
    for gene_id, span in sorted(gene_span.items()):
        tx_id = best_tx.get(gene_id, "")
        rows.append(
            {
                "gene_id": gene_id,
                "gene_span_bp": int(span),
                "cds_bp": int(best_cds.get(gene_id, 0)),
                "n_exons": int(n_exons.get(tx_id, 0)),
                "subgenome": subgenome_from_gene_id(gene_id),
            }
        )
    if not rows:
        raise ValueError(f"No gene features in {path}")
    return pd.DataFrame(rows)
