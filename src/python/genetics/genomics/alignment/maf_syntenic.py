"""Extract syntenic alignment columns from multi-species MAF files."""

from __future__ import annotations

import gzip
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterator, TextIO

from infra.utils.io import save_df_to_tsv, save_thresholds

_SCORE_RE = re.compile(r"score=([-\d.]+)")


@dataclass(frozen=True)
class MafRow:
    src: str
    species: str
    chrom: str
    start: int
    size: int
    strand: str
    src_size: int
    seq: str


@dataclass
class BlockStats:
    blocks_total: int = 0
    blocks_skipped_score: int = 0
    blocks_skipped_no_ref: int = 0
    blocks_processed: int = 0
    rows_dropped_cross_chr: int = 0
    columns_total: int = 0
    columns_emitted: int = 0


def chr_label_to_ref(chr_label: str) -> tuple[str, str]:
    """Map wheat chr label (e.g. 1A) to (ref_species, ref_chrom)."""
    subgenome = chr_label[-1].upper()
    ref_num = chr_label[:-1]
    ref_species = {"A": "traesA", "B": "traesB", "D": "traesD"}.get(subgenome)
    if ref_species is None:
        raise ValueError(f"Unsupported subgenome in chr label: {chr_label}")
    return ref_species, f"chr{ref_num}"


def load_species_prefixes(species_list_path: str | None) -> list[str] | None:
    if not species_list_path:
        return None
    prefixes: list[str] = []
    with open(species_list_path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            prefixes.append(line.split()[0])
    return prefixes


def _open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def _open_binary_write(path: str) -> BinaryIO:
    if path.endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8")
    return open(path, "w", encoding="utf-8")


def _parse_row(line: str) -> MafRow:
    parts = line.rstrip("\n").split()
    if len(parts) < 7 or parts[0] != "s":
        raise ValueError(f"Invalid MAF row: {line[:80]}")
    src = parts[1]
    species, _, chrom = src.partition(".")
    return MafRow(
        src=src,
        species=species,
        chrom=chrom,
        start=int(parts[2]),
        size=int(parts[3]),
        strand=parts[4],
        src_size=int(parts[5]),
        seq=parts[6],
    )


def _iter_blocks(handle: TextIO) -> Iterator[tuple[float | None, list[MafRow]]]:
    score: float | None = None
    rows: list[MafRow] = []
    for line in handle:
        if line.startswith("a "):
            if rows:
                yield score, rows
            match = _SCORE_RE.search(line)
            score = float(match.group(1)) if match else None
            rows = []
        elif line.startswith("s "):
            rows.append(_parse_row(line))
    if rows:
        yield score, rows


def _keep_row(row: MafRow, ref_chrom: str, filter_cross_chr: bool) -> bool:
    if not filter_cross_chr:
        return True
    if not row.chrom.startswith("chr"):
        return True
    return row.chrom == ref_chrom


def _build_species_columns(
    species_order: list[str] | None,
    ref_species: str,
) -> list[str] | None:
    """Return fixed output columns when a species list is provided."""
    if not species_order:
        return None
    cols: list[str] = []
    seen: set[str] = set()
    for sp in [ref_species, *species_order]:
        if sp not in seen:
            cols.append(sp)
            seen.add(sp)
    return cols


def _ordered_species(
    species_order: list[str] | None,
    ref_species: str,
    seen: set[str],
) -> list[str]:
    fixed = _build_species_columns(species_order, ref_species)
    if fixed is not None:
        return fixed
    return sorted(seen)


def extract_syntenic_sites_from_maf(
    maf_path: str,
    output_prefix: str,
    chr_label: str,
    *,
    species_list_path: str | None = None,
    min_species: int = 2,
    min_score: float = 0.0,
    filter_cross_chr: bool = True,
) -> dict[str, int | float | str]:
    """Extract k-of-n syntenic alignment columns from a subgenome MAF.

    Writes:
      - ``{output_prefix}.matrix.tsv.gz``: wide matrix (ref position + per-species base)
      - ``{output_prefix}.summary.tsv``: run counters
      - ``{output_prefix}.th.tsv``: compact threshold-style stats dict
    """
    if min_species < 2:
        raise ValueError("min_species must be >= 2")

    ref_species, ref_chrom = chr_label_to_ref(chr_label)
    species_order = load_species_prefixes(species_list_path)
    stats = BlockStats()
    seen_species: set[str] = set()
    species_cols = _build_species_columns(species_order, ref_species)

    matrix_path = f"{output_prefix}.matrix.tsv.gz"
    Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)

    with _open_text(maf_path) as handle, _open_binary_write(matrix_path) as out:
        if species_cols is not None:
            header = ["ref_chr", "ref_pos", "block_score", "n_species", *species_cols]
            out.write("\t".join(header) + "\n")

        for score, rows in _iter_blocks(handle):
            stats.blocks_total += 1
            if score is not None and score < min_score:
                stats.blocks_skipped_score += 1
                continue

            ref_rows = [r for r in rows if r.species == ref_species]
            if len(ref_rows) != 1:
                stats.blocks_skipped_no_ref += 1
                continue
            ref_row = ref_rows[0]
            if ref_row.chrom != ref_chrom:
                stats.blocks_skipped_no_ref += 1
                continue

            kept_rows = [ref_row]
            for row in rows:
                if row.species == ref_species:
                    continue
                if _keep_row(row, ref_chrom, filter_cross_chr):
                    kept_rows.append(row)
                else:
                    stats.rows_dropped_cross_chr += 1

            if len(kept_rows) < min_species:
                continue

            stats.blocks_processed += 1
            for row in kept_rows:
                seen_species.add(row.species)

            if species_cols is None:
                species_cols = _ordered_species(species_order, ref_species, seen_species)
                header = ["ref_chr", "ref_pos", "block_score", "n_species", *species_cols]
                out.write("\t".join(header) + "\n")

            max_len = max(len(r.seq) for r in kept_rows)
            positions = [r.start for r in kept_rows]
            ref_index = 0
            block_score = "" if score is None else f"{score:.6f}"

            for col in range(max_len):
                stats.columns_total += 1
                present: list[tuple[int, str]] = []
                ref_base = "-"
                for idx, row in enumerate(kept_rows):
                    base = row.seq[col] if col < len(row.seq) else "-"
                    if base != "-":
                        present.append((idx, base))
                        if idx == ref_index:
                            ref_base = base
                            ref_pos = positions[idx]
                    if base != "-":
                        positions[idx] += 1 if row.strand == "+" else -1

                if ref_base == "-":
                    continue
                if len(present) < min_species:
                    continue

                stats.columns_emitted += 1
                base_map = {kept_rows[idx].species: base for idx, base in present}
                row_vals = [
                    ref_chrom,
                    str(ref_pos),
                    block_score,
                    str(len(present)),
                    *(base_map.get(sp, "-") for sp in species_cols),
                ]
                out.write("\t".join(row_vals) + "\n")

    if species_cols is None:
        species_cols = _ordered_species(species_order, ref_species, seen_species)
        if stats.columns_emitted == 0:
            with _open_binary_write(matrix_path) as out:
                header = ["ref_chr", "ref_pos", "block_score", "n_species", *species_cols]
                out.write("\t".join(header) + "\n")

    summary = {
        "chr_label": chr_label,
        "ref_species": ref_species,
        "ref_chrom": ref_chrom,
        "maf_path": maf_path,
        "min_species": min_species,
        "min_score": min_score,
        "filter_cross_chr": filter_cross_chr,
        "blocks_total": stats.blocks_total,
        "blocks_processed": stats.blocks_processed,
        "blocks_skipped_score": stats.blocks_skipped_score,
        "blocks_skipped_no_ref": stats.blocks_skipped_no_ref,
        "rows_dropped_cross_chr": stats.rows_dropped_cross_chr,
        "columns_total": stats.columns_total,
        "columns_emitted": stats.columns_emitted,
        "species_in_output": len(species_cols),
    }
    save_df_to_tsv(
        __import__("pandas").DataFrame([summary]),
        f"{output_prefix}.summary.tsv",
    )
    save_thresholds(summary, f"{output_prefix}.th.tsv")
    return summary
