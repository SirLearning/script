"""Exclusive genomic-context classes for popdep site panels.

Classes (first match wins): gap, repeat, gc, cnv, else normal.

IWGSC RefSeq v1.0 ships a CLARITE TE GFF (repeat) but not GC or CNV tracks.
GC windows are computed from the v1.0 FASTA; CNV intervals come from the
published 1 Mb mapq_bin_popdep grid (cohort CNV fraction).
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.gene.gff_bed import _merge_half_open, collect_merged_gff_intervals
from infra.utils.io import save_df_to_tsv
from infra.wheat.ref_v1 import (
    _ref_segment_intervals_for_chromosome,
    get_chr_id_length,
)

REGION_CLASSES: tuple[str, ...] = ("gap", "repeat", "gc", "cnv", "normal")
CLASS_PRIORITY: tuple[str, ...] = ("gap", "repeat", "gc", "cnv")

REGION_CLASS_COLORS: dict[str, str] = {
    "normal": "#bdbdbd",
    "repeat": "#1b9e77",
    "gc": "#7570b3",
    "cnv": "#d95f02",
    "gap": "#252525",
}
REGION_CLASS_ALPHA: dict[str, float] = {
    "normal": 0.12,
    "repeat": 0.28,
    "gc": 0.65,
    "cnv": 0.80,
    "gap": 0.90,
}
# Last hue_order level is on top of scatter and stacked marginals.
REGION_CLASS_HUE_ORDER: tuple[str, ...] = ("normal", "repeat", "gc", "cnv", "gap")
REGION_CLASS_LEGEND_LABELS: dict[str, str] = {
    "normal": "normal",
    "repeat": "repeat",
    "gc": "GC",
    "cnv": "CNV",
    "gap": "gap",
}

REPEAT_GFF_TYPES: tuple[str, ...] = (
    "repeat_region",
    "nested_repeat",
    "repeat_fragment",
)

GC_WINDOW_BP = 1000
GC_HIGH = 0.55
GC_LOW = 0.35
CNV_FRAC_MIN = 0.20

BED_CLASSES: tuple[str, ...] = ("gap", "repeat", "gc", "cnv")


class IntervalSet:
    """Merged 0-based half-open intervals keyed by PLINK segment chr_id."""

    def __init__(self, arrays: dict[int, tuple[np.ndarray, np.ndarray]]):
        self.arrays = arrays

    @classmethod
    def from_interval_map(
        cls, by_chr: dict[int, list[tuple[int, int]]]
    ) -> IntervalSet:
        arrays: dict[int, tuple[np.ndarray, np.ndarray]] = {}
        for chr_id, raw in by_chr.items():
            merged = _merge_half_open(list(raw))
            if not merged:
                continue
            starts = np.fromiter((s for s, _ in merged), dtype=np.int64, count=len(merged))
            ends = np.fromiter((e for _, e in merged), dtype=np.int64, count=len(merged))
            arrays[int(chr_id)] = (starts, ends)
        return cls(arrays)

    def contains_many(self, chr_id: int, pos0: np.ndarray) -> np.ndarray:
        pos0 = np.asarray(pos0, dtype=np.int64)
        pair = self.arrays.get(int(chr_id))
        out = np.zeros(pos0.shape[0], dtype=bool)
        if pair is None:
            return out
        starts, ends = pair
        if starts.size == 0:
            return out
        idx = np.searchsorted(starts, pos0, side="right") - 1
        valid = idx >= 0
        if not valid.any():
            return out
        ii = idx[valid]
        out[valid] = pos0[valid] < ends[ii]
        return out

    def n_intervals(self) -> int:
        return int(sum(starts.size for starts, _ in self.arrays.values()))


def clip_linear_interval_to_segments(
    ref_chrom: str,
    start0: int,
    end0: int,
) -> list[tuple[int, int, int]]:
    """Clip a linear 0-based interval to PLINK segments ``(chr_id, start0, end0)``."""
    if end0 <= start0:
        return []
    try:
        intervals = _ref_segment_intervals_for_chromosome(ref_chrom)
    except (ValueError, KeyError):
        return []
    if not intervals:
        return []
    out: list[tuple[int, int, int]] = []
    for plink_chr, meta in intervals:
        seg_start = int(meta["ref_start"])
        seg_end = int(meta["ref_end"])
        overlap_s = max(int(start0), seg_start)
        overlap_e = min(int(end0), seg_end)
        if overlap_e > overlap_s:
            out.append((int(plink_chr), overlap_s - seg_start, overlap_e - seg_start))
    return out


def _add_linear_interval(
    by_chr: dict[int, list[tuple[int, int]]],
    ref_chrom: str,
    start0: int,
    end0: int,
) -> None:
    for chr_id, local_s, local_e in clip_linear_interval_to_segments(
        ref_chrom, start0, end0
    ):
        by_chr[int(chr_id)].append((int(local_s), int(local_e)))


def load_linear_bed_as_segments(path: str | Path) -> dict[int, list[tuple[int, int]]]:
    """Load a chr1A-style BED (0-based) and clip to PLINK segments."""
    by_chr: dict[int, list[tuple[int, int]]] = defaultdict(list)
    path = Path(path)
    opener = _open_text(path)
    with opener as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start0 = int(fields[1])
            end0 = int(fields[2])
            _add_linear_interval(by_chr, chrom, start0, end0)
    return dict(by_chr)


def load_segment_bed(path: str | Path) -> dict[int, list[tuple[int, int]]]:
    """Load a BED whose first column is a PLINK chr_id (int or chrNNN)."""
    by_chr: dict[int, list[tuple[int, int]]] = defaultdict(list)
    path = Path(path)
    with _open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chr_id = _parse_segment_chrom(fields[0])
            if chr_id is None:
                continue
            by_chr[chr_id].append((int(fields[1]), int(fields[2])))
    return dict(by_chr)


def write_segment_bed(
    by_chr: dict[int, list[tuple[int, int]]],
    path: str | Path,
) -> int:
    """Write merged segment intervals; returns interval count."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with path.open("w", encoding="utf-8") as handle:
        handle.write("#chr_id\tstart\tend\n")
        for chr_id in sorted(by_chr):
            for start, end in _merge_half_open(list(by_chr[chr_id])):
                handle.write(f"{chr_id}\t{start}\t{end}\n")
                n_written += 1
    return n_written


def te_gff_to_segment_intervals(
    gff_path: str | Path,
    *,
    feature_types: tuple[str, ...] = REPEAT_GFF_TYPES,
) -> dict[int, list[tuple[int, int]]]:
    """Merge IWGSC CLARITE TE features and clip to PLINK segments."""
    wanted = frozenset(feature_types)
    linear = collect_merged_gff_intervals(gff_path, feature_types=tuple(wanted))
    by_chr: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for chrom, intervals in linear.items():
        for start0, end0 in intervals:
            _add_linear_interval(by_chr, chrom, start0, end0)
    return dict(by_chr)


def gc_windows_from_sequence(
    seq: str,
    *,
    window_bp: int = GC_WINDOW_BP,
    gc_high: float = GC_HIGH,
    gc_low: float = GC_LOW,
    min_called_frac: float = 0.5,
) -> list[tuple[int, int]]:
    """0-based windows whose GC fraction is ≥ high or ≤ low (N ignored)."""
    if window_bp <= 0:
        raise ValueError("window_bp must be positive")
    if not seq:
        return []
    arr = np.frombuffer(seq.encode("ascii", "ignore"), dtype=np.uint8).copy()
    lower = (arr >= 97) & (arr <= 122)
    arr[lower] -= 32
    n = int(arr.size)
    out: list[tuple[int, int]] = []
    for start in range(0, n, window_bp):
        end = min(start + window_bp, n)
        window = arr[start:end]
        is_a = window == 65
        is_c = window == 67
        is_g = window == 71
        is_t = window == 84
        acgt = int(is_a.sum() + is_c.sum() + is_g.sum() + is_t.sum())
        width = int(end - start)
        if width == 0 or acgt == 0 or acgt / width < min_called_frac:
            continue
        gc = (int(is_g.sum()) + int(is_c.sum())) / acgt
        if gc >= gc_high or gc <= gc_low:
            out.append((start, end))
    return out


def gc_fraction(seq: str) -> float | None:
    """GC / (A+C+G+T); None when the window has no ACGT bases."""
    if not seq:
        return None
    arr = np.frombuffer(seq.encode("ascii", "ignore"), dtype=np.uint8).copy()
    lower = (arr >= 97) & (arr <= 122)
    arr[lower] -= 32
    acgt = int(((arr == 65) | (arr == 67) | (arr == 71) | (arr == 84)).sum())
    if acgt == 0:
        return None
    return int(((arr == 71) | (arr == 67)).sum()) / acgt


def gc_windows_from_fasta(
    fasta_path: str | Path,
    *,
    window_bp: int = GC_WINDOW_BP,
    gc_high: float = GC_HIGH,
    gc_low: float = GC_LOW,
    min_called_frac: float = 0.5,
) -> dict[int, list[tuple[int, int]]]:
    """Extreme-GC windows on FASTA contigs named 0..44 (vmap4 ABD)."""
    by_chr: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for chrom, seq in iter_fasta_records(fasta_path):
        chr_id = _parse_segment_chrom(chrom)
        if chr_id is None:
            continue
        length = get_chr_id_length(chr_id)
        if length is not None and len(seq) != int(length):
            print(
                f"[Warning] FASTA contig {chrom} length {len(seq)} "
                f"!= ref_v1 {length}"
            )
        by_chr[chr_id].extend(
            gc_windows_from_sequence(
                seq,
                window_bp=window_bp,
                gc_high=gc_high,
                gc_low=gc_low,
                min_called_frac=min_called_frac,
            )
        )
    return dict(by_chr)


def cnv_intervals_from_mapq_table(
    table: pd.DataFrame,
    *,
    cnv_frac_min: float = CNV_FRAC_MIN,
) -> dict[int, list[tuple[int, int]]]:
    """1 Mb bins with cohort CNV fraction ≥ threshold, clipped to segments.

    ``bin_start``/``bin_end`` are 1-based inclusive linear coordinates.
    """
    need = {"logical_chr", "bin_start", "bin_end", "cnv_frac"}
    missing = need - set(table.columns)
    if missing:
        raise ValueError(f"mapq popdep grid missing columns: {sorted(missing)}")
    by_chr: dict[int, list[tuple[int, int]]] = defaultdict(list)
    frac = pd.to_numeric(table["cnv_frac"], errors="coerce")
    keep = table.loc[frac >= float(cnv_frac_min)]
    for row in keep.itertuples(index=False):
        start1 = int(getattr(row, "bin_start"))
        end1 = int(getattr(row, "bin_end"))
        start0 = start1 - 1
        end0 = end1
        _add_linear_interval(by_chr, str(getattr(row, "logical_chr")), start0, end0)
    return dict(by_chr)


def cnv_intervals_from_mapq_dir(
    mapq_dir: str | Path,
    *,
    prefix: str = "allmapq_1mb",
    cnv_frac_min: float = CNV_FRAC_MIN,
) -> dict[int, list[tuple[int, int]]]:
    """Union CNV bins from ``{prefix}.chr*.popdep_bins.tsv`` under ``mapq_dir``."""
    root = Path(mapq_dir)
    paths = sorted(root.glob(f"{prefix}.chr*.popdep_bins.tsv"))
    if not paths:
        raise FileNotFoundError(
            f"No {prefix}.chr*.popdep_bins.tsv files in {root}"
        )
    by_chr: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for path in paths:
        table = pd.read_csv(path, sep="\t")
        part = cnv_intervals_from_mapq_table(table, cnv_frac_min=cnv_frac_min)
        for chr_id, intervals in part.items():
            by_chr[chr_id].extend(intervals)
    return dict(by_chr)


def classify_positions(
    chr_id: int,
    positions_1based,
    class_sets: dict[str, IntervalSet],
    *,
    priority: tuple[str, ...] = CLASS_PRIORITY,
) -> np.ndarray:
    """Exclusive class labels for 1-based positions on one PLINK segment."""
    pos0 = np.asarray(positions_1based, dtype=np.int64) - 1
    labels = np.full(pos0.shape[0], "normal", dtype=object)
    assigned = np.zeros(pos0.shape[0], dtype=bool)
    for cls in priority:
        interval_set = class_sets.get(cls)
        if interval_set is None:
            continue
        hit = interval_set.contains_many(int(chr_id), pos0)
        take = hit & ~assigned
        labels[take] = cls
        assigned |= take
    return labels


def load_class_interval_sets(annotation_dir: str | Path) -> dict[str, IntervalSet]:
    """Load ``{class}.seg.bed`` files written by ``build_region_class_beds``."""
    root = Path(annotation_dir)
    sets: dict[str, IntervalSet] = {}
    for cls in BED_CLASSES:
        path = root / f"{cls}.seg.bed"
        if not path.exists():
            raise FileNotFoundError(f"Missing region-class BED: {path}")
        sets[cls] = IntervalSet.from_interval_map(load_segment_bed(path))
    return sets


def build_region_class_beds(
    *,
    gap_bed: str | Path,
    te_gff: str | Path,
    fasta_path: str | Path,
    mapq_popdep_dir: str | Path,
    output_dir: str | Path,
    window_bp: int = GC_WINDOW_BP,
    gc_high: float = GC_HIGH,
    gc_low: float = GC_LOW,
    cnv_frac_min: float = CNV_FRAC_MIN,
    mapq_prefix: str = "allmapq_1mb",
) -> dict[str, int]:
    """Write four segment BEDs plus a one-row summary TSV."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    gap = load_linear_bed_as_segments(gap_bed)
    repeat = te_gff_to_segment_intervals(te_gff)
    gc = gc_windows_from_fasta(
        fasta_path, window_bp=window_bp, gc_high=gc_high, gc_low=gc_low
    )
    cnv = cnv_intervals_from_mapq_dir(
        mapq_popdep_dir, prefix=mapq_prefix, cnv_frac_min=cnv_frac_min
    )

    counts = {
        "gap": write_segment_bed(gap, out / "gap.seg.bed"),
        "repeat": write_segment_bed(repeat, out / "repeat.seg.bed"),
        "gc": write_segment_bed(gc, out / "gc.seg.bed"),
        "cnv": write_segment_bed(cnv, out / "cnv.seg.bed"),
    }
    save_df_to_tsv(pd.DataFrame([counts]), str(out / "region_class_bed.summary.tsv"))
    return counts


def iter_fasta_records(path: str | Path):
    """Yield ``(header_token, sequence)`` from a possibly gzipped FASTA."""
    seq_chunks: list[str] = []
    name: str | None = None
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_chunks)
                name = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if name is not None:
            yield name, "".join(seq_chunks)


def _open_text(path: Path):
    path = Path(path)
    if str(path).endswith(".gz"):
        import gzip

        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _parse_segment_chrom(raw: str) -> int | None:
    token = str(raw).strip()
    if token.startswith("chr") and token[3:].isdigit():
        return int(token[3:])
    if token.isdigit():
        return int(token)
    return None
