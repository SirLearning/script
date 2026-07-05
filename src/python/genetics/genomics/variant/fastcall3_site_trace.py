"""Read-only trace of MAC=0 sites through FastCall3 intermediates (.pos.txt, .ing.gz, .iac.gz)."""

from __future__ import annotations

import gzip
import logging
import struct
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import pandas as pd

from infra.utils.io import save_df_to_tsv

logger = logging.getLogger(__name__)

# Production run under 04runScreens/all uses 5 Mb bins for disc and scan file names.
DISC_BIN_SIZE = 5_000_000
ALLELE_NAMES: dict[int, str] = {0: "A", 1: "C", 2: "G", 3: "T", 4: "D", 5: "I"}
ING_TERMINATOR = 2_147_483_647
DEFAULT_CS_TAXA = ("CS_mp_2018_8X", "CS_sg_2014_3X", "CS_sg_2017_60X")


@dataclass(frozen=True)
class BinRange:
    chrom: int
    bin_start: int
    bin_end: int

    @property
    def filename_stem(self) -> str:
        return f"{self.chrom}_{self.bin_start}_{self.bin_end}"


@dataclass
class LibSiteRef:
    chrom: int
    position: int
    global_line: int
    bin_index: int
    bin: BinRange


@dataclass
class IacSiteRecord:
    taxon: str
    status: str
    ref_count: int | None = None
    alt_counts: list[int] = field(default_factory=list)
    allele_num: int | None = None


@dataclass
class IngSiteRecord:
    taxon: str
    found: bool
    allele_base: str | None = None
    indel_length: int | None = None
    allele_coding: int | None = None


def calc_bin_range(chrom: int, position: int, bin_size: int = DISC_BIN_SIZE) -> BinRange:
    """Return the disc/scan bin that contains ``position`` (1-based inclusive lower bound)."""
    bin_start = (position - 1) // bin_size * bin_size + 1
    bin_end = bin_start + bin_size
    return BinRange(chrom=chrom, bin_start=bin_start, bin_end=bin_end)


def decode_allele_pack(pack: int, bin_start: int) -> tuple[int, str, int]:
    rel_pos = pack >> 9
    coding = (pack >> 6) & 0x7
    indel_len = pack & 0x3F
    position = rel_pos + bin_start
    base = ALLELE_NAMES.get(coding, f"?{coding}")
    return position, base, indel_len


def pack_to_alt_string(pack: int) -> str:
    coding = (pack >> 6) & 0x7
    indel_len = pack & 0x3F
    base = ALLELE_NAMES.get(coding, f"?{coding}")
    if base == "I" and indel_len:
        return f"I{indel_len}"
    if base == "D" and indel_len:
        return f"-{indel_len}"
    return base


def _read_java_utf(handle) -> str:
    length = struct.unpack(">H", handle.read(2))[0]
    if length == 0:
        return ""
    return handle.read(length).decode("utf-8")


def resolve_pos_file(gen_dir: Path, chrom: int) -> Path:
    matches = sorted(gen_dir.glob(f"{chrom}_1_*.pos.txt"))
    if not matches:
        raise FileNotFoundError(f"No .pos.txt for chrom {chrom} under {gen_dir}")
    if len(matches) > 1:
        logger.warning("Multiple .pos.txt for chrom %s; using %s", chrom, matches[0])
    return matches[0]


def lookup_lib_site(pos_file: Path, chrom: int, position: int) -> LibSiteRef | None:
    """Locate ``position`` in a chromosome ``.pos.txt`` and derive its scan-bin index."""
    bin_range = calc_bin_range(chrom, position)
    global_line = 0
    bin_index = 0
    with pos_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            global_line += 1
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            pos = int(parts[1])
            if pos < bin_range.bin_start:
                continue
            if pos >= bin_range.bin_end:
                break
            if pos == position:
                return LibSiteRef(
                    chrom=chrom,
                    position=position,
                    global_line=global_line,
                    bin_index=bin_index,
                    bin=bin_range,
                )
            bin_index += 1
    return None


def read_ing_site(ing_path: Path, bin_start: int, target_pos: int) -> list[tuple[str, int]]:
    """Return decoded (base, indel_len) alt observations at ``target_pos`` in one .ing.gz."""
    hits: list[tuple[str, int]] = []
    with gzip.open(ing_path, "rb") as handle:
        _read_java_utf(handle)
        struct.unpack(">h", handle.read(2))
        file_bin_start = struct.unpack(">i", handle.read(4))[0]
        struct.unpack(">i", handle.read(4))
        if file_bin_start != bin_start:
            logger.warning("ing bin_start mismatch: %s expected %s", ing_path, bin_start)
        while True:
            pack = struct.unpack(">i", handle.read(4))[0]
            if pack == ING_TERMINATOR:
                break
            pos, base, indel_len = decode_allele_pack(pack, file_bin_start)
            if pos == target_pos:
                hits.append((base, indel_len))
    return hits


def read_iac_site(iac_path: Path, site_index: int) -> IacSiteRecord:
    taxon = ""
    with gzip.open(iac_path, "rb") as handle:
        taxon = _read_java_utf(handle)
        struct.unpack(">h", handle.read(2))
        struct.unpack(">i", handle.read(4))
        struct.unpack(">i", handle.read(4))
        position_num = struct.unpack(">i", handle.read(4))[0]
        if site_index >= position_num:
            return IacSiteRecord(taxon=taxon, status="index_out_of_range")
        for idx in range(position_num):
            allele_num = struct.unpack(">b", handle.read(1))[0]
            if idx != site_index:
                if allele_num >= 0:
                    handle.read(2 * allele_num)
                continue
            if allele_num < 0:
                return IacSiteRecord(taxon=taxon, status="missing", allele_num=allele_num)
            counts = list(struct.unpack(">" + "h" * allele_num, handle.read(2 * allele_num)))
            ref_count = counts[0] if counts else 0
            alt_counts = counts[1:] if len(counts) > 1 else []
            return IacSiteRecord(
                taxon=taxon,
                status="ok",
                ref_count=ref_count,
                alt_counts=alt_counts,
                allele_num=allele_num,
            )
    return IacSiteRecord(taxon=taxon, status="parse_error")


def classify_iac(record: IacSiteRecord) -> str:
    if record.status != "ok":
        return record.status
    ref_count = record.ref_count or 0
    alt_sum = sum(record.alt_counts)
    if ref_count > 0 and alt_sum == 0:
        return "ref_only"
    if ref_count == 0 and alt_sum > 0:
        return "alt_only"
    if ref_count > 0 and alt_sum > 0:
        return "mixed"
    return "zero_depth"


def _ing_worker(args: tuple[str, str, int, int]) -> dict:
    taxon, ing_path, bin_start, target_pos = args
    hits = read_ing_site(Path(ing_path), bin_start, target_pos)
    alleles = ",".join(f"{base}{indel_len}" if indel_len else base for base, indel_len in hits)
    return {"taxon": taxon, "ing_found": bool(hits), "ing_alleles": alleles}


def scan_ing_cohort(
    ing_root: Path,
    bin_range: BinRange,
    target_pos: int,
    max_workers: int = 32,
) -> pd.DataFrame:
    """Find all taxa whose disc ``.ing.gz`` records ``target_pos`` (read-only)."""
    bin_name = f"{bin_range.filename_stem}.ing.gz"
    tasks: list[tuple[str, str, int, int]] = []
    for taxon_dir in sorted(p for p in ing_root.iterdir() if p.is_dir()):
        ing_path = taxon_dir / bin_name
        if ing_path.is_file():
            tasks.append((taxon_dir.name, str(ing_path), bin_range.bin_start, target_pos))
    rows: list[dict] = []
    if not tasks:
        return pd.DataFrame(columns=["taxon", "ing_found", "ing_alleles"])
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(_ing_worker, task) for task in tasks]
        for fut in as_completed(futures):
            rows.append(fut.result())
    return pd.DataFrame(rows)


def read_iac_multi(iac_path: Path, site_indices: set[int]) -> dict[int, IacSiteRecord]:
    """Read several library indices from one ``.iac.gz`` in a single pass."""
    results: dict[int, IacSiteRecord] = {}
    taxon = ""
    with gzip.open(iac_path, "rb") as handle:
        taxon = _read_java_utf(handle)
        struct.unpack(">h", handle.read(2))
        struct.unpack(">i", handle.read(4))
        struct.unpack(">i", handle.read(4))
        position_num = struct.unpack(">i", handle.read(4))[0]
        wanted = {i for i in site_indices if 0 <= i < position_num}
        for idx in range(position_num):
            allele_num = struct.unpack(">b", handle.read(1))[0]
            if idx in wanted:
                if allele_num < 0:
                    results[idx] = IacSiteRecord(taxon=taxon, status="missing", allele_num=allele_num)
                else:
                    counts = list(struct.unpack(">" + "h" * allele_num, handle.read(2 * allele_num)))
                    results[idx] = IacSiteRecord(
                        taxon=taxon,
                        status="ok",
                        ref_count=counts[0] if counts else 0,
                        alt_counts=counts[1:] if len(counts) > 1 else [],
                        allele_num=allele_num,
                    )
            elif allele_num >= 0:
                handle.read(2 * allele_num)
    for idx in site_indices:
        if idx not in results:
            results[idx] = IacSiteRecord(taxon=taxon, status="index_out_of_range")
    return results


def _iac_multi_worker(args: tuple[str, str, list[int]]) -> list[dict]:
    taxon, iac_path, site_indices = args
    records = read_iac_multi(Path(iac_path), set(site_indices))
    out: list[dict] = []
    for idx, record in records.items():
        out.append(
            {
                "taxon": taxon,
                "iac_site_index": idx,
                "iac_status": classify_iac(record),
                "ref_count": record.ref_count,
                "alt_counts": ",".join(str(x) for x in record.alt_counts) if record.alt_counts else "",
                "allele_num": record.allele_num,
                "raw_status": record.status,
            }
        )
    return out


def scan_iac_cohort_bins(
    iac_root: Path,
    bin_range: BinRange,
    site_indices: list[int],
    max_workers: int = 32,
) -> pd.DataFrame:
    """Extract multiple library indices from the same scan bin across all taxa."""
    bin_name = f"{bin_range.filename_stem}.iac.gz"
    tasks: list[tuple[str, str, list[int]]] = []
    for taxon_dir in sorted(p for p in iac_root.iterdir() if p.is_dir()):
        iac_path = taxon_dir / bin_name
        if iac_path.is_file():
            tasks.append((taxon_dir.name, str(iac_path), site_indices))
    rows: list[dict] = []
    if not tasks:
        return pd.DataFrame(
            columns=[
                "taxon",
                "iac_site_index",
                "iac_status",
                "ref_count",
                "alt_counts",
                "allele_num",
                "raw_status",
            ]
        )
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(_iac_multi_worker, task) for task in tasks]
        for fut in as_completed(futures):
            rows.extend(fut.result())
    return pd.DataFrame(rows)


def scan_ing_for_taxa(
    ing_root: Path,
    bin_range: BinRange,
    target_pos: int,
    taxa: Iterable[str],
) -> pd.DataFrame:
    bin_name = f"{bin_range.filename_stem}.ing.gz"
    rows: list[dict] = []
    for taxon in taxa:
        ing_path = ing_root / taxon / bin_name
        if not ing_path.is_file():
            rows.append(
                {
                    "taxon": taxon,
                    "ing_file": "",
                    "ing_found": False,
                    "ing_alleles": "",
                }
            )
            continue
        hits = read_ing_site(ing_path, bin_range.bin_start, target_pos)
        alleles = ",".join(f"{base}{indel_len}" if indel_len else base for base, indel_len in hits)
        rows.append(
            {
                "taxon": taxon,
                "ing_file": str(ing_path),
                "ing_found": bool(hits),
                "ing_alleles": alleles,
            }
        )
    return pd.DataFrame(rows)


def _select_taxa_for_ing_detail(iac_df: pd.DataFrame, cs_taxa: tuple[str, ...], max_extra: int = 15) -> list[str]:
    selected = list(cs_taxa)
    if iac_df.empty:
        return selected
    for category in ("alt_only", "mixed", "ref_only"):
        subset = iac_df.loc[iac_df["iac_status"] == category, "taxon"].head(max_extra)
        for taxon in subset:
            if taxon not in selected:
                selected.append(taxon)
    return selected


def trace_mac_zero_sites(
    sites_table: pd.DataFrame,
    gen_dir: str | Path,
    ing_dir: str | Path,
    output_prefix: str,
    cs_taxa: tuple[str, ...] = DEFAULT_CS_TAXA,
    max_workers: int = 32,
) -> dict[str, pd.DataFrame]:
    """
    Trace MAC=0 sites through FastCall3 ``.pos.txt``, ``.ing.gz``, and ``.iac.gz``.

    Uses ``gen/*_1_*.pos.txt`` for library membership (same order as scan ``.iac.gz``).
    vlib is not required. All inputs are read-only.

    Writes:
      - ``{output_prefix}.summary.tsv`` — cohort counts per site
      - ``{output_prefix}.iac_by_taxa.tsv`` — scan read counts for every taxon
      - ``{output_prefix}.iac_cs.tsv`` — CS scan slice
      - ``{output_prefix}.ing_disc_by_taxa.tsv`` — disc discoveries across all taxa
      - ``{output_prefix}.ing_detail.tsv`` — disc check for CS + exemplar taxa
    """
    gen_dir = Path(gen_dir)
    ing_dir = Path(ing_dir)
    iac_root = gen_dir / "indiCounts"
    if not iac_root.is_dir():
        raise FileNotFoundError(f"Missing indiCounts under {gen_dir}")

    pos_cache: dict[int, Path] = {}
    site_specs: list[dict] = []
    summary_rows: list[dict] = []
    taxa_rows: list[dict] = []
    ing_disc_rows: list[dict] = []
    ing_rows: list[dict] = []

    for row_dict in sites_table.to_dict(orient="records"):
        chrom = int(row_dict.get("CHROM", row_dict.get("ref_id")))
        if "ref_pos" in row_dict and pd.notna(row_dict["ref_pos"]):
            position = int(row_dict["ref_pos"])
        else:
            position = int(str(row_dict["ID"]).split("-", 1)[1])

        if chrom not in pos_cache:
            pos_cache[chrom] = resolve_pos_file(gen_dir, chrom)
        lib_ref = lookup_lib_site(pos_cache[chrom], chrom, position)
        site_key = row_dict.get("ID", f"{chrom}-{position}")
        if lib_ref is None:
            summary_rows.append(
                {
                    **{k: row_dict.get(k) for k in ("subgenome", "ID", "REF", "ALT", "mac0_category")},
                    "CHROM": chrom,
                    "POS": position,
                    "in_lib_pos_txt": False,
                    "note": "not_in_pos_txt",
                }
            )
            continue
        site_specs.append(
            {
                "row": row_dict,
                "site_key": site_key,
                "chrom": chrom,
                "position": position,
                "lib_ref": lib_ref,
            }
        )

    # Batch scan .iac.gz by bin (one pass per taxon per bin for all site indices in that bin).
    bin_groups: dict[str, list[dict]] = {}
    for spec in site_specs:
        stem = spec["lib_ref"].bin.filename_stem
        bin_groups.setdefault(stem, []).append(spec)

    iac_by_site: dict[str, pd.DataFrame] = {}
    for stem, specs in bin_groups.items():
        bin_range = specs[0]["lib_ref"].bin
        indices = sorted({s["lib_ref"].bin_index for s in specs})
        iac_bin_df = scan_iac_cohort_bins(iac_root, bin_range, indices, max_workers=max_workers)
        for spec in specs:
            site_key = spec["site_key"]
            idx = spec["lib_ref"].bin_index
            site_df = iac_bin_df.loc[iac_bin_df["iac_site_index"] == idx].copy()
            site_df["site_id"] = site_key
            site_df["CHROM"] = spec["chrom"]
            site_df["POS"] = spec["position"]
            iac_by_site[site_key] = site_df
            taxa_rows.extend(site_df.to_dict(orient="records"))

    for spec in site_specs:
        row_dict = spec["row"]
        site_key = spec["site_key"]
        chrom = spec["chrom"]
        position = spec["position"]
        lib_ref = spec["lib_ref"]
        iac_df = iac_by_site.get(site_key, pd.DataFrame())
        status_counts = iac_df["iac_status"].value_counts().to_dict() if not iac_df.empty else {}

        ing_cohort_df = scan_ing_cohort(ing_dir, lib_ref.bin, position, max_workers=max_workers)
        disc_hits = ing_cohort_df.loc[ing_cohort_df["ing_found"]].copy() if not ing_cohort_df.empty else ing_cohort_df
        for _, ing_row in disc_hits.iterrows():
            ing_disc_rows.append(
                {"site_id": site_key, "CHROM": chrom, "POS": position, **ing_row.to_dict()}
            )

        ing_detail_taxa = _select_taxa_for_ing_detail(iac_df, cs_taxa)
        ing_df = scan_ing_for_taxa(ing_dir, lib_ref.bin, position, ing_detail_taxa)
        for _, ing_row in ing_df.iterrows():
            ing_rows.append({"site_id": site_key, "CHROM": chrom, "POS": position, **ing_row.to_dict()})

        summary_rows.append(
            {
                **{k: row_dict.get(k) for k in ("subgenome", "ID", "REF", "ALT", "mac0_category")},
                "CHROM": chrom,
                "POS": position,
                "in_lib_pos_txt": True,
                "pos_txt_line": lib_ref.global_line,
                "iac_bin": lib_ref.bin.filename_stem,
                "iac_site_index": lib_ref.bin_index,
                "iac_n_taxa": len(iac_df),
                "iac_missing": status_counts.get("missing", 0),
                "iac_ref_only": status_counts.get("ref_only", 0),
                "iac_alt_only": status_counts.get("alt_only", 0),
                "iac_mixed": status_counts.get("mixed", 0),
                "iac_zero_depth": status_counts.get("zero_depth", 0),
                "disc_ing_taxa_n": len(disc_hits),
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    taxa_df = pd.DataFrame(taxa_rows)
    ing_disc_df = pd.DataFrame(ing_disc_rows)
    ing_df_all = pd.DataFrame(ing_rows)

    save_df_to_tsv(summary_df, f"{output_prefix}.summary.tsv")
    save_df_to_tsv(taxa_df, f"{output_prefix}.iac_by_taxa.tsv")
    save_df_to_tsv(ing_disc_df, f"{output_prefix}.ing_disc_by_taxa.tsv")
    save_df_to_tsv(ing_df_all, f"{output_prefix}.ing_detail.tsv")

    cs_mask = taxa_df["taxon"].isin(cs_taxa) if not taxa_df.empty else pd.Series(dtype=bool)
    save_df_to_tsv(taxa_df.loc[cs_mask].copy(), f"{output_prefix}.iac_cs.tsv")

    logger.info("Wrote FastCall3 trace tables with prefix %s", output_prefix)
    return {
        "summary": summary_df,
        "iac_by_taxa": taxa_df,
        "ing_disc_by_taxa": ing_disc_df,
        "ing_detail": ing_df_all,
    }


def load_mac_zero_site_table(path: str | Path) -> pd.DataFrame:
    """Load a MAC=0 site table TSV (e.g. from ``test_thin.mac0_sites.all.tsv``)."""
    return pd.read_csv(path, sep="\t")
