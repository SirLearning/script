"""Tests for FastCall3 intermediate file decoders (read-only trace helpers)."""

from __future__ import annotations

import gzip
import struct
from pathlib import Path

import pandas as pd

from genetics.genomics.variant.fastcall3_site_trace import (
    calc_bin_range,
    decode_allele_pack,
    lookup_lib_site,
    read_iac_site,
    read_ing_site,
    trace_mac_zero_sites,
)

ING_TERMINATOR = 2_147_483_647


def _write_java_utf(handle, text: str) -> None:
    encoded = text.encode("utf-8")
    handle.write(struct.pack(">H", len(encoded)))
    handle.write(encoded)


def test_calc_bin_range():
    b = calc_bin_range(32, 105_721_177)
    assert b.bin_start == 105_000_001
    assert b.bin_end == 110_000_001
    assert b.filename_stem == "32_105000001_110000001"


def test_decode_allele_pack():
    bin_start = 105_000_001
    rel = 105_721_177 - bin_start
    pack = (rel << 9) + (5 << 6) + 1  # IUPAC I with indel len 1
    pos, base, indel = decode_allele_pack(pack, bin_start)
    assert pos == 105_721_177
    assert base == "I"
    assert indel == 1


def test_lookup_lib_site(tmp_path: Path):
    pos_file = tmp_path / "2_1_100.pos.txt"
    pos_file.write_text("2\t10\n2\t20\n2\t30\n2\t40\n", encoding="utf-8")
    ref = lookup_lib_site(pos_file, 2, 30)
    assert ref is not None
    assert ref.bin_index == 2
    assert ref.global_line == 3


def test_read_ing_and_iac_roundtrip(tmp_path: Path):
    bin_start = 100
    bin_end = 200
    target_pos = 150
    rel = target_pos - bin_start
    pack = (rel << 9) + (2 << 6) + 0  # G SNP

    ing_path = tmp_path / "taxon1.ing.gz"
    with gzip.open(ing_path, "wb") as handle:
        _write_java_utf(handle, "taxon1")
        handle.write(struct.pack(">h", 2))
        handle.write(struct.pack(">ii", bin_start, bin_end))
        handle.write(struct.pack(">i", pack))
        handle.write(struct.pack(">i", ING_TERMINATOR))

    hits = read_ing_site(ing_path, bin_start, target_pos)
    assert hits == [("G", 0)]

    iac_path = tmp_path / "taxon1.iac.gz"
    with gzip.open(iac_path, "wb") as handle:
        _write_java_utf(handle, "taxon1")
        handle.write(struct.pack(">h", 2))
        handle.write(struct.pack(">ii", bin_start, bin_end))
        handle.write(struct.pack(">i", 3))  # positionNum
        handle.write(struct.pack(">b", -1))  # missing
        handle.write(struct.pack(">b", 2))
        handle.write(struct.pack(">hh", 10, 0))  # ref_only at index 1
        handle.write(struct.pack(">b", 2))
        handle.write(struct.pack(">hh", 0, 8))  # alt_only at index 2

    missing = read_iac_site(iac_path, 0)
    assert missing.status == "missing"
    ref_only = read_iac_site(iac_path, 1)
    assert ref_only.ref_count == 10
    assert ref_only.alt_counts == [0]
    alt_only = read_iac_site(iac_path, 2)
    assert alt_only.ref_count == 0
    assert alt_only.alt_counts == [8]


def test_trace_mac_zero_sites_minimal(tmp_path: Path):
    gen_dir = tmp_path / "gen"
    ing_dir = tmp_path / "ing"
    iac_root = gen_dir / "indiCounts"
    iac_root.mkdir(parents=True)
    ing_dir.mkdir()

    pos_file = gen_dir / "2_1_1000.pos.txt"
    pos_file.write_text("2\t100\n2\t150\n2\t160\n", encoding="utf-8")

    bin_range = calc_bin_range(2, 150)
    bin_name = f"{bin_range.filename_stem}.iac.gz"
    ing_name = f"{bin_range.filename_stem}.ing.gz"

    rel = 150 - bin_range.bin_start
    pack = (rel << 9) + (3 << 6)

    for taxon in ("CS_sg_2017_60X", "sampleA"):
        taxon_iac = iac_root / taxon
        taxon_iac.mkdir()
        taxon_ing = ing_dir / taxon
        taxon_ing.mkdir()
        with gzip.open(taxon_iac / bin_name, "wb") as handle:
            _write_java_utf(handle, taxon)
            handle.write(struct.pack(">h", 2))
            handle.write(struct.pack(">iii", bin_range.bin_start, bin_range.bin_end, 3))
            handle.write(struct.pack(">b", -1))  # pos 100
            if taxon == "sampleA":
                handle.write(struct.pack(">b", 2))
                handle.write(struct.pack(">hh", 0, 5))  # pos 150 alt-only
            else:
                handle.write(struct.pack(">b", -1))  # pos 150 missing for CS
            handle.write(struct.pack(">b", -1))  # pos 160
        with gzip.open(taxon_ing / ing_name, "wb") as handle:
            _write_java_utf(handle, taxon)
            handle.write(struct.pack(">h", 2))
            handle.write(struct.pack(">ii", bin_range.bin_start, bin_range.bin_end))
            if taxon == "sampleA":
                handle.write(struct.pack(">i", pack))
            handle.write(struct.pack(">i", ING_TERMINATOR))

    sites = pd.DataFrame(
        [
            {
                "subgenome": "A",
                "ID": "2-150",
                "CHROM": 2,
                "ref_pos": 150,
                "REF": "G",
                "ALT": "T",
                "mac0_category": "ref_absent_hom_or_hap_alt_only",
            }
        ]
    )
    out_prefix = str(tmp_path / "trace")
    result = trace_mac_zero_sites(sites, gen_dir, ing_dir, out_prefix, max_workers=2)
    summary = result["summary"]
    assert bool(summary.loc[0, "in_lib_pos_txt"])
    assert summary.loc[0, "iac_alt_only"] == 1
    assert summary.loc[0, "iac_missing"] == 1
    assert summary.loc[0, "disc_ing_taxa_n"] == 1
