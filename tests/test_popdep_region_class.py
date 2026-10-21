"""Region-class intervals and exclusive popdep site labels."""

from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.variant.popdep_region_class import (
    IntervalSet,
    REGION_CLASS_ALPHA,
    REGION_CLASS_COLORS,
    REGION_CLASS_HUE_ORDER,
    REGION_CLASS_LEGEND_LABELS,
    classify_positions,
    clip_linear_interval_to_segments,
    cnv_intervals_from_mapq_table,
    gc_fraction,
    gc_windows_from_sequence,
    te_gff_to_segment_intervals,
    write_segment_bed,
)


def test_region_class_hue_order_puts_normal_at_bottom():
    assert REGION_CLASS_HUE_ORDER[0] == "normal"
    assert REGION_CLASS_HUE_ORDER[-1] == "gap"
    assert REGION_CLASS_HUE_ORDER == ("normal", "repeat", "gc", "cnv", "gap")


def test_region_class_colors_grey_normal_orange_cnv():
    assert REGION_CLASS_COLORS["normal"] == "#bdbdbd"
    assert REGION_CLASS_COLORS["repeat"] == "#1b9e77"
    assert REGION_CLASS_COLORS["cnv"] == "#d95f02"


def test_region_class_point_alpha_original():
    assert REGION_CLASS_ALPHA["normal"] == 0.12
    assert REGION_CLASS_ALPHA["repeat"] == 0.28
    assert REGION_CLASS_ALPHA["gc"] == 0.65
    assert REGION_CLASS_ALPHA["cnv"] == 0.80
    assert REGION_CLASS_ALPHA["gap"] == 0.90


def test_region_class_legend_labels_uppercase_gc_cnv():
    assert REGION_CLASS_LEGEND_LABELS["gc"] == "GC"
    assert REGION_CLASS_LEGEND_LABELS["cnv"] == "CNV"


def test_clip_linear_interval_stays_on_chr1a_part1():
    clipped = clip_linear_interval_to_segments("chr1A", 1_437_363, 1_437_463)
    assert clipped == [(1, 1_437_363, 1_437_463)]


def test_clip_linear_interval_splits_across_chr1a_parts():
    # part1 ends at 471304005; span 200 bp across the join.
    clipped = clip_linear_interval_to_segments("chr1A", 471_303_905, 471_304_105)
    assert clipped[0] == (1, 471_303_905, 471_304_005)
    assert clipped[1][0] == 2
    assert clipped[1][1] == 0
    assert clipped[1][2] == 100


def test_gc_fraction_and_extreme_windows():
    assert gc_fraction("ATATAT") == 0.0
    assert gc_fraction("GCGCGC") == 1.0
    assert gc_fraction("NNNN") is None
    high = gc_windows_from_sequence("GCGCGCGCGC", window_bp=10, gc_high=0.55, gc_low=0.35)
    low = gc_windows_from_sequence("ATATATATAT", window_bp=10, gc_high=0.55, gc_low=0.35)
    mid = gc_windows_from_sequence("ACGTACGTAC", window_bp=10, gc_high=0.55, gc_low=0.35)
    assert high == [(0, 10)]
    assert low == [(0, 10)]
    assert mid == []


def test_classify_positions_priority_gap_over_repeat():
    gap = IntervalSet.from_interval_map({1: [(10, 20)]})
    repeat = IntervalSet.from_interval_map({1: [(0, 50)]})
    gc = IntervalSet.from_interval_map({1: [(0, 50)]})
    cnv = IntervalSet.from_interval_map({1: [(0, 50)]})
    labels = classify_positions(
        1,
        [5, 15, 60],
        {"gap": gap, "repeat": repeat, "gc": gc, "cnv": cnv},
    )
    assert list(labels) == ["repeat", "gap", "normal"]


def test_cnv_intervals_from_mapq_table_threshold():
    table = pd.DataFrame(
        {
            "logical_chr": ["chr1A", "chr1A"],
            "bin_start": [1, 1_000_001],
            "bin_end": [1_000_000, 2_000_000],
            "cnv_frac": [0.25, 0.05],
        }
    )
    by_chr = cnv_intervals_from_mapq_table(table, cnv_frac_min=0.2)
    assert 1 in by_chr
    assert by_chr[1] == [(0, 1_000_000)]


def test_te_gff_to_segment_intervals_merges(tmp_path: Path):
    gff = tmp_path / "te.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1A\tCLARITE\trepeat_region\t11\t20\t.\t+\t.\tID=r1\n"
        "chr1A\tCLARITE\tnested_repeat\t18\t30\t.\t+\t.\tID=r2\n"
        "chr1A\tCLARITE\trepeat_fragment\t18\t22\t.\t+\t.\tParent=r2\n",
        encoding="utf-8",
    )
    by_chr = te_gff_to_segment_intervals(gff)
    merged = IntervalSet.from_interval_map(by_chr)
    hits = merged.contains_many(1, np.array([9, 10, 29, 30], dtype=np.int64))
    assert list(hits) == [False, True, True, False]


def test_write_and_contains_roundtrip(tmp_path: Path):
    bed = tmp_path / "gap.seg.bed"
    n = write_segment_bed({1: [(10, 20), (15, 25)]}, bed)
    assert n == 1
    text = bed.read_text(encoding="utf-8")
    assert "1\t10\t25\n" in text
