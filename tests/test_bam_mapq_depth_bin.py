"""Unit tests for bam_mapq_depth_bin helpers."""

from genetics.germplasm.sample.bam_mapq_depth_bin import (
    _cigar_ref_len,
    _logical_position,
    _mq_bin_index,
    logical_chromosome_lengths,
    TSV_COLUMNS,
)


def test_cigar_ref_len():
    assert _cigar_ref_len("100M") == 100
    assert _cigar_ref_len("50M30S20M") == 70
    assert _cigar_ref_len("10=5X2D80M") == 97


def test_mq_bin_index():
    assert _mq_bin_index(-1) == -1
    assert _mq_bin_index(0) == 0
    assert _mq_bin_index(9) == 0
    assert _mq_bin_index(10) == 1
    assert _mq_bin_index(59) == 5
    assert _mq_bin_index(60) == 5


def test_logical_position_chr1a():
    mapped = _logical_position(1, 1)
    assert mapped is not None
    logical, pos = mapped
    assert logical == "chr1A"
    assert pos >= 1


def test_logical_chromosome_lengths_cover_autosomes():
    lengths = logical_chromosome_lengths()
    assert "chr1A" in lengths
    assert "chr7D" in lengths
    assert all(v > 0 for v in lengths.values())


def test_tsv_columns_order():
    assert TSV_COLUMNS[-2:] == ("total_reads", "bin_mean_depth")
