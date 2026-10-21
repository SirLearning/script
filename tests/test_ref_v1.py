"""IWGSC RefSeq v1.0 numeric chr ids vs published parts names and split points."""

import pytest

from infra.wheat.ref_v1 import (
    REF_V1_CHR_ID_TO_PARTS_NAME,
    REF_V1_PARTS_NAME_TO_CHR_ID,
    get_ref_v1_chr_id_from_parts_name,
    get_ref_v1_parts_name,
    get_ref_v1_plink_chr_length,
    parts_name_position_to_ref,
    plink_segment_to_ref_position,
)

# Official IWGSC 2018 parts_to_chr.bed split points (0-based half-open on ori chr).
_IWGSC_PARTS_SPLITS = {
    "chr1A_part1": ("chr1A", 0, 471304005),
    "chr1A_part2": ("chr1A", 471304005, 594102056),
    "chr4D_part1": ("chr4D", 0, 451004620),
    "chr4D_part2": ("chr4D", 451004620, 509857067),
    "chr6D_part1": ("chr6D", 0, 450509124),
    "chr6D_part2": ("chr6D", 450509124, 473592718),
    "chrUn": ("chrUn", 0, 480980714),
}


def test_parts_names_are_unique_for_all_chr_ids():
    assert len(REF_V1_CHR_ID_TO_PARTS_NAME) == 45
    assert set(REF_V1_CHR_ID_TO_PARTS_NAME) == set(range(45))
    assert len(set(REF_V1_CHR_ID_TO_PARTS_NAME.values())) == 45
    assert REF_V1_PARTS_NAME_TO_CHR_ID["chr1A_part1"] == 1
    assert REF_V1_PARTS_NAME_TO_CHR_ID["chr1A_part2"] == 2


def test_numeric_id_roundtrips_with_published_parts_name():
    assert get_ref_v1_parts_name(0) == "chrUn"
    assert get_ref_v1_parts_name(1) == "chr1A_part1"
    assert get_ref_v1_parts_name("2") == "chr1A_part2"
    assert get_ref_v1_parts_name(43) == "Mit"
    assert get_ref_v1_parts_name(44) == "Chl"
    assert get_ref_v1_chr_id_from_parts_name("chr1A_part1") == 1
    assert get_ref_v1_chr_id_from_parts_name("1A_part1") == 1
    assert get_ref_v1_chr_id_from_parts_name("chrUn") == 0
    assert get_ref_v1_chr_id_from_parts_name("Mit") == 43
    for cid in range(45):
        assert get_ref_v1_chr_id_from_parts_name(get_ref_v1_parts_name(cid)) == cid


def test_unknown_parts_name_raises():
    with pytest.raises(ValueError, match="Unknown"):
        get_ref_v1_chr_id_from_parts_name("chr1A")
    with pytest.raises(ValueError, match="Unknown"):
        get_ref_v1_parts_name(45)


def test_split_lengths_match_iwgsc_2018_parts_bed():
    for parts_name, (ori_chr, ori_start, ori_end) in _IWGSC_PARTS_SPLITS.items():
        cid = get_ref_v1_chr_id_from_parts_name(parts_name)
        assert get_ref_v1_plink_chr_length(cid) == ori_end - ori_start
        local_end = ori_end - ori_start
        chrom, pos = parts_name_position_to_ref(parts_name, local_end)
        assert f"chr{chrom}" == ori_chr or chrom == ori_chr.replace("chr", "")
        assert pos == ori_end
        if ori_start == 0:
            assert parts_name_position_to_ref(parts_name, 1)[1] == 1
        else:
            assert parts_name_position_to_ref(parts_name, 1)[1] == ori_start + 1


def test_parts_local_pos_matches_plink_segment_mapping():
    assert parts_name_position_to_ref("chr1A_part2", 10) == plink_segment_to_ref_position(2, 10)
    assert parts_name_position_to_ref("chr1A_part2", 10) == ("1A", 471304005 + 10)
