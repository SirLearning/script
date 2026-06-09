"""Unit tests for infra.utils.io helpers."""

import os
import tempfile

import pytest

from infra.utils.io import detect_delimiter, load_df_generic, load_df_from_tsv


def test_detect_delimiter_tab():
    assert detect_delimiter("a\tb\tc") == "\t"


def test_detect_delimiter_comma():
    assert detect_delimiter("a,b,c") == ","


def test_detect_delimiter_whitespace():
    assert detect_delimiter("a  b  c") is None


def test_load_df_generic_tsv():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as handle:
        handle.write("col1\tcol2\n1\t2\n3\t4\n")
        path = handle.name
    try:
        df = load_df_generic(path)
        assert list(df.columns) == ["col1", "col2"]
        assert len(df) == 2
    finally:
        os.unlink(path)


def test_load_df_from_tsv_missing_file():
    with pytest.raises(FileNotFoundError):
        load_df_from_tsv("/nonexistent/path/file.tsv")


def test_load_df_generic_gz_tsv():
    import gzip

    with tempfile.NamedTemporaryFile(suffix=".tsv.gz", delete=False) as handle:
        path = handle.name
    try:
        with gzip.open(path, "wt") as handle:
            handle.write("x\ty\n10\t20\n")
        df = load_df_generic(path)
        assert list(df.columns) == ["x", "y"]
        assert df.iloc[0]["x"] == 10
    finally:
        os.unlink(path)
