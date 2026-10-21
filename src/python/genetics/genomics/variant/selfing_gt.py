"""Rewrite unphased VCF genotypes to a selfing-style phased encoding.

Wheat is highly inbred, so converting ``/`` to ``|`` in the GT field is a
haplotype approximation: homozygotes are unambiguous; heterozygotes keep the
listed allele order. Missing ``./.`` becomes ``.|.``. Other VCF columns are
unchanged.
"""

from __future__ import annotations

import gzip
from pathlib import Path


def _open_text(path: str | Path, mode: str):
    path = Path(path)
    if path.name.endswith(".gz"):
        return gzip.open(path, mode + "t", encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def _phase_gt_field(gt: str) -> str:
    """Replace unphased allele separators in one GT string."""
    if not gt or gt == ".":
        return gt
    return gt.replace("/", "|")


def _rewrite_sample_gt(sample: str, gt_index: int) -> str:
    parts = sample.split(":")
    if gt_index >= len(parts):
        return sample
    parts[gt_index] = _phase_gt_field(parts[gt_index])
    return ":".join(parts)


def rewrite_gt_slash_to_pipe(src: str | Path, dst: str | Path) -> int:
    """Copy a VCF, rewriting GT ``/`` to ``|``.

    Returns the number of variant rows written.
    """
    src = Path(src)
    dst = Path(dst)
    n_sites = 0
    with _open_text(src, "r") as incoming, _open_text(dst, "w") as outgoing:
        for line in incoming:
            if line.startswith("#"):
                outgoing.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                outgoing.write(line)
                n_sites += 1
                continue
            fmt = fields[8].split(":")
            try:
                gt_index = fmt.index("GT")
            except ValueError:
                outgoing.write(line)
                n_sites += 1
                continue
            fields[9:] = [_rewrite_sample_gt(sample, gt_index) for sample in fields[9:]]
            outgoing.write("\t".join(fields) + "\n")
            n_sites += 1
    return n_sites
