"""breachr — MAPQ- and depth-aware segment caller for 1 Mb bin profiles.

Chromosome breaks and dosage shifts (CNV, inter-homoeologue shadow mapping,
aneuploidy) leave joint signatures in per-bin read depth and MAPQ histograms.
This package consumes published ``*.mapq_depth_bin.tsv`` tables only; it does
not scan BAMs or build bins.
"""

from breachr.caller import call_mapq_depth_bin_tsv, call_mapq_depth_bins
from breachr.config import CallConfig
from breachr.io import load_mapq_depth_bin_tsv
from breachr.reconcile import (
    ReconcileConfig,
    parse_abnormal_chr_list,
    reconcile_flagged_samples,
    write_reconcile_outputs,
)

__all__ = [
    "CallConfig",
    "ReconcileConfig",
    "call_mapq_depth_bins",
    "call_mapq_depth_bin_tsv",
    "load_mapq_depth_bin_tsv",
    "parse_abnormal_chr_list",
    "reconcile_flagged_samples",
    "write_reconcile_outputs",
]
