"""Deprecated: use genetics.genomics.plink.results_io."""

from genetics.wheat._deprecation import warn_deprecated_shim

warn_deprecated_shim("plink_results", "genetics.genomics.plink.results_io")


from genetics.genomics.plink.results_io import *  # noqa: F403
