"""Deprecated: use genetics.genomics.variant.snp_qc_plot."""

from genetics.wheat._deprecation import warn_deprecated_shim

warn_deprecated_shim("snp_qc", "genetics.genomics.variant.snp_qc_plot")

from genetics.genomics.variant.snp_qc_plot import summarize_and_plot_snp_qc  # noqa: F401
