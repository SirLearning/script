"""Deprecated: use genetics.genomics.variant.cnv_plot."""

from genetics.wheat._deprecation import warn_deprecated_shim

warn_deprecated_shim("cnv", "genetics.genomics.variant.cnv_plot")


from genetics.genomics.variant.cnv_plot import plot_cnv_results  # noqa: F401
