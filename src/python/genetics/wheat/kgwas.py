"""Deprecated: use genetics.gwas.association_plot."""

from genetics.wheat._deprecation import warn_deprecated_shim

warn_deprecated_shim("kgwas", "genetics.gwas.association_plot")


from genetics.gwas.association_plot import plot_kgwas_results  # noqa: F401
