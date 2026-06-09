"""Deprecated: use genetics.genomics.variant.genetic_map_plot."""

from genetics.wheat._deprecation import warn_deprecated_shim

warn_deprecated_shim("genetic_map", "genetics.genomics.variant.genetic_map_plot")


from genetics.genomics.variant.genetic_map_plot import plot_genetic_map  # noqa: F401
