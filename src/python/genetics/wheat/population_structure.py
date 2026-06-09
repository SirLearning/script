"""Deprecated: use genetics.genomics.sample.pca_structure."""

from genetics.wheat._deprecation import warn_deprecated_shim

warn_deprecated_shim("population_structure", "genetics.genomics.sample.pca_structure")


from genetics.genomics.sample.pca_structure import plot_population_structure  # noqa: F401
