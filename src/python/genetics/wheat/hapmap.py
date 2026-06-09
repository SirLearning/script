"""Deprecated: use genetics.genomics.variant.hapmap_report."""

from genetics.wheat._deprecation import warn_deprecated_shim

warn_deprecated_shim("hapmap", "genetics.genomics.variant.hapmap_report")


from genetics.genomics.variant.hapmap_report import report_hapmap_table  # noqa: F401
