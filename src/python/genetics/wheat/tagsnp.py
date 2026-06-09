"""Deprecated: use genetics.genomics.variant.tagsnp_report."""

from genetics.wheat._deprecation import warn_deprecated_shim

warn_deprecated_shim("tagsnp", "genetics.genomics.variant.tagsnp_report")


from genetics.genomics.variant.tagsnp_report import report_tagsnp_selection  # noqa: F401
