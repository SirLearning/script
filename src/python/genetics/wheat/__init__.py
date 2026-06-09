"""Deprecated package: use genetics.genomics.* and genetics.gwas.association_plot."""

import warnings

from genetics.wheat._deprecation import warn_deprecated_shim

warn_deprecated_shim(
    "__init__",
    "genetics.genomics / genetics.gwas (see genetics.wheat.<module> docstrings)",
)
warnings.filterwarnings("default", category=DeprecationWarning, module=r"genetics\.wheat.*")

from genetics.genomics.plink.results_io import *  # noqa: F403
from genetics.genomics.sample.pca_structure import plot_population_structure
from genetics.genomics.variant.cnv_plot import plot_cnv_results
from genetics.genomics.variant.genetic_map_plot import plot_genetic_map
from genetics.genomics.variant.hapmap_report import report_hapmap_table
from genetics.genomics.variant.snp_qc_plot import summarize_and_plot_snp_qc
from genetics.genomics.variant.tagsnp_report import report_tagsnp_selection
from genetics.gwas.association_plot import plot_gwas_results, plot_kgwas_results
