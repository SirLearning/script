"""TagSNP reporting from PLINK2 --indep-pairwise .prune.in."""

from genetics.genomics.plink.results_io import summarize_tagsnp_from_prune


def report_tagsnp_selection(prune_in_path, output_prefix, max_tags=1000):
    summarize_tagsnp_from_prune(prune_in_path, output_prefix, max_tags=max_tags)
