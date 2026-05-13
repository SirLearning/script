nextflow.enable.dsl=2

/*
 * Minimal entry: re-run LD decay + cross-chromosome LD baseline plots only.
 * Expects existing PLINK2 .vcor files under:
 *   ${params.output_dir}/${params.job}/process/${params.mod}/variant/
 * for subgenomes A, B, D, Others.
 *
 * Config: pass -c /path/to/workflow/Genetics/nextflow.config (see tmp/README.md).
 */
include { variant_ld_decay_plot; variant_ld_crosschr_plot } from '../genotype/stats.nf'

workflow {
    if (!params.mod) {
        error "params.mod is required (e.g. test_thin, test_common_thin)."
    }
    if (!params.output_dir || !params.job) {
        error "params.output_dir and params.job are required."
    }

    subgenomes = ['A', 'B', 'D', 'Others']
    base = "${params.output_dir}/${params.job}/process/${params.mod}/variant"

    ch_ld = channel.from(
        subgenomes.collect { sg ->
            tuple(sg, "sub_${sg}", file("${base}/${sg}.info.ld.vcor", checkIfExists: true))
        }
    )
    ch_ld_cross = channel.from(
        subgenomes.collect { sg ->
            tuple(sg, "sub_${sg}", file("${base}/${sg}.info.ld.crosschr.vcor", checkIfExists: true))
        }
    )

    variant_ld_decay_plot(ch_ld)
    variant_ld_crosschr_plot(ch_ld_cross)
}
