nextflow.enable.dsl=2

/*
 * Classify a VEP-like table into project coding masks (plof_high / missense).
 * Does not run VEP. Use after a one-chromosome VEP/csq table exists.
 *
 *   nextflow run .../subworkflows/tmp/coding_mask_chr_pilot.nf \
 *     -c .../nextflow.config \
 *     --user_dir ... --output_dir ... --job test_plink --mod coding_mask_chr032 \
 *     --coding_vep_table /path/to/chr032.vep.tsv
 */

include { coding_mask_from_vep } from '../../modules/local/genotype/stats/stats_cohort.nf'

workflow {
    if (!params.output_dir || !params.job || !params.user_dir) {
        error 'coding_mask_chr_pilot: --output_dir, --job, and --user_dir are required.'
    }
    if (!params.mod) {
        error 'coding_mask_chr_pilot: --mod is required (publish folder).'
    }
    if (!params.coding_vep_table) {
        error 'coding_mask_chr_pilot: --coding_vep_table is required (VEP TSV).'
    }
    def prefix = params.coding_mask_prefix ?: "${params.mod}"
    coding_mask_from_vep(channel.value(tuple(prefix, file(params.coding_vep_table, checkIfExists: true))))
}
