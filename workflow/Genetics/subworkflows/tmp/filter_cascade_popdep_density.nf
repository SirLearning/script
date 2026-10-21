#!/usr/bin/env nextflow
/*
 * PopDep RelativeDepth mean-vs-SD cell-density filter on all n500 grids,
 * then het → MAC≥5 → WAP MAC≥2 on the test_thin variant library.
 *
 * Does not overwrite process/test_thin popdep.info.tsv (annotates in the work dir
 * from --popdep_density_dir). Publishes counts under stats/{mod}/.
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/filter_cascade_popdep_density.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_thin \
 *     --popdep_density_dir /data1/.../process/_popdep_n500_gam_credible_ntaxa500
 */

nextflow.enable.dsl=2

include {
    popdep_density_hist_chr
    popdep_density_merge
    variant_filter_cascade_popdep_density
} from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = 'test_thin'
params.filter_cascade_het_mod = 'test_common_only'
params.filter_cascade_wap_mod = 'test_wap'
params.filter_cascade_wap_keep = null
params.filter_cascade_het_alpha = '1e-6'

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'filter_cascade_popdep_density: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir) {
        error 'filter_cascade_popdep_density: --user_dir is required (conda stats)'
    }
    if (!params.src_dir) {
        error 'filter_cascade_popdep_density: --src_dir is required'
    }
    def popdep_root = params.popdep_density_dir ?: params.popdep_dir
    if (!popdep_root) {
        error 'filter_cascade_popdep_density: --popdep_density_dir is required'
    }

    def wap_mod = params.filter_cascade_wap_mod ?: 'test_wap'
    def keep_name = params.filter_cascade_wap_keep ?: "${wap_mod}.keep"
    def wap_keep = file("${params.output_dir}/${params.job}/process/${wap_mod}/${keep_name}", checkIfExists: true)
    def grids = Channel.fromPath("${popdep_root}/variant/chr*.popdep.txt.bgz", checkIfExists: true)

    log.info "${params.c_green}PopDep density cascade ${params.mod} from ${popdep_root}${params.c_reset}"
    popdep_density_hist_chr(grids)
    popdep_density_merge(popdep_density_hist_chr.out.npz.collect())
    variant_filter_cascade_popdep_density(
        popdep_density_merge.out.masks
            .combine(popdep_density_merge.out.thresholds)
            .combine(Channel.of(wap_keep))
    )
}
