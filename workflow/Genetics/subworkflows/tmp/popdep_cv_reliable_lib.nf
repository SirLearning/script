#!/usr/bin/env nextflow
/*
 * n500 RelMean vs RelCV reliable-site libraries.
 *
 * Stream all `_popdep_n500_gam_credible_ntaxa500` 1 bp grids into a
 * 500×500 linear histogram (RelMean [0, 6] × RelCV [0, 25] = 250_000
 * cells). Keep the densest non-empty cells until 60% and 80% of
 * in-range sites remain. Write two merged BED libraries (0-based
 * half-open, PopDep Chrom/Position).
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/popdep_cv_reliable_lib.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod _popdep_n500_cv_reliable \
 *     --popdep_density_dir /data1/.../process/_popdep_n500_gam_credible_ntaxa500
 */

nextflow.enable.dsl=2

include {
    popdep_cv_reliable_hist_chr
    popdep_cv_reliable_merge
    popdep_cv_reliable_lib_chr
    popdep_cv_reliable_concat
} from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = '_popdep_n500_cv_reliable'
params.popdep_density_dir = null

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'popdep_cv_reliable_lib: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir) {
        error 'popdep_cv_reliable_lib: --user_dir is required (conda stats)'
    }
    if (!params.src_dir) {
        error 'popdep_cv_reliable_lib: --src_dir is required'
    }
    def popdep_root = params.popdep_density_dir ?: params.popdep_dir
    if (!popdep_root) {
        error 'popdep_cv_reliable_lib: --popdep_density_dir is required'
    }

    def grid_glob = "${popdep_root}/variant/chr*.popdep.txt.bgz"
    log.info "${params.c_green}n500 RelMean-CV 500x500 keep 60/80 libraries from ${popdep_root}${params.c_reset}"

    popdep_cv_reliable_hist_chr(Channel.fromPath(grid_glob, checkIfExists: true))
    popdep_cv_reliable_merge(popdep_cv_reliable_hist_chr.out.npz.collect())
    popdep_cv_reliable_lib_chr(
        Channel.fromPath(grid_glob, checkIfExists: true)
            .combine(popdep_cv_reliable_merge.out.masks)
    )
    popdep_cv_reliable_concat(
        popdep_cv_reliable_lib_chr.out.bed60.collect(),
        popdep_cv_reliable_lib_chr.out.bed80.collect(),
        popdep_cv_reliable_lib_chr.out.info.collect(),
    )
}
