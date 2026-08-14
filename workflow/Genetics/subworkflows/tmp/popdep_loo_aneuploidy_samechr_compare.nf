#!/usr/bin/env nextflow
/*
 * Compare aneuploidy same-chromosome LOO influence vs chr1 baseline LOO.
 *
 * Prerequisite: popdep_n500_loo_aneuploidy_samechr.nf completed for manifest chromosomes
 *   and chr1 full LOO info TSV exists.
 *
 * Launch:
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/64run_popdep_n500_loo_aneuploidy_samechr
 *   nextflow run .../popdep_loo_aneuploidy_samechr_compare.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink
 */

nextflow.enable.dsl=2

include { popdep_loo_aneuploidy_samechr_compare } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'

workflow {
    if (!params.home_dir || !params.user_dir || !params.src_dir) {
        error 'popdep_loo_aneuploidy_samechr_compare: --home_dir, --user_dir, --src_dir are required'
    }

    def manifest = params.popdep_loo_aneuploidy_manifest
        ?: "${params.home_dir}/00data/popdep_bench/n500_aneuploidy_samechr_loo.manifest.tsv"
    def chr1_info = params.popdep_loo_chr1_info_tsv
        ?: "${params.output_dir}/${params.job}/process/_popdep_n500_loo_chr1_10mb/stats/n500_chr1_10000000_loo.info.tsv"

    if (!file(manifest).exists()) {
        error "popdep_loo_aneuploidy_samechr_compare: missing manifest ${manifest}"
    }
    if (!file(chr1_info).exists()) {
        error "popdep_loo_aneuploidy_samechr_compare: missing chr1 LOO summary ${chr1_info}"
    }

    popdep_loo_aneuploidy_samechr_compare(
        file(manifest, checkIfExists: true),
        file(chr1_info, checkIfExists: true),
    )
}
