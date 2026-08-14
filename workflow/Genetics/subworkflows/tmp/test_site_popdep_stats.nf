nextflow.enable.dsl=2

/*
 * Popdep-only stats for test_site mod (all 10 Mb sampled sites; variants highlighted).
 *
 * Prerequisite: test_site_build.nf has written
 *   {output_dir}/{job}/process/test_site/sites/panel.sites.tsv
 *
 * Publishes under {output_dir}/{job}/stats/test_site/{plots,info,logs}
 *
 * Example:
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/77run_test_site_popdep_stats
 *   nextflow run .../subworkflows/tmp/test_site_popdep_stats.nf -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink \
 *     --mod test_site
 */

include { test_site_popdep_plots } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = null

workflow {
    if (!params.mod) {
        error 'test_site_popdep_stats: --mod test_site is required.'
    }
    if (!params.output_dir || !params.job) {
        error 'test_site_popdep_stats: --output_dir and --job are required.'
    }

    def process_dir = params.process_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}"
    test_site_popdep_plots(channel.value(process_dir))
}
