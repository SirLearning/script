#!/usr/bin/env nextflow
/*
 * N500 leave-one-out PopDepFull benchmark: chr1 first 10 Mb.
 *
 * Prepares baseline + LOO taxaBamMaps, runs TIGER PopDepFull per case, BGZF+tabix,
 * and optional influence summary vs baseline.
 *
 * Publish layout (params.popdep_loo_publish_root):
 *   baseline/variant/chr001.popdep.txt.bgz (+ .tbi)
 *   baseline/logs/*.log
 *   loo_<taxa>/variant/chr001.popdep.txt.bgz
 *   loo_<taxa>/logs/*.log
 *   meta/cases.tsv
 *   stats/n500_chr1_10mb_loo.info.tsv
 *
 * Launch from a vmap4 run folder (conda run; screen for full 501-case sweep):
 *   mkdir -p /data/home/tusr1/01projects/vmap4/10stats.genome/58run_popdep_n500_loo_chr1_nf/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/58run_popdep_n500_loo_chr1_nf
 *   nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/popdep_n500_loo_chr1.nf \
 *     -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink \
 *     --popdep_loo_max_cases 2
 *
 * Full panel: omit --popdep_loo_max_cases or set --popdep_loo_max_cases 0
 */

nextflow.enable.dsl=2

include {
    prepare_popdep_loo_taxa_maps
    calc_population_depth_loo_case
    popdep_tiger_gz_to_bgzip_tabix_loo
} from '../../modules/local/genotype/processor/processor_depth.nf'

include { popdep_loo_influence_summary } from '../../modules/local/genotype/stats/stats_variant.nf'
include { getTigerJarConfig } from '../../modules/local/infra/infra_tiger.nf'

params.job = 'test_plink'

workflow {
    if (!params.home_dir || !params.user_dir || !params.src_dir) {
        error 'popdep_n500_loo_chr1: --home_dir, --user_dir, --src_dir are required'
    }
    if (!params.output_dir || !params.job) {
        error 'popdep_n500_loo_chr1: --output_dir and --job are required'
    }

    def source_map = params.popdep_loo_source_taxa_map
        ?: "${params.home_dir}/00data/05taxaBamMap/popdep/tb.N500_gam_credible_deep.taxaBamMap.txt"
    if (!file(source_map).exists()) {
        error "popdep_n500_loo_chr1: missing source taxa map ${source_map}"
    }

    def publish_root = params.popdep_loo_publish_root
        ?: "${params.output_dir}/${params.job}/process/${params.popdep_loo_publish_slug}"
    def chr_val = params.popdep_loo_chr.toString()
    def grid_id = String.format('chr%03d', chr_val.toInteger())
    def max_loo = params.popdep_loo_max_cases != null ? params.popdep_loo_max_cases.toInteger() : 0
    def stub_vcf = file("${params.home_dir}/00data/popdep_bench/stub.vcf.gz")
    def target_path = params.popdep_loo_target_taxa_file ?: ''

    log.info "popdep_n500_loo_chr1: source_map=${source_map}"
    log.info "popdep_n500_loo_chr1: publish_root=${publish_root}"
    log.info "popdep_n500_loo_chr1: chr=${chr_val} max_bp=${params.popdep_loo_bench_max_bp} max_loo_cases=${max_loo}"
    if (target_path) {
        log.info "popdep_n500_loo_chr1: target_taxa_file=${target_path}"
    }

    prepare_popdep_loo_taxa_maps(
        publish_root,
        file(source_map, checkIfExists: true),
        max_loo,
        target_path,
    )

    ch_prep_manifest = prepare_popdep_loo_taxa_maps.out.manifest
    ch_prep_maps = prepare_popdep_loo_taxa_maps.out.taxa_maps

    ch_cases = ch_prep_manifest
        .splitCsv(sep: '\t', header: true)
        .combine(ch_prep_maps)
        .map { row, maps_dir ->
            def rel = row.taxa_map.replaceFirst(/^taxa_maps\//, '')
            tuple(
                publish_root,
                row.case_id,
                row.chr.toString(),
                file("${maps_dir}/${rel}"),
                stub_vcf,
            )
        }
        .filter { pub_root, case_id, chr_id, taxa_map, vcf ->
            if (params.popdep_loo_force_rerun) {
                return true
            }
            !file("${pub_root}/${case_id}/variant/${grid_id}.popdep.txt.bgz").exists()
        }

    def pd_config = getTigerJarConfig(
        params.popdep_tiger_jar,
        params.home_dir,
        'PopDepFull',
    )
    log.info "popdep_n500_loo_chr1: TIGER jar=${pd_config.path} app=${pd_config.app_name}"
    ch_tiger_config = channel.value(tuple(pd_config.path, pd_config.app_name, pd_config.java_version))

    calc_population_depth_loo_case(ch_cases, ch_tiger_config)
    popdep_tiger_gz_to_bgzip_tabix_loo(
        calc_population_depth_loo_case.out.tiger_gz
            .map { case_id, chr_id, gz -> tuple(publish_root, case_id, chr_id, gz) },
    )

    if (params.popdep_loo_compare) {
        def baseline_grid = "${publish_root}/baseline/variant/${grid_id}.popdep.txt.bgz"
        def summary_prefix = params.popdep_loo_summary_prefix
            ?: "n500_chr${chr_val}_${params.popdep_loo_bench_max_bp}_loo"
        ch_popdep_done = popdep_tiger_gz_to_bgzip_tabix_loo.out.popdep
            .map { case_id, chr_id, bgz, tbi -> case_id }
            .collect()
            .ifEmpty( channel.value([]) )
        ch_summary_in = ch_prep_manifest
            .combine(ch_popdep_done)
            .map { entry ->
                def manifest = entry[0]
                def n_cases = manifest.readLines().size() - 1
                tuple(n_cases, publish_root, baseline_grid, grid_id, manifest, summary_prefix)
            }
        popdep_loo_influence_summary(ch_summary_in)
    }
}
