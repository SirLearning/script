#!/usr/bin/env nextflow
/*
 * N500 aneuploidy same-chromosome PopDep LOO (flagged samples on matched popdep_chr).
 *
 * Uses targeted LOO maps from n500_aneuploidy_samechr_chr<N>.taxa.txt (not full 500-case sweep).
 *
 * Publish: {output_dir}/{job}/process/_popdep_n500_loo_aneuploidy_samechr/chr<N>/
 *
 * Single chr (example chr4 — CRR875588 chr4B trisomic):
 *   mkdir -p /data/home/tusr1/01projects/vmap4/10stats.genome/64run_popdep_n500_loo_aneuploidy_samechr/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/64run_popdep_n500_loo_aneuploidy_samechr
 *   nextflow run .../popdep_n500_loo_aneuploidy_samechr.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink \
 *     --popdep_loo_chr 4
 *
 * All manifest chromosomes: run launch_aneuploidy_samechr_loo.sh from the run folder (sequential).
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
        error 'popdep_n500_loo_aneuploidy_samechr: --home_dir, --user_dir, --src_dir are required'
    }
    if (!params.output_dir || !params.job) {
        error 'popdep_n500_loo_aneuploidy_samechr: --output_dir and --job are required'
    }

    def chr_val = params.popdep_loo_chr.toString().toInteger()
    def publish_root = params.popdep_loo_publish_root
        ?: "${params.output_dir}/${params.job}/process/_popdep_n500_loo_aneuploidy_samechr/chr${chr_val}"
    def summary_prefix = params.popdep_loo_summary_prefix
        ?: "n500_aneuploidy_samechr_chr${chr_val}_${params.popdep_loo_bench_max_bp}_loo"
    def target_taxa_file = params.popdep_loo_target_taxa_file
        ?: "${params.home_dir}/00data/popdep_bench/n500_aneuploidy_samechr_chr${chr_val}.taxa.txt"

    params.popdep_loo_publish_root = publish_root
    params.popdep_loo_summary_prefix = summary_prefix
    params.popdep_loo_target_taxa_file = target_taxa_file

    def source_map = params.popdep_loo_source_taxa_map
        ?: "${params.home_dir}/00data/05taxaBamMap/popdep/tb.N500_gam_credible_deep.taxaBamMap.txt"
    if (!file(source_map).exists()) {
        error "popdep_n500_loo_aneuploidy_samechr: missing source taxa map ${source_map}"
    }
    if (!file(target_taxa_file).exists()) {
        error "popdep_n500_loo_aneuploidy_samechr: missing target taxa file ${target_taxa_file}"
    }

    def grid_id = String.format('chr%03d', chr_val)
    def stub_vcf = file("${params.home_dir}/00data/popdep_bench/stub.vcf.gz")

    log.info "popdep_n500_loo_aneuploidy_samechr: chr=${chr_val} publish_root=${publish_root}"
    log.info "popdep_n500_loo_aneuploidy_samechr: target_taxa_file=${target_taxa_file}"

    prepare_popdep_loo_taxa_maps(
        publish_root,
        file(source_map, checkIfExists: true),
        0,
        target_taxa_file,
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
    ch_tiger_config = channel.value(tuple(pd_config.path, pd_config.app_name, pd_config.java_version))

    calc_population_depth_loo_case(ch_cases, ch_tiger_config)
    popdep_tiger_gz_to_bgzip_tabix_loo(
        calc_population_depth_loo_case.out.tiger_gz
            .map { case_id, chr_id, gz -> tuple(publish_root, case_id, chr_id, gz) },
    )

    if (params.popdep_loo_compare) {
        def baseline_grid = "${publish_root}/baseline/variant/${grid_id}.popdep.txt.bgz"
        ch_popdep_done = popdep_tiger_gz_to_bgzip_tabix_loo.out.popdep
            .map { case_id, chr_id, bgz, tbi -> case_id }
            .collect()
            .ifEmpty(channel.value([]))
        ch_summary_in = ch_prep_manifest
            .combine(ch_popdep_done)
            .map { entry ->
                def manifest = entry[0]
                tuple(manifest.readLines().size() - 1, publish_root, baseline_grid, grid_id, manifest, summary_prefix)
            }
        popdep_loo_influence_summary(ch_summary_in)
    }
}
