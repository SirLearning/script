nextflow.enable.dsl=2

/*
 * Build per-chromosome TIGER population depth grids (full reference scan, not variant targets).
 *
 * Published once to frozen params.popdep_dir/variant (default main_raw under vmap4 test_plink/process).
 * Re-runs skip chromosomes that already have *.popdep.txt unless popdep_force_rerun=true.
 * Not part of main.nf — launch only via partial_router.nf --partial_task main_raw_popdepth.
 *
 * Launch (from vmap4 run dir, conda run, screen for full sweep):
 *   partial_router.nf --partial_task main_raw_popdepth --mod main_raw_popdepth --job test_plink ...
 */

include { calc_population_depth } from '../../../modules/local/genotype/processor/processor_depth.nf'
include { getRefV1SubChr } from '../../../modules/local/infra/infra_ref_v1.nf'
include { getTigerJarConfig } from '../../../modules/local/infra/infra_tiger.nf'
include {
    hasPopdepForChr
    countPopdeps
} from '../../../modules/local/infra/infra_popdep_reuse.nf'

workflow RUN_MAIN_RAW_POPDEPTH {
    main:
    if (params.mod != 'main_raw_popdepth') {
        error "RUN_MAIN_RAW_POPDEPTH: params.mod must be main_raw_popdepth (got: ${params.mod})."
    }
    if (!params.home_dir) {
        error 'RUN_MAIN_RAW_POPDEPTH: params.home_dir is required.'
    }
    if (!params.popdep_dir) {
        error 'RUN_MAIN_RAW_POPDEPTH: params.popdep_dir is required.'
    }

    def popdep_root = params.popdep_dir
    def chr_list = params.chr
        ? [params.chr.toString().replaceFirst(/^chr/, '')]
        : getRefV1SubChr('ALL')

    if (params.popdep_chr_exclude) {
        def exclude = params.popdep_chr_exclude.toString().split(/[,\\s]+/).collect { it.replaceFirst(/^chr/, '').trim() }.findAll { it }
        chr_list = chr_list.findAll { c -> !exclude.contains(c.toString()) }
        log.info "main_raw_popdepth: excluding chromosomes ${exclude} (popdep_chr_exclude)"
    }

    if (!params.popdep_force_rerun) {
        def before = chr_list.size()
        chr_list = chr_list.findAll { chr ->
            def id = String.format('chr%03d', chr.toString().toInteger())
            !hasPopdepForChr(popdep_root, id)
        }
        def skipped = before - chr_list.size()
        if (skipped > 0) {
            log.info "${params.c_green}Reuse frozen popdepth from:${params.c_reset} ${popdep_root}/variant (${skipped} chr skipped)"
        }
    }

    log.info "main_raw_popdepth: frozen publish root ${params.popdep_publish_dir ?: popdep_root}/variant"
    log.info "main_raw_popdepth: input VCFs ${popdep_root}/chrNNN.vcf.gz"
    log.info "main_raw_popdepth: chromosomes to compute ${chr_list}"

    def existing = countPopdeps(popdep_root)
    if (chr_list.isEmpty()) {
        log.info "${params.c_green}All popdepth references present (${existing} chr under ${popdep_root}/variant); nothing to run.${params.c_reset}"
    } else {
        def pd_config = getTigerJarConfig(params.popdep_tiger_jar, params.home_dir, params.popdep_tiger_app)
        log.info "main_raw_popdepth: TIGER ${params.popdep_tiger_jar} app=${pd_config.app_name}"
        ch_tiger_config = channel.value(tuple(pd_config.path, pd_config.app_name, pd_config.java_version))
        ch_vcf = channel.from(chr_list).map { chr ->
            def chr_int = chr.toString().toInteger()
            def id = String.format('chr%03d', chr_int)
            def vcf = file("${popdep_root}/${id}.vcf.gz", checkIfExists: true)
            tuple(id, chr_int.toString(), vcf)
        }
        calc_population_depth(ch_vcf, ch_tiger_config)
    }
}
