nextflow.enable.dsl=2

/*
 * Build per-chromosome TIGER population depth grids (full reference scan, not variant targets).
 *
 * PopDepCrossChr: one JVM pass over all BAMs (tb.ALL); length file Chr/Length/nTaxa per PopDepFull tb.A/B/D/ALL.
 * PopDepFull / PopDepCrossChr TIGER steps emit per-chr *.popdep.txt.gz; shared popdep_tiger_gz_to_bgzip_tabix
 * converts gzip → BGZF + tabix (label popdep_bgz_tabix, parallel per chr).
 *
 * Published once to frozen params.popdep_dir/variant (default main_raw under vmap4 test_plink/process).
 * Re-runs skip chromosomes that already have *.popdep.txt.bgz (or legacy *.popdep.txt) unless popdep_force_rerun=true.
 * Not part of main.nf — launch only via partial_router.nf --partial_task main_raw_popdepth.
 *
 * Launch (from vmap4 run dir, conda run, screen for full sweep):
 *   partial_router.nf --partial_task main_raw_popdepth --mod main_raw_popdepth --job test_plink ...
 */

include { calc_population_depth } from '../../../modules/local/genotype/processor/processor_depth.nf'
include { calc_population_depth_crosschr } from '../../../modules/local/genotype/processor/processor_depth.nf'
include { popdep_tiger_gz_to_bgzip_tabix } from '../../../modules/local/genotype/processor/processor_depth.nf'
include { getRefV1SubChr } from '../../../modules/local/infra/infra_ref_v1.nf'
include { countPopDepTaxaFile_v1 } from '../../../modules/local/infra/infra_ref_v1.nf'
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
    def chr_list = params.popdep_chr_list
        ? params.popdep_chr_list.toString().split(/[,\\s]+/).collect { it.replaceFirst(/^chr/, '').trim() }.findAll { it }
        : (params.chr
            ? [params.chr.toString().replaceFirst(/^chr/, '')]
            : getRefV1SubChr('ALL'))

    if (params.popdep_chr_list) {
        log.info "main_raw_popdepth: bench chr list ${chr_list} (popdep_chr_list)"
    }
    if (params.popdep_taxa_bam_file) {
        log.info "main_raw_popdepth: taxa-BAM override ${params.popdep_taxa_bam_file}"
    }
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

    def missing_vcf = params.popdep_skip_vcf_check ? [] : chr_list.findAll { chr ->
        def id = String.format('chr%03d', chr.toString().toInteger())
        !file("${popdep_root}/${id}.vcf.gz").exists()
    }
    if (missing_vcf) {
        log.warn "main_raw_popdepth: skipping chromosomes with no input VCF under ${popdep_root}: ${missing_vcf}"
        chr_list = chr_list.findAll { chr -> !missing_vcf.contains(chr.toString()) }
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
        if (pd_config.app_name == 'PopDepCrossChr') {
            def tb_path = params.popdep_taxa_bam_file ?: "${params.home_dir}/00data/05taxaBamMap/vmap4_v1/tb.ALL.txt"
            log.info "PopDepCrossChr: scan ${tb_path}; per-chr nTaxa A=${countPopDepTaxaFile_v1(params.home_dir, 'A')} B=${countPopDepTaxaFile_v1(params.home_dir, 'B')} D=${countPopDepTaxaFile_v1(params.home_dir, 'D')} Others/ALL=${countPopDepTaxaFile_v1(params.home_dir, 'ALL')}; threads=${params.popdep_crosschr_threads} memory=${params.popdep_crosschr_memory_gb}G"
            calc_population_depth_crosschr(channel.value(chr_list), ch_tiger_config)
            ch_crosschr_gz = calc_population_depth_crosschr.out.tiger_gz
                .flatMap { gzs ->
                    (gzs instanceof List ? gzs : [gzs]).collect { gz ->
                        def chr = gz.getName().replaceAll(/\.popdep\.txt\.gz$/, '')
                        def id = String.format('chr%03d', chr.toInteger())
                        tuple(id, chr.toString(), gz)
                    }
                }
            popdep_tiger_gz_to_bgzip_tabix(ch_crosschr_gz)
        } else {
            def full_threads = params.popdep_tiger_threads != null
                ? "${params.popdep_tiger_threads}".toInteger()
                : 32
            def full_mem_gb = params.popdep_tiger_memory_gb != null
                ? "${params.popdep_tiger_memory_gb}".toInteger()
                : 128
            log.info "PopDepFull: threads=${full_threads} memory=${full_mem_gb}G maxForks=${params.popdep_tiger_max_forks}"
            ch_vcf = channel.from(chr_list).map { chr ->
                def chr_int = chr.toString().toInteger()
                def id = String.format('chr%03d', chr_int)
                def vcf = params.popdep_skip_vcf_check
                    ? file("${params.home_dir}/00data/popdep_bench/stub.vcf.gz")
                    : file("${popdep_root}/${id}.vcf.gz", checkIfExists: true)
                tuple(id, chr_int.toString(), vcf)
            }
            calc_population_depth(ch_vcf, ch_tiger_config)
            popdep_tiger_gz_to_bgzip_tabix(calc_population_depth.out.tiger_gz)
        }
    }
}
