nextflow.enable.dsl=2

/*
 * Build per-chromosome site MQ reference grids (full reference scan, not variant targets).
 *
 * Uses bcftools mpileup over -r CHR:1-LEN for a fixed BAM subset; mean MAPQ from INFO/I16
 * (no bcftools call). Pads to one row per reference coordinate when mq_pad_all_positions=true.
 *
 * Published once to frozen params.mq_dir (default under vmap4 test_plink/process).
 * Re-runs skip chromosomes that already have *.site_mq.ref.* unless mq_force_rerun=true.
 * Not part of main.nf — launch only via partial_router.nf --partial_task abstract_mq_50_bams.
 *
 * Launch:
 *   partial_router.nf --partial_task abstract_mq_50_bams --mod abstract_mq_50_bams --job test_plink ...
 */

include { calc_site_mq_bcftools } from '../../../modules/local/genotype/processor/processor_mq.nf'
include {
    getRefV1SubChr
    getRefV1ChrLength
    getRefFastaForChr_v1
} from '../../../modules/local/infra/infra_ref_v1.nf'
include {
    hasSiteMqRef
    countSiteMqRefs
} from '../../../modules/local/infra/infra_mq_reuse.nf'

workflow RUN_ABSTRACT_MQ_50_BAMS {
    main:
    if (params.mod != 'abstract_mq_50_bams') {
        error "RUN_ABSTRACT_MQ_50_BAMS: params.mod must be abstract_mq_50_bams (got: ${params.mod})."
    }
    if (!params.home_dir) {
        error 'RUN_ABSTRACT_MQ_50_BAMS: params.home_dir is required.'
    }
    if (!params.mq_dir) {
        error 'RUN_ABSTRACT_MQ_50_BAMS: params.mq_dir is required.'
    }
    if (!params.mq_bam_list_file) {
        error 'RUN_ABSTRACT_MQ_50_BAMS: params.mq_bam_list_file is required.'
    }

    def mq_root = params.mq_dir
    def bam_list = file(params.mq_bam_list_file, checkIfExists: true)
    def chr_list = params.chr
        ? [params.chr.toString().replaceFirst(/^chr/, '')]
        : getRefV1SubChr('ALL')

    if (params.mq_chr_exclude) {
        def exclude = params.mq_chr_exclude.toString().split(/[,\\s]+/).collect { it.replaceFirst(/^chr/, '').trim() }.findAll { it }
        chr_list = chr_list.findAll { c -> !exclude.contains(c.toString()) }
        log.info "abstract_mq_50_bams: excluding chromosomes ${exclude} (mq_chr_exclude)"
    }

    if (!params.mq_force_rerun) {
        def before = chr_list.size()
        chr_list = chr_list.findAll { chr ->
            def id = String.format('chr%03d', chr.toString().toInteger())
            !hasSiteMqRef(mq_root, id)
        }
        def skipped = before - chr_list.size()
        if (skipped > 0) {
            log.info "${params.c_green}Reuse frozen site MQ from:${params.c_reset} ${mq_root}/reference (${skipped} chr skipped)"
        }
    }

    log.info "abstract_mq_50_bams: frozen publish root ${mq_root} (reference/, logs/)"
    log.info "abstract_mq_50_bams: full-chromosome mpileup reference (NOT variant-target list)"
    log.info "abstract_mq_50_bams: BAM list ${bam_list}"
    log.info "abstract_mq_50_bams: chromosomes to compute ${chr_list}"
    log.info "abstract_mq_50_bams: publish calls *.site_mq.calls.tsv.gz (CHROM POS REF ALT MQ float from I16) + ref grid *.site_mq.ref.*"

    def existing = countSiteMqRefs(mq_root)
    if (chr_list.isEmpty()) {
        log.info "${params.c_green}All site MQ references present (${existing} chr under ${mq_root}/reference); nothing to run.${params.c_reset}"
    } else {
        def ch = channel.from(chr_list).map { chr ->
            def chr_int = chr.toString().toInteger()
            def id = String.format('chr%03d', chr_int)
            def chr_len = getRefV1ChrLength(chr)
            def reference = file(getRefFastaForChr_v1(chr, params.home_dir), checkIfExists: true)
            tuple(id, chr, chr_len, reference, bam_list)
        }

        calc_site_mq_bcftools(ch)
    }
}
