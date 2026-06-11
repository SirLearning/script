nextflow.enable.dsl=2

/*
 * Build per-chromosome site MQ reference grids (full reference scan, not variant targets).
 *
 * Uses bcftools mpileup over -r CHR:1-LEN for a fixed BAM subset; mean MAPQ from INFO/I16
 * (no bcftools call). Pads to one row per reference coordinate when mq_pad_all_positions=true.
 *
 * Legacy 01test runs used 2_1_122798052.pos.txt (variant coordinates, ~29.8M/chr);
 * this workflow intentionally does NOT read that file.
 *
 * Launch:
 *   partial_router.nf --partial_task abstract_mq_50_bams --mod abstract_mq_50_bams ...
 */

include { calc_site_mq_bcftools } from '../../../modules/local/genotype/processor/processor_mq.nf'
include {
    getRefV1SubChr
    getRefV1ChrLength
    getRefFastaForChr_v1
} from '../../../modules/local/infra/infra_ref_v1.nf'

workflow RUN_ABSTRACT_MQ_50_BAMS {
    main:
    if (params.mod != 'abstract_mq_50_bams') {
        error "RUN_ABSTRACT_MQ_50_BAMS: params.mod must be abstract_mq_50_bams (got: ${params.mod})."
    }
    if (!params.home_dir) {
        error 'RUN_ABSTRACT_MQ_50_BAMS: params.home_dir is required.'
    }
    if (!params.output_dir || !params.job) {
        error 'RUN_ABSTRACT_MQ_50_BAMS: params.output_dir and params.job are required.'
    }
    if (!params.mq_bam_list_file) {
        error 'RUN_ABSTRACT_MQ_50_BAMS: params.mq_bam_list_file is required.'
    }

    def bam_list = file(params.mq_bam_list_file, checkIfExists: true)
    def chr_list = params.chr
        ? [params.chr.toString().replaceFirst(/^chr/, '')]
        : getRefV1SubChr('ALL')

    if (params.mq_chr_exclude) {
        def exclude = params.mq_chr_exclude.toString().split(/[,\\s]+/).collect { it.replaceFirst(/^chr/, '').trim() }.findAll { it }
        chr_list = chr_list.findAll { c -> !exclude.contains(c.toString()) }
        log.info "abstract_mq_50_bams: excluding chromosomes ${exclude} (mq_chr_exclude)"
    }

    if (chr_list.isEmpty()) {
        error 'RUN_ABSTRACT_MQ_50_BAMS: no chromosomes left after mq_chr_exclude filter.'
    }

    log.info "abstract_mq_50_bams: full-chromosome mpileup reference (NOT variant-target list)"
    log.info "abstract_mq_50_bams: BAM list ${bam_list}"
    log.info "abstract_mq_50_bams: chromosomes ${chr_list}"
    log.info "abstract_mq_50_bams: publish calls *.site_mq.calls.tsv.gz (CHROM POS REF ALT MQ float from I16) + ref grid *.site_mq.ref.*"

    def ch = channel.from(chr_list).map { chr ->
        def chr_int = chr.toString().toInteger()
        def id = String.format('chr%03d', chr_int)
        def chr_len = getRefV1ChrLength(chr)
        def reference = file(getRefFastaForChr_v1(chr, params.home_dir), checkIfExists: true)
        tuple(id, chr, chr_len, reference, bam_list)
    }

    calc_site_mq_bcftools(ch)
}
