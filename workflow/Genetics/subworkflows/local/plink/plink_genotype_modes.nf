nextflow.enable.dsl=2

/*
 * Composed DSL2 workflows for PLINK / PLINK2 genotype tracks (non-wheat modes).
 * Process definitions live under workflow/Genetics/modules/local/{genotype,dynamic,static,integrated}/.
 */

include {
    plink_processor as PLINK_PROCESSOR
    test_plink_processor as TEST_PLINK_PROCESSOR
    test_plink_camp as TEST_PLINK_CAMP_PROCESSOR
    test_common_thin_processor as TEST_COMMON_THIN_PROCESSOR
    test_common_only_processor as TEST_COMMON_ONLY_PROCESSOR
    test_common_inter_processor as TEST_COMMON_INTER_PROCESSOR
    test_common_inter_refchr_processor as TEST_COMMON_INTER_REFCHR_PROCESSOR
    test_wap_refchr_processor as TEST_WAP_REFCHR_PROCESSOR
    test_wap_thin_refchr_processor as TEST_WAP_THIN_REFCHR_PROCESSOR
    test_hex_thin_processor as TEST_HEX_THIN_PROCESSOR
    test_wap_common_refchr_processor as TEST_WAP_COMMON_REFCHR_PROCESSOR
    test_rare_only_processor as TEST_RARE_ONLY_PROCESSOR
    test_hex_processor as TEST_HEX_PROCESSOR
    test_wap_processor as TEST_WAP_PROCESSOR
    test_wap_exon_processor as TEST_WAP_EXON_PROCESSOR
    test_collinear_processor as TEST_COLLINEAR_PROCESSOR
} from '../../../modules/local/genotype/processor/processor.nf'
include {
    plink_stats as PLINK_STATS
    test_plink_stats as TEST_PLINK_STATS
} from '../../../modules/local/genotype/stats/stats.nf'
include {
    variant_ld_r2_decay_along_subgenome_plot
    variant_ld_decay_plot
    variant_ld_decay_rare_common_plot
} from '../../../modules/local/genotype/stats/stats_variant.nf'

workflow RUN_V1_PLINK {
    take:
    ch_vcf

    main:
    def processor_out = PLINK_PROCESSOR(ch_vcf)
    PLINK_STATS(
        processor_out.smiss,
        processor_out.scount,
        processor_out.vmiss,
        processor_out.gcount,
        processor_out.afreq,
        processor_out.hardy,
        processor_out.popdep)
}

workflow RUN_TEST_PLINK {
    take:
    ch_vcf

    main:
    def processor_out = TEST_PLINK_PROCESSOR(ch_vcf)
    TEST_PLINK_STATS(
        processor_out.ld,
        processor_out.ld_cross,
        processor_out.smiss,
        processor_out.scount,
        processor_out.vmiss,
        processor_out.gcount,
        processor_out.afreq,
        processor_out.hardy,
        processor_out.mq,
        processor_out.popdep,
        processor_out.merged_bfile,
        processor_out.merged_pfile)
}

workflow RUN_TEST_PLINK_CAMP {
    take:
    ch_vcf
    camp_vmap4_map_tsv

    main:
    def processor_out = TEST_PLINK_CAMP_PROCESSOR(ch_vcf, camp_vmap4_map_tsv)
    TEST_PLINK_STATS(
        processor_out.ld,
        processor_out.ld_cross,
        processor_out.smiss,
        processor_out.scount,
        processor_out.vmiss,
        processor_out.gcount,
        processor_out.afreq,
        processor_out.hardy,
        processor_out.mq,
        processor_out.popdep,
        processor_out.merged_bfile,
        processor_out.merged_pfile)
}
workflow RUN_TEST_COMMON_THIN {
    take:
    ch_vcf

    main:
    def processor_out = TEST_COMMON_THIN_PROCESSOR(ch_vcf)
    TEST_PLINK_STATS(
        processor_out.ld,
        processor_out.ld_cross,
        processor_out.smiss,
        processor_out.scount,
        processor_out.vmiss,
        processor_out.gcount,
        processor_out.afreq,
        processor_out.hardy,
        processor_out.mq,
        processor_out.popdep,
        processor_out.merged_bfile,
        processor_out.merged_pfile)
}

workflow RUN_TEST_COMMON_ONLY {
    take:
    ch_vcf

    main:
    def processor_out = TEST_COMMON_ONLY_PROCESSOR(ch_vcf)
    TEST_PLINK_STATS(
        processor_out.ld,
        processor_out.ld_cross,
        processor_out.smiss,
        processor_out.scount,
        processor_out.vmiss,
        processor_out.gcount,
        processor_out.afreq,
        processor_out.hardy,
        processor_out.mq,
        processor_out.popdep,
        processor_out.merged_bfile,
        processor_out.merged_pfile)
}

workflow RUN_TEST_COMMON_INTER {
    take:
    ch_vcf

    main:
    def processor_out = TEST_COMMON_INTER_PROCESSOR(ch_vcf)
    TEST_PLINK_STATS(
        processor_out.ld,
        processor_out.ld_cross,
        processor_out.smiss,
        processor_out.scount,
        processor_out.vmiss,
        processor_out.gcount,
        processor_out.afreq,
        processor_out.hardy,
        processor_out.mq,
        processor_out.popdep,
        processor_out.merged_bfile,
        processor_out.merged_pfile)
}

workflow RUN_TEST_COMMON_INTER_REFCHR {
    main:
    def processor_out = TEST_COMMON_INTER_REFCHR_PROCESSOR()
    TEST_PLINK_STATS(
        processor_out.ld,
        processor_out.ld_cross,
        processor_out.smiss,
        processor_out.scount,
        processor_out.vmiss,
        processor_out.gcount,
        processor_out.afreq,
        processor_out.hardy,
        processor_out.mq,
        processor_out.popdep,
        processor_out.merged_bfile,
        processor_out.merged_pfile)
    processor_out.ld_along
        .filter { id, _chrom, _vcor, _focals -> id.toString() in ['A', 'B', 'D'] }
        .groupTuple(by: 0)
        .set { ch_sg_refchr }
    variant_ld_r2_decay_along_subgenome_plot(ch_sg_refchr)
}

workflow RUN_TEST_WAP_REFCHR {
    main:
    def processor_out = TEST_WAP_REFCHR_PROCESSOR()
    variant_ld_decay_plot(processor_out.ld)
    processor_out.ld
        .join(processor_out.afreq, by: [0, 1])
        .set { ch_ld_afreq }
    variant_ld_decay_rare_common_plot(ch_ld_afreq)
}

workflow RUN_TEST_WAP_THIN_REFCHR {
    main:
    def processor_out = TEST_WAP_THIN_REFCHR_PROCESSOR()
    variant_ld_decay_plot(processor_out.ld)
    processor_out.ld
        .join(processor_out.afreq, by: [0, 1])
        .set { ch_ld_afreq }
    variant_ld_decay_rare_common_plot(ch_ld_afreq)
}

workflow RUN_TEST_HEX_THIN {
    main:
    def processor_out = TEST_HEX_THIN_PROCESSOR()
    variant_ld_decay_plot(processor_out.ld)
    processor_out.ld
        .join(processor_out.afreq, by: [0, 1])
        .set { ch_ld_afreq }
    variant_ld_decay_rare_common_plot(ch_ld_afreq)
}

workflow RUN_TEST_WAP_COMMON_REFCHR {
    take:
    ch_vcf

    main:
    def processor_out = TEST_WAP_COMMON_REFCHR_PROCESSOR(ch_vcf)
    processor_out.ld_along
        .filter { id, _chrom, _vcor, _focals -> id.toString() in ['A', 'B', 'D'] }
        .groupTuple(by: 0)
        .set { ch_sg_wap }
    variant_ld_r2_decay_along_subgenome_plot(ch_sg_wap)
}

workflow RUN_TEST_RARE_ONLY {
    take:
    ch_vcf

    main:
    def processor_out = TEST_RARE_ONLY_PROCESSOR(ch_vcf)
    TEST_PLINK_STATS(
        processor_out.ld,
        processor_out.ld_cross,
        processor_out.smiss,
        processor_out.scount,
        processor_out.vmiss,
        processor_out.gcount,
        processor_out.afreq,
        processor_out.hardy,
        processor_out.mq,
        processor_out.popdep,
        processor_out.merged_bfile,
        processor_out.merged_pfile)
}

workflow RUN_TEST_HEX {
    take:
    ch_vcf

    main:
    def processor_out = TEST_HEX_PROCESSOR(ch_vcf)
    TEST_PLINK_STATS(
        processor_out.ld,
        processor_out.ld_cross,
        processor_out.smiss,
        processor_out.scount,
        processor_out.vmiss,
        processor_out.gcount,
        processor_out.afreq,
        processor_out.hardy,
        processor_out.mq,
        processor_out.popdep,
        processor_out.merged_bfile,
        processor_out.merged_pfile)
}

workflow RUN_TEST_WAP_EXON {
    take:
    ch_vcf

    main:
    // Exon extract only: no windowed LD / MQ / popdep. Those tracks do not
    // scale to an unthinned coding callset and are not inputs to VEP or burden.
    TEST_WAP_EXON_PROCESSOR(ch_vcf)
}

workflow RUN_TEST_WAP {
    take:
    ch_vcf

    main:
    def processor_out = TEST_WAP_PROCESSOR(ch_vcf)
    TEST_PLINK_STATS(
        processor_out.ld,
        processor_out.ld_cross,
        processor_out.smiss,
        processor_out.scount,
        processor_out.vmiss,
        processor_out.gcount,
        processor_out.afreq,
        processor_out.hardy,
        processor_out.mq,
        processor_out.popdep,
        processor_out.merged_bfile,
        processor_out.merged_pfile)
}

workflow RUN_TEST_COLLINEAR {
    take:
    ch_vcf
    collinear_root

    main:
    def processor_out = TEST_COLLINEAR_PROCESSOR(ch_vcf, collinear_root)
    TEST_PLINK_STATS(
        processor_out.ld,
        processor_out.ld_cross,
        processor_out.smiss,
        processor_out.scount,
        processor_out.vmiss,
        processor_out.gcount,
        processor_out.afreq,
        processor_out.hardy,
        processor_out.mq,
        processor_out.popdep,
        processor_out.merged_bfile,
        processor_out.merged_pfile)
}