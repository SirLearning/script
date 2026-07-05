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
} from '../../../modules/local/genotype/processor/processor.nf'
include {
    plink_stats as PLINK_STATS
    test_plink_stats as TEST_PLINK_STATS
} from '../../../modules/local/genotype/stats/stats.nf'

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
        processor_out.popdep)
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
        processor_out.popdep)
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
        processor_out.popdep)
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
        processor_out.popdep)
}

