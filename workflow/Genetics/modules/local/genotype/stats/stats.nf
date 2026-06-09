nextflow.enable.dsl=2

include { getTaxaBamMapFile_v1 } from '../../infra/infra_ref_v1.nf'

include {
    variant_missing_stats
    variant_maf_stats
    variant_mac_stats
    variant_mac_dist_log_redraw
    variant_popdep_mahalanobis
    variant_ld_decay_plot
    variant_ld_crosschr_plot
    variant_popdep_stats
} from './stats_variant.nf'

include {
    sample_missing_stats
    sample_coverage_stats
    sample_heterozygosity_stats
    sample_mapping_rate_stats
    sample_ref_ibs_stats
    sample_king_stats
    sample_ibs_stats
    sample_germplasm_dedup
    plink_pca
    PLOT_PLINK_PCA
} from './stats_sample.nf'

include { vcftools_vcf_qc_r; assess_plink_debug_plots; dumpnice_vcf_qc_assess } from './stats_assess.nf'

include {
    plot_plink2_population_structure
    report_plink2_tagsnp
    plot_plink2_snp_qc
    plot_cnv_calls
    plot_genetic_map
    report_hapmap_table
} from './stats_integrated.nf'

include {
    report_plink_chr_variant_counts
    plot_thin_common_chr_variant_compare
} from './stats_chr_report.nf'


// --- Genotype statistics workflows ---

workflow test_plink_stats {
    take:
    ld
    ld_cross
    // Expect a channel: [ val(id), path(smiss) ]
    smiss
    scount
    vmiss
    gcount
    afreq
    hardy

    main:
    // prepare input
    def idxstats = file("${params.output_dir}/vmap4_v1_idxstat_summary.txt")
    def group_file = file("${params.output_dir}/sample_group.txt")
    def tbm_dir = file("${params.home_dir}/00data/05taxaBamMap")
    def dbone_dir = file("${params.user_dir}/git/DBone/Service/src/main/resources/raw/20251208")
    def tbm = smiss.map { id, chr, _smiss_path ->
        def subgenome_tbm = ''
        if (chr.startsWith('sub')) {
            if (id == "Others") {
                subgenome_tbm = "${params.home_dir}/00data/05taxaBamMap/all.ALL.taxaBamMap.txt"
            } else {
                subgenome_tbm = "${params.home_dir}/00data/05taxaBamMap/all.${id}.taxaBamMap.txt"
            }
        } else {
            subgenome_tbm = getTaxaBamMapFile_v1(chr, params.home_dir)
        }
        subgenome_tbm = file(subgenome_tbm)
        tuple(id, chr, subgenome_tbm)
    }
    // 1 sample stats
    def miss_out = sample_missing_stats(smiss)
    def cov_out = sample_coverage_stats(smiss, tbm)
    def het_out = sample_heterozygosity_stats(scount, smiss)
    def mr_out = sample_mapping_rate_stats(idxstats, group_file, smiss)
    // 1.1 kinship
    def ref_ibs_out = sample_ref_ibs_stats(group_file, idxstats, scount, smiss)
    // 1.2 germplasm dedup
    def dedup_out = sample_germplasm_dedup(smiss, tbm_dir, dbone_dir)

    // collect sample info and thresholds
    def sample_info = miss_out.info.combine(het_out.info).combine(mr_out.info).combine(ref_ibs_out.info).combine(cov_out.info)
    def sample_th = miss_out.th.combine(het_out.th).combine(dedup_out.th).combine(cov_out.th)
    
    // 2 variant stats
    // 2.1 Basic Variant Stats (Missing, MAF)
    def vmiss_out = variant_missing_stats(vmiss)
    def maf_out = variant_maf_stats(afreq)
    variant_mac_stats(gcount)
    variant_ld_decay_plot(ld)
    variant_ld_crosschr_plot(ld_cross)

    // emit:
}

workflow plink_stats {
    take:
    // Expect a channel: [ val(id), path(smiss) ]
    smiss
    scount
    vmiss
    gcount
    afreq
    hardy
    popdep

    main:
    // def stats_out = sample_missing_stats(smiss)
    def vmiss_out = variant_missing_stats(vmiss)
    def maf_out = variant_maf_stats(afreq)
    // def popdep_stats_out = variant_popdep_stats(popdep, vmiss)
    def popdep_maha_out = variant_popdep_mahalanobis(popdep)

    // emit:
}
