nextflow.enable.dsl=2

/*
 * Wheat integrated study design: compose processor (PLINK2/awk), GWAS, and stats plot modules.
 * Invoked from main.nf via RUN_WHEAT_INTEGRATED or tmp/wheat_integrated_from_plink.nf.
 */

include {
    listMergedSubgenomeTestPfileTuples
    hasMergedSubgenomeTestPfiles
    hasPlinkBasicInfoForMergedTests
    listMergedSubgenomeSnpQcPlotTuples
} from '../../modules/local/genotype/utils.nf'
include {
    plink2_pca
    plink2_tagsnp_prune
    mk_plink_basic_info
    awk_depth_cnv_call
} from '../../modules/local/genotype/processor.nf'
include {
    plot_plink2_population_structure
    report_plink2_tagsnp
    plot_plink2_snp_qc
    plot_cnv_calls
    plot_genetic_map
    report_hapmap_table
} from '../../modules/local/genotype/stats.nf'
include { plink2_gwas_glm; gcta_gwas; plot_gwas_association } from '../../modules/local/static/gwas/gwas.nf'

workflow WHEAT_STUDY_FROM_PLINK {
    take:
    source_mod
    wheat_mod

    main:
    def proc_dir = params.process_dir ?: "${params.output_dir}/${params.job}/process/${source_mod}"
    if (!hasMergedSubgenomeTestPfiles(proc_dir)) {
        error "WHEAT_STUDY_FROM_PLINK: missing merged *_test.plink2 under ${proc_dir}"
    }
    log.info "${params.c_green}Wheat study ${wheat_mod} from PLINK2 (${source_mod}) in ${proc_dir}${params.c_reset}"
    params.wheat_plink_source_mod = source_mod
    params.wheat_integrated_mod = wheat_mod
    WHEAT_STUDY_ON_PFILES(channel.from(listMergedSubgenomeTestPfileTuples(proc_dir)), wheat_mod)
}

workflow WHEAT_STUDY_ON_PFILES {
    take:
    merged_pfile
    wheat_mod

    main:
    if (wheat_mod == 'wheat_pca_tsne') {
        def pca_out = plink2_pca(merged_pfile)
        pca_out.pca
            .map { id, chr, eigenvec, eigenval ->
                tuple(id, chr, eigenvec, eigenval, id, wheat_mod)
            }
            .set { ch_plot }
        plot_plink2_population_structure(
            ch_plot.map { id, chr, ev, ea, prefix, pm -> tuple(id, chr, ev, ea) },
            ch_plot.map { id, chr, ev, ea, prefix, pm -> prefix },
            ch_plot.map { id, chr, ev, ea, prefix, pm -> pm },
        )
    } else if (wheat_mod == 'wheat_tagsnp') {
        def tag_out = plink2_tagsnp_prune(merged_pfile)
        tag_out.prune
            .map { id, chr, prune_in -> tuple(id, chr, prune_in, id, wheat_mod) }
            .set { ch_tag }
        report_plink2_tagsnp(
            ch_tag.map { id, chr, pi, prefix, pm -> tuple(id, chr, pi) },
            ch_tag.map { id, chr, pi, prefix, pm -> prefix },
            ch_tag.map { id, chr, pi, prefix, pm -> pm },
        )
    } else if (wheat_mod == 'wheat_snp_qc') {
        def proc_dir = params.process_dir ?: "${params.output_dir}/${params.job}/process/${params.wheat_plink_source_mod ?: params.mod}"
        def ch_qc
        if (hasPlinkBasicInfoForMergedTests(proc_dir)) {
            log.info "${params.c_green}SNP QC plots from existing mk_plink_basic_info outputs in ${proc_dir}/variant${params.c_reset}"
            ch_qc = channel.from(listMergedSubgenomeSnpQcPlotTuples(proc_dir))
        } else {
            log.info "${params.c_yellow}Running mk_plink_basic_info before SNP QC plots (no variant/*.info.afreq yet)${params.c_reset}"
            def basic = mk_plink_basic_info(merged_pfile)
            ch_qc = basic.afreq
                .join(basic.vmiss)
                .map { id, chr, af, id2, chr2, vm ->
                    tuple(af, vm, id)
                }
        }
        plot_plink2_snp_qc(
            ch_qc.map { af, vm, prefix -> af },
            ch_qc.map { af, vm, prefix -> vm },
            ch_qc.map { af, vm, prefix -> prefix },
        )
    } else {
        error "WHEAT_STUDY_ON_PFILES: unsupported wheat_mod '${wheat_mod}'."
    }
}

workflow RUN_WHEAT_INTEGRATED {
    main:
    if (!params.output_dir) {
        error "RUN_WHEAT_INTEGRATED: params.output_dir is required."
    }
    if (!params.job) {
        error "RUN_WHEAT_INTEGRATED: params.job is required."
    }
    def prefix = "${params.mod}"

    if (params.mod in ['wheat_pca_tsne', 'wheat_tagsnp', 'wheat_snp_qc'] && params.wheat_plink_source_mod) {
        WHEAT_STUDY_FROM_PLINK(params.wheat_plink_source_mod, params.mod)
    } else if (params.mod == 'wheat_pca_tsne' && params.wheat_plink_eigenvec && params.wheat_plink_eigenval) {
        plot_plink2_population_structure(
            channel.of(tuple('custom', 'custom', file(params.wheat_plink_eigenvec), file(params.wheat_plink_eigenval))),
            channel.of(prefix),
            channel.of(params.mod),
        )
    } else if (params.mod == 'wheat_snp_qc') {
        if (!params.wheat_afreq_input || !params.wheat_vmiss_input) {
            error "wheat_snp_qc requires wheat_afreq_input and wheat_vmiss_input, or wheat_plink_source_mod."
        }
        plot_plink2_snp_qc(file(params.wheat_afreq_input), file(params.wheat_vmiss_input), prefix)
    } else if (params.mod == 'wheat_tagsnp' && params.wheat_prune_in) {
        report_plink2_tagsnp(
            channel.of(tuple('custom', 'custom', file(params.wheat_prune_in))),
            channel.of(prefix),
            channel.of(params.mod),
        )
    } else if (params.mod == 'wheat_cnv') {
        if (!params.wheat_table_input) {
            error "wheat_cnv requires params.wheat_table_input (depth matrix)."
        }
        def cnv_out = awk_depth_cnv_call(file(params.wheat_table_input), prefix)
        plot_cnv_calls(cnv_out.cnv, prefix)
    } else if (params.mod == 'wheat_genetic_map') {
        if (!params.wheat_map_input) {
            error "wheat_genetic_map requires params.wheat_map_input (precomputed map)."
        }
        plot_genetic_map(file(params.wheat_map_input), prefix)
    } else if (params.mod == 'wheat_hapmap') {
        if (!params.wheat_hapmap_input) {
            error "wheat_hapmap requires params.wheat_hapmap_input."
        }
        report_hapmap_table(file(params.wheat_hapmap_input), prefix)
    } else if (params.mod == 'wheat_gwas') {
        if (!params.wheat_gwas_phenotype) {
            error "wheat_gwas requires params.wheat_gwas_phenotype."
        }
        if (params.wheat_gwas_association) {
            plot_gwas_association(file(params.wheat_gwas_association), params.wheat_gwas_tool, prefix)
        } else if (params.wheat_gwas_tool == 'gcta' && params.wheat_plink_source_mod) {
            def ch_p = channel.from(listMergedSubgenomeTestPfileTuples(
                params.process_dir ?: "${params.output_dir}/${params.job}/process/${params.wheat_plink_source_mod}"
            ))
            def gcta = gcta_gwas(ch_p, file(params.wheat_gwas_phenotype), params.wheat_gwas_trait)
            plot_gwas_association(gcta.glm, gcta.source, prefix)
        } else if (params.wheat_plink_source_mod) {
            def proc_dir = params.process_dir ?: "${params.output_dir}/${params.job}/process/${params.wheat_plink_source_mod}"
            def ch_p = channel.from(listMergedSubgenomeTestPfileTuples(proc_dir))
            def glm = plink2_gwas_glm(ch_p, file(params.wheat_gwas_phenotype), params.wheat_gwas_trait)
            plot_gwas_association(glm.glm, glm.source, prefix)
        } else {
            error "wheat_gwas requires wheat_plink_source_mod or wheat_gwas_association."
        }
    } else if (params.mod == 'wheat_kgwas') {
        if (!params.wheat_kgwas_glm) {
            error "wheat_kgwas requires wheat_kgwas_glm (plink2 --glm output)."
        }
        plot_gwas_association(file(params.wheat_kgwas_glm), 'plink2', prefix)
    } else {
        error "Unknown wheat mod '${params.mod}'."
    }
}
