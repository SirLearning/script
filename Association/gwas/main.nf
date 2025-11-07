#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- Include GWAS models ---
include { RUN_GWAS_PLINK } from './plink'
include { RUN_GWAS_GAPIT } from './gapit'
// include { RUN_GWAS_RMVP } from './rmvp' // Example for another model

// --- GWAS Execution Workflow ---
workflow GWAS_RUN {
    take:
        ch_input // channel: [ val(meta), path(vcf), path(pca), path(pheno) ]

    main:
        ch_input.branch {
            plink: params.gwas_models.split(',').contains('plink')
            gapit: params.gwas_models.split(',').contains('gapit')
            // rmvp: params.gwas_models.split(',').contains('rmvp')
            other: true
        }
        .set { gwas_ch }

        def plink_res = Channel.empty()
        def gapit_res = Channel.empty()

        if (gwas_ch.plink) {
            RUN_GWAS_PLINK(gwas_ch.plink)
            plink_res = RUN_GWAS_PLINK.out.plink_linear
        }
        if (gwas_ch.gapit) {
            RUN_GWAS_GAPIT(gwas_ch.gapit)
            gapit_res = RUN_GWAS_GAPIT.out.gapit_outputs
        }
        // if (gwas_ch.rmvp) {
        //     RUN_GWAS_RMVP(gwas_ch.rmvp)
        // }

        gwas_ch.other.view { "Skipping GWAS for models not in: ${params.gwas_models}" }

        plink_res.mix(gapit_res).set { all_results }

    emit:
        gwas_results = all_results
}
