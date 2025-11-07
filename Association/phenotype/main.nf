#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- Include modules ---
include { PHENOTYPE_PROCESS } from './process'
include { PHENOTYPE_STATS } from './stats'

// --- Phenotype Processing Workflow ---
workflow PHENOTYPE_PIPELINE {
    take:
        ch_pheno // channel: [ val(meta), path(pheno_file) ]

    main:
        // 1. Process phenotype data (e.g., from database)
        PHENOTYPE_PROCESS(ch_pheno)

        // 2. Perform statistical analysis (BLUP/BLUE)
        PHENOTYPE_STATS(PHENOTYPE_PROCESS.out.processed_pheno)

    emit:
        processed_pheno = PHENOTYPE_STATS.out.stats_pheno
}
