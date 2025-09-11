#!/usr/bin/env nextflow

/*
 * Main Vmap4 Data Processing Pipeline
 * Entry point for all workflows based on 00seq_data_process.ipynb
 * Author: Generated from Vmap4 workflow
 * Date: 2025-08-11
 */

nextflow.enable.dsl=2

// Import workflow modules
include { data_processing_workflow } from './DataProcess/overall/data_processing.nf'
include { alignment_workflow } from './DataProcess/alignment/sequence_alignment.nf' 
include { depth_qc_workflow } from './DataProcess/overall/depth_qc.nf'
include { fastcall2_workflow } from './DataProcess/calling/run/runFastCall2.nf'
include { fastcall3_workflow } from './DataProcess/calling/run/runFastCall3.nf'
include { performance_workflow } from './DataProcess/calling/perf/fastcall.nf'

// Parameters
params.workflow_type = "alignment"  // data_processing, alignment, depth_qc, fastcall2, fastcall3, performance
params.help = false

def helpMessage() {
    log.info """
    ========================================
    Vmap4 Data Processing Pipeline
    ========================================
    
    Usage:
        nextflow run main.nf --workflow_type <workflow> [options]
    
    Workflow Types:
        data_processing     Download and process SRA data
        alignment          Align reads to reference genome
        depth_qc           Calculate depth and run quality control
        fastcall2          Run FastCall2 variant calling
        fastcall3          Run FastCall3 variant calling
        performance        Run performance analysis
    
    Global Options:
        --workflow_type     Type of workflow to run (required)
        --help              Show this help message
    
    Workflow-specific options:
        See individual workflow help messages for details
    
    Examples:
        # Run alignment workflow
        nextflow run main.nf --workflow_type alignment --fastq_dir ./fastq --reference genome.fa
        
        # Run performance analysis
        nextflow run main.nf --workflow_type performance --benchmark_mode true
    """.stripIndent()
}

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    log.info """\
    ========================================
    Vmap4 Data Processing Pipeline
    ========================================
    Workflow type: ${params.workflow_type}
    Nextflow version: ${workflow.nextflow.version}
    Pipeline version: ${workflow.manifest.version ?: 'N/A'}
    ========================================
    """.stripIndent()

    // Route to appropriate workflow based on type
    if (params.workflow_type == 'data_processing') {
        data_processing_workflow()
    } else if (params.workflow_type == 'alignment') {
        alignment_workflow()
    } else if (params.workflow_type == 'depth_qc') {
        depth_qc_workflow()
    } else if (params.workflow_type == 'fastcall2') {
        fastcall2_workflow()
    } else if (params.workflow_type == 'fastcall3') {
        fastcall3_workflow()
    } else if (params.workflow_type == 'performance') {
        performance_workflow()
    } else {
        log.error "Unknown workflow type: ${params.workflow_type}"
        log.error "Valid options: data_processing, alignment, depth_qc, fastcall2, fastcall3, performance"
        exit 1
    }
}

