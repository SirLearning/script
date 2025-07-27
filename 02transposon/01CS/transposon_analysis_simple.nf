#!/usr/bin/env nextflow

/*
 * Transposon Analysis Pipeline
 * 
 * This pipeline processes transposon data and generates composition statistics
 */

nextflow.enable.dsl = 2

// Parameters
params.input_data = "newplot.txt"
params.gff_file = "chr1A.gff3"
params.lib_file = "*.lib"
params.output_dir = "results"
params.help = false

/*
 * Process 1: Unified Transposon Analysis
 */
process UNIFIED_TRANSPOSON_ANALYSIS {
    tag "unified_analysis"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path input_file
    path script_file
    path gff_file
    path lib_file
    
    output:
    path "composition_stats.txt", optional: true
    path "composition_summary.txt", optional: true
    path "gff_category_counts.txt", optional: true
    path "gff_analysis_report.txt", optional: true
    path "lib_stats_*.txt", optional: true
    path "qc_report.txt", optional: true
    path "data_validation.txt", optional: true
    path "transposon_analysis_report.html"
    path "transposon_analysis_report.txt"
    path "analysis_results.json"
    
    script:
    """
    # Prepare arguments for the unified script
    ARGS="--mode all --input ${input_file} --output ."
    
    # Add GFF file if it exists and is not empty
    if [ -f "${gff_file}" ] && [ -s "${gff_file}" ]; then
        ARGS="\$ARGS --gff ${gff_file}"
    fi
    
    # Add library file if it exists and is not empty
    if [ -f "${lib_file}" ] && [ -s "${lib_file}" ]; then
        ARGS="\$ARGS --lib ${lib_file}"
    fi
    
    # Run the unified analysis script
    python3 ${script_file} \$ARGS
    
    # Ensure all expected output files exist (create empty ones if missing)
    touch composition_stats.txt composition_summary.txt
    touch gff_category_counts.txt gff_analysis_report.txt
    touch qc_report.txt data_validation.txt
    touch transposon_analysis_report.html transposon_analysis_report.txt
    touch analysis_results.json
    """
}

/*
 * Workflow definition
 */
workflow {
    // Show help if requested
    if (params.help) {
        println """
        Transposon Analysis Pipeline (Unified Version)
        
        Usage:
        nextflow run transposon_analysis_simple.nf [options]
        
        Options:
        --input_data      Input data file (default: newplot.txt)
        --gff_file        GFF3 annotation file (default: chr1A.gff3) [optional]
        --lib_file        Library file (default: *.lib) [optional]
        --output_dir      Output directory (default: results)
        --help            Show this help message
        
        Features:
        - Unified analysis using a single Python script
        - Transposon composition analysis
        - GFF3 category counting (if file exists)
        - Library statistics (if file exists)
        - Quality control and validation
        - Comprehensive HTML and text reports
        """
        exit 0
    }
    
    // Print pipeline info
    println """
    T R A N S P O S O N   A N A L Y S I S   P I P E L I N E   ( U N I F I E D )
    =======================================================================
    Input data        : ${params.input_data}
    GFF file          : ${params.gff_file}
    Library file      : ${params.lib_file}
    Output directory  : ${params.output_dir}
    """
    
    // Define input channels
    input_data_ch = Channel.fromPath(params.input_data, checkIfExists: true)
    unified_script_ch = Channel.fromPath("transposon_analyzer.py", checkIfExists: true)
    gff_file_ch = Channel.fromPath(params.gff_file, checkIfExists: false).ifEmpty(file("empty.gff"))
    lib_file_ch = Channel.fromPath(params.lib_file, checkIfExists: false).ifEmpty(file("empty.lib"))
    
    // Run unified analysis
    UNIFIED_TRANSPOSON_ANALYSIS(input_data_ch, unified_script_ch, gff_file_ch, lib_file_ch)
}
