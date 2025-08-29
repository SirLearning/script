#!/usr/bin/env nextflow

/*
 * Sequence Data Processing Pipeline
 * Based on methods from 00seq_data_process.ipynb
 * Includes: data download, quality control, subsampling, and preparation
 * Author: Generated from Vmap4 workflow
 * Date: 2025-08-11
 */

nextflow.enable.dsl=2

// Parameters
params.sra_list = null
params.output_dir = "data_processing_output"
params.reference = null
params.threads = 20
params.memory = "8G"
params.help = false

// SRA download parameters
params.max_size = "200G"
params.use_prefetch = true

// Quality control parameters
params.run_fastqc = true

// Subsampling parameters
params.target_depth = 1.0
params.subsample = false

// Data preparation parameters
params.compress_output = true

def helpMessage() {
    log.info """
    ========================================
    Sequence Data Processing Pipeline
    ========================================
    
    Usage:
        nextflow run data_processing.nf --sra_list <sra_ids.txt> [options]
    
    Required parameters:
        --sra_list          File containing SRA accession IDs (one per line)
    
    Optional parameters:
        --output_dir        Output directory (default: data_processing_output)
        --threads           Number of threads (default: 20)
        --memory            Memory allocation (default: 8G)
        --reference         Reference genome for depth calculation
        --max_size          Maximum download size (default: 200G)
        --use_prefetch      Use prefetch before fasterq-dump (default: true)
        --run_fastqc        Run FastQC quality control (default: true)
        --subsample         Enable subsampling (default: false)
        --target_depth      Target sequencing depth for subsampling (default: 1.0)
        --compress_output   Compress FASTQ output (default: true)
    
    Example:
        nextflow run data_processing.nf \\
            --sra_list sra_accessions.txt \\
            --output_dir ./processed_data \\
            --reference reference.fa \\
            --subsample \\
            --target_depth 5.0
    """.stripIndent()
}

workflow data_processing_workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    // Parameter validation
    if (!params.sra_list) {
        log.error "SRA list file is required. Use --sra_list"
        exit 1
    }

    log.info """\
    ========================================
    Sequence Data Processing Pipeline
    ========================================
    SRA list         : ${params.sra_list}
    Output directory : ${params.output_dir}
    Threads          : ${params.threads}
    Memory           : ${params.memory}
    Max download size: ${params.max_size}
    Use prefetch     : ${params.use_prefetch}
    Run FastQC       : ${params.run_fastqc}
    Subsample        : ${params.subsample}
    Target depth     : ${params.target_depth}
    Compress output  : ${params.compress_output}
    ========================================
    """.stripIndent()

    // Read SRA accession list
    sra_ch = Channel
        .fromPath(params.sra_list)
        .splitText()
        .map { it.trim() }
        .filter { it != "" }

    // Step 1: Download SRA data
    downloaded_data = download_sra_data(sra_ch)
    
    // Step 2: Extract FASTQ files
    fastq_files = extract_fastq(downloaded_data)
    
    // Step 3: Quality control (optional)
    if (params.run_fastqc) {
        qc_results = run_fastqc(fastq_files)
    }
    
    // Step 4: Subsampling (optional)
    if (params.subsample) {
        if (!params.reference) {
            log.warn "Reference genome required for subsampling. Skipping subsampling step."
            final_fastq = fastq_files
        } else {
            final_fastq = subsample_reads(fastq_files)
        }
    } else {
        final_fastq = fastq_files
    }
    
    // Step 5: Generate MD5 checksums
    checksums = generate_md5_checksums(final_fastq)
    
    // Step 6: Generate processing summary
    summary = generate_processing_summary(final_fastq.collect(), checksums.collect())
}

process download_sra_data {
    tag "${sra_id}"
    memory params.memory
    cpus params.threads
    time '12.h'
    errorStrategy 'retry'
    maxRetries 3
    publishDir "${params.output_dir}/raw_data", mode: 'copy', pattern: "*.sra"
    publishDir "${params.output_dir}/logs", mode: 'copy', pattern: "*_download.log"
    
    input:
    val sra_id
    
    output:
    tuple val(sra_id), path("${sra_id}/${sra_id}.sra"), emit: sra_files
    path "${sra_id}_download.log", emit: log
    
    script:
    def prefetch_cmd = params.use_prefetch ? "prefetch -X ${params.max_size} -O ./ ${sra_id} &&" : ""
    """
    echo "Starting download for SRA: ${sra_id}" > ${sra_id}_download.log
    echo "Date: \$(date)" >> ${sra_id}_download.log
    echo "Max size: ${params.max_size}" >> ${sra_id}_download.log
    echo "Use prefetch: ${params.use_prefetch}" >> ${sra_id}_download.log
    echo "" >> ${sra_id}_download.log
    
    # Create directory
    mkdir -p ${sra_id}
    
    # Download SRA file
    if [ "${params.use_prefetch}" = "true" ]; then
        echo "Step 1: Downloading with prefetch..." >> ${sra_id}_download.log
        prefetch -X ${params.max_size} -O ./ ${sra_id} >> ${sra_id}_download.log 2>&1
        
        if [ ! -f "${sra_id}/${sra_id}.sra" ]; then
            echo "Error: SRA file not found after prefetch" >> ${sra_id}_download.log
            exit 1
        fi
    else
        echo "Skipping prefetch step" >> ${sra_id}_download.log
        # Create placeholder SRA file for direct fasterq-dump
        touch ${sra_id}/${sra_id}.sra
    fi
    
    echo "Download completed successfully" >> ${sra_id}_download.log
    echo "File size: \$(ls -lh ${sra_id}/${sra_id}.sra | awk '{print \$5}')" >> ${sra_id}_download.log
    """
}

process extract_fastq {
    tag "${sra_id}"
    memory params.memory
    cpus params.threads
    time '8.h'
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${params.output_dir}/fastq", mode: 'copy', pattern: "*.fastq.gz"
    publishDir "${params.output_dir}/logs", mode: 'copy', pattern: "*_extract.log"
    
    input:
    tuple val(sra_id), path(sra_file)
    
    output:
    tuple val(sra_id), path("${sra_id}_1.fastq.gz"), path("${sra_id}_2.fastq.gz"), emit: fastq_files
    path "${sra_id}_extract.log", emit: log
    
    script:
    def sra_path = params.use_prefetch ? sra_file : sra_id
    """
    echo "Starting FASTQ extraction for: ${sra_id}" > ${sra_id}_extract.log
    echo "Date: \$(date)" >> ${sra_id}_extract.log
    echo "Threads: ${params.threads}" >> ${sra_id}_extract.log
    echo "SRA file: ${sra_path}" >> ${sra_id}_extract.log
    echo "" >> ${sra_id}_extract.log
    
    # Extract FASTQ files
    echo "Running fasterq-dump..." >> ${sra_id}_extract.log
    fasterq-dump -e ${params.threads} --split-3 \\
        --defline-seq '@\$sn[_\$rn]/\$ri' \\
        --defline-qual '+' \\
        -O ./ ${sra_path} >> ${sra_id}_extract.log 2>&1
    
    # Check if extraction was successful
    if [ ! -f "${sra_id}_1.fastq" ] || [ ! -f "${sra_id}_2.fastq" ]; then
        echo "Error: FASTQ extraction failed" >> ${sra_id}_extract.log
        exit 1
    fi
    
    # Compress FASTQ files
    if [ "${params.compress_output}" = "true" ]; then
        echo "Compressing FASTQ files..." >> ${sra_id}_extract.log
        pigz -3 -p ${params.threads} ${sra_id}_1.fastq ${sra_id}_2.fastq
    else
        # Rename for consistency
        mv ${sra_id}_1.fastq ${sra_id}_1.fastq.gz
        mv ${sra_id}_2.fastq ${sra_id}_2.fastq.gz
    fi
    
    # Clean up SRA file if prefetch was used
    if [ "${params.use_prefetch}" = "true" ] && [ -f "${sra_file}" ]; then
        echo "Cleaning up SRA file..." >> ${sra_id}_extract.log
        rm -f ${sra_file}
        rmdir ${sra_id} 2>/dev/null || true
    fi
    
    echo "FASTQ extraction completed" >> ${sra_id}_extract.log
    echo "Forward reads: \$(ls -lh ${sra_id}_1.fastq.gz | awk '{print \$5}')" >> ${sra_id}_extract.log
    echo "Reverse reads: \$(ls -lh ${sra_id}_2.fastq.gz | awk '{print \$5}')" >> ${sra_id}_extract.log
    """
}

process run_fastqc {
    tag "${sra_id}"
    memory '4g'
    cpus 4
    time '2.h'
    publishDir "${params.output_dir}/quality_control", mode: 'copy'
    
    input:
    tuple val(sra_id), path(fastq1), path(fastq2)
    
    output:
    tuple val(sra_id), path("*_fastqc.html"), path("*_fastqc.zip"), emit: qc_files
    
    script:
    """
    fastqc -t 4 -o ./ ${fastq1} ${fastq2}
    """
}

process subsample_reads {
    tag "${sra_id}"
    memory params.memory
    cpus 4
    time '4.h'
    publishDir "${params.output_dir}/subsampled", mode: 'copy', pattern: "*.fastq.gz"
    publishDir "${params.output_dir}/logs", mode: 'copy', pattern: "*_subsample.log"
    
    input:
    tuple val(sra_id), path(fastq1), path(fastq2)
    
    output:
    tuple val(sra_id), path("${sra_id}_sub_1.fastq.gz"), path("${sra_id}_sub_2.fastq.gz"), emit: subsampled_fastq
    path "${sra_id}_subsample.log", emit: log
    
    script:
    """
    echo "Starting subsampling for: ${sra_id}" > ${sra_id}_subsample.log
    echo "Target depth: ${params.target_depth}" >> ${sra_id}_subsample.log
    echo "Reference: ${params.reference}" >> ${sra_id}_subsample.log
    echo "" >> ${sra_id}_subsample.log
    
    # Calculate total reads
    total_reads=\$(zcat ${fastq1} | wc -l | awk '{print \$1/4}')
    echo "Total reads: \$total_reads" >> ${sra_id}_subsample.log
    
    # Estimate genome size (for wheat, approximately 17Gb)
    genome_size=17000000000
    
    # Calculate read length (estimate from first read)
    read_length=\$(zcat ${fastq1} | head -n 2 | tail -n 1 | wc -c)
    echo "Estimated read length: \$read_length" >> ${sra_id}_subsample.log
    
    # Calculate current depth
    current_depth=\$(echo "scale=2; \$total_reads * \$read_length * 2 / \$genome_size" | bc -l)
    echo "Current depth: \$current_depth" >> ${sra_id}_subsample.log
    
    # Calculate sampling ratio
    if [ \$(echo "\$current_depth > ${params.target_depth}" | bc -l) -eq 1 ]; then
        sampling_ratio=\$(echo "scale=6; ${params.target_depth} / \$current_depth" | bc -l)
        target_reads=\$(echo "\$total_reads * \$sampling_ratio" | bc -l | cut -d'.' -f1)
        
        echo "Sampling ratio: \$sampling_ratio" >> ${sra_id}_subsample.log
        echo "Target reads: \$target_reads" >> ${sra_id}_subsample.log
        
        # Subsample using seqtk
        seqtk sample -s100 ${fastq1} \$target_reads > ${sra_id}_sub_1.fastq
        seqtk sample -s100 ${fastq2} \$target_reads > ${sra_id}_sub_2.fastq
        
        # Compress
        pigz -3 ${sra_id}_sub_1.fastq ${sra_id}_sub_2.fastq
        
        echo "Subsampling completed" >> ${sra_id}_subsample.log
    else
        echo "Current depth is lower than target. No subsampling needed." >> ${sra_id}_subsample.log
        # Just copy the original files
        cp ${fastq1} ${sra_id}_sub_1.fastq.gz
        cp ${fastq2} ${sra_id}_sub_2.fastq.gz
    fi
    """
}

process generate_md5_checksums {
    tag "${sra_id}"
    memory '2g'
    cpus 2
    time '1.h'
    publishDir "${params.output_dir}/checksums", mode: 'copy'
    
    input:
    tuple val(sra_id), path(fastq1), path(fastq2)
    
    output:
    tuple val(sra_id), path("${sra_id}_checksums.md5"), emit: checksums
    
    script:
    """
    md5sum ${fastq1} > ${sra_id}_checksums.md5
    md5sum ${fastq2} >> ${sra_id}_checksums.md5
    """
}

process generate_processing_summary {
    tag "processing_summary"
    memory '4g'
    cpus 2
    time '1.h'
    publishDir "${params.output_dir}/summary", mode: 'copy'
    
    input:
    path(fastq_files)
    path(checksum_files)
    
    output:
    path "processing_summary.txt"
    path "sample_info.csv"
    
    script:
    """
    echo "Sequence Data Processing Summary" > processing_summary.txt
    echo "===============================" >> processing_summary.txt
    echo "Date: \$(date)" >> processing_summary.txt
    echo "Pipeline version: ${workflow.manifest.version ?: 'N/A'}" >> processing_summary.txt
    echo "Nextflow version: ${workflow.nextflow.version}" >> processing_summary.txt
    echo "" >> processing_summary.txt
    echo "Parameters used:" >> processing_summary.txt
    echo "- Threads: ${params.threads}" >> processing_summary.txt
    echo "- Memory: ${params.memory}" >> processing_summary.txt
    echo "- Max download size: ${params.max_size}" >> processing_summary.txt
    echo "- Use prefetch: ${params.use_prefetch}" >> processing_summary.txt
    echo "- Run FastQC: ${params.run_fastqc}" >> processing_summary.txt
    echo "- Subsample: ${params.subsample}" >> processing_summary.txt
    echo "- Target depth: ${params.target_depth}" >> processing_summary.txt
    echo "- Compress output: ${params.compress_output}" >> processing_summary.txt
    echo "" >> processing_summary.txt
    
    # Count processed samples
    sample_count=\$(ls -1 *_1.fastq.gz 2>/dev/null | wc -l)
    echo "Total samples processed: \$sample_count" >> processing_summary.txt
    echo "" >> processing_summary.txt
    
    # Create CSV header
    echo "Sample,Forward_reads_file,Reverse_reads_file,Forward_size_MB,Reverse_size_MB,MD5_forward,MD5_reverse" > sample_info.csv
    
    # Process each sample
    for fq1 in *_1.fastq.gz; do
        if [ -f "\$fq1" ]; then
            # Extract sample name
            sample=\$(basename \$fq1 _1.fastq.gz)
            fq2="\${sample}_2.fastq.gz"
            
            if [ -f "\$fq2" ]; then
                # Get file sizes
                size1=\$(ls -l \$fq1 | awk '{printf "%.2f", \$5/1024/1024}')
                size2=\$(ls -l \$fq2 | awk '{printf "%.2f", \$5/1024/1024}')
                
                # Get MD5 checksums
                md5_file="\${sample}_checksums.md5"
                if [ -f "\$md5_file" ]; then
                    md5_1=\$(grep \$fq1 \$md5_file | awk '{print \$1}')
                    md5_2=\$(grep \$fq2 \$md5_file | awk '{print \$1}')
                else
                    md5_1="N/A"
                    md5_2="N/A"
                fi
                
                # Write to CSV
                echo "\$sample,\$fq1,\$fq2,\$size1,\$size2,\$md5_1,\$md5_2" >> sample_info.csv
                
                # Write to summary
                echo "Sample: \$sample" >> processing_summary.txt
                echo "  Forward reads: \$fq1 (\$size1 MB)" >> processing_summary.txt
                echo "  Reverse reads: \$fq2 (\$size2 MB)" >> processing_summary.txt
                echo "  MD5 forward: \$md5_1" >> processing_summary.txt
                echo "  MD5 reverse: \$md5_2" >> processing_summary.txt
                echo "" >> processing_summary.txt
            fi
        fi
    done
    
    # Summary statistics
    echo "=== OVERALL STATISTICS ===" >> processing_summary.txt
    if [ -f sample_info.csv ] && [ \$(wc -l < sample_info.csv) -gt 1 ]; then
        total_size=\$(tail -n +2 sample_info.csv | awk -F',' '{sum+=\$4+\$5} END {printf "%.2f", sum}')
        echo "Total data size: \$total_size MB" >> processing_summary.txt
        avg_size=\$(tail -n +2 sample_info.csv | awk -F',' '{sum+=\$4+\$5; count++} END {if(count>0) printf "%.2f", sum/count; else print "0"}')
        echo "Average sample size: \$avg_size MB" >> processing_summary.txt
    fi
    echo "" >> processing_summary.txt
    echo "Processing completed on: \$(date)" >> processing_summary.txt
    """
}
