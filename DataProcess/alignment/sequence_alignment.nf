#!/usr/bin/env nextflow

/*
 * Sequence Alignment Pipeline for wheat genome analysis
 * Based on the methods from 00seq_data_process.ipynb
 * Author: Generated from Vmap4 workflow
 * Date: 2025-08-11
 */

nextflow.enable.dsl=2

// Parameters
params.fastq_dir = null
params.reference = null
params.output_dir = "alignment_output"
params.threads = 20
params.memory = "4G"
params.help = false

// BWA alignment parameters
params.use_bwa_mem2 = true  // Use bwa-mem2 by default (faster)
params.read_group_platform = "illumina"
params.read_group_id = "Triticum"

// Sample parameters
params.sample_list = null
params.fq_list = null

def helpMessage() {
    log.info """
    ========================================
    Sequence Alignment Pipeline
    ========================================
    
    Usage:
        nextflow run sequence_alignment.nf --fastq_dir <fastq_directory> --reference <ref.fa> --sample_list <samples.txt> --fq_list <fqlist.txt>
    
    Required parameters:
        --fastq_dir         Directory containing FASTQ files
        --reference         Reference genome fasta file (should be indexed)
        --sample_list       File containing sample names (one per line)
        --fq_list           File containing FASTQ file prefixes (one per line)
    
    Optional parameters:
        --output_dir        Output directory (default: alignment_output)
        --threads           Number of threads (default: 20)
        --memory            Memory allocation per sort process (default: 4G)
        --use_bwa_mem2      Use bwa-mem2 instead of bwa (default: true)
        --read_group_platform   Read group platform (default: illumina)
        --read_group_id     Read group ID (default: Triticum)
    
    Example:
        nextflow run sequence_alignment.nf \\
            --fastq_dir /path/to/fastq \\
            --reference /path/to/reference.fa.gz \\
            --sample_list samples.txt \\
            --fq_list fqlist.txt \\
            --threads 20
    """.stripIndent()
}

workflow alignment_workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    // Parameter validation
    if (!params.fastq_dir) {
        log.error "FASTQ directory is required. Use --fastq_dir"
        exit 1
    }
    if (!params.reference) {
        log.error "Reference genome file is required. Use --reference"
        exit 1
    }
    if (!params.sample_list) {
        log.error "Sample list file is required. Use --sample_list"
        exit 1
    }
    if (!params.fq_list) {
        log.error "FASTQ list file is required. Use --fq_list"
        exit 1
    }

    log.info """\
    ========================================
    Sequence Alignment Pipeline
    ========================================
    FASTQ directory  : ${params.fastq_dir}
    Reference genome : ${params.reference}
    Sample list      : ${params.sample_list}
    FASTQ list       : ${params.fq_list}
    Output directory : ${params.output_dir}
    Threads          : ${params.threads}
    Memory per sort  : ${params.memory}
    Use BWA-MEM2     : ${params.use_bwa_mem2}
    ========================================
    """.stripIndent()

    // Read sample and FASTQ lists
    sample_ch = Channel
        .fromPath(params.sample_list)
        .splitText()
        .map { it.trim() }
        .filter { it != "" }

    fq_ch = Channel
        .fromPath(params.fq_list)
        .splitText()
        .map { it.trim() }
        .filter { it != "" }

    // Combine samples and FASTQ files
    sample_fq_ch = sample_ch.combine(fq_ch)

    // Run alignment workflow
    alignment_results = align_reads(sample_fq_ch)
    
    // Generate alignment statistics
    alignment_stats = generate_alignment_stats(alignment_results.collect())
}

process align_reads {
    tag "${sample}_${fq_prefix}"
    memory { task.attempt < 3 ? params.memory : (task.attempt * 2).toString() + 'G' }
    cpus params.threads
    time '24.h'
    errorStrategy 'retry'
    maxRetries 3
    publishDir "${params.output_dir}/bam", mode: 'copy', pattern: "*.rmdup.bam*"
    publishDir "${params.output_dir}/logs", mode: 'copy', pattern: "*.log"
    
    input:
    tuple val(sample), val(fq_prefix)
    
    output:
    tuple val(sample), path("${fq_prefix}.rmdup.bam"), path("${fq_prefix}.rmdup.bam.bai"), emit: bam_files
    path "${fq_prefix}_alignment.log", emit: log
    
    script:
    def aligner = params.use_bwa_mem2 ? "bwa-mem2" : "bwa"
    def read_group = "@RG\\tID:${params.read_group_id}\\tPL:${params.read_group_platform}\\tSM:${sample}"
    """
    echo "Starting alignment for sample: ${sample}, FASTQ prefix: ${fq_prefix}" > ${fq_prefix}_alignment.log
    echo "Date: \$(date)" >> ${fq_prefix}_alignment.log
    echo "Aligner: ${aligner}" >> ${fq_prefix}_alignment.log
    echo "Reference: ${params.reference}" >> ${fq_prefix}_alignment.log
    echo "Threads: ${params.threads}" >> ${fq_prefix}_alignment.log
    echo "Memory: ${params.memory}" >> ${fq_prefix}_alignment.log
    echo "" >> ${fq_prefix}_alignment.log
    
    # Check if FASTQ files exist
    if [ ! -f "${params.fastq_dir}/${fq_prefix}/${fq_prefix}_f1.fastq.gz" ]; then
        echo "Error: Forward FASTQ file not found: ${params.fastq_dir}/${fq_prefix}/${fq_prefix}_f1.fastq.gz" >> ${fq_prefix}_alignment.log
        exit 1
    fi
    
    if [ ! -f "${params.fastq_dir}/${fq_prefix}/${fq_prefix}_r2.fastq.gz" ]; then
        echo "Error: Reverse FASTQ file not found: ${params.fastq_dir}/${fq_prefix}/${fq_prefix}_r2.fastq.gz" >> ${fq_prefix}_alignment.log
        exit 1
    fi
    
    echo "FASTQ files verified" >> ${fq_prefix}_alignment.log
    echo "Forward: ${params.fastq_dir}/${fq_prefix}/${fq_prefix}_f1.fastq.gz" >> ${fq_prefix}_alignment.log
    echo "Reverse: ${params.fastq_dir}/${fq_prefix}/${fq_prefix}_r2.fastq.gz" >> ${fq_prefix}_alignment.log
    echo "" >> ${fq_prefix}_alignment.log
    
    # Step 1: Alignment
    echo "Step 1: Running alignment..." >> ${fq_prefix}_alignment.log
    start_time=\$(date +%s)
    
    ${aligner} mem -t ${params.threads} \\
        -R "${read_group}" \\
        ${params.reference} \\
        ${params.fastq_dir}/${fq_prefix}/${fq_prefix}_f1.fastq.gz \\
        ${params.fastq_dir}/${fq_prefix}/${fq_prefix}_r2.fastq.gz \\
        | samtools view -S -b -> ${fq_prefix}.bam
    
    alignment_time=\$(date +%s)
    echo "Alignment completed in \$((alignment_time - start_time)) seconds" >> ${fq_prefix}_alignment.log
    
    # Step 2: Sort by name for fixmate
    echo "Step 2: Sorting by name for fixmate..." >> ${fq_prefix}_alignment.log
    samtools sort -n -m ${params.memory} -@ ${params.threads} \\
        -o ${fq_prefix}.namesort.bam -O bam ${fq_prefix}.bam
    
    namesort_time=\$(date +%s)
    echo "Name sorting completed in \$((namesort_time - alignment_time)) seconds" >> ${fq_prefix}_alignment.log
    
    # Step 3: Fix mate information
    echo "Step 3: Fixing mate information..." >> ${fq_prefix}_alignment.log
    samtools fixmate -@ ${params.threads} -m ${fq_prefix}.namesort.bam ${fq_prefix}.fixmate.bam
    rm -f ${fq_prefix}.namesort.bam
    
    fixmate_time=\$(date +%s)
    echo "Fixmate completed in \$((fixmate_time - namesort_time)) seconds" >> ${fq_prefix}_alignment.log
    
    # Step 4: Sort by position
    echo "Step 4: Sorting by position..." >> ${fq_prefix}_alignment.log
    samtools sort -m ${params.memory} -@ ${params.threads} \\
        -o ${fq_prefix}.fixmate.pos.bam -O bam ${fq_prefix}.fixmate.bam
    rm -f ${fq_prefix}.fixmate.bam
    
    possort_time=\$(date +%s)
    echo "Position sorting completed in \$((possort_time - fixmate_time)) seconds" >> ${fq_prefix}_alignment.log
    
    # Step 5: Mark and remove duplicates
    echo "Step 5: Marking and removing duplicates..." >> ${fq_prefix}_alignment.log
    samtools markdup -@ ${params.threads} -r ${fq_prefix}.fixmate.pos.bam ${fq_prefix}.rmdup.bam
    rm -f ${fq_prefix}.fixmate.pos.bam
    rm -f ${fq_prefix}.bam
    
    markdup_time=\$(date +%s)
    echo "Duplicate removal completed in \$((markdup_time - possort_time)) seconds" >> ${fq_prefix}_alignment.log
    
    # Step 6: Index final BAM
    echo "Step 6: Indexing final BAM..." >> ${fq_prefix}_alignment.log
    samtools index -@ ${params.threads} ${fq_prefix}.rmdup.bam
    
    index_time=\$(date +%s)
    echo "Indexing completed in \$((index_time - markdup_time)) seconds" >> ${fq_prefix}_alignment.log
    
    # Generate final statistics
    echo "" >> ${fq_prefix}_alignment.log
    echo "=== FINAL STATISTICS ===" >> ${fq_prefix}_alignment.log
    echo "Total processing time: \$((index_time - start_time)) seconds" >> ${fq_prefix}_alignment.log
    echo "Final BAM file: ${fq_prefix}.rmdup.bam" >> ${fq_prefix}_alignment.log
    echo "BAM file size: \$(ls -lh ${fq_prefix}.rmdup.bam | awk '{print \$5}')" >> ${fq_prefix}_alignment.log
    
    # Get alignment statistics
    samtools flagstat ${fq_prefix}.rmdup.bam >> ${fq_prefix}_alignment.log
    
    echo "Alignment pipeline completed successfully!" >> ${fq_prefix}_alignment.log
    """
}

process generate_alignment_stats {
    tag "alignment_stats"
    memory '8g'
    cpus 4
    time '2.h'
    publishDir "${params.output_dir}/stats", mode: 'copy'
    
    input:
    path(bam_files)
    
    output:
    path "alignment_summary.txt"
    path "sample_stats.csv"
    
    script:
    """
    echo "Alignment Pipeline Summary Report" > alignment_summary.txt
    echo "=================================" >> alignment_summary.txt
    echo "Date: \$(date)" >> alignment_summary.txt
    echo "Reference genome: ${params.reference}" >> alignment_summary.txt
    echo "Aligner used: ${params.use_bwa_mem2 ? 'bwa-mem2' : 'bwa'}" >> alignment_summary.txt
    echo "Threads per job: ${params.threads}" >> alignment_summary.txt
    echo "Memory per sort: ${params.memory}" >> alignment_summary.txt
    echo "" >> alignment_summary.txt
    
    # Count processed samples
    bam_count=\$(ls -1 *.rmdup.bam 2>/dev/null | wc -l)
    echo "Total samples processed: \$bam_count" >> alignment_summary.txt
    echo "" >> alignment_summary.txt
    
    # Create CSV header
    echo "Sample,Total_reads,Mapped_reads,Properly_paired,Duplicates,Mapping_rate,File_size_MB" > sample_stats.csv
    
    # Process each BAM file
    for bam in *.rmdup.bam; do
        if [ -f "\$bam" ]; then
            sample=\$(basename \$bam .rmdup.bam)
            echo "Processing \$sample..." >> alignment_summary.txt
            
            # Get flagstat information
            flagstat_output=\$(samtools flagstat \$bam)
            
            # Extract statistics
            total_reads=\$(echo "\$flagstat_output" | head -n1 | awk '{print \$1}')
            mapped_reads=\$(echo "\$flagstat_output" | grep "mapped (" | head -n1 | awk '{print \$1}')
            properly_paired=\$(echo "\$flagstat_output" | grep "properly paired" | awk '{print \$1}')
            duplicates=\$(echo "\$flagstat_output" | grep "duplicates" | awk '{print \$1}')
            
            # Calculate mapping rate
            if [ "\$total_reads" -gt 0 ]; then
                mapping_rate=\$(echo "scale=2; \$mapped_reads * 100 / \$total_reads" | bc -l)
            else
                mapping_rate="0"
            fi
            
            # Get file size in MB
            file_size=\$(ls -l \$bam | awk '{printf "%.2f", \$5/1024/1024}')
            
            # Write to CSV
            echo "\$sample,\$total_reads,\$mapped_reads,\$properly_paired,\$duplicates,\$mapping_rate,\$file_size" >> sample_stats.csv
            
            # Write to summary
            echo "  Total reads: \$total_reads" >> alignment_summary.txt
            echo "  Mapped reads: \$mapped_reads (\$mapping_rate%)" >> alignment_summary.txt
            echo "  Properly paired: \$properly_paired" >> alignment_summary.txt
            echo "  Duplicates: \$duplicates" >> alignment_summary.txt
            echo "  File size: \$file_size MB" >> alignment_summary.txt
            echo "" >> alignment_summary.txt
        fi
    done
    
    # Summary statistics
    echo "=== OVERALL STATISTICS ===" >> alignment_summary.txt
    if [ -f sample_stats.csv ] && [ \$(wc -l < sample_stats.csv) -gt 1 ]; then
        echo "Average mapping rate: \$(tail -n +2 sample_stats.csv | awk -F',' '{sum+=\$6; count++} END {if(count>0) printf "%.2f%%", sum/count; else print "N/A"}')" >> alignment_summary.txt
        echo "Total disk space used: \$(tail -n +2 sample_stats.csv | awk -F',' '{sum+=\$7} END {printf "%.2f GB", sum/1024}')" >> alignment_summary.txt
    fi
    echo "" >> alignment_summary.txt
    echo "Report generated on: \$(date)" >> alignment_summary.txt
    """
}
