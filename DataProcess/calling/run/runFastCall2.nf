#!/usr/bin/env nextflow

/*
 * High-throughput FastCall2 pipeline for variant calling
 * Author: SirLearning
 * Date: 2025-07-25
 */

nextflow.enable.dsl=2

// Parameters
params.reference = null
params.taxaBamMap = null
params.tiger_jar = null
params.samtools_path = null
params.output_dir = "fastcall2_output"
params.threads = 32
params.memory = "100g"
params.help = false

// FastCall2 disc parameters
params.disc_min_depth = 30
params.disc_min_qual = 20
params.disc_min_snp_count = 2
params.disc_min_allele_freq = 0.2
params.disc_min_coverage = 3
params.disc_max_missing = 0.8
params.disc_min_het_freq = 0.35
params.disc_max_het_freq = 0.2

// FastCall2 scan parameters
params.scan_min_depth = 30
params.scan_min_qual = 20
params.scan_p_value = 0.05

// Chromosome list (can be modified based on your reference)
params.chromosomes = [
    "1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D",
    "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D",
    "7A", "7B", "7D"
]

def helpMessage() {
    log.info """
    ========================================
    FastCall2 High-throughput Pipeline
    ========================================
    
    Usage:
        nextflow run runFastCall2.nf --reference <ref.fa> --taxaBamMap <map.txt> --tiger_jar <TIGER.jar> --samtools_path <samtools>
    
    Required parameters:
        --reference         Reference genome fasta file
        --taxaBamMap        Taxa-BAM mapping file
        --tiger_jar         Path to TIGER jar file
        --samtools_path     Path to samtools executable
    
    Optional parameters:
        --output_dir        Output directory (default: fastcall2_output)
        --threads           Number of threads (default: 32)
        --memory            Memory allocation (default: 100g)
        --chromosomes       List of chromosomes to process (default: wheat chromosomes)
    
    Disc parameters:
        --disc_min_depth    Minimum depth for disc (default: 30)
        --disc_min_qual     Minimum quality for disc (default: 20)
        --disc_min_snp_count    Minimum SNP count (default: 2)
        --disc_min_allele_freq  Minimum allele frequency (default: 0.2)
        --disc_min_coverage     Minimum coverage (default: 3)
        --disc_max_missing      Maximum missing rate (default: 0.8)
        --disc_min_het_freq     Minimum heterozygous frequency (default: 0.35)
        --disc_max_het_freq     Maximum heterozygous frequency (default: 0.2)
    
    Scan parameters:
        --scan_min_depth    Minimum depth for scan (default: 30)
        --scan_min_qual     Minimum quality for scan (default: 20)
        --scan_p_value      P-value threshold (default: 0.05)
    """.stripIndent()
}

workflow fastcall2_workflow {
    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    // Parameter validation
    if (!params.reference) {
        log.error "Reference genome file is required. Use --reference"
        exit 1
    }
    if (!params.taxaBamMap) {
        log.error "Taxa-BAM mapping file is required. Use --taxaBamMap"
        exit 1
    }
    if (!params.tiger_jar) {
        log.error "TIGER jar file is required. Use --tiger_jar"
        exit 1
    }
    if (!params.samtools_path) {
        log.error "Samtools path is required. Use --samtools_path"
        exit 1
    }

    log.info """\
    ========================================
    FastCall2 High-throughput Pipeline
    ========================================
    Reference genome : ${params.reference}
    Taxa-BAM map     : ${params.taxaBamMap}
    TIGER jar        : ${params.tiger_jar}
    Samtools path    : ${params.samtools_path}
    Output directory : ${params.output_dir}
    Threads          : ${params.threads}
    Memory           : ${params.memory}
    Chromosomes      : ${params.chromosomes.join(', ')}
    ========================================
    """.stripIndent()

    // Create output directories
    def output_dir = file(params.output_dir)
    def pipeline_info_dir = file("./pipeline_info")
    output_dir.mkdirs()
    pipeline_info_dir.mkdirs()

    // Create chromosome channel
    chromosome_ch = Channel.fromList(params.chromosomes)
    
    // Step 1: Run disc for each chromosome
    disc_results = fastcall2_disc(chromosome_ch)
    
    // Step 2: Run blib for each chromosome
    blib_results = fastcall2_blib(disc_results)
    
    // Step 3: Run scan for each chromosome
    scan_results = fastcall2_scan(blib_results)
    
    // Step 4: Collect all results
    collect_results(scan_results.collect())
}

process fastcall2_disc {
    tag "disc_${chromosome}"
    memory params.memory
    cpus params.threads
    time '24.h'
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${params.output_dir}/disc", mode: 'copy', pattern: "*.ing"
    
    input:
    val chromosome
    
    output:
    tuple val(chromosome), path("*.ing"), emit: disc_files
    path "disc_${chromosome}.log", emit: log
    
    script:
    """
    java -Xmx${params.memory} -jar ${params.tiger_jar} \\
        -app FastCall2 \\
        -mod disc \\
        -a ${params.reference} \\
        -b ${params.taxaBamMap} \\
        -c 0 \\
        -d ${params.disc_min_depth} \\
        -e ${params.disc_min_qual} \\
        -f ${params.disc_min_snp_count} \\
        -g ${params.disc_min_allele_freq} \\
        -h ${params.disc_min_coverage} \\
        -i ${params.disc_max_missing} \\
        -j ${params.disc_min_het_freq} \\
        -k ${params.disc_max_het_freq} \\
        -l ${chromosome} \\
        -m ${params.threads} \\
        -n ./ \\
        -o ${params.samtools_path} \\
        > disc_${chromosome}.log 2>&1
    """
}

process fastcall2_blib {
    tag "blib_${chromosome}"
    memory params.memory
    cpus params.threads
    time '24.h'
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${params.output_dir}/blib", mode: 'copy', pattern: "*.lib.gz"
    
    input:
    tuple val(chromosome), path(disc_files)
    
    output:
    tuple val(chromosome), path("*.lib.gz"), emit: blib_files
    path "blib_${chromosome}.log", emit: log
    
    script:
    """
    # Check if disc files exist
    if [ ! -f *.ing ]; then
        echo "Error: No .ing files found from disc step"
        exit 1
    fi
    
    java -Xmx${params.memory} -jar ${params.tiger_jar} \\
        -app FastCall2 \\
        -mod blib \\
        -a ${params.reference} \\
        -b 1 \\
        -c 2 \\
        -d ${params.threads} \\
        -e ./ \\
        -f ./ \\
        > blib_${chromosome}.log 2>&1
    """
}

process fastcall2_scan {
    tag "scan_${chromosome}"
    memory params.memory
    cpus params.threads
    time '24.h'
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${params.output_dir}/scan", mode: 'copy', pattern: "*.{vcf,vcf.gz}"
    
    input:
    tuple val(chromosome), path(blib_files)
    
    output:
    tuple val(chromosome), path("*.{vcf,vcf.gz}"), emit: vcf_files
    path "scan_${chromosome}.log", emit: log
    
    script:
    def lib_file = blib_files.find { it.name.endsWith('.lib.gz') }
    """
    # Verify that lib file exists
    if [ ! -f "${lib_file}" ]; then
        echo "Error: Library file ${lib_file} not found"
        exit 1
    fi
    
    java -Xmx${params.memory} -jar ${params.tiger_jar} \\
        -app FastCall2 \\
        -mod scan \\
        -a ${params.reference} \\
        -b ${params.taxaBamMap} \\
        -c ${lib_file} \\
        -d 1 \\
        -e 0 \\
        -f ${params.scan_min_depth} \\
        -g ${params.scan_min_qual} \\
        -h ${params.scan_p_value} \\
        -i ${params.samtools_path} \\
        -j ${params.threads} \\
        -k ./ \\
        > scan_${chromosome}.log 2>&1
    """
}

process collect_results {
    tag "collect_results"
    memory '16g'
    cpus 4
    time '4.h'
    publishDir "${params.output_dir}/final", mode: 'copy'
    
    input:
    path(vcf_files)
    
    output:
    path "merged_variants.vcf.gz", optional: true
    path "summary_stats.txt"
    
    script:
    """
    # Create output directory if it doesn't exist
    mkdir -p final_output
    
    # Check if VCF files exist
    vcf_count=\$(ls -1 *.vcf *.vcf.gz 2>/dev/null | wc -l)
    if [ \$vcf_count -eq 0 ]; then
        echo "Warning: No VCF files found to merge"
        touch merged_variants.vcf.gz
    else
        # Merge all VCF files
        if command -v bcftools &> /dev/null; then
            echo "Merging \$vcf_count VCF files using bcftools..."
            bcftools concat ${vcf_files} | bcftools sort -Oz -o merged_variants.vcf.gz
            bcftools index merged_variants.vcf.gz
            echo "VCF files merged successfully"
        else
            echo "Warning: bcftools not found. VCF files are available separately in scan directory."
            touch merged_variants.vcf.gz
        fi
    fi
    
    # Generate summary statistics
    echo "FastCall2 Pipeline Summary" > summary_stats.txt
    echo "=========================" >> summary_stats.txt
    echo "Date: \$(date)" >> summary_stats.txt
    echo "Number of chromosomes processed: ${params.chromosomes.size()}" >> summary_stats.txt
    echo "Chromosomes: ${params.chromosomes.join(', ')}" >> summary_stats.txt
    echo "Reference genome: ${params.reference}" >> summary_stats.txt
    echo "Taxa-BAM mapping: ${params.taxaBamMap}" >> summary_stats.txt
    echo "" >> summary_stats.txt
    echo "Parameters used:" >> summary_stats.txt
    echo "- Min depth: ${params.disc_min_depth}" >> summary_stats.txt
    echo "- Min quality: ${params.disc_min_qual}" >> summary_stats.txt
    echo "- Min allele frequency: ${params.disc_min_allele_freq}" >> summary_stats.txt
    echo "- P-value threshold: ${params.scan_p_value}" >> summary_stats.txt
    echo "" >> summary_stats.txt
    echo "VCF files processed: \$vcf_count" >> summary_stats.txt
    echo "" >> summary_stats.txt
    echo "Output files:" >> summary_stats.txt
    ls -la *.vcf* >> summary_stats.txt 2>/dev/null || echo "No VCF files found" >> summary_stats.txt
    """
}
