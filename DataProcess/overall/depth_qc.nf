#!/usr/bin/env nextflow

/*
 * Depth Calculation and Quality Control Pipeline
 * Based on methods from 00seq_data_process.ipynb
 * Includes: depth counting, indexing, quality assessment
 * Author: Generated from Vmap4 workflow
 * Date: 2025-08-11
 */

nextflow.enable.dsl=2

// Parameters
params.bam_dir = null
params.bam_list = null
params.output_dir = "depth_qc_output"
params.threads = 10
params.memory = "4G"
params.help = false

// Depth calculation parameters
params.use_mosdepth = true
params.window_size = 500
params.min_mapping_quality = 20

// Quality control parameters
params.run_flagstat = true
params.generate_coverage_plots = false

def helpMessage() {
    log.info """
    ========================================
    Depth Calculation and Quality Control Pipeline
    ========================================
    
    Usage:
        nextflow run depth_qc.nf --bam_dir <bam_directory> --bam_list <bam_files.txt>
    
    Required parameters:
        --bam_dir           Directory containing BAM files
        --bam_list          File containing BAM file names (one per line)
    
    Optional parameters:
        --output_dir        Output directory (default: depth_qc_output)
        --threads           Number of threads (default: 10)
        --memory            Memory allocation (default: 4G)
        --use_mosdepth      Use mosdepth for depth calculation (default: true)
        --window_size       Window size for depth calculation (default: 500)
        --min_mapping_quality  Minimum mapping quality (default: 20)
        --run_flagstat      Run samtools flagstat (default: true)
        --generate_coverage_plots  Generate coverage plots (default: false)
    
    Example:
        nextflow run depth_qc.nf \\
            --bam_dir /path/to/bam \\
            --bam_list bam_files.txt \\
            --threads 20 \\
            --generate_coverage_plots
    """.stripIndent()
}

workflow depth_qc_workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    // Parameter validation
    if (!params.bam_dir) {
        log.error "BAM directory is required. Use --bam_dir"
        exit 1
    }
    if (!params.bam_list) {
        log.error "BAM list file is required. Use --bam_list"
        exit 1
    }

    log.info """\
    ========================================
    Depth Calculation and Quality Control Pipeline
    ========================================
    BAM directory    : ${params.bam_dir}
    BAM list         : ${params.bam_list}
    Output directory : ${params.output_dir}
    Threads          : ${params.threads}
    Memory           : ${params.memory}
    Use mosdepth     : ${params.use_mosdepth}
    Window size      : ${params.window_size}
    Min mapping qual : ${params.min_mapping_quality}
    Run flagstat     : ${params.run_flagstat}
    Generate plots   : ${params.generate_coverage_plots}
    ========================================
    """.stripIndent()

    // Read BAM file list
    bam_ch = Channel
        .fromPath(params.bam_list)
        .splitText()
        .map { it.trim() }
        .filter { it != "" }

    // Step 1: Index BAM files if needed
    indexed_bams = index_bam_files(bam_ch)
    
    // Step 2: Calculate depth
    depth_results = calculate_depth(indexed_bams)
    
    // Step 3: Quality control (optional)
    if (params.run_flagstat) {
        qc_results = run_quality_control(indexed_bams)
    }
    
    // Step 4: Extract depth statistics
    depth_stats = extract_depth_stats(depth_results)
    
    // Step 5: Generate summary report
    summary_report = generate_depth_summary(
        depth_stats.collect(),
        params.run_flagstat ? qc_results.collect() : Channel.empty().collect()
    )
    
    // Step 6: Generate coverage plots (optional)
    if (params.generate_coverage_plots) {
        coverage_plots = generate_coverage_plots(depth_results.collect())
    }
}

process index_bam_files {
    tag "${bam_name}"
    memory params.memory
    cpus params.threads
    time '4.h'
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${params.output_dir}/indexed_bams", mode: 'copy', pattern: "*.bai"
    publishDir "${params.output_dir}/logs", mode: 'copy', pattern: "*_index.log"
    
    input:
    val bam_name
    
    output:
    tuple val(bam_name), path("${bam_name}"), path("${bam_name}.bai"), emit: indexed_bam
    path "${bam_name}_index.log", emit: log
    
    script:
    """
    echo "Processing BAM file: ${bam_name}" > ${bam_name}_index.log
    echo "Date: \$(date)" >> ${bam_name}_index.log
    echo "BAM directory: ${params.bam_dir}" >> ${bam_name}_index.log
    echo "" >> ${bam_name}_index.log
    
    # Check if BAM file exists
    if [ ! -f "${params.bam_dir}/${bam_name}" ]; then
        echo "Error: BAM file not found: ${params.bam_dir}/${bam_name}" >> ${bam_name}_index.log
        exit 1
    fi
    
    # Copy BAM file to working directory
    echo "Copying BAM file..." >> ${bam_name}_index.log
    cp "${params.bam_dir}/${bam_name}" ./
    
    # Check if index already exists
    if [ -f "${params.bam_dir}/${bam_name}.bai" ]; then
        echo "Index file already exists, copying..." >> ${bam_name}_index.log
        cp "${params.bam_dir}/${bam_name}.bai" ./
    else
        echo "Creating index file..." >> ${bam_name}_index.log
        samtools index -@ ${params.threads} ${bam_name}
    fi
    
    # Verify index was created
    if [ ! -f "${bam_name}.bai" ]; then
        echo "Error: Failed to create index file" >> ${bam_name}_index.log
        exit 1
    fi
    
    echo "Indexing completed successfully" >> ${bam_name}_index.log
    echo "BAM size: \$(ls -lh ${bam_name} | awk '{print \$5}')" >> ${bam_name}_index.log
    echo "Index size: \$(ls -lh ${bam_name}.bai | awk '{print \$5}')" >> ${bam_name}_index.log
    """
}

process calculate_depth {
    tag "${bam_name}"
    memory params.memory
    cpus params.threads
    time '6.h'
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${params.output_dir}/depth", mode: 'copy', pattern: "*.{summary.txt,mosdepth.*}"
    publishDir "${params.output_dir}/logs", mode: 'copy', pattern: "*_depth.log"
    
    input:
    tuple val(bam_name), path(bam_file), path(bai_file)
    
    output:
    tuple val(bam_name), path("${bam_name}.mosdepth.summary.txt"), emit: depth_summary
    path "${bam_name}.mosdepth.*", emit: depth_files
    path "${bam_name}_depth.log", emit: log
    
    script:
    def base_name = bam_name.replaceAll(/\.bam$/, "")
    """
    echo "Starting depth calculation for: ${bam_name}" > ${bam_name}_depth.log
    echo "Date: \$(date)" >> ${bam_name}_depth.log
    echo "Method: ${params.use_mosdepth ? 'mosdepth' : 'samtools depth'}" >> ${bam_name}_depth.log
    echo "Threads: ${params.threads}" >> ${bam_name}_depth.log
    echo "Window size: ${params.window_size}" >> ${bam_name}_depth.log
    echo "Min mapping quality: ${params.min_mapping_quality}" >> ${bam_name}_depth.log
    echo "" >> ${bam_name}_depth.log
    
    start_time=\$(date +%s)
    
    if [ "${params.use_mosdepth}" = "true" ]; then
        echo "Running mosdepth..." >> ${bam_name}_depth.log
        mosdepth -t ${params.threads} \\
                 -n \\
                 -x \\
                 -Q ${params.min_mapping_quality} \\
                 ${base_name} \\
                 ${bam_file} >> ${bam_name}_depth.log 2>&1
        
        # Check if mosdepth succeeded
        if [ ! -f "${base_name}.mosdepth.summary.txt" ]; then
            echo "Error: mosdepth failed to generate summary file" >> ${bam_name}_depth.log
            exit 1
        fi
        
        # Rename files to include .bam in the name for consistency
        for file in ${base_name}.mosdepth.*; do
            if [ -f "\$file" ]; then
                new_name=\$(echo \$file | sed 's/${base_name}/${bam_name}/')
                mv \$file \$new_name
            fi
        done
        
    else
        echo "Running samtools depth..." >> ${bam_name}_depth.log
        samtools depth -@ ${params.threads} \\
                      -q ${params.min_mapping_quality} \\
                      ${bam_file} | \\
        awk '{sum+=\$3; count++} END {
            if(count > 0) {
                mean_depth = sum/count;
                print "chrom\\tlength\\tbases\\tmean\\tmin\\tmax";
                print "total\\t" count "\\t" sum "\\t" mean_depth "\\t0\\t" max_depth;
            }
        }' > ${bam_name}.mosdepth.summary.txt
    fi
    
    end_time=\$(date +%s)
    elapsed_time=\$((end_time - start_time))
    
    echo "Depth calculation completed in \$elapsed_time seconds" >> ${bam_name}_depth.log
    
    # Extract key statistics
    if [ -f "${bam_name}.mosdepth.summary.txt" ]; then
        mean_depth=\$(tail -n 1 ${bam_name}.mosdepth.summary.txt | awk '{print \$4}')
        echo "Mean depth: \$mean_depth" >> ${bam_name}_depth.log
    fi
    
    echo "Output files:" >> ${bam_name}_depth.log
    ls -la ${bam_name}.mosdepth.* >> ${bam_name}_depth.log 2>/dev/null || echo "No mosdepth files found" >> ${bam_name}_depth.log
    """
}

process run_quality_control {
    tag "${bam_name}"
    memory '4g'
    cpus 4
    time '2.h'
    publishDir "${params.output_dir}/qc", mode: 'copy', pattern: "*.{flagstat,stats}"
    publishDir "${params.output_dir}/logs", mode: 'copy', pattern: "*_qc.log"
    
    input:
    tuple val(bam_name), path(bam_file), path(bai_file)
    
    output:
    tuple val(bam_name), path("${bam_name}.flagstat"), path("${bam_name}.stats"), emit: qc_files
    path "${bam_name}_qc.log", emit: log
    
    script:
    """
    echo "Starting quality control for: ${bam_name}" > ${bam_name}_qc.log
    echo "Date: \$(date)" >> ${bam_name}_qc.log
    echo "" >> ${bam_name}_qc.log
    
    # Run samtools flagstat
    echo "Running samtools flagstat..." >> ${bam_name}_qc.log
    samtools flagstat ${bam_file} > ${bam_name}.flagstat
    
    # Run samtools stats
    echo "Running samtools stats..." >> ${bam_name}_qc.log
    samtools stats ${bam_file} > ${bam_name}.stats
    
    # Extract key QC metrics
    echo "=== QC SUMMARY ===" >> ${bam_name}_qc.log
    
    # From flagstat
    total_reads=\$(grep "total" ${bam_name}.flagstat | head -1 | awk '{print \$1}')
    mapped_reads=\$(grep "mapped (" ${bam_name}.flagstat | head -1 | awk '{print \$1}')
    duplicates=\$(grep "duplicates" ${bam_name}.flagstat | awk '{print \$1}')
    
    if [ "\$total_reads" -gt 0 ]; then
        mapping_rate=\$(echo "scale=2; \$mapped_reads * 100 / \$total_reads" | bc -l)
        duplicate_rate=\$(echo "scale=2; \$duplicates * 100 / \$total_reads" | bc -l)
    else
        mapping_rate="0"
        duplicate_rate="0"
    fi
    
    echo "Total reads: \$total_reads" >> ${bam_name}_qc.log
    echo "Mapped reads: \$mapped_reads (\$mapping_rate%)" >> ${bam_name}_qc.log
    echo "Duplicates: \$duplicates (\$duplicate_rate%)" >> ${bam_name}_qc.log
    
    # From stats
    insert_size=\$(grep "^SN" ${bam_name}.stats | grep "insert size average" | cut -f3)
    echo "Average insert size: \$insert_size" >> ${bam_name}_qc.log
    
    echo "Quality control completed" >> ${bam_name}_qc.log
    """
}

process extract_depth_stats {
    tag "${bam_name}"
    memory '2g'
    cpus 2
    time '1.h'
    publishDir "${params.output_dir}/depth_stats", mode: 'copy'
    
    input:
    tuple val(bam_name), path(depth_summary)
    
    output:
    tuple val(bam_name), path("${bam_name}_depth_stats.txt"), emit: depth_stats
    
    script:
    """
    echo "Extracting depth statistics for: ${bam_name}" > ${bam_name}_depth_stats.txt
    echo "Date: \$(date)" >> ${bam_name}_depth_stats.txt
    echo "" >> ${bam_name}_depth_stats.txt
    
    # Extract statistics from mosdepth summary
    if [ -f "${depth_summary}" ]; then
        echo "=== DEPTH STATISTICS ===" >> ${bam_name}_depth_stats.txt
        
        # Get overall statistics (usually the last line for total)
        total_line=\$(tail -n 1 ${depth_summary})
        
        if [ -n "\$total_line" ]; then
            chrom=\$(echo "\$total_line" | awk '{print \$1}')
            length=\$(echo "\$total_line" | awk '{print \$2}')
            bases=\$(echo "\$total_line" | awk '{print \$3}')
            mean_depth=\$(echo "\$total_line" | awk '{print \$4}')
            
            echo "Sample: ${bam_name}" >> ${bam_name}_depth_stats.txt
            echo "Chromosome/Region: \$chrom" >> ${bam_name}_depth_stats.txt
            echo "Total length: \$length bp" >> ${bam_name}_depth_stats.txt
            echo "Total bases covered: \$bases" >> ${bam_name}_depth_stats.txt
            echo "Mean depth: \$mean_depth" >> ${bam_name}_depth_stats.txt
            
            # Calculate coverage percentage
            if [ "\$length" -gt 0 ] && [ "\$bases" -gt 0 ]; then
                coverage_pct=\$(echo "scale=2; \$bases * 100 / \$length" | bc -l)
                echo "Coverage percentage: \$coverage_pct%" >> ${bam_name}_depth_stats.txt
            fi
        else
            echo "Warning: Could not extract depth statistics" >> ${bam_name}_depth_stats.txt
        fi
        
        echo "" >> ${bam_name}_depth_stats.txt
        echo "=== RAW MOSDEPTH SUMMARY ===" >> ${bam_name}_depth_stats.txt
        cat ${depth_summary} >> ${bam_name}_depth_stats.txt
    else
        echo "Error: Depth summary file not found" >> ${bam_name}_depth_stats.txt
    fi
    """
}

process generate_depth_summary {
    tag "depth_summary"
    memory '8g'
    cpus 4
    time '2.h'
    publishDir "${params.output_dir}/summary", mode: 'copy'
    
    input:
    path(depth_stats_files)
    path(qc_files)
    
    output:
    path "depth_summary_report.txt"
    path "sample_depth_summary.csv"
    path "depth_qc_combined.csv"
    
    script:
    """
    echo "Depth Calculation and Quality Control Summary" > depth_summary_report.txt
    echo "=============================================" >> depth_summary_report.txt
    echo "Date: \$(date)" >> depth_summary_report.txt
    echo "Total samples processed: \$(ls *_depth_stats.txt | wc -l)" >> depth_summary_report.txt
    echo "" >> depth_summary_report.txt
    echo "Parameters used:" >> depth_summary_report.txt
    echo "- Use mosdepth: ${params.use_mosdepth}" >> depth_summary_report.txt
    echo "- Window size: ${params.window_size}" >> depth_summary_report.txt
    echo "- Min mapping quality: ${params.min_mapping_quality}" >> depth_summary_report.txt
    echo "- Threads: ${params.threads}" >> depth_summary_report.txt
    echo "" >> depth_summary_report.txt
    
    # Create CSV header for depth summary
    echo "Sample,Mean_Depth,Total_Length,Total_Bases,Coverage_Percentage" > sample_depth_summary.csv
    
    # Process each depth stats file
    for stats_file in *_depth_stats.txt; do
        if [ -f "\$stats_file" ]; then
            sample=\$(basename \$stats_file _depth_stats.txt)
            
            # Extract depth information
            mean_depth=\$(grep "Mean depth:" \$stats_file | awk '{print \$3}' | head -1)
            total_length=\$(grep "Total length:" \$stats_file | awk '{print \$3}' | head -1)
            total_bases=\$(grep "Total bases covered:" \$stats_file | awk '{print \$4}' | head -1)
            coverage_pct=\$(grep "Coverage percentage:" \$stats_file | awk '{print \$3}' | sed 's/%//' | head -1)
            
            # Handle missing values
            mean_depth=\${mean_depth:-"N/A"}
            total_length=\${total_length:-"N/A"}
            total_bases=\${total_bases:-"N/A"}
            coverage_pct=\${coverage_pct:-"N/A"}
            
            echo "\$sample,\$mean_depth,\$total_length,\$total_bases,\$coverage_pct" >> sample_depth_summary.csv
            
            # Add to summary report
            echo "Sample: \$sample" >> depth_summary_report.txt
            echo "  Mean depth: \$mean_depth" >> depth_summary_report.txt
            echo "  Coverage: \$coverage_pct%" >> depth_summary_report.txt
            echo "" >> depth_summary_report.txt
        fi
    done
    
    # Create combined depth and QC CSV if QC files are available
    echo "Sample,Mean_Depth,Coverage_Percentage,Total_Reads,Mapped_Reads,Mapping_Rate,Duplicates,Duplicate_Rate" > depth_qc_combined.csv
    
    # Process QC files if available
    qc_count=\$(ls *.flagstat 2>/dev/null | wc -l)
    if [ \$qc_count -gt 0 ]; then
        echo "=== QUALITY CONTROL SUMMARY ===" >> depth_summary_report.txt
        
        for flagstat_file in *.flagstat; do
            if [ -f "\$flagstat_file" ]; then
                sample=\$(basename \$flagstat_file .flagstat)
                
                # Extract QC metrics
                total_reads=\$(grep "total" \$flagstat_file | head -1 | awk '{print \$1}')
                mapped_reads=\$(grep "mapped (" \$flagstat_file | head -1 | awk '{print \$1}')
                duplicates=\$(grep "duplicates" \$flagstat_file | awk '{print \$1}')
                
                # Calculate rates
                if [ "\$total_reads" -gt 0 ]; then
                    mapping_rate=\$(echo "scale=2; \$mapped_reads * 100 / \$total_reads" | bc -l)
                    duplicate_rate=\$(echo "scale=2; \$duplicates * 100 / \$total_reads" | bc -l)
                else
                    mapping_rate="0"
                    duplicate_rate="0"
                fi
                
                # Get depth info for this sample
                depth_line=\$(grep "^\$sample," sample_depth_summary.csv)
                if [ -n "\$depth_line" ]; then
                    mean_depth=\$(echo "\$depth_line" | cut -d',' -f2)
                    coverage_pct=\$(echo "\$depth_line" | cut -d',' -f5)
                else
                    mean_depth="N/A"
                    coverage_pct="N/A"
                fi
                
                echo "\$sample,\$mean_depth,\$coverage_pct,\$total_reads,\$mapped_reads,\$mapping_rate,\$duplicates,\$duplicate_rate" >> depth_qc_combined.csv
                
                echo "Sample: \$sample" >> depth_summary_report.txt
                echo "  Total reads: \$total_reads" >> depth_summary_report.txt
                echo "  Mapped reads: \$mapped_reads (\$mapping_rate%)" >> depth_summary_report.txt
                echo "  Duplicates: \$duplicates (\$duplicate_rate%)" >> depth_summary_report.txt
                echo "" >> depth_summary_report.txt
            fi
        done
    else
        echo "No QC files found. QC analysis was not performed." >> depth_summary_report.txt
    fi
    
    # Generate summary statistics
    echo "=== OVERALL STATISTICS ===" >> depth_summary_report.txt
    if [ -f sample_depth_summary.csv ] && [ \$(wc -l < sample_depth_summary.csv) -gt 1 ]; then
        avg_depth=\$(tail -n +2 sample_depth_summary.csv | awk -F',' '{if(\$2!="N/A") {sum+=\$2; count++}} END {if(count>0) printf "%.2f", sum/count; else print "N/A"}')
        avg_coverage=\$(tail -n +2 sample_depth_summary.csv | awk -F',' '{if(\$5!="N/A") {sum+=\$5; count++}} END {if(count>0) printf "%.2f", sum/count; else print "N/A"}')
        echo "Average depth across samples: \$avg_depth" >> depth_summary_report.txt
        echo "Average coverage across samples: \$avg_coverage%" >> depth_summary_report.txt
    fi
    
    if [ -f depth_qc_combined.csv ] && [ \$(wc -l < depth_qc_combined.csv) -gt 1 ]; then
        avg_mapping=\$(tail -n +2 depth_qc_combined.csv | awk -F',' '{if(\$6!="N/A") {sum+=\$6; count++}} END {if(count>0) printf "%.2f", sum/count; else print "N/A"}')
        echo "Average mapping rate: \$avg_mapping%" >> depth_summary_report.txt
    fi
    
    echo "" >> depth_summary_report.txt
    echo "Report generated on: \$(date)" >> depth_summary_report.txt
    """
}

process generate_coverage_plots {
    tag "coverage_plots"
    memory '16g'
    cpus 4
    time '4.h'
    publishDir "${params.output_dir}/plots", mode: 'copy'
    
    input:
    path(depth_files)
    
    output:
    path "coverage_distribution.png"
    path "depth_comparison.png"
    path "coverage_summary.txt"
    
    script:
    """
    echo "Generating coverage plots..." > coverage_summary.txt
    echo "Date: \$(date)" >> coverage_summary.txt
    echo "" >> coverage_summary.txt
    
    # Count available depth files
    bed_count=\$(ls *.mosdepth.regions.bed.gz 2>/dev/null | wc -l)
    summary_count=\$(ls *.mosdepth.summary.txt 2>/dev/null | wc -l)
    
    echo "Mosdepth regions files: \$bed_count" >> coverage_summary.txt
    echo "Mosdepth summary files: \$summary_count" >> coverage_summary.txt
    echo "" >> coverage_summary.txt
    
    if [ \$summary_count -gt 0 ]; then
        # Create R script for plotting
        cat > plot_coverage.R << 'EOF'
library(ggplot2)
library(dplyr)
library(readr)

# Read all summary files
summary_files <- list.files(pattern = "\\\\.mosdepth\\\\.summary\\\\.txt$")
coverage_data <- data.frame()

for (file in summary_files) {
    sample_name <- gsub("\\\\.mosdepth\\\\.summary\\\\.txt$", "", file)
    data <- read_tsv(file, col_names = c("chrom", "length", "bases", "mean", "min", "max"))
    data\$sample <- sample_name
    coverage_data <- rbind(coverage_data, data)
}

# Coverage distribution plot
if (nrow(coverage_data) > 0) {
    p1 <- ggplot(coverage_data, aes(x = sample, y = mean, fill = sample)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Mean Coverage by Sample", 
             x = "Sample", 
             y = "Mean Depth") +
        guides(fill = FALSE)
    
    ggsave("coverage_distribution.png", p1, width = 12, height = 8, dpi = 300)
    
    # Depth comparison plot
    p2 <- ggplot(coverage_data, aes(x = reorder(sample, mean), y = mean)) +
        geom_point(size = 3, color = "steelblue") +
        geom_hline(yintercept = median(coverage_data\$mean, na.rm = TRUE), 
                   linetype = "dashed", color = "red") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Depth Comparison Across Samples",
             subtitle = paste("Red line: median depth =", 
                             round(median(coverage_data\$mean, na.rm = TRUE), 2)),
             x = "Sample (ordered by depth)", 
             y = "Mean Depth")
    
    ggsave("depth_comparison.png", p2, width = 12, height = 8, dpi = 300)
    
    cat("Coverage plots generated successfully\\n")
} else {
    cat("No coverage data found for plotting\\n")
}
EOF
        
        # Run R script
        if command -v Rscript &> /dev/null; then
            Rscript plot_coverage.R >> coverage_summary.txt 2>&1
            echo "R plotting completed" >> coverage_summary.txt
        else
            echo "R not available. Creating placeholder plots..." >> coverage_summary.txt
            # Create placeholder images
            convert -size 800x600 xc:white -pointsize 20 -draw "text 50,300 'R not available for plotting'" coverage_distribution.png 2>/dev/null || touch coverage_distribution.png
            convert -size 800x600 xc:white -pointsize 20 -draw "text 50,300 'R not available for plotting'" depth_comparison.png 2>/dev/null || touch depth_comparison.png
        fi
    else
        echo "No summary files found for plotting" >> coverage_summary.txt
        touch coverage_distribution.png
        touch depth_comparison.png
    fi
    
    echo "Plot generation completed" >> coverage_summary.txt
    """
}
