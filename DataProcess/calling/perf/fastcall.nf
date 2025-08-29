#!/usr/bin/env nextflow

/*
 * FastCall2 Performance Analysis Pipeline
 * Based on methods from 00seq_data_process.ipynb
 * Includes performance monitoring for disc, blib, and scan steps
 * Author: Generated from Vmap4 workflow
 * Date: 2025-08-11
 */

nextflow.enable.dsl=2

// Parameters
params.reference = null
params.taxaBamMap = null
params.tiger_jar = null
params.samtools_path = null
params.output_dir = "performance_output"
params.threads = 32
params.memory = "100g"
params.help = false

// Performance monitoring parameters
params.monitor_cpu = true
params.monitor_memory = true
params.monitor_io = true
params.benchmark_mode = false
params.test_data_size = "1M"  // For test data generation

// FastCall2 disc parameters
params.disc_min_depth = 30
params.disc_min_qual = 20
params.disc_min_snp_count = 2
params.disc_min_allele_freq = 0.2
params.disc_min_coverage = 3
params.disc_max_missing = 0.8
params.disc_min_het_freq = 0.35
params.disc_max_het_freq = 0.2

// Test chromosome for performance analysis
params.test_chromosome = "1"

def helpMessage() {
    log.info """
    ========================================
    FastCall2 Performance Analysis Pipeline
    ========================================
    
    Usage:
        nextflow run fastcall.nf --reference <ref.fa> --taxaBamMap <map.txt> --tiger_jar <TIGER.jar> --samtools_path <samtools>
    
    Required parameters:
        --reference         Reference genome fasta file
        --taxaBamMap        Taxa-BAM mapping file
        --tiger_jar         Path to TIGER jar file
        --samtools_path     Path to samtools executable
    
    Optional parameters:
        --output_dir        Output directory (default: performance_output)
        --threads           Number of threads (default: 32)
        --memory            Memory allocation (default: 100g)
        --test_chromosome   Chromosome for testing (default: 1)
        --benchmark_mode    Run in benchmark mode (default: false)
        --test_data_size    Size of test data (1M, 10M, 100M) (default: 1M)
    
    Performance monitoring:
        --monitor_cpu       Monitor CPU usage (default: true)
        --monitor_memory    Monitor memory usage (default: true)
        --monitor_io        Monitor I/O usage (default: true)
    """.stripIndent()
}

workflow performance_workflow {
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
    FastCall2 Performance Analysis Pipeline
    ========================================
    Reference genome  : ${params.reference}
    Taxa-BAM map      : ${params.taxaBamMap}
    TIGER jar         : ${params.tiger_jar}
    Samtools path     : ${params.samtools_path}
    Output directory  : ${params.output_dir}
    Test chromosome   : ${params.test_chromosome}
    Benchmark mode    : ${params.benchmark_mode}
    Test data size    : ${params.test_data_size}
    ========================================
    """.stripIndent()

    // Create test data if in benchmark mode
    if (params.benchmark_mode) {
        test_data = create_test_data()
        test_reference = test_data[0]
        test_taxa_map = test_data[1]
    } else {
        test_reference = Channel.value(params.reference)
        test_taxa_map = Channel.value(params.taxaBamMap)
    }

    // Performance analysis workflow
    perf_disc = performance_analysis_disc(test_reference, test_taxa_map)
    perf_blib = performance_analysis_blib(perf_disc[0], test_reference)
    perf_scan = performance_analysis_scan(perf_blib[0], test_reference, test_taxa_map)
    
    // Collect and analyze all performance data
    performance_report = generate_performance_report(
        perf_disc[1].collect(),
        perf_blib[1].collect(), 
        perf_scan[1].collect()
    )
}

process create_test_data {
    tag "create_test_data"
    memory '16g'
    cpus 4
    time '2.h'
    publishDir "${params.output_dir}/test_data", mode: 'copy'
    
    output:
    tuple path("test_reference.fa"), path("test_taxaBamMap.txt"), emit: test_data
    
    script:
    def size_param = params.test_data_size == "1M" ? "1000000" : params.test_data_size == "10M" ? "10000000" : "100000000"
    """
    echo "Creating test reference genome..." > test_data_creation.log
    
    # Create test reference (subset of original)
    if [ "${params.test_data_size}" = "1M" ]; then
        # Extract 1Mb region
        seqkit subseq -r 1:1000000 ${params.reference} -o test_reference.fa
    elif [ "${params.test_data_size}" = "10M" ]; then
        # Extract 10Mb region  
        seqkit subseq -r 1:10000000 ${params.reference} -o test_reference.fa
    else
        # Extract 100Mb region
        seqkit subseq -r 1:100000000 ${params.reference} -o test_reference.fa
    fi
    
    # Create subset of taxa mapping (first 10 samples)
    head -n 10 ${params.taxaBamMap} > test_taxaBamMap.txt
    
    echo "Test data created successfully" >> test_data_creation.log
    echo "Reference size: \$(ls -lh test_reference.fa | awk '{print \$5}')" >> test_data_creation.log
    echo "Taxa samples: \$(wc -l test_taxaBamMap.txt | awk '{print \$1}')" >> test_data_creation.log
    """
}

process performance_analysis_disc {
    tag "perf_disc"
    memory params.memory
    cpus params.threads
    time '24.h'
    publishDir "${params.output_dir}/performance/disc", mode: 'copy'
    
    input:
    path reference
    path taxaBamMap
    
    output:
    path "*.ing", emit: disc_files
    tuple path("disc_performance.json"), path("disc_system_stats.txt"), emit: perf_data
    
    script:
    """
    # Start performance monitoring
    echo "Starting performance monitoring for disc step" > disc_system_stats.txt
    echo "Date: \$(date)" >> disc_system_stats.txt
    echo "Hostname: \$(hostname)" >> disc_system_stats.txt
    echo "CPU cores: \$(nproc)" >> disc_system_stats.txt
    echo "Memory: \$(free -h | grep Mem | awk '{print \$2}')" >> disc_system_stats.txt
    echo "Disk space: \$(df -h . | tail -1 | awk '{print \$4}')" >> disc_system_stats.txt
    echo "" >> disc_system_stats.txt
    
    # Record start time and system state
    start_time=\$(date +%s.%N)
    start_cpu=\$(cat /proc/stat | grep ^cpu | head -1 | awk '{print \$2+\$3+\$4+\$5+\$6+\$7+\$8}')
    start_mem=\$(free | grep Mem | awk '{print \$3}')
    
    # Start system monitoring in background
    (
        while [ -f disc_running.flag ]; do
            echo "\$(date),\$(cat /proc/loadavg | awk '{print \$1}'),\$(free | grep Mem | awk '{print \$3/\$2*100}'),\$(iostat -d 1 1 | tail -n +4 | awk '{sum+=\$4} END {print sum}')" >> disc_monitor.csv
            sleep 5
        done
    ) &
    monitor_pid=\$!
    
    # Create flag file
    touch disc_running.flag
    
    # Run FastCall2 disc with performance monitoring
    /usr/bin/time -v java -Xmx${params.memory} -jar ${params.tiger_jar} \\
        -app FastCall2 \\
        -mod disc \\
        -a \$reference \\
        -b \$taxaBamMap \\
        -c 0 \\
        -d ${params.disc_min_depth} \\
        -e ${params.disc_min_qual} \\
        -f ${params.disc_min_snp_count} \\
        -g ${params.disc_min_allele_freq} \\
        -h ${params.disc_min_coverage} \\
        -i ${params.disc_max_missing} \\
        -j ${params.disc_min_het_freq} \\
        -k ${params.disc_max_het_freq} \\
        -l ${params.test_chromosome} \\
        -m ${params.threads} \\
        -n ./ \\
        -o ${params.samtools_path} \\
        > disc_output.log 2> disc_time_output.txt
    
    # Stop monitoring
    rm -f disc_running.flag
    kill \$monitor_pid 2>/dev/null || true
    
    # Record end time and system state
    end_time=\$(date +%s.%N)
    end_cpu=\$(cat /proc/stat | grep ^cpu | head -1 | awk '{print \$2+\$3+\$4+\$5+\$6+\$7+\$8}')
    end_mem=\$(free | grep Mem | awk '{print \$3}')
    
    # Calculate metrics
    elapsed_time=\$(echo "\$end_time - \$start_time" | bc -l)
    
    # Extract memory usage from time output
    max_memory=\$(grep "Maximum resident set size" disc_time_output.txt | awk '{print \$NF}')
    cpu_time=\$(grep "User time" disc_time_output.txt | awk '{print \$NF}')
    sys_time=\$(grep "System time" disc_time_output.txt | awk '{print \$NF}')
    
    # Create JSON performance report
    cat > disc_performance.json << EOF
{
    "step": "disc",
    "timestamp": "\$(date -Iseconds)",
    "parameters": {
        "threads": ${params.threads},
        "memory": "${params.memory}",
        "chromosome": "${params.test_chromosome}",
        "min_depth": ${params.disc_min_depth},
        "min_qual": ${params.disc_min_qual}
    },
    "performance": {
        "elapsed_time_seconds": \$elapsed_time,
        "cpu_time_seconds": \$cpu_time,
        "system_time_seconds": \$sys_time,
        "max_memory_kb": \$max_memory,
        "threads_used": ${params.threads}
    },
    "input_files": {
        "reference": "\$reference",
        "taxa_bam_map": "\$taxaBamMap"
    },
    "output_files": \$(ls *.ing | wc -l)
}
EOF
    
    # Add performance summary to stats
    echo "=== DISC PERFORMANCE SUMMARY ===" >> disc_system_stats.txt
    echo "Elapsed time: \$elapsed_time seconds" >> disc_system_stats.txt
    echo "CPU time: \$cpu_time seconds" >> disc_system_stats.txt
    echo "System time: \$sys_time seconds" >> disc_system_stats.txt
    echo "Max memory: \$max_memory KB" >> disc_system_stats.txt
    echo "Output files: \$(ls *.ing | wc -l)" >> disc_system_stats.txt
    """
}

process performance_analysis_blib {
    tag "perf_blib"
    memory params.memory
    cpus params.threads
    time '24.h'
    publishDir "${params.output_dir}/performance/blib", mode: 'copy'
    
    input:
    path disc_files
    path reference
    
    output:
    path "*.lib.gz", emit: blib_files
    tuple path("blib_performance.json"), path("blib_system_stats.txt"), emit: perf_data
    
    script:
    """
    # Performance monitoring setup
    echo "Starting performance monitoring for blib step" > blib_system_stats.txt
    echo "Date: \$(date)" >> blib_system_stats.txt
    echo "Input files: \$(ls *.ing | wc -l)" >> blib_system_stats.txt
    echo "" >> blib_system_stats.txt
    
    start_time=\$(date +%s.%N)
    
    # Start monitoring
    touch blib_running.flag
    (
        while [ -f blib_running.flag ]; do
            echo "\$(date),\$(cat /proc/loadavg | awk '{print \$1}'),\$(free | grep Mem | awk '{print \$3/\$2*100}')" >> blib_monitor.csv
            sleep 5
        done
    ) &
    monitor_pid=\$!
    
    # Run FastCall2 blib
    /usr/bin/time -v java -Xmx${params.memory} -jar ${params.tiger_jar} \\
        -app FastCall2 \\
        -mod blib \\
        -a \$reference \\
        -b 1 \\
        -c 2 \\
        -d ${params.threads} \\
        -e ./ \\
        -f ./ \\
        > blib_output.log 2> blib_time_output.txt
    
    # Stop monitoring
    rm -f blib_running.flag
    kill \$monitor_pid 2>/dev/null || true
    
    end_time=\$(date +%s.%N)
    elapsed_time=\$(echo "\$end_time - \$start_time" | bc -l)
    
    # Extract performance metrics
    max_memory=\$(grep "Maximum resident set size" blib_time_output.txt | awk '{print \$NF}')
    cpu_time=\$(grep "User time" blib_time_output.txt | awk '{print \$NF}')
    
    # Create JSON report
    cat > blib_performance.json << EOF
{
    "step": "blib",
    "timestamp": "\$(date -Iseconds)",
    "parameters": {
        "threads": ${params.threads},
        "memory": "${params.memory}"
    },
    "performance": {
        "elapsed_time_seconds": \$elapsed_time,
        "cpu_time_seconds": \$cpu_time,
        "max_memory_kb": \$max_memory
    },
    "input_files": \$(ls *.ing | wc -l),
    "output_files": \$(ls *.lib.gz | wc -l)
}
EOF
    
    echo "=== BLIB PERFORMANCE SUMMARY ===" >> blib_system_stats.txt
    echo "Elapsed time: \$elapsed_time seconds" >> blib_system_stats.txt
    echo "CPU time: \$cpu_time seconds" >> blib_system_stats.txt
    echo "Max memory: \$max_memory KB" >> blib_system_stats.txt
    """
}

process performance_analysis_scan {
    tag "perf_scan"
    memory params.memory
    cpus params.threads
    time '24.h'
    publishDir "${params.output_dir}/performance/scan", mode: 'copy'
    
    input:
    path blib_files
    path reference
    path taxaBamMap
    
    output:
    path "*.vcf*", emit: vcf_files
    tuple path("scan_performance.json"), path("scan_system_stats.txt"), emit: perf_data
    
    script:
    def lib_file = "*.lib.gz"
    """
    echo "Starting performance monitoring for scan step" > scan_system_stats.txt
    echo "Date: \$(date)" >> scan_system_stats.txt
    
    start_time=\$(date +%s.%N)
    
    # Find the library file
    lib_file=\$(ls *.lib.gz | head -1)
    
    # Start monitoring
    touch scan_running.flag
    (
        while [ -f scan_running.flag ]; do
            echo "\$(date),\$(cat /proc/loadavg | awk '{print \$1}'),\$(free | grep Mem | awk '{print \$3/\$2*100}')" >> scan_monitor.csv
            sleep 5
        done
    ) &
    monitor_pid=\$!
    
    # Run FastCall2 scan
    /usr/bin/time -v java -Xmx${params.memory} -jar ${params.tiger_jar} \\
        -app FastCall2 \\
        -mod scan \\
        -a \$reference \\
        -b \$taxaBamMap \\
        -c \$lib_file \\
        -d 1 \\
        -e 0 \\
        -f 30 \\
        -g 20 \\
        -h 0.05 \\
        -i ${params.samtools_path} \\
        -j ${params.threads} \\
        -k ./ \\
        > scan_output.log 2> scan_time_output.txt
    
    # Stop monitoring
    rm -f scan_running.flag
    kill \$monitor_pid 2>/dev/null || true
    
    end_time=\$(date +%s.%N)
    elapsed_time=\$(echo "\$end_time - \$start_time" | bc -l)
    
    # Extract performance metrics
    max_memory=\$(grep "Maximum resident set size" scan_time_output.txt | awk '{print \$NF}')
    cpu_time=\$(grep "User time" scan_time_output.txt | awk '{print \$NF}')
    
    # Create JSON report
    cat > scan_performance.json << EOF
{
    "step": "scan",
    "timestamp": "\$(date -Iseconds)",
    "parameters": {
        "threads": ${params.threads},
        "memory": "${params.memory}"
    },
    "performance": {
        "elapsed_time_seconds": \$elapsed_time,
        "cpu_time_seconds": \$cpu_time,
        "max_memory_kb": \$max_memory
    },
    "input_files": 1,
    "output_files": \$(ls *.vcf* | wc -l)
}
EOF
    
    echo "=== SCAN PERFORMANCE SUMMARY ===" >> scan_system_stats.txt
    echo "Elapsed time: \$elapsed_time seconds" >> scan_system_stats.txt
    echo "CPU time: \$cpu_time seconds" >> scan_system_stats.txt
    echo "Max memory: \$max_memory KB" >> scan_system_stats.txt
    """
}

process generate_performance_report {
    tag "performance_report"
    memory '8g'
    cpus 4
    time '2.h'
    publishDir "${params.output_dir}/reports", mode: 'copy'
    
    input:
    path disc_perf
    path blib_perf
    path scan_perf
    
    output:
    path "fastcall2_performance_report.html"
    path "performance_summary.json"
    path "performance_comparison.csv"
    
    script:
    """
    # Combine all performance data
    echo '{"pipeline_performance": [' > performance_summary.json
    
    # Add disc performance
    cat \$(ls *disc_performance.json | head -1) >> performance_summary.json
    echo ',' >> performance_summary.json
    
    # Add blib performance  
    cat \$(ls *blib_performance.json | head -1) >> performance_summary.json
    echo ',' >> performance_summary.json
    
    # Add scan performance
    cat \$(ls *scan_performance.json | head -1) >> performance_summary.json
    echo ']}' >> performance_summary.json
    
    # Create CSV comparison
    echo "Step,Elapsed_Time_Seconds,CPU_Time_Seconds,Max_Memory_KB,Threads" > performance_comparison.csv
    
    # Extract data from JSON files
    for json_file in *_performance.json; do
        step=\$(grep '"step"' \$json_file | cut -d'"' -f4)
        elapsed=\$(grep '"elapsed_time_seconds"' \$json_file | cut -d':' -f2 | tr -d ' ,')
        cpu_time=\$(grep '"cpu_time_seconds"' \$json_file | cut -d':' -f2 | tr -d ' ,')
        memory=\$(grep '"max_memory_kb"' \$json_file | cut -d':' -f2 | tr -d ' ,')
        threads=\$(grep '"threads_used"' \$json_file | cut -d':' -f2 | tr -d ' ,}' | head -1)
        
        echo "\$step,\$elapsed,\$cpu_time,\$memory,\$threads" >> performance_comparison.csv
    done
    
    # Generate HTML report
    cat > fastcall2_performance_report.html << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>FastCall2 Performance Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        .metric { background-color: #f9f9f9; padding: 10px; margin: 10px 0; }
        .summary { background-color: #e7f3ff; padding: 15px; margin: 20px 0; }
    </style>
</head>
<body>
    <h1>FastCall2 Performance Analysis Report</h1>
    
    <div class="summary">
        <h2>Executive Summary</h2>
        <p>This report contains performance analysis results for the FastCall2 variant calling pipeline.</p>
        <p>Generated on: EOF
    date >> fastcall2_performance_report.html
    cat >> fastcall2_performance_report.html << 'EOF'
</p>
    </div>
    
    <h2>Performance Metrics by Step</h2>
    <table>
        <tr>
            <th>Step</th>
            <th>Elapsed Time (s)</th>
            <th>CPU Time (s)</th>
            <th>Max Memory (KB)</th>
            <th>Efficiency (%)</th>
        </tr>
EOF
    
    # Add performance data to HTML table
    tail -n +2 performance_comparison.csv | while IFS=, read step elapsed cpu_time memory threads; do
        efficiency=\$(echo "scale=2; \$cpu_time / \$elapsed * 100" | bc -l 2>/dev/null || echo "N/A")
        echo "        <tr><td>\$step</td><td>\$elapsed</td><td>\$cpu_time</td><td>\$memory</td><td>\$efficiency</td></tr>" >> fastcall2_performance_report.html
    done
    
    cat >> fastcall2_performance_report.html << 'EOF'
    </table>
    
    <h2>Resource Utilization</h2>
    <div class="metric">
        <h3>Memory Usage</h3>
        <p>Peak memory usage across all steps and recommendations for optimization.</p>
    </div>
    
    <div class="metric">
        <h3>CPU Efficiency</h3>
        <p>CPU utilization efficiency indicates how well the parallel processing is working.</p>
    </div>
    
    <h2>Recommendations</h2>
    <ul>
        <li>Monitor memory usage patterns for optimal allocation</li>
        <li>Consider thread optimization based on CPU efficiency metrics</li>
        <li>Analyze I/O patterns for storage optimization</li>
    </ul>
    
    <footer>
        <p>Report generated by FastCall2 Performance Analysis Pipeline</p>
    </footer>
</body>
</html>
EOF
    """
}