#!/usr/bin/env nextflow

/*
 * FastCall3 Disc Performance Analysis Pipeline
 * Author: Performance tuning for disc step
 * Date: 2025-09-18
 * Purpose: Test different parallel configurations to find optimal CPU/Memory utilization
 */

nextflow.enable.dsl=2

// Performance testing parameters
params.home_dir = null
params.pop = "test"
params.job = "perf_test"
params.tiger_jar = "TIGER_F3_20250915.jar"
params.java_lib = "/data/dazheng/lib/jvm"
params.samtools_path = "samtools"
params.output_dir = "perf_output"
params.help = false

// Test configuration parameters
params.test_chromosomes = ["1", "2"]  // Use multiple chromosomes for testing
params.parallel_configs = [
    [parallel: 2, cpus: 16, memory: "64g"],
    [parallel: 4, cpus: 8, memory: "32g"],
    [parallel: 6, cpus: 8, memory: "24g"],
    [parallel: 8, cpus: 6, memory: "20g"],
    [parallel: 10, cpus: 6, memory: "16g"],
    [parallel: 12, cpus: 4, memory: "16g"],
    [parallel: 16, cpus: 4, memory: "12g"],
    [parallel: 20, cpus: 4, memory: "10g"]
]

// Monitor interval in seconds
params.monitor_interval = 30
params.test_duration_limit = 7200  // 2 hours max per test

// Import functions from runFastCall3.nf
def getPopulationConfig(pop, home_dir) {
    def popConfigs = [
        "test": [
            bam_dir: "${home_dir}/01testData/02bam/01test",
            depth_dir: "${home_dir}/01testData/04depth/01test",
            reference: "${home_dir}/01testData/03ref/abd_1M.fa",
            description: "Test dataset"
        ],
        "A": [
            bam_dir: "${home_dir}/00data/02bam/bam1/A",
            depth_dir: "${home_dir}/00data/04depth/01A",
            reference: "${home_dir}/00data/03ref/01A/a_iwgscV1.fa.gz",
            description: "A genome population"
        ]
    ]
    
    if (!popConfigs.containsKey(pop)) {
        def validPops = popConfigs.keySet().join(", ")
        throw new Exception("Unknown population: ${pop}. Valid options: ${validPops}")
    }
    
    return popConfigs[pop]
}

def getJavaSetupScript(javaVersion, javaLibDir) {
    def javaVersionMap = [
        "java8": "jdk-8", 
        "java17": "jdk-17",
        "java21": "jdk-21"
    ]
    
    def javaDir = javaVersionMap[javaVersion]
    def javaHome = "${javaLibDir}/${javaDir}"
    
    return """
    export JAVA_HOME=${javaHome}
    export PATH=\$JAVA_HOME/bin:\$PATH
    java -version 2>&1 | head -3
    """.stripIndent()
}

def getTigerJarConfig(tigerJarPath) {
    def jarFile = file(tigerJarPath)
    def jarName = jarFile.name
    
    def config = [
        java_version: "java17",
        fastcall_version: "FastCall3",
        app_name: "FastCall3"
    ]
    
    return [
        path: tigerJarPath,
        name: jarName,
        java_version: config.java_version,
        fastcall_version: config.fastcall_version,
        app_name: config.app_name
    ]
}

def resolveTigerJarPath(tigerJar, homeDir) {
    def candidatePaths = [
        "${homeDir}/lib/${tigerJar}",
        "${homeDir}/software/${tigerJar}",
        "${homeDir}/tools/${tigerJar}",
        "${homeDir}/${tigerJar}",
        "./${tigerJar}"
    ]
    
    for (path in candidatePaths) {
        if (file(path).exists()) {
            log.info "Found TIGER jar at: ${path}"
            return path
        }
    }
    
    throw new Exception("TIGER jar file not found: ${tigerJar}")
}

// Help message
def helpMessage() {
    log.info """
    ========================================
    FastCall3 Disc Performance Analysis
    ========================================
    
    Usage:
        nextflow run perfDiscCpuMem.nf --home_dir <home> [options]
    
    Required parameters:
        --home_dir          Home directory containing data and libraries
    
    Optional parameters:
        --pop               Population name (default: test)
        --job               Job name (default: perf_test)
        --tiger_jar         TIGER jar filename (default: TIGER_F3_20250915.jar)
        --java_lib          Java installation directory (default: /data/dazheng/lib/jvm)
        --output_dir        Output directory (default: perf_output)
        --test_chromosomes  Chromosomes to test (default: ["1A", "1B", "1D"])
        --monitor_interval  Resource monitoring interval in seconds (default: 30)
        
    Performance test configurations:
        Each test runs with different parallel/cpus/memory combinations:
        - parallel=2,  cpus=16, memory=64g
        - parallel=4,  cpus=8,  memory=32g
        - parallel=6,  cpus=8,  memory=24g
        - parallel=8,  cpus=6,  memory=20g
        - parallel=10, cpus=6,  memory=16g
        - parallel=12, cpus=4,  memory=16g
        - parallel=16, cpus=4,  memory=12g
        - parallel=20, cpus=4,  memory=10g
        
    Output:
        - Individual performance logs for each configuration
        - Resource utilization charts
        - Summary report with recommendations
        
    Examples:
        # Basic test with test dataset
        nextflow run perfDiscCpuMem.nf --home_dir /data/project
        
        # Test with A population
        nextflow run perfDiscCpuMem.nf --home_dir /data/project --pop A
        
        # Custom test chromosomes
        nextflow run perfDiscCpuMem.nf --home_dir /data/project \\
            --test_chromosomes '["1", "2"]'
    """.stripIndent()
}

// Main workflow
workflow {
    perfDiscCpuMem_workflow()
}

workflow perfDiscCpuMem_workflow {
    // Show help if requested
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    // Parameter validation
    if (!params.home_dir) {
        log.error "Home directory is required. Use --home_dir"
        exit 1
    }
    
    // Setup paths
    def popConfig = getPopulationConfig(params.pop, params.home_dir)
    def bam_dir = popConfig.bam_dir
    def depth_dir = popConfig.depth_dir
    def reference = popConfig.reference
    def tiger_jar_path = resolveTigerJarPath(params.tiger_jar, params.home_dir)
    def tiger_jar_config = getTigerJarConfig(tiger_jar_path)
    def workshop_jar = "${params.home_dir}/lib/Workshop.jar"
    
    log.info """
    ========================================
    FastCall3 Disc Performance Analysis
    ========================================
    Home directory   : ${params.home_dir}
    Population       : ${params.pop}
    BAM directory    : ${bam_dir}
    Reference genome : ${reference}
    TIGER jar        : ${tiger_jar_config.path}
    Test chromosomes : ${params.test_chromosomes.join(', ')}
    Output directory : ${params.output_dir}
    ========================================
    """.stripIndent()
    
    // Step 1: Prepare taxa-BAM mapping if needed
    if (!file("${params.home_dir}/00data/05taxaBamMap/${params.job}.taxaBamMap.txt").exists()) {
        prep_results = prepareTaxaBamMap(
            file(workshop_jar),
            depth_dir,
            bam_dir,
            params.home_dir,
            params.job
        )
        taxa_bam_map = prep_results.taxa_bam_map
    } else {
        taxa_bam_map = Channel.fromPath("${params.home_dir}/00data/05taxaBamMap/${params.job}.taxaBamMap.txt")
        log.info "Using existing taxa-BAM mapping file"
    }
    
    // Step 2: Create test configuration channel
    config_ch = Channel.fromList(params.parallel_configs)
        .map { config ->
            [config.parallel, config.cpus, config.memory, config]
        }
    
    // Step 3: Run performance tests
    perf_results = runPerfTest(
        config_ch,
        Channel.fromList(params.test_chromosomes),
        file(reference),
        taxa_bam_map,
        file(tiger_jar_config.path),
        params.samtools_path,
        tiger_jar_config,
        tiger_jar_path,
        reference
    )
    
    // Step 4: Collect and analyze results
    final_report = analyzeResults(
        perf_results.logs.collect(),
        perf_results.resource_logs.collect()
    )
}

process prepareTaxaBamMap {
    tag "prepare_taxa_bam_map"
    memory '16g'
    cpus 4
    publishDir "${params.home_dir}/00data/05taxaBamMap", mode: 'copy'
    
    input:
    path workshop_jar
    val depth_dir
    val bam_dir
    val home_dir
    val job
    
    output:
    path "${job}.taxaBamMap.txt", emit: taxa_bam_map
    path "${job}.taxaRunMap.txt", emit: taxa_run_map, optional: true
    
    script:
    def javaSetup = getJavaSetupScript("java21", params.java_lib)
    """
    echo "Preparing taxa-BAM mapping file..."
    ${javaSetup}
    
    java -Xmx16g -jar "${workshop_jar}" \\
        -d ${depth_dir} \\
        -b ${bam_dir} \\
        -o ${job}.taxaBamMap.txt \\
        -t ${job}.taxaRunMap.txt
    """
}

process runPerfTest {
    tag "perf_test_p${parallel}_c${cpus}_m${memory}"
    memory '8g'
    cpus 2
    publishDir "${params.output_dir}/perf_tests", mode: 'copy'
    time '4h'
    
    input:
    tuple val(parallel), val(cpus), val(memory), val(config)
    each chromosome
    path reference
    path taxa_bam_map
    path tiger_jar
    val samtools_path
    val tiger_jar_config
    val tiger_jar_path
    val reference_path
    
    output:
    path "perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log", emit: logs
    path "resource_monitor_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log", emit: resource_logs
    
    script:
    def javaSetup = getJavaSetupScript(tiger_jar_config.java_version, params.java_lib)
    """
    echo "Starting performance test: parallel=${parallel}, cpus=${cpus}, memory=${memory}, chromosome=${chromosome}" > perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
    echo "Test started at: \$(date)" >> perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
    
    ${javaSetup} >> perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log 2>&1
    
    # Start resource monitoring in background
    (
        echo "timestamp,cpu_percent,memory_used_gb,memory_available_gb,load_avg" > resource_monitor_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
        while true; do
            timestamp=\$(date +"%Y-%m-%d %H:%M:%S")
            cpu_percent=\$(top -bn1 | grep "Cpu(s)" | awk '{print \$2}' | cut -d'%' -f1)
            memory_info=\$(free -g | grep "Mem:")
            memory_used=\$(echo \$memory_info | awk '{print \$3}')
            memory_available=\$(echo \$memory_info | awk '{print \$7}')
            load_avg=\$(uptime | awk -F'load average:' '{print \$2}' | cut -d',' -f1 | xargs)
            echo "\$timestamp,\$cpu_percent,\$memory_used,\$memory_available,\$load_avg" >> resource_monitor_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
            sleep ${params.monitor_interval}
        done
    ) &
    monitor_pid=\$!
    
    # Simulate parallel disc processes
    start_time=\$(date +%s)
    pids=()
    
    for i in \$(seq 1 ${parallel}); do
        (
            echo "Starting disc process \$i for chromosome ${chromosome}" >> perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
            
            # Create separate work directory for each process
            mkdir -p work_\$i
            cd work_\$i
            
            # Create symbolic links to required files
            ln -s ../perf_test.taxaBamMap.txt ./
            ln -s ../a_iwgscV1.fa.gz ./
            
            /usr/bin/time -v java -Xmx${memory} -jar ${tiger_jar_path} \\
                -app ${tiger_jar_config.app_name} \\
                -mod disc \\
                -a ${reference_path} \\
                -b ${taxa_bam_map} \\
                -c 0 \\
                -d 30 \\
                -e 20 \\
                -f 2 \\
                -g 0.2 \\
                -h 3 \\
                -i 0.8 \\
                -j 0.35 \\
                -k 0.2 \\
                -l ${chromosome} \\
                -m ${cpus} \\
                -n ./ \\
                -o ${samtools_path} \\
                > disc_\$i.log 2>&1
            
            echo "Disc process \$i completed at: \$(date)" >> ../perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
        ) &
        pids+=(\$!)
    done
    
    # Wait for all processes to complete or timeout
    timeout ${params.test_duration_limit} bash -c "
        for pid in \${pids[@]}; do
            wait \$pid
        done
    "
    
    # Stop resource monitoring
    kill \$monitor_pid 2>/dev/null || true
    
    end_time=\$(date +%s)
    duration=\$((end_time - start_time))
    
    echo "Test completed at: \$(date)" >> perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
    echo "Total duration: \$duration seconds" >> perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
    
    # Collect resource usage statistics
    echo "=== Resource Usage Summary ===" >> perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
    echo "Configuration: parallel=${parallel}, cpus=${cpus}, memory=${memory}" >> perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
    echo "Chromosome: ${chromosome}" >> perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
    echo "Duration: \$duration seconds" >> perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
    
    # Extract peak memory usage from time output
    for i in \$(seq 1 ${parallel}); do
        if [ -f work_\$i/disc_\$i.log ]; then
            echo "Process \$i resource usage:" >> perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
            grep -E "(Maximum resident set size|User time|System time|Percent of CPU)" work_\$i/disc_\$i.log >> perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log || true
        fi
    done
    
    # Calculate average resource utilization
    if [ -f resource_monitor_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log ]; then
        echo "=== Average Resource Utilization ===" >> perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
        awk -F',' 'NR>1 {cpu_sum+=\$2; mem_sum+=\$3; count++} END {
            if(count>0) {
                printf "Average CPU: %.2f%%\\n", cpu_sum/count;
                printf "Average Memory Used: %.2f GB\\n", mem_sum/count;
                printf "Samples: %d\\n", count;
            }
        }' resource_monitor_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log >> perf_test_p${parallel}_c${cpus}_m${memory}_chr${chromosome}.log
    fi
    """
}

process analyzeResults {
    tag "analyze_results"
    memory '8g'
    cpus 4
    publishDir "${params.output_dir}/final_report", mode: 'copy'
    
    input:
    path(perf_logs)
    path(resource_logs)
    
    output:
    path "performance_analysis_report.txt"
    path "resource_utilization_summary.csv"
    path "recommendations.txt"
    
    script:
    """
    echo "FastCall3 Disc Performance Analysis Report" > performance_analysis_report.txt
    echo "==========================================" >> performance_analysis_report.txt
    echo "Analysis Date: \$(date)" >> performance_analysis_report.txt
    echo "Test Configurations: ${params.parallel_configs.size()}" >> performance_analysis_report.txt
    echo "Test Chromosomes: ${params.test_chromosomes.join(', ')}" >> performance_analysis_report.txt
    echo "" >> performance_analysis_report.txt
    
    # Create CSV header
    echo "config,parallel,cpus,memory,chromosome,duration_sec,avg_cpu_percent,avg_memory_gb,peak_memory_kb,efficiency_score" > resource_utilization_summary.csv
    
    # Process each performance log
    for log_file in ${perf_logs}; do
        echo "Processing: \$log_file" >> performance_analysis_report.txt
        
        # Extract configuration from filename
        config=\$(basename \$log_file .log | sed 's/perf_test_//')
        parallel=\$(echo \$config | grep -o 'p[0-9]*' | sed 's/p//')
        cpus=\$(echo \$config | grep -o 'c[0-9]*' | sed 's/c//')
        memory=\$(echo \$config | grep -o 'm[0-9]*g' | sed 's/m//;s/g//')
        chromosome=\$(echo \$config | grep -o 'chr[^_]*' | sed 's/chr//')
        
        # Extract duration
        duration=\$(grep "Total duration:" \$log_file | awk '{print \$3}' || echo "0")
        
        # Extract average CPU and memory from the log
        avg_cpu=\$(grep "Average CPU:" \$log_file | awk '{print \$3}' | sed 's/%//' || echo "0")
        avg_memory=\$(grep "Average Memory Used:" \$log_file | awk '{print \$4}' || echo "0")
        
        # Extract peak memory usage
        peak_memory=\$(grep "Maximum resident set size" \$log_file | awk '{print \$6}' | sort -n | tail -1 || echo "0")
        
        # Calculate efficiency score (CPU * Memory utilization / duration)
        efficiency=\$(awk "BEGIN {if(\$duration>0) print (\$avg_cpu * \$avg_memory)/\$duration; else print 0}")
        
        # Add to CSV
        echo "\$config,\$parallel,\$cpus,\$memory,\$chromosome,\$duration,\$avg_cpu,\$avg_memory,\$peak_memory,\$efficiency" >> resource_utilization_summary.csv
        
        # Add detailed analysis to report
        echo "Configuration: parallel=\$parallel, cpus=\$cpus, memory=\${memory}g, chromosome=\$chromosome" >> performance_analysis_report.txt
        echo "  Duration: \$duration seconds" >> performance_analysis_report.txt
        echo "  Average CPU: \$avg_cpu%" >> performance_analysis_report.txt
        echo "  Average Memory: \$avg_memory GB" >> performance_analysis_report.txt
        echo "  Peak Memory: \$peak_memory KB" >> performance_analysis_report.txt
        echo "  Efficiency Score: \$efficiency" >> performance_analysis_report.txt
        echo "" >> performance_analysis_report.txt
    done
    
    # Generate recommendations
    echo "FastCall3 Disc Performance Recommendations" > recommendations.txt
    echo "==========================================" >> recommendations.txt
    echo "Generated: \$(date)" >> recommendations.txt
    echo "" >> recommendations.txt
    
    # Find best configurations
    echo "=== Best Configurations Analysis ===" >> recommendations.txt
    echo "" >> recommendations.txt
    
    # Best efficiency
    best_efficiency=\$(awk -F',' 'NR>1 {print \$10, \$0}' resource_utilization_summary.csv | sort -nr | head -1)
    echo "Highest Efficiency Configuration:" >> recommendations.txt
    echo "\$best_efficiency" | awk '{print \$2}' | awk -F',' '{printf "  Parallel: %s, CPUs: %s, Memory: %sg, Efficiency: %s\\n", \$2, \$3, \$4, \$10}' >> recommendations.txt
    echo "" >> recommendations.txt
    
    # Fastest completion
    fastest=\$(awk -F',' 'NR>1 {print \$6, \$0}' resource_utilization_summary.csv | sort -n | head -1)
    echo "Fastest Completion Configuration:" >> recommendations.txt
    echo "\$fastest" | awk '{print \$2}' | awk -F',' '{printf "  Parallel: %s, CPUs: %s, Memory: %sg, Duration: %s seconds\\n", \$2, \$3, \$4, \$6}' >> recommendations.txt
    echo "" >> recommendations.txt
    
    # Best CPU utilization
    best_cpu=\$(awk -F',' 'NR>1 {print \$7, \$0}' resource_utilization_summary.csv | sort -nr | head -1)
    echo "Best CPU Utilization Configuration:" >> recommendations.txt
    echo "\$best_cpu" | awk '{print \$2}' | awk -F',' '{printf "  Parallel: %s, CPUs: %s, Memory: %sg, CPU: %s%%\\n", \$2, \$3, \$4, \$7}' >> recommendations.txt
    echo "" >> recommendations.txt
    
    # Memory efficiency
    best_memory=\$(awk -F',' 'NR>1 && \$8>0 {print \$8/\$4, \$0}' resource_utilization_summary.csv | sort -nr | head -1)
    echo "Best Memory Efficiency Configuration:" >> recommendations.txt
    echo "\$best_memory" | awk '{print \$2}' | awk -F',' '{printf "  Parallel: %s, CPUs: %s, Memory: %sg, Memory Usage: %s GB\\n", \$2, \$3, \$4, \$8}' >> recommendations.txt
    echo "" >> recommendations.txt
    
    echo "=== General Recommendations ===" >> recommendations.txt
    echo "" >> recommendations.txt
    echo "1. Review the efficiency scores to find the optimal balance" >> recommendations.txt
    echo "2. Consider your system's total resources when choosing parallel levels" >> recommendations.txt
    echo "3. Monitor actual resource usage during production runs" >> recommendations.txt
    echo "4. Adjust based on your specific data characteristics" >> recommendations.txt
    echo "" >> recommendations.txt
    echo "Note: Higher efficiency scores indicate better resource utilization" >> recommendations.txt
    echo "per unit time. Consider both efficiency and absolute performance." >> recommendations.txt
    
    echo "Analysis completed successfully" >> performance_analysis_report.txt
    """
}
