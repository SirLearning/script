#!/usr/bin/env nextflow

/*
 * Advanced FastCall3 pipeline for variant calling with enhanced features
 * Author: SirLearning
 * Date: 2025-11-19
 */

nextflow.enable.dsl=2

// --- Required Parameters ---
// 1. directory
params.home_dir = null
params.user_dir = null
params.output_dir = null
// 2. work specific
params.job = null
params.mod = null   // includes "full", "disc_only", "blib_only", "scan_only", "from_disc", "from_blib"
params.server = null  // server name

// --- Optional parameters ---
// 1. directory
params.ing_dir = null
params.vlib_dir = null
params.gen_dir = null
// 3. software settings
params.tiger_jar = "TIGER_F3_20251121.jar"
params.workshop_jar = "Workshop.jar"
params.samtools = "samtools"
// 4. single chromosome selector
params.chr = null
// 5. generate taxa-BAM mapping file on the server name input
params.gen_tbm = null  // server name for taxaBamMap generation

// --- software specific parameters ---
// FastCall3 disc parameters
params.disc_min_MQ = 30
params.disc_min_BQ = 20
params.disc_min_DC = 2
params.disc_min_DR = 0.2
params.disc_max_DR = 3
params.disc_HoR = 0.8
params.disc_HeR = 0.35
params.disc_TDR = 0.2
// FastCall3 blib parameters
params.blib_MAO = 2
// FastCall3 scan parameters
params.scan_min_MQ = 30
params.scan_min_BQ = 20
params.scan_error_rate = 0.05

// --- help info ---
params.help = false

// --- 107 command line ---
// nextflow run /data/dazheng/git/script/DataProcess/calling/run/runFastCall3.nf --home_dir /data/dazheng/01projects/vmap4 --user_dir /data/dazheng --output_dir /data/dazheng/01projects/vmap4/04test.FC3/07a.new.wf --job test_chr12 --mod scan2 --server s107 --gen_tbm 1

// --- 243 command line ---
// screen -dmS all bash -c "cd /data/home/tusr1/01projects/vmap4/04runScreens/02run && source ~/.bashrc && conda activate run && nextflow run /data/home/tusr1/git/script/DataProcess/calling/run/runFastCall3.nf --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 --output_dir /data/home/tusr1/01projects/vmap4/04runScreens --job all --mod from_blib --server s243 --gen_tbm 1"
// nextflow run /data/home/tusr1/git/script/DataProcess/calling/run/runFastCall3.nf --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 --output_dir /data/home/tusr1/01projects/vmap4/04runScreens --job all --mod from_blib --server s243 --gen_tbm 1

// --- old command line cache ---
// nextflow run /data/dazheng/git/script/DataProcess/calling/run/runFastCall3.nf --home_dir /data/dazheng/01projects/vmap4 --java_lib /data/dazheng/lib/jvm --pop chr1 --job test_ABD --workflow_mode disc_only --tiger_jar TIGER_F3_20250915.jar
// nextflow run /data/dazheng/git/script/DataProcess/calling/run/runFastCall3.nf --home_dir /data/dazheng/01projects/vmap4 --java_lib /data/dazheng/lib/jvm --pop chr1 --job test_ABD --workflow_mode from_disc --tiger_jar TIGER_F3_20251013.jar --ing_dir /data/dazheng/01projects/vmap4/04testFastCall3/01chr1/ing --vlib_dir /data/dazheng/01projects/vmap4/04testFastCall3/01chr1/vLib --gen_dir /data/dazheng/01projects/vmap4/04testFastCall3/01chr1/gen
// nextflow run /data/dazheng/git/script/DataProcess/calling/run/runFastCall3.nf --home_dir /data/dazheng/01projects/vmap4 --java_lib /data/dazheng/lib/jvm --pop chr1 --job test_ABD --workflow_mode from_disc --tiger_jar TIGER_F3_20251016.jar --ing_dir /data/dazheng/01projects/vmap4/04testFastCall3/01chr1/ing --vlib_dir /data/dazheng/01projects/vmap4/04testFastCall3/01chr1/vLib --gen_dir /data/dazheng/01projects/vmap4/04testFastCall3/01chr1/gen
// screen -dmS run_A_disc bash -c "cd /data/home/tusr1/01projects/runScreens/01A/disc && source ~/.bashrc && conda activate run && nextflow run /data/home/tusr1/01projects/DataProcess/calling/run/runFastCall3.nf --home_dir /data/home/tusr1/01projects/vmap4 --java_lib /data/home/tusr1/lib/jvm --pop A --job run_A_disc --workflow_mode disc_only --tiger_jar TIGER_F3_20250915.jar"


workflow {
    runFastCall3_workflow()
}

workflow runFastCall3_workflow {
    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    // Parameter validation
    if (!params.home_dir) {
        log.error "Home directory is required. Use --home_dir"
        exit 1
    }
    if (!params.user_dir) {
        log.error "User directory is required. Use --user_dir"
        exit 1
    }
    if (!params.job) {
        log.error "Job name is required. Use --job"
        exit 1
    }
    if (!params.mod) {
        log.error "Workflow mode is required. Use --mod"
        exit 1
    }

    // Step 1: prepare input information

    // --- set default variables ---
    def wheat_chrs = (0..44).collect { it.toString() }
    def output_path = "${params.output_dir}/${params.job}"
    def taxaBamMap_dir = file("${params.home_dir}/00data/05taxaBamMap")

    // --- retrieve parameters ---
    // 1. directory
    def ing_dir_resolved  = (params.ing_dir  ?: "${output_path}/ing").toString().replaceAll(/\/+/, "/")
    def vlib_dir_resolved = (params.vlib_dir ?: "${output_path}/vLib").toString().replaceAll(/\/+/, "/")
    def gen_dir_resolved  = (params.gen_dir  ?: "${output_path}/gen").toString().replaceAll(/\/+/, "/")
    // 2. work specific config
    def job_config = getJobConfig(params.job, params.home_dir)
    // Resolve and configure TIGER jar automatically
    def software_config = getSoftwareConfig(params.home_dir, params.user_dir, params.tiger_jar, params.workshop_jar, params.samtools)
    // create software input tuple
    def tiger_jar_input = tuple(
        software_config.tiger_jar_config.path,
        software_config.tiger_jar_config.app_name,
        software_config.tiger_jar_config.java_version,
        software_config.tiger_jar_config.java_lib_dir
    )
    def workshop_jar_input = tuple(
        software_config.workshop_jar_config.path,
        software_config.workshop_jar_config.java_version,
        software_config.workshop_jar_config.java_lib_dir
    )
    def samtools_input = software_config.samtools_config.path
    // --- validate directories ---
    // ensure essential directories exist
    def pathValidation = validatePaths([
        "Home directory": params.home_dir,
        "TIGER jar": software_config.tiger_jar_config.path,
        "Workshop jar": software_config.workshop_jar_config.path
    ])
    if (!pathValidation.isValid) {
        log.error "Path validation failed:"
        pathValidation.errors.each { log.error "  - ${it}" }
        exit 1
    }
    // create output directories
    def output_dir = file(output_path)
    output_dir.mkdirs()
    file(ing_dir_resolved).mkdirs()
    file(vlib_dir_resolved).mkdirs()
    file(gen_dir_resolved).mkdirs()
    // create pipeline info directory
    def pipeline_info_dir = file("./pipeline_info/${params.job}")
    pipeline_info_dir.mkdirs()

    // --- set taxaBamMap.txt and reference genome for different populations ---
    def server_pop = getServerPopulations(params.server)
    log.info "Server populations for taxaBamMap generation: ${server_pop}"
    def tbm_files = []
    if (params.gen_tbm) {
        log.info "Mode: Generating new taxaBamMap files."
        def pop_ch = Channel.fromList(server_pop)
        def pop_input_ch = pop_ch.map { pop ->
            def pop_config = getPopulationConfig(pop, params.home_dir)
            return pop_config
        }
        prepared_tbm_ch = prepareTaxaBamMap(workshop_jar_input, pop_input_ch).taxa_bam_map
    } else {
        tbm_files = Channel.fromList(server_pop)
            .map { pop -> tuple(pop, file("${params.home_dir}/00data/05taxaBamMap/${pop}.taxaBamMap.txt")) }

        log.info "Mode: Using existing taxaBamMap files."
        def existing_files = server_pop
            .collect { pop -> file("${params.home_dir}/00data/05taxaBamMap/${pop}.taxaBamMap.txt") }
            .findAll { it.exists() }
        if (existing_files.size() != server_pop.size()) {
            log.warn "Missing some existing taxaBamMap files. Found ${existing_files.size()}/${server_pop.size()}."
        }
        prepared_tbm_ch = Channel.fromPath(existing_files)
    }
    def pop_tbm_files_list_ch = prepared_tbm_ch.collect()
        .map { tbm_file_list ->
            def all_files = tbm_file_list
            def tuple_list = []
            if (job_config.a_pop.any()) {
                def a_pop_files = tbm_file_list
                    .findAll { f -> job_config.a_pop.any { pop_name -> f.name == "${pop_name}.taxaBamMap.txt" } }
                tuple_list << ["A", a_pop_files]
            }
            if (job_config.b_pop.any()) {
                def b_pop_files = tbm_file_list
                    .findAll { f -> job_config.b_pop.any { pop_name -> f.name == "${pop_name}.taxaBamMap.txt" } }
                tuple_list << ["B", b_pop_files]
            }
            if (job_config.d_pop.any()) {
                def d_pop_files = tbm_file_list
                    .findAll { f -> job_config.d_pop.any { pop_name -> f.name == "${pop_name}.taxaBamMap.txt" } }
                tuple_list << ["D", d_pop_files]
            }
            tuple_list << ["ALL", all_files]
            return tuple_list
        }
        .flatMap { it }
        .map { item ->
            log.info "Processing sub-genome: ${item[0]}"
            log.info "  Included taxaBamMap files: ${item[1].collect { it.name }.toList()}"
            def sub_genome = item[0]
            def pop_tbm_files = item[1]
            return tuple(pop_tbm_files, job_config.name, sub_genome)
        }
    def sub_genome_tbm_ch = concatTaxaBamMap(pop_tbm_files_list_ch).tbm
        .collect()
        .map { flat_list -> 
            log.info "show tuples: ${flat_list}"
            flat_list.collate(2).collectEntries()
        }

    sub_genome_tbm_ch
        .view { sub_genome_tbm -> 
            log.info """\
            ========================================
            FastCall3 Advanced Pipeline
            ========================================
            Home directory          : ${params.home_dir}
            User directory          : ${params.user_dir}
            Job name                : ${params.job}
            TIGER jar               : ${software_config.tiger_jar_config.path}
            FastCall version        : ${software_config.tiger_jar_config.app_name}
            TIGER java version      : ${software_config.tiger_jar_config.java_version}
            Workshop jar            : ${software_config.workshop_jar_config.path}
            Workshop java version   : ${software_config.workshop_jar_config.java_version}
            Samtools path           : ${software_config.samtools_config.path}
            All taxaBamMap          : ${sub_genome_tbm["ALL"]}
            Output directory        : ${output_path}
            ING directory           : ${ing_dir_resolved}
            vLib directory          : ${vlib_dir_resolved}
            GEN directory           : ${gen_dir_resolved}
            ========================================
            """.stripIndent()
        }
    
    // Step 2: run workflow in chromosome channel

    // Determine chromosomes from fai file if default, otherwise use provided chromosomes.
    // If --chr is specified, it overrides --chromosomes and runs only that chromosome.
    def chromosomes
    if (job_config.chroms && job_config.chroms.size() > 0) {
        chromosomes = job_config.chroms
        log.info "Using chromosomes from job config: ${chromosomes}"
    } else {
        // chromosomes = ["0","43","44"]
        // chromosomes = ["43"]
        chromosomes = []
        if (job_config.a_pop && job_config.a_pop.size() > 0) {
            // chromosomes += ["1","2","7","8","13","14","19","20","25","26","31","32","37","38"]
            // chromosomes += ["1","2","7","13","14","19","20","25","26","31","38"]  // temporary remove some chr for finished job
            // chromosomes += ["7","13","14","19","25","26","38"]  // temporary remove some chr for finished job
            // chromosomes += ["7","13","38"]  // temporary remove some chr for finished job
            // chromosomes += ["38"]  // temporary remove some chr for finished job
            chromosomes += ["20"]  // temporary remove some chr for finished job
        }
        if (job_config.b_pop && job_config.b_pop.size() > 0) {
            // chromosomes += ["3","4","9","10","15","16","21","22","27","28","33","34","39","40"]
            // chromosomes += ["4","9","10","15","16","21","22","27","28","33","34","39","40"]
            // chromosomes += ["9","10","15","16","21","22","27","33","34","39","40"]
            // chromosomes += ["9","10","15","21","22","27","33","39","40"]
            // chromosomes += ["10","27","40"]
            chromosomes += ["3","27","33"]
        }
        if (job_config.d_pop && job_config.d_pop.size() > 0) {
            // chromosomes += ["5","6","11","12","17","18","23","24","29","30","35","36","41","42"]
            // chromosomes += ["6","11","12","17","18","23","24","29","30","35","36","41","42"]
            // chromosomes += ["6","11","12","17","18","23","24","29","35","36","41","42"]
            // chromosomes += ["6","11","12","17","23","24","29","35","36","41","42"]
            // chromosomes += ["11","17","23","24","29","35","36","41","42"]
            // chromosomes += ["29","41"]
            // chromosomes += ["29","41"]
            // chromosomes += ["5","17"]
            chromosomes += ["5"]
        }
    }
    if (params.chr) {
        chromosomes = params.chr.split("|").collect { it.trim() }
        log.info "Using single chromosome from --chr: ${chromosomes[0]}"
    }
    log.info "Pipeline setup completed. Processing ${chromosomes.size()} chromosomes: ${chromosomes.size() > 0 ? chromosomes[0] + '...' + chromosomes[-1] : 'none'}"
    // Create chromosome channel
    chr_config_ch = sub_genome_tbm_ch
        .map { sub_genome_tbm ->
            def chr_config_list = []
            for (chromosome in chromosomes) {
                chr_config_list << getChrConfig(chromosome, params.home_dir, sub_genome_tbm)
            }
            return chr_config_list
        }
        .flatMap { it }

    log.info "Running workflow in mode: ${params.mod}"

    switch (params.mod) {
        case "full":
            // Complete workflow: disc -> blib -> scan -> collect
            disc_results = fastcall3_disc(
                chr_config_ch, 
                tiger_jar_input,
                samtools_input,
                ing_dir_resolved
            )
            blib_results = fastcall3_blib(
                disc_results.chr_config,
                tiger_jar_input,
                ing_dir_resolved,
                vlib_dir_resolved
            )
            scan_results = fastcall3_scan(
                blib_results.chr_config_lib,
                tiger_jar_input,
                samtools_input,
                gen_dir_resolved
            )
            break
            
        case "disc_only":
            // Only run disc analysis
            disc_results = fastcall3_disc(
                chr_config_ch, 
                tiger_jar_input,
                samtools_input,
                ing_dir_resolved
            )
            log.info "Disc analysis completed. Use 'from_disc' mode to continue."
            break
            
        case "from_disc":
            blib_results = fastcall3_blib(
                chr_config_ch,
                tiger_jar_input,
                ing_dir_resolved,
                vlib_dir_resolved
            )

            scan_results = fastcall3_scan(
                blib_results.chr_config_lib,
                tiger_jar_input,
                samtools_input,
                gen_dir_resolved
            )

            break
            
        case "blib_only":
            // Only run blib generation from existing disc (ing) results
            blib_results = fastcall3_blib(
                chr_config_ch,
                tiger_jar_input,
                ing_dir_resolved,
                vlib_dir_resolved
            )
            log.info "Blib generation completed. Use 'from_blib' mode to continue."
            break
            
        case "from_blib":
        case "scan_only":
            // Start from existing blib (vLib) results (only scan)
            def blib_dir = vlib_dir_resolved
            if (!file(blib_dir).exists()) {
                log.error "Blib directory not found: ${blib_dir}"
                log.error "Run workflow with 'blib_only', 'from_disc', or 'full' mode first"
                exit 1
            }
            def blib_chr_config_lib_ch = chr_config_ch
                .map { chromo, ref, fai, gzi, tbm ->
                    def lib_file = file(blib_dir).listFiles()
                        .find { f ->
                            if (!f.name.endsWith(".lib.gz")) return false
                            def fname = f.name
                            def m = (fname =~ /^([^_]+)_.*\.lib\.gz$/)
                            def lib_chr = m.find() ? m.group(1) : fname.tokenize('.')[0]
                            return lib_chr == chromo
                        }
                    if (lib_file) {
                        return tuple(chromo, ref, fai, gzi, tbm, lib_file)
                    } else {
                        log.warn "No .lib.gz file found for chromosome: ${chromo} in ${blib_dir}"
                        return null
                    }
                }
                .filter { it != null }
            if (!blib_chr_config_lib_ch) {
                log.error "No .lib.gz file found for chromosome: ${chr_config_ch.chr} in ${blib_dir}"
                exit 1
            }
            scan_results = fastcall3_scan(
                blib_chr_config_lib_ch,
                tiger_jar_input,
                samtools_input,
                gen_dir_resolved
            )
            break
        
        case "scan2":
            // Start from existing blib (vLib) results (only scan)
            def blib_dir = vlib_dir_resolved
            if (!file(blib_dir).exists()) {
                log.error "Blib directory not found: ${blib_dir}"
                log.error "Run workflow with 'blib_only', 'from_disc', or 'full' mode first"
                exit 1
            }
            def blib_chr_config_lib_ch = chr_config_ch
                .map { chromo, ref, fai, gzi, tbm ->
                    def lib_file = file(blib_dir).listFiles()
                        .find { f ->
                            if (!f.name.endsWith(".lib.gz")) return false
                            def fname = f.name
                            def m = (fname =~ /^([^_]+)_.*\.lib\.gz$/)
                            def lib_chr = m.find() ? m.group(1) : fname.tokenize('.')[0]
                            return lib_chr == chromo
                        }
                    if (lib_file) {
                        return tuple(chromo, ref, fai, gzi, tbm, lib_file)
                    } else {
                        log.warn "No .lib.gz file found for chromosome: ${chromo} in ${blib_dir}"
                        return null
                    }
                }
                .filter { it != null }
            if (!blib_chr_config_lib_ch) {
                log.error "No .lib.gz file found for chromosome: ${chr_config_ch.chr} in ${blib_dir}"
                exit 1
            }
            scan_results = fastcall3_scan2(
                blib_chr_config_lib_ch,
                tiger_jar_input,
                samtools_input,
                gen_dir_resolved
            )
            break

        default:
            log.error "Unknown workflow mode: ${workflow_mode}"
            log.error "Available modes: full, disc_only, blib_only, scan_only, from_disc, from_blib"
            exit 1
    }
    
}

process check_existing_files {
    tag "check_existing_files"
    memory '1g'
    cpus 1
    
    input:
    val output_dir
    val file_pattern
    val stage_name
    
    output:
    stdout
    
    script:
    """
    echo "Checking for existing ${stage_name} files..."
    find ${output_dir}/${stage_name} -name "${file_pattern}" 2>/dev/null | head -5
    count=\$(find ${output_dir}/${stage_name} -name "${file_pattern}" 2>/dev/null | wc -l)
    echo "Found \$count ${stage_name} files"
    """
}

process prepareTaxaBamMap {
    tag "prepare_taxa_bam_map_${pop}"
    memory '16g'
    cpus 4
    publishDir "${params.home_dir}/00data/05taxaBamMap", mode: 'copy'
    
    input:
    tuple path(workshop_jar), val(java_version), val(java_lib_dir)
    tuple val(pop), val(depth_dir), val(bam_dir)
    
    output:
    path "${pop}.taxaBamMap.txt", emit: taxa_bam_map
    path "${pop}.taxaRunMap.txt", emit: taxa_run_map
    path "prepare_tbm_${pop}.log", emit: log
    
    script:
    def output_bam_file = "${pop}.taxaBamMap.txt"
    def output_run_file = "${pop}.taxaRunMap.txt"
    def javaSetup = getJavaSetupScript(java_version, java_lib_dir)
    """
    echo "Preparing taxa-BAM mapping file for population ${pop}..." > prepare_taxa_bam_map_${pop}.log
    echo "Using Java 21 for prepareTaxaBamMap process" >> prepare_taxa_bam_map_${pop}.log
    
    ${javaSetup} >> prepare_taxa_bam_map_${pop}.log 2>&1
    
    # Use TaxaBamMap.java only
    echo "Running TaxaBamMap from Workshop jar..." >> prepare_taxa_bam_map_${pop}.log
    
    java -Xmx16g -jar "${workshop_jar}" \\
        -d ${depth_dir} \\
        -b ${bam_dir} \\
        -o ${output_bam_file} \\
        -t ${output_run_file} \\
        >> prepare_tbm_${pop}.log 2>&1
    
    echo "TaxaBamMap execution completed" >> prepare_taxa_bam_map_${pop}.log
    """
}

process concatTaxaBamMap {
    tag "concat_taxa_bam_map_${sub_genome}"
    memory '16g'
    cpus 4
    publishDir "${params.home_dir}/00data/05taxaBamMap", mode: 'copy'
    
    input:
    tuple path(taxa_bam_maps), val(job_name), val(sub_genome)
    
    output:
    tuple val(sub_genome), path("${job_name}.${sub_genome}.taxaBamMap.txt"), emit: tbm
    path "${sub_genome}.concat.log", emit: log
    
    script:
    """
    echo "Concatenating taxa-BAM mapping files (${taxa_bam_maps.size()})..." > ${sub_genome}.concat.log

    cat ${taxa_bam_maps.join(' ')} > ${job_name}.${sub_genome}.taxaBamMap.txt
    echo -e "Taxa\tCoverage-Of-All-Bams\tBams(A list of bams of the taxon, seperated by the delimiter of Tab)" > ${job_name}.${sub_genome}.taxaBamMap.txt.header
    cat ${job_name}.${sub_genome}.taxaBamMap.txt >> ${job_name}.${sub_genome}.taxaBamMap.txt.header
    mv ${job_name}.${sub_genome}.taxaBamMap.txt.header ${job_name}.${sub_genome}.taxaBamMap.txt

    echo "Concatenation completed. Output: ${job_name}.${sub_genome}.taxaBamMap.txt" >> ${sub_genome}.concat.log
    """
}

process fastcall3_disc {
    tag "disc_${chromosome}"
    publishDir("${params.output_dir}/${params.job}/ing", mode: 'copy', overwrite: true)
    
    input:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map)
    tuple val(tiger_jar), val(app_name), val(java_version), val(java_lib_dir)
    val samtools_path
    val ing_path
    
    output:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map), emit: chr_config
    path "disc_${chromosome}.log", emit: log
    
    script:
    def javaSetup = getJavaSetupScript(java_version, java_lib_dir)
    """
    echo "Starting disc analysis for chromosome ${chromosome}..." > disc_${chromosome}.log
    echo "Using ${java_version} for ${app_name} disc process" >> disc_${chromosome}.log
    
    ${javaSetup} >> disc_${chromosome}.log 2>&1
    
    # Monitor system resources
    echo "System resources before TIGER execution:" >> disc_${chromosome}.log
    echo "Memory: \$(free -h | grep Mem)" >> disc_${chromosome}.log
    echo "CPU cores (allocated): ${task.cpus}" >> disc_${chromosome}.log
    echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g" >> disc_${chromosome}.log
    echo "TIGER jar: ${tiger_jar}" >> disc_${chromosome}.log
    echo "FastCall version: ${app_name}" >> disc_${chromosome}.log
    
    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod disc \\
        -a ${reference} \\
        -b ${taxa_bam_map} \\
        -c 0 \\
        -d ${params.disc_min_MQ} \\
        -e ${params.disc_min_BQ} \\
        -f ${params.disc_min_DC} \\
        -g ${params.disc_min_DR} \\
        -h ${params.disc_max_DR} \\
        -i ${params.disc_HoR} \\
        -j ${params.disc_HeR} \\
        -k ${params.disc_TDR} \\
        -l ${chromosome} \\
        -m ${task.cpus} \\
        -n ${ing_path} \\
        -o ${samtools_path} \\
        >> disc_${chromosome}.log 2>&1

    echo "Disc analysis for chromosome ${chromosome} completed" >> disc_${chromosome}.log
    """
}

process fastcall3_blib {
    tag "blib_${chromosome}"
    publishDir("${params.output_dir}/${params.job}/vLib", mode: 'copy', overwrite: true)
    
    input:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map)
    tuple val(tiger_jar), val(app_name), val(java_version), val(java_lib_dir)
    val ing_path
    val vLib_path
    
    output:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map), path("${chromosome}_*.lib.gz"), emit: chr_config_lib
    path "blib_${chromosome}.log", emit: log
    
    script:
    def javaSetup = getJavaSetupScript(java_version, java_lib_dir)
    """
    echo "Starting blib generation for chromosome ${chromosome}..." > blib_${chromosome}.log
    echo "Using ${java_version} for ${app_name} blib process" >> blib_${chromosome}.log
    
    ${javaSetup} >> blib_${chromosome}.log 2>&1
    
    # Check if disc files exist
    # Check if disc files exist under ing_path for this chromosome (search recursively)
    if ! find ${ing_path} -type f -name "${chromosome}_*.ing.gz" | grep -q . ; then
        echo "Error: No .ing.gz files found for chromosome ${chromosome} under ${ing_path}" >> blib_${chromosome}.log
        exit 1
    fi
    echo "Counting input disc files under ${ing_path}..." >> blib_${chromosome}.log
    cnt_ing=\$(find ${ing_path} -type f -name "${chromosome}_*.ing.gz" | wc -l)
    echo "Found \$cnt_ing .ing.gz chunks for chromosome ${chromosome}" >> blib_${chromosome}.log
    
    echo "Input disc files (showing up to 5 examples):" >> blib_${chromosome}.log
    find ${ing_path} -type f -name "${chromosome}_*.ing.gz" | head -n 5 >> blib_${chromosome}.log 2>/dev/null || true
    
    # Monitor system resources
    echo "System resources before TIGER execution:" >> blib_${chromosome}.log
    echo "Memory: \$(free -h | grep Mem)" >> blib_${chromosome}.log
    echo "CPU cores (allocated): ${task.cpus}" >> blib_${chromosome}.log
    echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g" >> blib_${chromosome}.log
    echo "TIGER jar: ${tiger_jar}" >> blib_${chromosome}.log
    echo "FastCall version: ${app_name}" >> blib_${chromosome}.log
    
    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod blib \\
        -a ${reference} \\
        -b ${chromosome} \\
        -c ${params.blib_MAO} \\
        -d ${task.cpus} \\
        -e ${ing_path} \\
        -f ${vLib_path} \\
        >> blib_${chromosome}.log 2>&1
    
    echo "Blib generation for chromosome ${chromosome} completed" >> blib_${chromosome}.log

    # Symlink produced lib for this chromosome back into workdir so Nextflow can capture it
    if ls ${vLib_path}/${chromosome}_*.lib.gz >/dev/null 2>&1; then
        ln -sf ${vLib_path}/${chromosome}_*.lib.gz . 2>> blib_${chromosome}.log || true
    else
        # If naming differs, just link the latest .lib.gz
        latest_lib=\$(ls -1t ${vLib_path}/*.lib.gz 2>/dev/null | head -n1)
        [ -n "\$latest_lib" ] && ln -sf "\$latest_lib" . 2>> blib_${chromosome}.log || true
    fi
    """
}

process fastcall3_scan {
    tag "scan_${chromosome}"
    publishDir("${params.output_dir}/${params.job}/gen", mode: 'move', overwrite: false)
    
    input:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map), path(blib_file)
    tuple val(tiger_jar), val(app_name), val(java_version), val(java_lib_dir)
    val samtools_path
    val gen_path
    
    output:
    // tuple val(chromosome), path("*.vcf*"), emit: vcf_files
    path "scan_${chromosome}.log", emit: log
    
    script:
    def javaSetup = getJavaSetupScript(java_version, java_lib_dir)
    // FIX: Typo 'param.home_dir' -> 'params.home_dir'
    // def chr_config = getChrConfig(chromosome, params.home_dir)
    def chr_int = chromosome.toString().toInteger()
    def vcf_name = (chr_int < 10) ? "chr00${chr_int}.vcf" : "chr0${chr_int}.vcf"
    """
    echo "Starting scan analysis for chromosome ${chromosome}..." > scan_${chromosome}.log
    echo "Using ${java_version} for ${app_name} scan process" >> scan_${chromosome}.log
    
    ${javaSetup} >> scan_${chromosome}.log 2>&1
    
    # Resolve the lib file (staged by Nextflow) for this chromosome
    lib_file_path="${blib_file}"
    if [ ! -f "\$lib_file_path" ]; then
        echo "Error: Library file \$lib_file_path not found" >> scan_${chromosome}.log
        exit 1
    fi
    echo "Input library file: \$lib_file_path" >> scan_${chromosome}.log
    echo "Library file size: \$(stat -c%s \"\$lib_file_path\") bytes" >> scan_${chromosome}.log
    
    # Monitor system resources
    echo "System resources before TIGER execution:" >> scan_${chromosome}.log
    echo "Memory: \$(free -h | grep Mem)" >> scan_${chromosome}.log
    echo "CPU cores (allocated): ${task.cpus}" >> scan_${chromosome}.log
    echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g" >> scan_${chromosome}.log
    echo "TIGER jar: ${tiger_jar}" >> scan_${chromosome}.log
    echo "FastCall version: ${app_name}" >> scan_${chromosome}.log

    java -Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod scan \\
        -a ${reference} \\
        -b ${taxa_bam_map} \\
        -c "\$lib_file_path" \\
        -d ${chromosome} \\
        -e 0 \\
        -f ${params.scan_min_MQ} \\
        -g ${params.scan_min_BQ} \\
        -h ${params.scan_error_rate} \\
        -i ${samtools_path} \\
        -j ${task.cpus} \\
        -k ${gen_path} \\
        >> scan_${chromosome}.log 2>&1
    
    bgzip -f -@ ${task.cpus} ${gen_path}/VCF/${vcf_name}

    echo "Scan analysis for chromosome ${chromosome} completed" >> scan_${chromosome}.log
    """
}

process fastcall3_scan2 {
    tag "scan2_${chromosome}"
    publishDir("${params.output_dir}/${params.job}/gen", mode: 'move', overwrite: false)
    
    input:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map), path(blib_file)
    tuple val(tiger_jar), val(app_name), val(java_version), val(java_lib_dir)
    val samtools_path
    val gen_path
    
    output:
    // tuple val(chromosome), path("*.vcf*"), emit: vcf_files
    path "scan_${chromosome}.log", emit: log
    
    script:
    def javaSetup = getJavaSetupScript(java_version, java_lib_dir)
    // FIX: Typo 'param.home_dir' -> 'params.home_dir'
    // def chr_config = getChrConfig(chromosome, params.home_dir)
    def chr_int = chromosome.toString().toInteger()
    def vcf_name = (chr_int < 10) ? "chr00${chr_int}.vcf" : "chr0${chr_int}.vcf"
    """
    echo "Starting scan analysis for chromosome ${chromosome}..." > scan_${chromosome}.log
    echo "Using ${java_version} for ${app_name} scan process" >> scan_${chromosome}.log
    
    ${javaSetup} >> scan_${chromosome}.log 2>&1
    
    # Resolve the lib file (staged by Nextflow) for this chromosome
    lib_file_path="${blib_file}"
    if [ ! -f "\$lib_file_path" ]; then
        echo "Error: Library file \$lib_file_path not found" >> scan_${chromosome}.log
        exit 1
    fi
    echo "Input library file: \$lib_file_path" >> scan_${chromosome}.log
    echo "Library file size: \$(stat -c%s \"\$lib_file_path\") bytes" >> scan_${chromosome}.log
    
    # Monitor system resources
    echo "System resources before TIGER execution:" >> scan_${chromosome}.log
    echo "Memory: \$(free -h | grep Mem)" >> scan_${chromosome}.log
    echo "CPU cores (allocated): ${task.cpus}" >> scan_${chromosome}.log
    echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g" >> scan_${chromosome}.log
    echo "TIGER jar: ${tiger_jar}" >> scan_${chromosome}.log
    echo "FastCall version: ${app_name}" >> scan_${chromosome}.log

    java -Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod scan2 \\
        -a ${reference} \\
        -b ${taxa_bam_map} \\
        -c "\$lib_file_path" \\
        -d ${chromosome} \\
        -e 0 \\
        -f ${params.scan_min_MQ} \\
        -g ${params.scan_min_BQ} \\
        -h ${params.scan_error_rate} \\
        -i ${samtools_path} \\
        -j ${task.cpus} \\
        -k ${gen_path} \\
        >> scan_${chromosome}.log 2>&1
    
    bgzip -f -@ ${task.cpus} ${gen_path}/VCF/${vcf_name}

    echo "Scan analysis for chromosome ${chromosome} completed" >> scan_${chromosome}.log
    """
}

def getJobConfig(job, home_dir) {
    def jobConfigs = [
        "chr1": [
            name: "chr1",
            a_pop: ["A"],
            b_pop: [],
            d_pop: [],
            chroms: ["1"]
        ],
        "test": [
            name: "test",
            a_pop: ["test"],
            b_pop: ["test"],
            d_pop: ["test"],
            chroms: []
        ],
        "test_chr1": [
            name: "test_chr1",
            a_pop: ["test_chr1"],
            b_pop: ["test_chr1"],
            d_pop: ["test_chr1"],
            chroms: []
        ],
        "test_chr12": [
            name: "test_chr12",
            a_pop: ["test_chr12"],
            b_pop: [],
            d_pop: [],
            chroms: ["1", "2"]
        ],
        "all": [
            name: "all",
            a_pop: ["A", "AB", "ABD", "WAP", "HZNU", "Nature", "w115", "w66", "w203", "w204", "w243"],
            b_pop: ["S", "AB", "ABD", "WAP", "HZNU", "Nature", "w115", "w66", "w203", "w204", "w243"],
            d_pop: ["D",       "ABD", "WAP", "HZNU", "Nature", "w115", "w66", "w203", "w204", "w243"],
            chroms: []
        ]
    ]
    
    if (jobConfigs.containsKey(job)) {
        return jobConfigs[job]
    } else {
        log.error "Unknown job configuration: ${job}"
        log.error "Available jobs: ${jobConfigs.keySet().join(', ')}"
        exit 1
    }
}

def getServerPopulations(tbm_gen_server) {
    def serverPopulations = [
        "s115": ["S", "D", "A", "AB", "ABD", "WAP", "HZNU", "Nature"],
        "s107": ["S", "D", "A", "AB", "ABD", "WAP", "HZNU", "Nature", "test_chr12"],
        "s66": ["S", "D", "A", "AB", "ABD", "WAP", "HZNU", "Nature"],
        "s203": ["S", "D", "A", "AB", "ABD", "WAP", "HZNU", "Nature"],
        "s204": ["S", "D", "A", "AB", "ABD", "WAP", "HZNU", "Nature"],
        "s243": ["S", "D", "A", "AB", "ABD", "WAP", "HZNU", "Nature", "w115", "w66", "w203", "w204", "w243"]
    ]

    return serverPopulations[tbm_gen_server]
}

// Population configuration function
def getPopulationConfig(pop, home_dir) {
    def popConfigs = [
        "chr1": [
            bam_dir: "${home_dir}/00data/02bam/bam1/A",
            depth_dir: "${home_dir}/00data/04depth/01A"
        ],
        "test": [
            bam_dir: "${home_dir}/01testData/02bam/01test",
            depth_dir: "${home_dir}/01testData/04depth/01test"
        ],
        "test_chr1": [
            bam_dir: "${home_dir}/01testData/02bam/02test1chr",
            depth_dir: "${home_dir}/01testData/04depth/01test"
        ],
        "test_chr12": [
            bam_dir: "${home_dir}/00data/02bam/bam1/A",
            depth_dir: "${home_dir}/00data/04depth/09test12"
        ],
        "A": [
            bam_dir: "${home_dir}/00data/02bam/bam1/A",
            depth_dir: "${home_dir}/00data/04depth/01A"
        ],
        "AB": [
            bam_dir: "${home_dir}/00data/02bam/bam1/AB",
            depth_dir: "${home_dir}/00data/04depth/02AB"
        ],
        "ABD": [
            bam_dir: "${home_dir}/00data/02bam/bam1/ABD",
            depth_dir: "${home_dir}/00data/04depth/03ABD"
        ],
        "D": [
            bam_dir: "${home_dir}/00data/02bam/bam1/D",
            depth_dir: "${home_dir}/00data/04depth/04D"
        ],
        "HZNU": [
            bam_dir: "${home_dir}/00data/02bam/bam1/HZNU",
            depth_dir: "${home_dir}/00data/04depth/05HZNU"
        ],
        "Nature": [
            bam_dir: "${home_dir}/00data/02bam/bam1/Nature",
            depth_dir: "${home_dir}/00data/04depth/06Nature"
        ],
        "S": [
            bam_dir: "${home_dir}/00data/02bam/bam1/S",
            depth_dir: "${home_dir}/00data/04depth/07S"
        ],
        "WAP": [
            bam_dir: "${home_dir}/00data/02bam/bam1/ABD",
            depth_dir: "${home_dir}/00data/04depth/08WAP"
        ],
        "w115": [
            bam_dir: "${home_dir}/00data/02bam/bam2/115",
            depth_dir: "${home_dir}/00data/02bam/bam2/115"
        ],
        "w203": [
            bam_dir: "${home_dir}/00data/02bam/bam2/203",
            depth_dir: "${home_dir}/00data/02bam/bam2/203"
        ],
        "w204": [
            bam_dir: "${home_dir}/00data/02bam/bam2/204",
            depth_dir: "${home_dir}/00data/02bam/bam2/204"
        ],
        "w243": [
            bam_dir: "${home_dir}/00data/02bam/bam2/243",
            depth_dir: "${home_dir}/00data/02bam/bam2/243"
        ],
        "w66": [
            bam_dir: "${home_dir}/00data/02bam/bam2/66",
            depth_dir: "${home_dir}/00data/02bam/bam2/66"
        ],
    ]
    
    if (!popConfigs.containsKey(pop)) {
        // log.error "Unknown population configuration: ${pop}"
        def validPops = popConfigs.keySet().join(", ")
        throw new Exception("Unknown population: ${pop}. Valid options: ${validPops}")
    }

    def config = popConfigs[pop]
    
    return tuple(pop, config.depth_dir, config.bam_dir)
}

def getChrConfig(chr, home_dir, sub_genome_tbm) {
    // Normalize possible prefixes like 'chr1' -> '1'
    def normalized = chr.toString().replaceFirst(/^chr/, '')

    def subGenomeConfigs = [
        "A": [
            chromosome: chr,
            reference: "${home_dir}/00data/03ref/01A/a_iwgscV1.fa.gz",
            fai_idx: "${home_dir}/00data/03ref/01A/a_iwgscV1.fa.gz.fai",
            gzi_idx: "${home_dir}/00data/03ref/01A/a_iwgscV1.fa.gz.gzi"
        ],
        "B": [
            chromosome: chr,
            reference: "${home_dir}/00data/03ref/05B/b_iwgscV1.fa.gz",
            fai_idx: "${home_dir}/00data/03ref/05B/b_iwgscV1.fa.gz.fai",
            gzi_idx: "${home_dir}/00data/03ref/05B/b_iwgscV1.fa.gz.gzi"
        ],
        "D": [
            chromosome: chr,
            reference: "${home_dir}/00data/03ref/04D/d_iwgscV1.fa.gz",
            fai_idx: "${home_dir}/00data/03ref/04D/d_iwgscV1.fa.gz.fai",
            gzi_idx: "${home_dir}/00data/03ref/04D/d_iwgscV1.fa.gz.gzi"
        ],
        "ALL": [
            chromosome: chr,
            reference: "${home_dir}/00data/03ref/01A/a_iwgscV1.fa.gz",
            fai_idx: "${home_dir}/00data/03ref/01A/a_iwgscV1.fa.gz.fai",
            gzi_idx: "${home_dir}/00data/03ref/01A/a_iwgscV1.fa.gz.gzi"
        ]
    ]

    def groupA = ["1","2","7","8","13","14","19","20","25","26","31","32","37","38"]
    def groupB = ["3","4","9","10","15","16","21","22","27","28","33","34","39","40"]
    def groupD = ["5","6","11","12","17","18","23","24","29","30","35","36","41","42"]
    def groupOthers = ["0","43","44"]

    def sub_genome
    if (groupA.contains(normalized)) {
        sub_genome = "A"
    } else if (groupB.contains(normalized)) {
        sub_genome = "B"
    } else if (groupD.contains(normalized)) {
        sub_genome = "D"
    } else if (groupOthers.contains(normalized)) {
        sub_genome = "ALL"
    } else {
        throw new IllegalArgumentException("Unknown chromosome: ${chr} (normalized: ${normalized})")
    }
    def config = subGenomeConfigs[sub_genome]
    def ref = file(config.reference)
    def fai = file(config.fai_idx)
    def gzi = file(config.gzi_idx)
    def tbm_file = file(sub_genome_tbm[sub_genome])

    return tuple(config.chromosome, ref, fai, gzi, tbm_file)
}

// Java version management function
def getJavaSetupScript(javaVersion, javaLibDir) {
    def javaVersionMap = [
        "java8": "jdk-8",
        "java11": "jdk-11", 
        "java17": "jdk-17",
        "java21": "jdk-21"
    ]
    
    def javaDir = javaVersionMap[javaVersion]
    if (!javaDir) {
        def validVersions = javaVersionMap.keySet().join(", ")
        throw new Exception("Unknown Java version: ${javaVersion}. Valid options: ${validVersions}")
    }
    
    def javaHome = "${javaLibDir}/${javaDir}"
    
    return """
    # Switch to ${javaVersion} (${javaDir})
    export JAVA_HOME=${javaHome}
    export PATH=\$JAVA_HOME/bin:\$PATH
    
    # Verify Java version
    if [ ! -d "\$JAVA_HOME" ]; then
        echo "Error: Java directory not found: \$JAVA_HOME"
        echo "Available Java versions in ${javaLibDir}:"
        ls -la ${javaLibDir}/ || echo "Java lib directory not accessible"
        exit 1
    fi
    
    java -version 2>&1 | head -3
    """.stripIndent()
}

def getSoftwareConfig(home_dir, user_dir, tiger_jar, workshop_jar, samtools) {

    def java_lib_dir = "${user_dir}/lib/jvm"

    def resolved_tiger_jar = "${home_dir}/lib/${tiger_jar}"
    def resolved_workshop_jar = "${home_dir}/lib/${workshop_jar}"

    def tiger_jar_config = getTigerJarConfig(resolved_tiger_jar, java_lib_dir)
    def workshop_jar_config = getWorkshopJarConfig(resolved_workshop_jar, java_lib_dir)
    def samtools_config = getSamtoolsConfig(samtools)
    
    return [
        tiger_jar_config: tiger_jar_config,
        workshop_jar_config: workshop_jar_config,
        samtools_config: samtools_config
    ]
}

// TIGER jar configuration function
def getTigerJarConfig(tigerJarPath, java_lib_dir) {
    def tiger_jar_versions = [
        "TIGER_F3_20251121.jar": [
            java_version: "java17",
            app_name: "FastCall3"
        ],
        "TIGER_F3_20251118.jar": [
            java_version: "java17",
            app_name: "FastCall3"
        ],
        "TIGER_F3_20251016.jar": [
            java_version: "java17",
            app_name: "FastCall3"
        ],
        "TIGER_F3_20250915.jar": [
            java_version: "java17",
            app_name: "FastCall3"
        ],
        "TIGER_20250526.jar": [
            java_version: "java8", 
            app_name: "FastCall2"
        ]
    ]

    def jarFile = file(tigerJarPath)
    def jarName = jarFile.name
    
    // Get configuration from version mapping
    def config = tiger_jar_versions[jarName]
    
    if (!config) {
        // Try to infer from filename patterns
        if (jarName.contains("F3") || jarName.contains("FastCall3")) {
            log.warn "Unknown TIGER jar: ${jarName}. Assuming FastCall3 configuration."
            config = [
                java_version: "java17",
                app_name: "FastCall3"
            ]
        } else if (jarName.contains("2023") || jarName.contains("FastCall2")) {
            log.warn "Unknown TIGER jar: ${jarName}. Assuming FastCall2 configuration."
            config = [
                java_version: "java8",
                app_name: "FastCall2" 
            ]
        } else {
            log.warn "Unknown TIGER jar: ${jarName}. Using default FastCall3 configuration."
            config = [
                java_version: "java17",
                app_name: "FastCall3"
            ]
        }
    }
    
    return [
        path: tigerJarPath,
        app_name: config.app_name,
        java_version: config.java_version,
        java_lib_dir: java_lib_dir,
    ]
}

def getWorkshopJarConfig(workshopJarPath, java_lib_dir) {
    return [
        path: workshopJarPath,
        java_version: "java21",
        java_lib_dir: java_lib_dir
    ]
}

def getSamtoolsConfig(samtoolsPath) {
    return [
        path: samtoolsPath
    ]
}


// Path validation function
def validatePaths(pathMap) {
    def errors = []
    def isValid = true
    
    pathMap.each { name, path ->
        if (!path) {
            errors << "${name} is not specified"
            isValid = false
        } else {
            def pathFile = file(path)
            if (!pathFile.exists()) {
                errors << "${name} does not exist: ${path}"
                isValid = false
            }
        }
    }
    
    return [isValid: isValid, errors: errors]
}

def helpMessage() {
    log.info """
    ========================================
    FastCall3 Advanced Pipeline
    ========================================
    
    Usage:
        nextflow run runFastCall3.nf \
          --home_dir <HOME> \
          --user_dir <USER_BASE> \
          --job <JOB_NAME> \
          --mod <full|disc_only|blib_only|scan_only|from_disc|from_blib>

    Required parameters:
        --home_dir       Project home with data/libs (e.g. /data/.../vmap4)
        --user_dir       User base containing lib/jvm (e.g. /data/dazheng)
        --job            Job name (used for outputs under output/<job>)
        --mod            Workflow mode: full | disc_only | blib_only | scan_only | from_disc | from_blib

    Optional directories:
        --ing_dir        Disc outputs (.ing.gz). Default: output/<job>/ing
        --vlib_dir       Blib outputs (.lib.gz). Default: output/<job>/vLib
        --gen_dir        Scan outputs (VCF). Default: output/<job>/gen

    Software settings:
        --tiger_jar      TIGER jar filename under <home_dir>/lib (default: auto-detect TIGER_F3_*.jar)
        --workshop_jar   Workshop jar filename under <home_dir>/lib (default: Workshop.jar)
        --samtools       Samtools executable (default: samtools)

    Chromosome selection:
        --chr            Chromosome(s) to run. Multiple can be separated by '|', e.g. "1|2|7".
                         If omitted, chromosomes are derived from job configuration A/B/D groups + others (0,43,44).

    Taxa-BAM mapping (TBM):
        --tbm_gen_server Optional server key to generate TBM (e.g. s243). When provided, TBMs are prepared
                         per server populations, then concatenated for ALL/A/B/D.

    Java/runtime management:
        - Java installs are expected in <user_dir>/lib/jvm: jdk-8, jdk-11, jdk-17, jdk-21
        - prepareTaxaBamMap uses Java 21
        - FastCall3 (disc/blib/scan) uses the Java version inferred from --tiger_jar (typically Java 17)

    Disc parameters:
        --disc_min_MQ    Minimum mapping quality (default: 30)
        --disc_min_BQ    Minimum base quality (default: 20)
        --disc_min_DC    Minimum read depth count (default: 2)
        --disc_min_DR    Minimum read depth ratio (default: 0.2)
        --disc_max_DR    Maximum read depth ratio (default: 3)
        --disc_HoR       Homozygous ratio (default: 0.8)
        --disc_HeR       Heterozygous ratio (default: 0.35)
        --disc_TDR       Third allele depth ratio (default: 0.2)

    Blib parameters:
        --blib_MAO       Minor allele occurrence threshold (default: 2)

    Scan parameters:
        --scan_min_MQ    Minimum mapping quality (default: 30)
        --scan_min_BQ    Minimum base quality (default: 20)
        --scan_error_rate Combined error rate of sequencing/misalignment (default: 0.05)

    Execution semantics:
        - full: blib starts per chromosome only after disc emits .ing.gz for that chromosome;
                scan starts per chromosome only after blib emits the corresponding .lib.gz.
        - from_disc: chromosomes are derived from existing --ing_dir and filtered by reference .fai.
        - scan_only/from_blib: consumes existing .lib.gz under --vlib_dir for selected chromosomes.
        - Staging avoids input filename collisions; blib discovers .ing.gz files at runtime.

    Examples:
        # Full pipeline with defaults
        nextflow run runFastCall3.nf \
          --home_dir /data/project/vmap4 \
          --user_dir /data/dazheng \
          --job test_full \
          --mod full

        # From existing disc (provide ing/vlib/gen roots)
        nextflow run runFastCall3.nf \
          --home_dir /data/project/vmap4 \
          --user_dir /data/dazheng \
          --job test_from_disc \
          --mod from_disc \
          --tiger_jar TIGER_F3_20251016.jar \
          --ing_dir /path/to/ing \
          --vlib_dir /path/to/vLib \
          --gen_dir /path/to/gen \
          --chr 1|2

        # Only scan from existing blib
        nextflow run runFastCall3.nf \
          --home_dir /data/project/vmap4 \
          --user_dir /data/dazheng \
          --job test_scan \
          --mod scan_only \
          --vlib_dir /path/to/vLib \
          --gen_dir /path/to/gen \
          --chr 1
    """.stripIndent()
}

