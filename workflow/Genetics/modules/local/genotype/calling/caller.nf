#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { getCallingJobConfig; getPopulationConfig; getChrConfig } from '../../infra/infra_job_config.nf'
include { getSoftwareConfig; validatePaths; getJavaSetupScript } from '../../infra/infra_tools.nf'
include { getServerPopulations } from '../../infra/infra_tiger.nf'

include {
    check_existing_files
    prepareTaxaBamMap
    concatTaxaBamMap
} from './caller_prep.nf'

include {
    fastcall3_disc
    fastcall3_blib
    fastcall3_scan
    fastcall3_scan2
} from './caller_fastcall.nf'

/*
 * Advanced FastCall3 pipeline for variant calling with enhanced features
 * Author: SirLearning
 * Date: 2026-02-02
 * version: 1.0
 */

// --- help info ---
params.help = false

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
        // version 0.1
        nextflow run /data/home/tusr1/01projects/DataProcess/calling/run/runFastCall3.nf \
            --home_dir /data/home/tusr1/01projects/vmap4 \
            --java_lib /data/home/tusr1/lib/jvm \
            --pop A \
            --job run_A_disc \
            --workflow_mode disc_only \
            --tiger_jar TIGER_F3_20250915.jar"
        
        // version 1.0
        nextflow run /data/home/tusr1/git/script/workflow/Genetics/modules/local/genotype/calling/caller.nf \
            --home_dir /data/home/tusr1/01projects/vmap4 \
            --user_dir /data/home/tusr1 \
            --output_dir /data/home/tusr1/01projects/vmap4/04runScreens \
            --vlib_dir {/path/to/vLib} \
            --gen_dir {/path/to/gen} \
            --server {s243, s107} \
            --job {all, test_chr12} \
            --mod {full, from_blib, scan_only, scan2} \
            {--gen_tbm 1} \
            {--chr 1|2} \
            {--tiger_jar TIGER_F3_20251016.jar}
        // blib MAO=1 rebuild
        nextflow run /data/home/tusr1/git/script/workflow/Genetics/modules/local/genotype/calling/caller.nf \
            --home_dir /data/home/tusr1/01projects/vmap4 \
            --user_dir /data/home/tusr1 \
            --output_dir /data/home/tusr1/01projects/vmap4/04runScreens \
            --ing_dir /data/home/tusr1/01projects/vmap4/04runScreens/all/ing \
            --server s243 \
            --gen_tbm 1 \
            --job rebuild \
            --mod from_disc \
            --chr 2 \
            --blib_MAO 1 \
            --tiger_jar TIGER_F3_2M_scan_20260129.jar
        nextflow run /data/home/tusr1/git/script/workflow/Genetics/modules/local/genotype/calling/caller.nf \
            --home_dir /data/home/tusr1/01projects/vmap4 \
            --user_dir /data/home/tusr1 \
            --output_dir /data/home/tusr1/01projects/vmap4/04runScreens \
            --server s243 \
            --job all \
            --mod scan2 \
            --chr 32
    
    screen prefix commands:
        // run_A_disc
        screen -dmS run_A_disc bash -c "\
        cd /data/home/tusr1/01projects/vmap4/04runScreens/01A/disc && \
        source ~/.bashrc && conda activate run && \
        "
        // all
        screen -dmS all bash -c "\
        cd /data/home/tusr1/01projects/vmap4/04runScreens/02run && \
        source ~/.bashrc && conda activate run && \
        "
        screen -dmS all bash -c "\
        cd /data/home/tusr1/01projects/vmap4/04runScreens/06run.32 && \
        source ~/.bashrc && conda activate run && \
        "
    screen -dmS all bash -c "\
        cd /data/home/tusr1/01projects/vmap4/04runScreens/06run.32 && \
        source ~/.bashrc && conda activate run && \
        nextflow run /data/home/tusr1/git/script/workflow/Genetics/modules/local/genotype/calling/caller.nf \
            --home_dir /data/home/tusr1/01projects/vmap4 \
            --user_dir /data/home/tusr1 \
            --output_dir /data/home/tusr1/01projects/vmap4/04runScreens \
            --server s243 \
            --job all \
            --mod from_blib \
            --chr 32"
    """.stripIndent()
}

workflow {
    run_FastCall3()
}

workflow run_FastCall3 {
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
    def output_path = "${params.output_dir}/${params.job}"

    // --- retrieve parameters ---
    // 1. directory
    def ing_dir_resolved  = (params.ing_dir  ?: "${output_path}/ing").toString().replaceAll(/\/+/, "/")
    def vlib_dir_resolved = (params.vlib_dir ?: "${output_path}/vLib").toString().replaceAll(/\/+/, "/")
    def gen_dir_resolved  = (params.gen_dir  ?: "${output_path}/gen").toString().replaceAll(/\/+/, "/")
    // 2. work specific config
    def job_config = getCallingJobConfig(params.job)
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
        pathValidation.errors.each { err -> log.error "  - ${err}" }
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
        def pop_ch = channel.fromList(server_pop)
        def pop_input_ch = pop_ch.map { pop ->
            def pop_config = getPopulationConfig(pop, params.home_dir)
            return pop_config
        }
        prepared_tbm_ch = prepareTaxaBamMap(workshop_jar_input, pop_input_ch).taxa_bam_map
    } else {
        tbm_files = channel.fromList(server_pop)
            .map { pop -> tuple(pop, file("${params.home_dir}/00data/05taxaBamMap/${pop}.taxaBamMap.txt")) }

        log.info "Mode: Using existing taxaBamMap files."
        def existing_files = server_pop
            .collect { pop -> file("${params.home_dir}/00data/05taxaBamMap/${pop}.taxaBamMap.txt") }
            .findAll { f -> f.exists() }
        if (existing_files.size() != server_pop.size()) {
            log.warn "Missing some existing taxaBamMap files. Found ${existing_files.size()}/${server_pop.size()}."
        }
        prepared_tbm_ch = channel.fromPath(existing_files)
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
        .flatMap { list -> list }
        .map { item ->
            log.info "Processing sub-genome: ${item[0]}"
            log.info "  Included taxaBamMap files: ${item[1].collect { f -> f.name }.toList()}"
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
            // chromosomes += ["20"]  // temporary remove some chr for finished job
            chromosomes += ["14"]  // temporary remove some chr for finished job
        }
        if (job_config.b_pop && job_config.b_pop.size() > 0) {
            // chromosomes += ["3","4","9","10","15","16","21","22","27","28","33","34","39","40"]
            // chromosomes += ["4","9","10","15","16","21","22","27","28","33","34","39","40"]
            // chromosomes += ["9","10","15","16","21","22","27","33","34","39","40"]
            // chromosomes += ["9","10","15","21","22","27","33","39","40"]
            // chromosomes += ["10","27","40"]
            // chromosomes += ["3","27","33"]
            chromosomes += []
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
            // chromosomes += ["5"]
            chromosomes += []
        }
    }
    if (params.chr) {
        chromosomes = params.chr.toString().tokenize("|").collect { chr -> chr.trim() }
        log.info "Using single chromosome from --chr: ${chromosomes[0]}"
    }
    log.info "Pipeline setup completed. Processing ${chromosomes.size()} chromosomes: ${chromosomes.size() > 0 ? chromosomes[0] + '...' + chromosomes[-1] : 'none'}"
    // Create chromosome channel
    chr_config_ch = sub_genome_tbm_ch
        .map { sub_genome_tbm ->
            return chromosomes.collect { chromosome ->
                getChrConfig(chromosome, params.home_dir, sub_genome_tbm)
            }
        }
        .flatMap { list -> list }

    log.info "Running workflow in mode: ${params.mod}"
    if (params.mod == "full") {
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
    } else if (params.mod == "disc_only") {
        // Only run disc generation
        disc_results = fastcall3_disc(
            chr_config_ch, 
            tiger_jar_input,
            samtools_input,
            ing_dir_resolved
        )
        log.info "DISC generation completed. Use 'from_disc' mode to continue."
    } else if (params.mod == "from_disc") {
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
    } else if (params.mod == "blib_only") {
        // Only run blib generation from existing disc (ing) results
        blib_results = fastcall3_blib(
            chr_config_ch,
            tiger_jar_input,
            ing_dir_resolved,
            vlib_dir_resolved
        )
        log.info "Blib generation completed. Use 'from_blib' mode to continue."
    } else if (params.mod == "from_blib" || params.mod == "scan_only") {
        // Start from existing blib (vLib) results (only scan)
        def blib_chr_config_lib_ch = load_lib_files(
            vlib_dir_resolved,
            chr_config_ch
        ).ch_lib
        scan_results = fastcall3_scan(
            blib_chr_config_lib_ch,
            tiger_jar_input,
            samtools_input,
            gen_dir_resolved
        )
    } else if (params.mod == "scan2") {
        // Start from existing blib (vLib) results (only scan)
        def blib_chr_config_lib_ch = load_lib_files(
            vlib_dir_resolved,
            chr_config_ch
        ).ch_lib
        scan_results = fastcall3_scan2(
            blib_chr_config_lib_ch,
            tiger_jar_input,
            samtools_input,
            gen_dir_resolved
        )
    } else {
        log.error "Unknown workflow mode: ${params.mod}"
        log.error "Available modes: full, disc_only, blib_only, scan_only, from_disc, from_blib"
        exit 1
    }    
}

workflow load_lib_files {
    /*
     * Load library files from a given directory matching a pattern
     */
    take:
    vlib_dir_resolved
    chr_config_ch
    
    main:
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
        .filter { f -> f != null }
    
    emit:
    ch_lib = blib_chr_config_lib_ch
}
