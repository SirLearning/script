#!/usr/bin/env nextflow

/*
 * Advanced FastCall3 pipeline for variant calling with enhanced features
 * Author: Based on FastCall2 workflow with improvements
 * Date: 2025-09-11
 */

nextflow.enable.dsl=2

// Parameters
params.home_dir = null
params.bam_dir = null
params.depth_dir = null  // Add depth directory parameter
params.reference = null
params.pop = null
params.job = null
params.tiger_jar = "TIGER_F3_20250911.jar"
params.workshop_jar = null
params.samtools_path = null
params.output_dir = "output"
params.threads = 16
params.memory = "64g"
params.help = false
params.java_lib = "/data/dazheng/lib/jvm"  // Java installation base directory
// Single chromosome selector (optional). If set, only this chromosome will be processed.
params.chr = null
// Optional: user-provided taxa-BAM mapping file to bypass generation
params.taxa_bam_map = null
// Output sub-folders (can be overridden at runtime). ing_dir, vlib_dir, gen_dir determine where .ing.gz, .lib.gz and VCFs go
// Note: Do not assign defaults here to avoid Nextflow "defined multiple times" warnings when users pass CLI values.

// TIGER jar version compatibility and configuration
params.tiger_jar_versions = [
    "TIGER_F3_20251016.jar": [
        java_version: "java17",
        fastcall_version: "FastCall3",
        app_name: "FastCall3"
    ],
    "TIGER_F3_20250915.jar": [
        java_version: "java17",
        fastcall_version: "FastCall3",
        app_name: "FastCall3"
    ],
    "TIGER_20250526.jar": [
        java_version: "java8", 
        fastcall_version: "FastCall2",
        app_name: "FastCall2"
    ]
]

// Workflow control parameters
params.workflow_mode = "full"  // Options: "full", "disc_only", "blib_only", "scan_only", "from_disc", "from_blib"
params.skip_taxa_map = false   // Skip taxa-BAM map generation if already exists
params.resume_from_checkpoint = false  // Resume from existing intermediate files

// FastCall3 disc parameters
params.disc_min_depth = 30
params.disc_min_qual = 20
params.disc_min_snp_count = 2
params.disc_min_allele_freq = 0.2
params.disc_min_coverage = 3
params.disc_max_missing = 0.8
params.disc_min_het_freq = 0.35
params.disc_max_het_freq = 0.2

// FastCall3 scan parameters
params.scan_min_depth = 30
params.scan_min_qual = 20
params.scan_p_value = 0.05

// Chromosome list (can be modified based on your reference)
params.chromosomes = (0..44).collect { it.toString() }

// execute by screen command line (on 243):
// screen -dmS run_A_disc bash -c "cd /data/home/tusr1/01projects/runScreens/01A/disc && source ~/.bashrc && conda activate run && nextflow run /data/home/tusr1/01projects/DataProcess/calling/run/runFastCall3.nf --home_dir /data/home/tusr1/01projects/vmap4 --java_lib /data/home/tusr1/lib/jvm --pop A --job run_A_disc --workflow_mode disc_only --tiger_jar TIGER_F3_20250915.jar"
// nextflow run /data/dazheng/git/script/DataProcess/calling/run/runFastCall3.nf --home_dir /data/dazheng/01projects/vmap4 --java_lib /data/dazheng/lib/jvm --pop ABD --job test_ABD --workflow_mode disc_only --tiger_jar TIGER_F3_20250915.jar
// nextflow run /data/dazheng/git/script/DataProcess/calling/run/runFastCall3.nf --home_dir /data/dazheng/01projects/vmap4 --java_lib /data/dazheng/lib/jvm --pop chr1 --job test_ABD --workflow_mode disc_only --tiger_jar TIGER_F3_20250915.jar

// nextflow run /data/dazheng/git/script/DataProcess/calling/run/runFastCall3.nf --home_dir /data/dazheng/01projects/vmap4 --java_lib /data/dazheng/lib/jvm --pop chr1 --job test_ABD --workflow_mode from_disc --tiger_jar TIGER_F3_20251013.jar --ing_dir /data/dazheng/01projects/vmap4/04testFastCall3/01chr1/ing --vlib_dir /data/dazheng/01projects/vmap4/04testFastCall3/01chr1/vLib --gen_dir /data/dazheng/01projects/vmap4/04testFastCall3/01chr1/gen

// nextflow run /data/dazheng/git/script/DataProcess/calling/run/runFastCall3.nf --home_dir /data/dazheng/01projects/vmap4 --java_lib /data/dazheng/lib/jvm --pop chr1 --job test_ABD --workflow_mode from_disc --tiger_jar TIGER_F3_20251016.jar --ing_dir /data/dazheng/01projects/vmap4/04testFastCall3/01chr1/ing --vlib_dir /data/dazheng/01projects/vmap4/04testFastCall3/01chr1/vLib --gen_dir /data/dazheng/01projects/vmap4/04testFastCall3/01chr1/gen

// Population configuration function
def getPopulationConfig(pop, home_dir) {
    def popConfigs = [
        "chr1": [
            bam_dir: "${home_dir}/00data/02bam/bam1/A",
            depth_dir: "${home_dir}/00data/04depth/01A",
            reference: "${home_dir}/00data/03ref/05chr1/chr1.fa.gz",
            description: "Test dataset chromosome 1"
        ],
        "test": [
            bam_dir: "${home_dir}/01testData/02bam/01test",
            depth_dir: "${home_dir}/01testData/04depth/01test",
            reference: "${home_dir}/01testData/03ref/abd_1M.fa",
            description: "Test dataset"
        ],
        "test_chr1": [
            bam_dir: "${home_dir}/01testData/02bam/02test1chr",
            depth_dir: "${home_dir}/01testData/04depth/01test",
            reference: "${home_dir}/01testData/03ref/abd_1_1M.fa",
            description: "Test dataset chromosome 1"
        ],
        "A": [
            bam_dir: "${home_dir}/00data/02bam/bam1/A",
            depth_dir: "${home_dir}/00data/04depth/01A",
            reference: "${home_dir}/00data/03ref/01A/a_iwgscV1.fa.gz",
            description: "A genome population"
        ],
        "AB": [
            bam_dir: "${home_dir}/00data/02bam/bam1/AB",
            depth_dir: "${home_dir}/00data/04depth/02AB",
            reference: "${home_dir}/00data/03ref/02AB/ab_iwgscV1.fa.gz",
            description: "AB genome population"
        ],
        "ABD": [
            bam_dir: "${home_dir}/00data/02bam/bam1/ABD",
            depth_dir: "${home_dir}/00data/04depth/03ABD",
            reference: "${home_dir}/00data/03ref/03ABD/abd_iwgscV1.fa.gz",
            description: "ABD genome population"
        ],
        "D": [
            bam_dir: "${home_dir}/00data/02bam/bam1/D",
            depth_dir: "${home_dir}/00data/04depth/04D",
            reference: "${home_dir}/00data/03ref/04D/d_iwgscV1.fa.gz",
            description: "D genome population"
        ],
        "HZNU": [
            bam_dir: "${home_dir}/00data/02bam/bam1/HZNU",
            depth_dir: "${home_dir}/00data/04depth/05HZNU",
            reference: "${home_dir}/00data/03ref/03ABD/abd_iwgscV1.fa.gz",
            description: "HZNU collection"
        ],
        "Nature": [
            bam_dir: "${home_dir}/00data/02bam/bam1/Nature",
            depth_dir: "${home_dir}/00data/04depth/06Nature",
            reference: "${home_dir}/00data/03ref/03ABD/abd_iwgscV1.fa.gz",
            description: "Nature publication samples"
        ],
        "S": [
            bam_dir: "${home_dir}/00data/02bam/bam1/S",
            depth_dir: "${home_dir}/00data/04depth/07S",
            reference: "${home_dir}/00data/03ref/06B/b_iwgscV1.fa.gz",
            description: "S genome population"
        ],
        "WAP": [
            bam_dir: "${home_dir}/00data/02bam/bam1/ABD",
            depth_dir: "${home_dir}/00data/04depth/08WAP",
            reference: "${home_dir}/00data/03ref/03ABD/abd_iwgscV1.fa.gz",
            description: "WAP collection"
        ],
        "WAP2": [
            bam_dir: "${home_dir}/00data/02bam/bam1/ABD",
            depth_dir: "${home_dir}/00data/04depth/08WAP",
            reference: "${home_dir}/00data/03ref/05WAP2/wap_batch2.fa.gz",
            description: "WAP collection"
        ],
        "w115": [
            bam_dir: "${home_dir}/00data/02bam/bam2/115",
            depth_dir: "${home_dir}/00data/02bam/bam2/115",
            reference: "${home_dir}/00data/03ref/03ABD/abd_iwgscV1.fa.gz",
            description: "115 collection"
        ],
        "w203": [
            bam_dir: "${home_dir}/00data/02bam/bam2/203",
            depth_dir: "${home_dir}/00data/02bam/bam2/203",
            reference: "${home_dir}/00data/03ref/03ABD/abd_iwgscV1.fa.gz",
            description: "203 collection"
        ],
        "w204": [
            bam_dir: "${home_dir}/00data/02bam/bam2/204",
            depth_dir: "${home_dir}/00data/02bam/bam2/204",
            reference: "${home_dir}/00data/03ref/03ABD/abd_iwgscV1.fa.gz",
            description: "204 collection"
        ],
        "w243": [
            bam_dir: "${home_dir}/00data/02bam/bam2/243",
            depth_dir: "${home_dir}/00data/02bam/bam2/243",
            reference: "${home_dir}/00data/03ref/03ABD/abd_iwgscV1.fa.gz",
            description: "243 collection"
        ],
        "w66": [
            bam_dir: "${home_dir}/00data/02bam/bam2/66",
            depth_dir: "${home_dir}/00data/02bam/bam2/66",
            reference: "${home_dir}/00data/03ref/03ABD/abd_iwgscV1.fa.gz",
            description: "66 collection"
        ],
    ]
    
    if (!popConfigs.containsKey(pop)) {
        def validPops = popConfigs.keySet().join(", ")
        throw new Exception("Unknown population: ${pop}. Valid options: ${validPops}")
    }
    
    return popConfigs[pop]
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

// TIGER jar configuration function
def getTigerJarConfig(tigerJarPath) {
    def jarFile = file(tigerJarPath)
    def jarName = jarFile.name
    
    // Get configuration from version mapping
    def config = params.tiger_jar_versions[jarName]
    
    if (!config) {
        // Try to infer from filename patterns
        if (jarName.contains("F3") || jarName.contains("FastCall3")) {
            log.warn "Unknown TIGER jar: ${jarName}. Assuming FastCall3 configuration."
            config = [
                java_version: "java17",
                fastcall_version: "FastCall3", 
                app_name: "FastCall3"
            ]
        } else if (jarName.contains("2023") || jarName.contains("FastCall2")) {
            log.warn "Unknown TIGER jar: ${jarName}. Assuming FastCall2 configuration."
            config = [
                java_version: "java8",
                fastcall_version: "FastCall2",
                app_name: "FastCall2" 
            ]
        } else {
            log.warn "Unknown TIGER jar: ${jarName}. Using default FastCall3 configuration."
            config = [
                java_version: "java17",
                fastcall_version: "FastCall3",
                app_name: "FastCall3"
            ]
        }
    }
    
    return [
        path: tigerJarPath,
        name: jarName,
        java_version: config.java_version,
        fastcall_version: config.fastcall_version,
        app_name: config.app_name
    ]
}

// TIGER jar path resolution function
def resolveTigerJarPath(tigerJar, homeDir) {
    def candidatePaths = []
    
    // If tiger_jar is provided and is absolute path
    if (tigerJar && file(tigerJar).isAbsolute()) {
        candidatePaths << tigerJar
    }
    
    // If tiger_jar is provided as relative path or filename
    if (tigerJar && !file(tigerJar).isAbsolute()) {
        candidatePaths << "${homeDir}/${tigerJar}"
        candidatePaths << "${homeDir}/lib/${tigerJar}"
    }
    
    // Default locations based on jar name
    def jarName = tigerJar ?: params.tiger_jar
    if (homeDir) {
        candidatePaths << "${homeDir}/lib/${jarName}"
        candidatePaths << "${homeDir}/software/${jarName}"
        candidatePaths << "${homeDir}/tools/${jarName}"
        candidatePaths << "${homeDir}/bin/${jarName}"
        candidatePaths << "${homeDir}/${jarName}"
    }
    
    // System-wide locations
    candidatePaths << "/usr/local/lib/${jarName}"
    candidatePaths << "/opt/tiger/${jarName}"
    candidatePaths << "./${jarName}"
    
    // Find the first existing file
    for (path in candidatePaths) {
        if (file(path).exists()) {
            log.info "Found TIGER jar at: ${path}"
            return path
        }
    }
    
    // If not found, provide helpful error message
    def searchedPaths = candidatePaths.join("\n  - ")
    throw new Exception("""
TIGER jar file not found. Searched locations:
  - ${searchedPaths}

To fix this issue:
1. Specify the full path: --tiger_jar /path/to/your/TIGER.jar
2. Place the jar file in one of the standard locations
3. Check the filename is correct
""".stripIndent())
}

// WORKSHOP jar path resolution function
def resolveWorkshopJarPath(workshopJar, homeDir) {
    def candidatePaths = []

    if (workshopJar) {
        def wf = file(workshopJar)
        if (wf.isAbsolute()) {
            candidatePaths << workshopJar
        } else {
            // relative path provided, try relative to CWD and homeDir/lib
            candidatePaths << workshopJar
            if (homeDir) candidatePaths << "${homeDir}/lib/${workshopJar}"
        }
    }

    // Standard locations
    if (homeDir) candidatePaths << "${homeDir}/lib/Workshop.jar"
    if (homeDir) candidatePaths << "${homeDir}/software/Workshop.jar"
    if (homeDir) candidatePaths << "${homeDir}/tools/Workshop.jar"
    candidatePaths << "./Workshop.jar"
    candidatePaths << "/usr/local/lib/Workshop.jar"

    for (path in candidatePaths) {
        if (file(path).exists()) {
            log.info "Found Workshop jar at: ${path}"
            return path
        }
    }

    def searched = candidatePaths.unique().join("\n  - ")
    throw new Exception("""
Workshop jar file not found. Searched locations:
  - ${searched}

To fix this issue:
1. Specify the full path: --workshop_jar /path/to/Workshop.jar
2. Place the jar file under ${homeDir}/lib/Workshop.jar
3. Check the filename is correct
""".stripIndent())
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

// Function to read chromosomes from fai file
def getChromosomesFromFai(referencePath) {
    def faiPath = "${referencePath}.fai"
    def faiFile = file(faiPath)
    
    if (!faiFile.exists()) {
        log.warn "FAI file not found: ${faiPath}. Using default chromosome list."
        return (0..44).collect { it.toString() }
    }
    
    try {
        def chromosomes = []
        faiFile.eachLine { line ->
            if (line.trim()) {
                def parts = line.split('\t')
                if (parts.length >= 1) {
                    chromosomes << parts[0].trim()
                }
            }
        }
        
        if (chromosomes.isEmpty()) {
            log.warn "No chromosomes found in FAI file: ${faiPath}. Using default chromosome list."
            return (0..44).collect { it.toString() }
        }
        
        log.info "Found ${chromosomes.size()} chromosomes in FAI file: ${chromosomes.join(', ')}"
        return chromosomes
        
    } catch (Exception e) {
        log.warn "Error reading FAI file ${faiPath}: ${e.message}. Using default chromosome list."
        return (0..44).collect { it.toString() }
    }
}

def helpMessage() {
    log.info """
    ========================================
    FastCall3 Advanced Pipeline
    ========================================
    
    Usage:
        nextflow run runFastCall3.nf --home_dir <home> --pop <population> --job <job_name>
    
    Required parameters:
        --home_dir          Home directory containing data and libraries
        --pop               Population/Project name
                           Options: A, AB, ABD, D, HZNU, Nature, S, WAP, Watkins
        --job               Job name for output identification
    
    TIGER jar configuration:
        --tiger_jar    TIGER jar filename (default: TIGER_F3_20250911.jar)
                           Common options:
                           - TIGER_F3_20250911.jar (FastCall3 latest)
                           - TIGER_F3_20250910.jar (FastCall3 previous)  
                           - TIGER_20250526.jar (FastCall2)
                           - Tiger.jar (legacy)
        
        TIGER jar search locations (in order):
        1. Absolute path specified by --tiger_jar
        2. \${home_dir}/lib/\${jar_name}
        3. \${home_dir}/software/\${jar_name}
        4. \${home_dir}/tools/\${jar_name}
        5. \${home_dir}/bin/\${jar_name}
        6. \${home_dir}/\${jar_name}
        7. /usr/local/lib/\${jar_name}
        8. /opt/tiger/\${jar_name}
        9. ./\${jar_name}
        
    Other required tools:
        --workshop_jar      Path to Workshop jar file (default: \${home_dir}/lib/Workshop.jar)
    
    Optional parameters:
        --bam_dir           BAM files directory (auto-configured based on population)
        --depth_dir         Depth files directory (auto-configured based on population)
        --reference         Reference genome fasta file
        --samtools_path     Path to samtools executable
        --output_dir        Output directory (default: output)
        --ing_dir           Directory for disc outputs (.ing.gz). Default: <output_dir>/<job>/ing
        --vlib_dir          Directory for blib outputs (.lib.gz). Default: <output_dir>/<job>/vLib
        --gen_dir           Directory for scan outputs (VCF). Default: <output_dir>/<job>/gen
        --taxa_bam_map      Path to an existing taxaBamMap.txt; if provided, the pipeline will use it directly and skip generation
        --threads           Number of threads (default: 32)
        --memory            Memory allocation (default: 100g)
        --chr               Single chromosome to process (overrides --chromosomes)
        --chromosomes       List of chromosomes to process (default: 0-44)
        --java_lib          Java installation base directory (default: /data/dazheng/lib/jvm)
        
    Java version management:
        The pipeline automatically manages Java versions using the --java_lib parameter.
        Java versions are expected to be in subdirectories: jdk-8, jdk-11, jdk-17, jdk-21
        - prepareTaxaBamMap process uses Java 21 (jdk-21)
        - FastCall3 processes (disc, blib, scan) use Java 17 (jdk-17)
        
        TIGER jar version compatibility:
        - TIGER_F3_*.jar: Requires Java 17 or 21
        - TIGER_202*.jar: Requires Java 11 or 17  
        - Tiger.jar: Requires Java 8 or 11
        
        Example directory structure:
        /data/dazheng/lib/jvm/
        ├── jdk-8/
        ├── jdk-11/
        ├── jdk-17/
        └── jdk-21/
        
    Workflow control:
        --workflow_mode     Workflow execution mode (default: full)
                           Options: full, disc_only, blib_only, scan_only, from_disc, from_blib
        --skip_taxa_map     Skip taxa-BAM map generation (default: false)
        --resume_from_checkpoint Resume from existing intermediate files (default: false)
    
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
        
    Examples:
        # Basic usage with auto-detected TIGER jar
        nextflow run runFastCall3.nf --home_dir /data/project --pop A --job test_run
        
        # Specify custom TIGER jar
        nextflow run runFastCall3.nf --home_dir /data/project --pop A --job test_run \\
            --tiger_jar /path/to/custom/TIGER.jar
            
        # Use different TIGER jar version
        nextflow run runFastCall3.nf --home_dir /data/project --pop A --job test_run \\
            --tiger_jar_name TIGER_20250526.jar
            
        # Resume from disc results
        nextflow run runFastCall3.nf --home_dir /data/project --pop A --job test_run \\
            --workflow_mode from_disc
    """.stripIndent()
}

// Main workflow entry point
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
    if (!params.pop) {
        log.error "Population name is required. Use --pop"
        exit 1
    }
    if (!params.job) {
        log.error "Job name is required. Use --job"
        exit 1
    }

    // Set default paths based on population
    def popConfig = getPopulationConfig(params.pop, params.home_dir)
    def bam_dir = params.bam_dir ?: popConfig.bam_dir
    def depth_dir = params.depth_dir ?: popConfig.depth_dir
    
    // Resolve and configure TIGER jar automatically
    def tiger_jar_path = resolveTigerJarPath(params.tiger_jar, params.home_dir)
    def tiger_jar_config = getTigerJarConfig(tiger_jar_path)
    
    def workshop_jar_path = resolveWorkshopJarPath(params.workshop_jar, params.home_dir)
    def reference = params.reference ?: popConfig.reference
    def samtools_path = params.samtools_path ?: "samtools"

    def output_path = params.output_dir ?: "${params.home_dir}/output/${params.job}"

    // Resolve sub output dirs from params or derive sensible defaults (avoid reassigning params to prevent warnings)
    def ing_dir_resolved  = (params.ing_dir  ?: "${output_path}/ing").toString().replaceAll(/\/+/, "/")
    def vlib_dir_resolved = (params.vlib_dir ?: "${output_path}/vLib").toString().replaceAll(/\/+/, "/")
    def gen_dir_resolved  = (params.gen_dir  ?: "${output_path}/gen").toString().replaceAll(/\/+/, "/")

    // Validate essential paths
    def pathValidation = validatePaths([
        "Home directory": params.home_dir,
        "BAM directory": bam_dir,
        "Reference genome": reference,
        "TIGER jar": tiger_jar_config.path,
        "Workshop jar": workshop_jar_path
    ])
    
    if (!pathValidation.isValid) {
        log.error "Path validation failed:"
        pathValidation.errors.each { log.error "  - ${it}" }
        exit 1
    }

    // Pre-resolve reference sidecar files (if exist) so tools can reuse them without regenerating
    def ref_fai_path  = file("${reference}.fai").exists() ? "${reference}.fai" : ""
    def ref_gzi_path  = file("${reference}.gzi").exists() ? "${reference}.gzi" : ""

    log.info """\
    ========================================
    FastCall3 Advanced Pipeline
    ========================================
    Home directory   : ${params.home_dir}
    Population       : ${params.pop} (${popConfig.description})
    BAM directory    : ${bam_dir}
    Depth directory  : ${depth_dir}
    Job name         : ${params.job}
    Reference genome : ${reference}
    TIGER jar        : ${tiger_jar_config.path}
    TIGER jar name   : ${tiger_jar_config.name}
    FastCall version : ${tiger_jar_config.fastcall_version}
    Java version     : ${tiger_jar_config.java_version}
    App name         : ${tiger_jar_config.app_name}
    Workshop jar     : ${workshop_jar_path}
    Samtools path    : ${samtools_path}
    Output directory : ${output_path}
    ING directory    : ${ing_dir_resolved}
    vLib directory   : ${vlib_dir_resolved}
    GEN directory    : ${gen_dir_resolved}
    Java library dir : ${params.java_lib}
    Workflow mode    : ${params.workflow_mode}
    Skip taxa map    : ${params.skip_taxa_map}
    Resume checkpoint: ${params.resume_from_checkpoint}
    Threads          : ${params.threads}
    Memory           : ${params.memory}
    ========================================
    """.stripIndent()

    // Create output directories
    def output_dir = file(output_path)
    def pipeline_info_dir = file("./pipeline_info/${params.job}")
    output_dir.mkdirs()
    pipeline_info_dir.mkdirs()
    // Ensure sub-dirs exist
    file(ing_dir_resolved).mkdirs()
    file(vlib_dir_resolved).mkdirs()
    file(gen_dir_resolved).mkdirs()

    def taxaBamMap_dir = file("${params.home_dir}/00data/05taxaBamMap")

    // Step 1: Prepare or load taxa-BAM mapping
    def taxa_bam_map
    if (params.taxa_bam_map) {
        // User provided a mapping file; validate and use it directly
        def provided = file(params.taxa_bam_map)
        if (!provided.exists()) {
            log.error "Provided --taxa_bam_map does not exist: ${params.taxa_bam_map}"
            exit 1
        }
        taxa_bam_map = Channel.fromPath(provided)
        log.info "Using user-provided taxa-BAM mapping: ${provided}"
    } else if (!params.skip_taxa_map) {
        // Generate mapping
        prep_results = prepareTaxaBamMap(
            file(workshop_jar_path), 
            depth_dir,
            bam_dir,
            params.home_dir,
            params.job
        )
        taxa_bam_map = prep_results.taxa_bam_map
        log.info "Generated taxa-BAM mapping for job: ${params.job}"
    } else {
        // Use existing taxa-BAM mapping file
        def existing_taxa_map = "${params.output_dir}/${params.job}/taxa_bam_map/${params.job}.taxaBamMap.txt"
        def fallback_taxa_map = "${taxaBamMap_dir}/${params.job}.taxaBamMap.txt"
        if (file(existing_taxa_map).exists()) {
            taxa_bam_map = Channel.fromPath(existing_taxa_map)
            log.info "Using existing taxa-BAM mapping: ${existing_taxa_map}"
        } else if (file(fallback_taxa_map).exists()) {
            taxa_bam_map = Channel.fromPath(fallback_taxa_map)
            log.info "Using fallback taxa-BAM mapping: ${fallback_taxa_map}"
        } else {
            log.error "Taxa-BAM mapping file not found in either:"
            log.error "  - ${existing_taxa_map}"
            log.error "  - ${fallback_taxa_map}"
            log.error "Provide --taxa_bam_map, or set --skip_taxa_map false to generate it"
            exit 1
        }
    }
    
    // Determine chromosomes from fai file if default, otherwise use provided chromosomes.
    // If --chr is specified, it overrides --chromosomes and runs only that chromosome.
    def chromosomes
    if (params.chr) {
        chromosomes = [params.chr.toString()]
        log.info "Using single chromosome from --chr: ${chromosomes[0]}"
    } else {
        def chromParam = params.chromosomes
        if (chromParam == (0..44).collect { it.toString() }) {
            // Use default params.chromosomes, get from fai file
            chromosomes = getChromosomesFromFai(reference)
            log.info "Using chromosomes from fai file: ${chromosomes.size()} chromosomes found"
        } else if (chromParam instanceof String || chromParam instanceof GString) {
            // Single chromosome passed as string
            chromosomes = [chromParam.toString()]
            log.info "Using single chromosome: ${chromosomes[0]}"
        } else if (chromParam instanceof List) {
            // Multiple chromosomes provided
            chromosomes = chromParam.collect { it.toString() }
            log.info "Using provided chromosomes: ${chromosomes.join(', ')}"
        } else {
            // Fallback to fai file
            chromosomes = getChromosomesFromFai(reference)
            log.info "Fallback: using chromosomes from fai file: ${chromosomes.size()} chromosomes found"
        }
    }
    
    // Create chromosome channel
    chromosome_ch = Channel.fromList(chromosomes)
    
    log.info "Pipeline setup completed. Processing ${chromosomes.size()} chromosomes: ${chromosomes.size() > 0 ? chromosomes[0] + '...' + chromosomes[-1] : 'none'}"
    
    // Generate population statistics
    pop_stats = generate_population_stats(
        bam_dir,
        depth_dir,
        params.pop
    )
    
    // Workflow execution based on mode
    def workflow_mode = params.workflow_mode ?: "full"
    log.info "Running workflow in mode: ${workflow_mode}"
    
    switch (workflow_mode) {
        case "full":
            // Complete workflow: disc -> blib -> scan -> collect
            disc_results = fastcall3_disc(
                chromosome_ch, 
                file(reference),
                ref_fai_path,
                ref_gzi_path,
                taxa_bam_map,
                file(tiger_jar_config.path),
                ing_dir_resolved,
                samtools_path,
                tiger_jar_config
            )
            
            blib_results = fastcall3_blib(
                    chromosome_ch,
                file(reference),
                ref_fai_path,
                ref_gzi_path,
                file(tiger_jar_config.path),
                tiger_jar_config,
                ing_dir_resolved,
                vlib_dir_resolved
            )
            
            scan_results = fastcall3_scan(
                blib_results.blib_files,
                file(reference),
                ref_fai_path,
                ref_gzi_path,
                taxa_bam_map,
                file(tiger_jar_config.path),
                samtools_path,
                tiger_jar_config,
                gen_dir_resolved
            )
            
            collect_results(
                scan_results.vcf_files.map{ chr, f -> f }.collect(),
                chromosomes,
                tiger_jar_config.name,
                tiger_jar_config.fastcall_version,
                tiger_jar_config.java_version
            )
            break
            
        case "disc_only":
            // Only run disc analysis
            disc_results = fastcall3_disc(
                chromosome_ch, 
                file(reference),
                ref_fai_path,
                ref_gzi_path,
                taxa_bam_map,
                file(tiger_jar_config.path),
                params.ing_dir,
                samtools_path,
                tiger_jar_config
            )
            log.info "Disc analysis completed. Use 'from_disc' mode to continue."
            break
            
        case "from_disc":
            // Start from existing disc (ing) results
            def disc_dir = ing_dir_resolved
            if (!file(disc_dir).exists()) {
                log.error "Disc directory not found: ${disc_dir}"
                log.error "Run workflow with 'disc_only' or 'full' mode first"
                exit 1
            }
            
            disc_files_ch = Channel.fromPath("${disc_dir}/**/*.ing.gz")
                .map { f ->
                    def fname = f.name
                    def m = (fname =~ /^([^_]+)_.*\.ing\.gz$/)
                    def chr = m.find() ? m.group(1) : f.parent.name.tokenize('_')[0]
                    tuple(chr, f)
                }
                .groupTuple()
                .filter { chr, file -> params.chr ? chr == params.chr.toString() : true }
            if (params.chr) {
                log.info "Filtering disc inputs to chromosome: ${params.chr}"
            }
            
            // Check if we have disc files (recursively search subfolders like A001/A002)
            def disc_count = 0
            def disc_root = new File(disc_dir as String)
            if (disc_root.exists()) {
                disc_root.eachFileRecurse(groovy.io.FileType.FILES) { f ->
                    if (f.name.endsWith('.ing.gz')) disc_count++
                }
            }
            if (disc_count == 0) {
                log.error "No .ing.gz files found in ${disc_dir}"
                exit 1
            }
            log.info "Found ${disc_count} disc files to process"
            
            // Reduce to chromosome-only channel to avoid staging all .ing.gz
            def chr_only_from_disc = disc_files_ch.map { chr, files -> chr }.distinct()
            blib_results = fastcall3_blib(
                chr_only_from_disc,
                file(reference),
                ref_fai_path,
                ref_gzi_path,
                file(tiger_jar_config.path),
                tiger_jar_config,
                ing_dir_resolved,
                vlib_dir_resolved
            )
            
            scan_results = fastcall3_scan(
                blib_results.blib_files,
                file(reference),
                ref_fai_path,
                ref_gzi_path,
                taxa_bam_map,
                file(tiger_jar_config.path),
                samtools_path,
                tiger_jar_config,
                params.gen_dir
            )
            
            collect_results(
                scan_results.vcf_files.map{ chr, f -> f }.collect(),
                chromosomes,
                tiger_jar_config.name,
                tiger_jar_config.fastcall_version,
                tiger_jar_config.java_version
            )
            break
            
        case "blib_only":
            // Only run blib generation from existing disc (ing) results
            def disc_dir = params.ing_dir
            if (!file(disc_dir).exists()) {
                log.error "Disc directory not found: ${disc_dir}"
                exit 1
            }
            
            disc_files_ch = Channel.fromPath("${disc_dir}/**/*.ing.gz")
                .map { f ->
                    def fname = f.name
                    def m = (fname =~ /^([^_]+)_.*\.ing\.gz$/)
                    def chr = m.find() ? m.group(1) : f.parent.name.tokenize('_')[0]
                    tuple(chr, f)
                }
                .groupTuple()
                .filter { chr, file -> params.chr ? chr == params.chr.toString() : true }
            if (params.chr) {
                log.info "Filtering disc inputs to chromosome: ${params.chr}"
            }
            
            // Reduce to chromosome-only channel to avoid staging all .ing.gz
            def chr_only_blib = disc_files_ch.map { chr, files -> chr }.distinct()
            blib_results = fastcall3_blib(
                chr_only_blib,
                file(reference),
                ref_fai_path,
                ref_gzi_path,
                file(tiger_jar_config.path),
                tiger_jar_config,
                params.ing_dir,
                params.vlib_dir
            )
            log.info "Blib generation completed. Use 'from_blib' mode to continue."
            break
            
        case "from_blib":
            // Start from existing blib (vLib) results
            def blib_dir = vlib_dir_resolved
            if (!file(blib_dir).exists()) {
                log.error "Blib directory not found: ${blib_dir}"
                log.error "Run workflow with 'blib_only', 'from_disc', or 'full' mode first"
                exit 1
            }
            
            blib_files_ch = Channel.fromPath("${blib_dir}/*.lib.gz")
                .map { f ->
                    def fname = f.name
                    def m = (fname =~ /^([^_]+)_.*\.lib\.gz$/)
                    def chr = m.find() ? m.group(1) : fname.tokenize('.')[0]
                    tuple(chr, f)
                }
                .filter { chr, file -> params.chr ? chr == params.chr.toString() : true }
            if (params.chr) {
                log.info "Filtering blib inputs to chromosome: ${params.chr}"
            }
            
            // Check if we have blib files
            blib_count = file(blib_dir).exists() ? file(blib_dir).list().findAll { it.endsWith('.lib.gz') }.size() : 0
            if (blib_count == 0) {
                log.error "No .lib.gz files found in ${blib_dir}"
                exit 1
            }
            log.info "Found ${blib_count} blib files to process"
            
            scan_results = fastcall3_scan(
                blib_files_ch,
                file(reference),
                ref_fai_path,
                ref_gzi_path,
                taxa_bam_map,
                file(tiger_jar_config.path),
                samtools_path,
                tiger_jar_config,
                gen_dir_resolved
            )
            
            collect_results(
                scan_results.vcf_files.map{ chr, f -> f }.collect(),
                chromosomes,
                tiger_jar_config.name,
                tiger_jar_config.fastcall_version,
                tiger_jar_config.java_version
            )
            break
            
        case "scan_only":
            // Only run scan from existing blib (vLib) results
            def blib_dir = vlib_dir_resolved
            if (!file(blib_dir).exists()) {
                log.error "Blib directory not found: ${blib_dir}"
                exit 1
            }
            
            blib_files_ch = Channel.fromPath("${blib_dir}/*.lib.gz")
                .map { f ->
                    def fname = f.name
                    def m = (fname =~ /^([^_]+)_.*\.lib\.gz$/)
                    def chr = m.find() ? m.group(1) : fname.tokenize('.')[0]
                    tuple(chr, f)
                }
                .filter { chr, file -> params.chr ? chr == params.chr.toString() : true }
            if (params.chr) {
                log.info "Filtering blib inputs to chromosome: ${params.chr}"
            }
            
            scan_results = fastcall3_scan(
                blib_files_ch,
                file(reference),
                ref_fai_path,
                ref_gzi_path,
                taxa_bam_map,
                file(tiger_jar_config.path),
                samtools_path,
                tiger_jar_config,
                gen_dir_resolved
            )
            
            collect_results(
                scan_results.vcf_files.map{ chr, f -> f }.collect(),
                chromosomes,
                tiger_jar_config.name,
                tiger_jar_config.fastcall_version,
                tiger_jar_config.java_version
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
    tag "prepare_taxa_bam_map_${params.pop}"
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
    path "${job}.taxaRunMap.txt", emit: taxa_run_map
    path "prepare_taxa_bam_map_${params.pop}.log", emit: log
    
    script:
    def output_bam_file = "${job}.taxaBamMap.txt"
    def output_run_file = "${job}.taxaRunMap.txt"
    def javaSetup = getJavaSetupScript("java21", params.java_lib)
    """
    echo "Preparing taxa-BAM mapping file for population ${params.pop}..." > prepare_taxa_bam_map_${params.pop}.log
    echo "Using Java 21 for prepareTaxaBamMap process" >> prepare_taxa_bam_map_${params.pop}.log
    
    ${javaSetup} >> prepare_taxa_bam_map_${params.pop}.log 2>&1
    
    # Use TaxaBamMap.java only
    echo "Running TaxaBamMap from Workshop jar..." >> prepare_taxa_bam_map_${params.pop}.log
    
    java -Xmx16g -jar "${workshop_jar}" \\
        -d ${depth_dir} \\
        -b ${bam_dir} \\
        -o ${output_bam_file} \\
        -t ${output_run_file} \\
        >> prepare_taxa_bam_map_${params.pop}.log 2>&1
    
    echo "TaxaBamMap execution completed" >> prepare_taxa_bam_map_${params.pop}.log
    
    # Copy to standard location for reuse
    mkdir -p ${home_dir}/00data/05taxaBamMap
    cp ${output_bam_file} ${home_dir}/00data/05taxaBamMap/
    cp ${output_run_file} ${home_dir}/00data/05taxaBamMap/ 2>/dev/null || true
    """
}

process fastcall3_disc {
    tag "disc_${chromosome}"
    memory params.memory
    cpus params.threads
    errorStrategy 'retry'
    maxRetries 2
    publishDir({ (params.ing_dir ?: "${params.output_dir}/${params.job}/ing").toString().replaceAll(/\/+/, "/") }, mode: 'copy', pattern: "*.ing.gz")
    
    input:
    val chromosome
    path reference
    val ref_fai_path
    val ref_gzi_path
    path taxa_bam_map
    path tiger_jar
    val ing_path
    val samtools_path
    val tiger_jar_config
    
    output:
    tuple val(chromosome), path("*.ing.gz"), emit: disc_files
    path "disc_${chromosome}.log", emit: log
    
    script:
    def javaSetup = getJavaSetupScript(tiger_jar_config.java_version, params.java_lib)
    """
    echo "Starting disc analysis for chromosome ${chromosome}..." > disc_${chromosome}.log
    echo "Using ${tiger_jar_config.java_version} for ${tiger_jar_config.fastcall_version} disc process" >> disc_${chromosome}.log
    
    ${javaSetup} >> disc_${chromosome}.log 2>&1

    # Link reference sidecar files if provided to avoid regenerating indices
    ref_base=\$(basename "${reference}")
    # .fai
    if [ -n "${ref_fai_path}" ] && [ -f "${ref_fai_path}" ]; then
        ln -sf "${ref_fai_path}" "\${ref_base}.fai" 2>> disc_${chromosome}.log || true
        echo "Linked FAI: \${ref_base}.fai -> ${ref_fai_path}" >> disc_${chromosome}.log
    fi
    # .gzi for gzipped fasta (POSIX compatible)
    case "\${ref_base}" in
      *.gz)
        if [ -n "${ref_gzi_path}" ] && [ -f "${ref_gzi_path}" ]; then
            ln -sf "${ref_gzi_path}" "\${ref_base}.gzi" 2>> disc_${chromosome}.log || true
            echo "Linked GZI: \${ref_base}.gzi -> ${ref_gzi_path}" >> disc_${chromosome}.log
        fi
        ;;
    esac
    # no .dict required
    
    # Monitor system resources
    echo "System resources before TIGER execution:" >> disc_${chromosome}.log
    echo "Memory: \$(free -h | grep Mem)" >> disc_${chromosome}.log
    echo "CPU cores (allocated): ${task.cpus}" >> disc_${chromosome}.log
    echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g" >> disc_${chromosome}.log
    echo "TIGER jar: ${tiger_jar_config.name}" >> disc_${chromosome}.log
    echo "FastCall version: ${tiger_jar_config.fastcall_version}" >> disc_${chromosome}.log
    
    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${tiger_jar_config.app_name} \\
        -mod disc \\
        -a ${reference} \\
        -b ${taxa_bam_map} \\
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
        -m ${task.cpus} \\
        -n ${ing_path} \\
        -o ${samtools_path} \\
        >> disc_${chromosome}.log 2>&1
    
    echo "Disc analysis for chromosome ${chromosome} completed" >> disc_${chromosome}.log

    # Symlink produced .ing.gz back into workdir so Nextflow can capture as output
    if ls ${ing_path}/${chromosome}_*.ing.gz >/dev/null 2>&1; then
        ln -sf ${ing_path}/${chromosome}_*.ing.gz . 2>> disc_${chromosome}.log || true
    else
        # fallback: pull any .ing.gz for this run
        find ${ing_path} -maxdepth 2 -type f -name "${chromosome}*.ing.gz" -exec ln -sf {} . \\; 2>> disc_${chromosome}.log || true
    fi
    """
}

process fastcall3_blib {
    tag "blib_${chromosome}"
    memory params.memory
    cpus params.threads
    errorStrategy 'retry'
    maxRetries 2
    publishDir({ (params.vlib_dir ?: "${params.output_dir}/${params.job}/vLib").toString().replaceAll(/\/+/, "/") }, mode: 'copy', pattern: "*.lib.gz")
    
    input:
    val chromosome
    path reference
    val ref_fai_path
    val ref_gzi_path
    path tiger_jar
    val tiger_jar_config
    val ing_path
    val vLib_path
    
    output:
    tuple val(chromosome), path("*.lib.gz"), emit: blib_files
    path "blib_${chromosome}.log", emit: log
    
    script:
    def javaSetup = getJavaSetupScript(tiger_jar_config.java_version, params.java_lib)
    """
    echo "Starting blib generation for chromosome ${chromosome}..." > blib_${chromosome}.log
    echo "Using ${tiger_jar_config.java_version} for ${tiger_jar_config.fastcall_version} blib process" >> blib_${chromosome}.log
    
    ${javaSetup} >> blib_${chromosome}.log 2>&1

    # Link reference sidecar files if provided to avoid regenerating indices
    ref_base=\$(basename "${reference}")
    if [ -n "${ref_fai_path}" ] && [ -f "${ref_fai_path}" ]; then
        ln -sf "${ref_fai_path}" "\${ref_base}.fai" 2>> blib_${chromosome}.log || true
        echo "Linked FAI: \${ref_base}.fai -> ${ref_fai_path}" >> blib_${chromosome}.log
    fi
    case "\${ref_base}" in
      *.gz)
        if [ -n "${ref_gzi_path}" ] && [ -f "${ref_gzi_path}" ]; then
            ln -sf "${ref_gzi_path}" "\${ref_base}.gzi" 2>> blib_${chromosome}.log || true
            echo "Linked GZI: \${ref_base}.gzi -> ${ref_gzi_path}" >> blib_${chromosome}.log
        fi
        ;;
    esac
    # no .dict required
    
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
    echo "TIGER jar: ${tiger_jar_config.name}" >> blib_${chromosome}.log
    echo "FastCall version: ${tiger_jar_config.fastcall_version}" >> blib_${chromosome}.log
    
    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${tiger_jar_config.app_name} \\
        -mod blib \\
        -a ${reference} \\
        -b ${chromosome} \\
        -c 2 \\
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
    memory params.memory
    cpus params.threads
    errorStrategy 'retry'
    maxRetries 2
    // Do not publish/copy outputs here; TIGER already writes VCFs under gen_path/VCF
    
    input:
    tuple val(chromosome), path(blib_files)
    path reference
    val ref_fai_path
    val ref_gzi_path
    path taxa_bam_map
    path tiger_jar
    val samtools_path
    val tiger_jar_config
    val gen_path
    
    output:
    tuple val(chromosome), path("*.vcf*"), emit: vcf_files
    path "scan_${chromosome}.log", emit: log
    
    script:
    def javaSetup = getJavaSetupScript(tiger_jar_config.java_version, params.java_lib)
    """
    echo "Starting scan analysis for chromosome ${chromosome}..." > scan_${chromosome}.log
    echo "Using ${tiger_jar_config.java_version} for ${tiger_jar_config.fastcall_version} scan process" >> scan_${chromosome}.log
    
    ${javaSetup} >> scan_${chromosome}.log 2>&1

    # Link reference sidecar files if provided to avoid regenerating indices
    ref_base=\$(basename "${reference}")
    if [ -n "${ref_fai_path}" ] && [ -f "${ref_fai_path}" ]; then
        ln -sf "${ref_fai_path}" "\${ref_base}.fai" 2>> scan_${chromosome}.log || true
        echo "Linked FAI: \${ref_base}.fai -> ${ref_fai_path}" >> scan_${chromosome}.log
    fi
    case "\${ref_base}" in
      *.gz)
        if [ -n "${ref_gzi_path}" ] && [ -f "${ref_gzi_path}" ]; then
            ln -sf "${ref_gzi_path}" "\${ref_base}.gzi" 2>> scan_${chromosome}.log || true
            echo "Linked GZI: \${ref_base}.gzi -> ${ref_gzi_path}" >> scan_${chromosome}.log
        fi
        ;;
    esac
    # no .dict required
    
    # Resolve the lib file (staged by Nextflow) for this chromosome
    lib_file_path="${blib_files}"
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
    echo "TIGER jar: ${tiger_jar_config.name}" >> scan_${chromosome}.log
    echo "FastCall version: ${tiger_jar_config.fastcall_version}" >> scan_${chromosome}.log
    
    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${tiger_jar_config.app_name} \\
        -mod scan \\
        -a ${reference} \\
        -b ${taxa_bam_map} \\
        -c "\$lib_file_path" \\
        -d ${chromosome} \\
        -e 0 \\
        -f 30 \\
        -g 20 \\
        -h 0.05 \\
        -i ${samtools_path} \\
        -j ${task.cpus} \\
        -k ${gen_path} \\
        >> scan_${chromosome}.log 2>&1
    
    echo "Scan analysis for chromosome ${chromosome} completed" >> scan_${chromosome}.log

    # Symlink produced VCFs back into workdir for Nextflow capture only (no publishing)
    found_vcf=0
    if find ${gen_path} -type f -name "${chromosome}*.vcf*" | grep -q . ; then
        find ${gen_path} -type f -name "${chromosome}*.vcf*" -exec ln -sf {} . \\; 2>> scan_${chromosome}.log || true
        found_vcf=1
    fi
    if [ "\$found_vcf" -eq 0 ]; then
        # fallback: any VCFs
        find ${gen_path} -type f -name "*.vcf*" -exec ln -sf {} . \\; 2>> scan_${chromosome}.log || true
    fi
    """
}

process collect_results {
    tag "collect_results"
    memory '16g'
    cpus 4
    publishDir "${params.output_dir}/${params.job}/final", mode: 'move'
    
    input:
    path(vcf_files)
    val(chromosomes_list)
    val(tiger_jar_name)
    val(tiger_fastcall_version)
    val(tiger_java_version)
    
    output:
    path "merged_variants.vcf.gz", optional: true
    path "summary_stats.txt"
    path "pipeline_report.html", optional: true
    
    script:
    """
    echo "Collecting and processing results..." > collect_results.log
    
    # Create output directory if it doesn't exist
    mkdir -p final_output
    
    # Check if VCF files exist
    vcf_count=\$(ls -1 *.vcf *.vcf.gz 2>/dev/null | wc -l)
    if [ \$vcf_count -eq 0 ]; then
        echo "Warning: No VCF files found to merge" >> collect_results.log
        touch merged_variants.vcf.gz
    else
        # Merge all VCF files
        if command -v bcftools &> /dev/null; then
            echo "Merging \$vcf_count VCF files using bcftools..." >> collect_results.log
            bcftools concat ${vcf_files} | bcftools sort -Oz -o merged_variants.vcf.gz
            bcftools index merged_variants.vcf.gz
            echo "VCF files merged successfully" >> collect_results.log
        else
            echo "Warning: bcftools not found. VCF files are available separately in scan directory." >> collect_results.log
            touch merged_variants.vcf.gz
        fi
    fi
    
    # Generate summary statistics
    echo "FastCall3 Pipeline Summary" > summary_stats.txt
    echo "=========================" >> summary_stats.txt
    echo "Date: \$(date)" >> summary_stats.txt
    echo "Population: ${params.pop}" >> summary_stats.txt
    echo "Job: ${params.job}" >> summary_stats.txt
    echo "Number of chromosomes processed: ${chromosomes_list.size()}" >> summary_stats.txt
    echo "Chromosomes: ${chromosomes_list.join(', ')}" >> summary_stats.txt
    echo "Home directory: ${params.home_dir}" >> summary_stats.txt
    echo "" >> summary_stats.txt
    echo "TIGER Configuration:" >> summary_stats.txt
    echo "- TIGER jar: ${tiger_jar_name}" >> summary_stats.txt
    echo "- FastCall version: ${tiger_fastcall_version}" >> summary_stats.txt
    echo "- Java version: ${tiger_java_version}" >> summary_stats.txt
    echo "" >> summary_stats.txt
    echo "Parameters used:" >> summary_stats.txt
    echo "- Min depth: ${params.disc_min_depth}" >> summary_stats.txt
    echo "- Min quality: ${params.disc_min_qual}" >> summary_stats.txt
    echo "- Min allele frequency: ${params.disc_min_allele_freq}" >> summary_stats.txt
    echo "- P-value threshold: ${params.scan_p_value}" >> summary_stats.txt
    echo "- Threads: ${params.threads}" >> summary_stats.txt
    echo "- Memory: ${params.memory}" >> summary_stats.txt
    echo "" >> summary_stats.txt
    echo "VCF files processed: \$vcf_count" >> summary_stats.txt
    echo "" >> summary_stats.txt
    echo "Output files:" >> summary_stats.txt
    ls -la *.vcf* >> summary_stats.txt 2>/dev/null || echo "No VCF files found" >> summary_stats.txt
    
    # Generate simple HTML report if possible
    if command -v python3 &> /dev/null; then
        python3 -c "
import datetime
html = '''
<!DOCTYPE html>
<html>
<head><title>FastCall3 Pipeline Report</title></head>
<body>
<h1>FastCall3 Pipeline Report</h1>
<p><strong>Date:</strong> ''' + str(datetime.datetime.now()) + '''</p>
<p><strong>Population:</strong> ${params.pop}</p>
<p><strong>Job:</strong> ${params.job}</p>
<p><strong>Chromosomes processed:</strong> ${chromosomes_list.size()}</p>
<p><strong>VCF files generated:</strong> \$vcf_count</p>
</body>
</html>
'''
with open('pipeline_report.html', 'w') as f:
    f.write(html)
        " || touch pipeline_report.html
    else
        touch pipeline_report.html
    fi
    
    echo "Results collection completed" >> collect_results.log
    """
}

process generate_population_stats {
    tag "stats_${params.pop}"
    memory '8g'
    cpus 2
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'move'
    
    input:
    val bam_dir
    val depth_dir
    val pop
    
    output:
    path "population_stats_${pop}.txt"
    path "bam_file_list_${pop}.txt"
    
    script:
    """
    echo "Population Statistics for ${pop}" > population_stats_${pop}.txt
    echo "=================================" >> population_stats_${pop}.txt
    echo "Date: \$(date)" >> population_stats_${pop}.txt
    echo "BAM directory: ${bam_dir}" >> population_stats_${pop}.txt
    echo "Depth directory: ${depth_dir}" >> population_stats_${pop}.txt
    echo "" >> population_stats_${pop}.txt
    
    # Count BAM files
    bam_count=\$(find ${bam_dir} -name "*.bam" 2>/dev/null | wc -l)
    echo "Number of BAM files: \$bam_count" >> population_stats_${pop}.txt
    
    # List BAM files
    echo "BAM Files:" > bam_file_list_${pop}.txt
    find ${bam_dir} -name "*.bam" 2>/dev/null | sort >> bam_file_list_${pop}.txt
    
    # BAM file sizes
    echo "" >> population_stats_${pop}.txt
    echo "BAM file sizes:" >> population_stats_${pop}.txt
    find ${bam_dir} -name "*.bam" -exec ls -lh {} \\; 2>/dev/null | awk '{print \$9 " : " \$5}' | sort >> population_stats_${pop}.txt
    
    # Check depth files if directory exists
    if [ -d "${depth_dir}" ]; then
        depth_count=\$(find ${depth_dir} -name "*" -type f 2>/dev/null | wc -l)
        echo "" >> population_stats_${pop}.txt
        echo "Number of depth files: \$depth_count" >> population_stats_${pop}.txt
    else
        echo "" >> population_stats_${pop}.txt
        echo "Depth directory not found: ${depth_dir}" >> population_stats_${pop}.txt
    fi
    
    echo "" >> population_stats_${pop}.txt
    echo "Statistics generation completed at: \$(date)" >> population_stats_${pop}.txt
    """
}
