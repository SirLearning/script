nextflow.enable.dsl=2

// --- Utility functions for vcf calling workflows ---

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
    
	def libDir = javaLibDir ?: "${params.user_dir}/lib/jvm"
    def javaHome = "${libDir}/${javaDir}"
    
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

// Get job-specific configuration
def getJobConfig(job, home_dir) {
    def jobConfigs = [
        "chr1": [
            vcf_dir: "${home_dir}/00data/06vcf/01chr1"
        ],
        "test": [
            vcf_dir: "${home_dir}/00data/06vcf/02test"
        ],
        "test_first": [
            vcf_file: "${home_dir}/00data/06vcf/02test/chr001.f1M.vcf"
        ],
        "test_middle": [
            vcf_file: "${home_dir}/00data/06vcf/02test/chr001.m1M.vcf"
        ],
        "test_last": [
            vcf_file: "${home_dir}/00data/06vcf/02test/chr001.l1M.vcf"
        ],
        "vmap4": [
            vcf_dir: "${home_dir}/00data/06vcf/03vmap4"  // pending
        ]
    ]

    if (!jobConfigs.containsKey(job)) {
        log.error "Unknown job specified: ${job}"
        System.exit(1)
    }

    return jobConfigs[job]
}

def getVcfIdFromPath(vcf_path) {
    def file_name = vcf_path.getName()
    // The id is the chromosome number with modifications, like vcf for chromosome 2 is chr001.vcf, for chromosome 10 is chr010.vcf
    // Here we remove the .vcf or .vcf.gz extension to get the id
    def tmp_id = file_name.replaceAll(/\.vcf$|\.vcf\.gz$/, '')
    def chromosome = null
    // Check if it matches pattern chrNNN (e.g., chr001, chr012)
    def matcher = (tmp_id =~ /^chr(\d{3})$/)
    if (matcher.matches()) {
        chromosome = matcher[0][1].toInteger()
        return [tmp_id, chromosome]
    } else {
        log.error "VCF file name does not match expected pattern 'chrNNN.vcf' or 'chrNNN.vcf.gz': ${file_name}"
        System.exit(1)
    }
}

// used for FastCall taxaBamMap
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
        ],
        "TIGER_PD_20260112.jar": [
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
