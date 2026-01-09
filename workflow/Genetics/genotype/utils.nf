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
