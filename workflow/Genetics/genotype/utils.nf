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
        ],
        "test_plink": [
            vcf_dir: "/data1/dazheng_tusr1/vmap4.VCF.v1"
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
        // Convert to integer to strip leading zeros (001->1), then to string for type consistency
        chromosome = matcher[0][1].toInteger().toString()
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

def getTigerJarConfig(tiger_jar_name, home_dir) {
    def tiger_jar_versions = [
        "TIGER_F3_2M_scan_20260129.jar": [
            java_version: "java17",
            app_name: "FastCall3"
        ],
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
            java_version: "java17", 
            app_name: "PopDep"
        ],
        "TIGER_PD_20260129.jar": [
            java_version: "java17", 
            app_name: "PopDep"
        ]
    ]

    def config = tiger_jar_versions[tiger_jar_name]
    def tigerJar = file("${home_dir}/lib/${tiger_jar_name}")
    
    return [
        path: tigerJar,
        app_name: config.app_name,
        java_version: config.java_version
    ]
}

def getPopDepTaxaBamFile_v1(chr, home_dir) {
    def normalized = chr.toString().replaceFirst(/^chr/, '')

    def subGenomeConfigs = [
        "A": "${home_dir}/00data/05taxaBamMap/vmap4_v1/tb.A.txt",
        "B": "${home_dir}/00data/05taxaBamMap/vmap4_v1/tb.B.txt",
        "D": "${home_dir}/00data/05taxaBamMap/vmap4_v1/tb.D.txt",
        "ALL": "${home_dir}/00data/05taxaBamMap/vmap4_v1/tb.ALL.txt"
    ]

    def groupA = getRefV1SubChr("A")
    def groupB = getRefV1SubChr("B")
    def groupD = getRefV1SubChr("D")
    def groupOthers = getRefV1SubChr("Others")

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
    def tb_file = subGenomeConfigs[sub_genome]
    return tb_file
}

def getRefV1SubChr(sub_genome) {
    def subGenomeChrMap = [
        "A": ["1","2","7","8","13","14","19","20","25","26","31","32","37","38"],
        "B": ["3","4","9","10","15","16","21","22","27","28","33","34","39","40"],
        "D": ["5","6","11","12","17","18","23","24","29","30","35","36","41","42"],
        "Others": ["0", "43", "44"],
        "ALL": [
            "0","1","2","3","4","5","6","7","8","9","10",
            "11","12","13","14","15","16","17","18","19","20",
            "21","22","23","24","25","26","27","28","29","30",
            "31","32","33","34","35","36","37","38","39","40",
            "41","42","43","44"
            ]
    ]
    return subGenomeChrMap[sub_genome]
}

def getRefV1ChrLength(chr) {
    def chrLengths = [
        0: 480980714,
        1: 471304005,
        2: 122798051,
        3: 438720154,
        4: 251131716,
        5: 452179604,
        6: 43273582,
        7: 462376173,
        8: 318422384,
        9: 453218924,
        10: 348037791,
        11: 462216879,
        12: 189635730,
        13: 454103970,
        14: 296739669,
        15: 448155269,
        16: 382674495,
        17: 476235359,
        18: 139317064,
        19: 452555092,
        20: 292033065,
        21: 451014251,
        22: 222603248,
        23: 451004620,
        24: 58852447,
        25: 453230519,
        26: 256543224,
        27: 451372872,
        28: 261776885,
        29: 451901030,
        30: 114179647,
        31: 452440856,
        32: 165638404,
        33: 452077197,
        34: 268911281,
        35: 450509124,
        36: 23083594,
        37: 450046986,
        38: 286659250,
        39: 453822637,
        40: 296797748,
        41: 453812268,
        42: 184873787,
        43: 452528,
        44: 134545
    ]
    return chrLengths[chr]
}
