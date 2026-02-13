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

def getSoftwareConfig(home_dir, user_dir, tiger_jar, workshop_jar, samtools) {

    def java_lib_dir = "${user_dir}/lib/jvm"

    def resolved_workshop_jar = "${home_dir}/lib/${workshop_jar}"

    def tiger_jar_config = getTigerJarConfig(tiger_jar, home_dir)
    def workshop_jar_config = getWorkshopJarConfig(resolved_workshop_jar, java_lib_dir)
    def samtools_config = getSamtoolsConfig(samtools)
    
    return [
        tiger_jar_config: tiger_jar_config,
        workshop_jar_config: workshop_jar_config,
        samtools_config: samtools_config
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
        ],
        "rebuild": [
            vcf_file: "/data/home/tusr1/01projects/vmap4/04runScreens/rebuild/gen/VCF/chr002.vcf.gz"
        ]
    ]

    if (!jobConfigs.containsKey(job)) {
        log.error "Unknown job specified: ${job}"
        System.exit(1)
    }

    return jobConfigs[job]
}

def getCallingJobConfig(job, home_dir) {
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
        ],
        "rebuild": [
            name: "rebuild",
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
        "TIGER_PD_20260130.jar": [
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

def getTaxaBamMapFile_v1(chr, home_dir) {
    def normalized = chr.toString().replaceFirst(/^chr/, '')

    def subGenomeConfigs = [
        "A": "${home_dir}/00data/05taxaBamMap/all.A.taxaBamMap.txt",
        "B": "${home_dir}/00data/05taxaBamMap/all.B.taxaBamMap.txt",
        "D": "${home_dir}/00data/05taxaBamMap/all.D.taxaBamMap.txt",
        "ALL": "${home_dir}/00data/05taxaBamMap/all.ALL.taxaBamMap.txt"
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
    def tbm_file = subGenomeConfigs[sub_genome]
    return tbm_file
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

def getRefV1NameChr(name) {
    def nameChrMap = [
        "chrUn": "0",
        "chr1A": ["1", "2"],
        "chr1B": ["3", "4"],
        "chr1D": ["5", "6"],
        "chr2A": ["7", "8"],
        "chr2B": ["9", "10"],
        "chr2D": ["11", "12"],
        "chr3A": ["13", "14"],
        "chr3B": ["15", "16"],
        "chr3D": ["17", "18"],
        "chr4A": ["19", "20"],
        "chr4B": ["21", "22"],
        "chr4D": ["23", "24"],
        "chr5A": ["25", "26"],
        "chr5B": ["27", "28"],
        "chr5D": ["29", "30"],
        "chr6A": ["31", "32"],
        "chr6B": ["33", "34"],
        "chr6D": ["35", "36"],
        "chr7A": ["37", "38"],
        "chr7B": ["39", "40"],
        "chr7D": ["41", "42"],
        "Mit": "43",
        "Chl": "44"
    ]
    return nameChrMap[name]
}

def getRefV1ChrName(chr) {
    def chrNames = [
        "0": "chrUn",
        "1": "chr1A",
        "2": "chr1A",
        "3": "chr1B",
        "4": "chr1B",
        "5": "chr1D",
        "6": "chr1D",
        "7": "chr2A",
        "8": "chr2A",
        "9": "chr2B",
        "10": "chr2B",
        "11": "chr2D",
        "12": "chr2D",
        "13": "chr3A",
        "14": "chr3A",
        "15": "chr3B",
        "16": "chr3B",
        "17": "chr3D",
        "18": "chr3D",
        "19": "chr4A",
        "20": "chr4A",
        "21": "chr4B",
        "22": "chr4B",
        "23": "chr4D",
        "24": "chr4D",
        "25": "chr5A",
        "26": "chr5A",
        "27": "chr5B",
        "28": "chr5B",
        "29": "chr5D",
        "30": "chr5D",
        "31": "chr6A",
        "32": "chr6A",
        "33": "chr6B",
        "34": "chr6B",
        "35": "chr6D",
        "36": "chr6D",
        "37": "chr7A",
        "38": "chr7A",
        "39": "chr7B",
        "40": "chr7B",
        "41": "chr7D",
        "42": "chr7D",
        "43": "Mit",
        "44": "Chl"
    ]
    return chrNames[chr]
}

def getRefV1ChrOffset(chr) {
    def chrOffsets = [
        "0": "0",
        "1": "0",
        "2": "471304005",
        "3": "0",
        "4": "438720154",
        "5": "0",
        "6": "452179604",
        "7": "0",
        "8": "462376173",
        "9": "0",
        "10": "453218924",
        "11": "0",
        "12": "462216879",
        "13": "0",
        "14": "454103970",
        "15": "0",
        "16": "448155269",
        "17": "0",
        "18": "476235359",
        "19": "0",
        "20": "452555092",
        "21": "0",
        "22": "451014251",
        "23": "0",
        "24": "451004620",
        "25": "0",
        "26": "453230519",
        "27": "0",
        "28": "451372872",
        "29": "0",
        "30": "451901030",
        "31": "0",
        "32": "452440856",
        "33": "0",
        "34": "452077197",
        "35": "0",
        "36": "450509124",
        "37": "0",
        "38": "450046986",
        "39": "0",
        "40": "453822637",
        "41": "0",
        "42": "453812268",
        "43": "0",
        "44": "0"
    ]
    return chrOffsets[chr]
}

def getRefV1ChrLength(chr) {
    def chrLengths = [
        "0": "480980714",
        "1": "471304005",
        "2": "122798051",
        "3": "438720154",
        "4": "251131716",
        "5": "452179604",
        "6": "43273582",
        "7": "462376173",
        "8": "318422384",
        "9": "453218924",
        "10": "348037791",
        "11": "462216879",
        "12": "189635730",
        "13": "454103970",
        "14": "296739669",
        "15": "448155269",
        "16": "382674495",
        "17": "476235359",
        "18": "139317064",
        "19": "452555092",
        "20": "292033065",
        "21": "451014251",
        "22": "222603248",
        "23": "451004620",
        "24": "58852447",
        "25": "453230519",
        "26": "256543224",
        "27": "451372872",
        "28": "261776885",
        "29": "451901030",
        "30": "114179647",
        "31": "452440856",
        "32": "165638404",
        "33": "452077197",
        "34": "268911281",
        "35": "450509124",
        "36": "23083594",
        "37": "450046986",
        "38": "286659250",
        "39": "453822637",
        "40": "296797748",
        "41": "453812268",
        "42": "184873787",
        "43": "452528",
        "44": "134545"
    ]
    return chrLengths[chr]
}
