nextflow.enable.dsl=2

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

def getCallingJobConfig(job) {
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
