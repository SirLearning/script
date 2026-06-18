nextflow.enable.dsl=2

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

def getPopDepTaxaBamFileAll_v1(home_dir) {
    return "${home_dir}/00data/05taxaBamMap/vmap4_v1/tb.ALL.txt"
}

def getPopDepSubGenomeForChr_v1(chr) {
    def normalized = chr.toString().replaceFirst(/^chr/, '')

    def groupA = getRefV1SubChr("A")
    def groupB = getRefV1SubChr("B")
    def groupD = getRefV1SubChr("D")
    def groupOthers = getRefV1SubChr("Others")

    if (groupA.contains(normalized)) {
        return "A"
    }
    if (groupB.contains(normalized)) {
        return "B"
    }
    if (groupD.contains(normalized)) {
        return "D"
    }
    if (groupOthers.contains(normalized)) {
        return "ALL"
    }
    throw new IllegalArgumentException("Unknown chromosome: ${chr} (normalized: ${normalized})")
}

def countPopDepTaxaFile_v1(home_dir, sub_genome) {
    def tbPaths = [
        "A": "${home_dir}/00data/05taxaBamMap/vmap4_v1/tb.A.txt",
        "B": "${home_dir}/00data/05taxaBamMap/vmap4_v1/tb.B.txt",
        "D": "${home_dir}/00data/05taxaBamMap/vmap4_v1/tb.D.txt",
        "ALL": "${home_dir}/00data/05taxaBamMap/vmap4_v1/tb.ALL.txt",
    ]
    def path = tbPaths[sub_genome]
    if (!path) {
        throw new IllegalArgumentException("Unknown PopDep sub_genome: ${sub_genome}")
    }
    def tbFile = file(path)
    if (!tbFile.exists()) {
        throw new FileNotFoundException("PopDep taxa-BAM file not found: ${path}")
    }
    return tbFile.readLines().size() - 1
}

def getPopDepNTaxaForChr_v1(chr, home_dir) {
    return countPopDepTaxaFile_v1(home_dir, getPopDepSubGenomeForChr_v1(chr))
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

def getRefFastaForChr_v1(chr, home_dir) {
    def normalized = chr.toString().replaceFirst(/^chr/, '')

    def subGenomeRefs = [
        "A": "${home_dir}/00data/03ref/01A/a_iwgscV1.fa.gz",
        "B": "${home_dir}/00data/03ref/05B/b_iwgscV1.fa.gz",
        "D": "${home_dir}/00data/03ref/04D/d_iwgscV1.fa.gz",
        "ALL": "${home_dir}/00data/03ref/01A/a_iwgscV1.fa.gz"
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
    return subGenomeRefs[sub_genome]
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
