nextflow.enable.dsl=2

// Frozen per-chromosome TIGER popdepth grids under main_raw/variant (one-off partial build; reuse thereafter).

def popdepChrPath(popdep_dir, id) {
    return file("${popdep_dir}/variant/${id}.popdep.txt")
}

def popdepChrBgzPath(popdep_dir, id) {
    return file("${popdep_dir}/variant/${id}.popdep.txt.bgz")
}

def hasPopdepForChr(popdep_dir, id) {
    return popdepChrBgzPath(popdep_dir, id).exists() ||
        popdepChrPath(popdep_dir, id).exists()
}

def countPopdeps(popdep_dir) {
    def variant_dir = new File("${popdep_dir}/variant")
    if (!variant_dir.isDirectory()) {
        return 0
    }
    def bgz = variant_dir.list({ _dir, name -> name.endsWith('.popdep.txt.bgz') } as FilenameFilter)?.length ?: 0
    if (bgz > 0) {
        return bgz
    }
    return variant_dir.list({ _dir, name -> name.endsWith('.popdep.txt') && !name.endsWith('.popdep.txt.bgz') } as FilenameFilter)?.length ?: 0
}

// True when process_dir/variant already has per-subgenome popdepth annotation (*.popdep.info.tsv).
def hasVariantPopdepForMergedTests(process_dir) {
    if (!process_dir) {
        return false
    }
    def subgenomes = ['A', 'B', 'D', 'Others']
    return subgenomes.every { sg ->
        file("${process_dir}/variant/${sg}.popdep.info.tsv").exists()
    }
}
