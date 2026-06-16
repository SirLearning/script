nextflow.enable.dsl=2

// Frozen per-chromosome TIGER popdepth grids under main_raw/variant (one-off partial build; reuse thereafter).

def popdepChrPath(popdep_dir, id) {
    return file("${popdep_dir}/variant/${id}.popdep.txt")
}

def hasPopdepForChr(popdep_dir, id) {
    return popdepChrPath(popdep_dir, id).exists()
}

def countPopdeps(popdep_dir) {
    def variant_dir = new File("${popdep_dir}/variant")
    if (!variant_dir.isDirectory()) {
        return 0
    }
    variant_dir.list({ _dir, name -> name.endsWith('.popdep.txt') } as FilenameFilter)?.length ?: 0
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
