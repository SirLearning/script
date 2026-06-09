nextflow.enable.dsl=2

// True when process_dir already holds merged per-subgenome plink2 outputs (A/B/D/Others _test.plink2).
def hasMergedSubgenomeTestPfiles(process_dir) {
    if (!process_dir) {
        return false
    }
    def subgenomes = ['A', 'B', 'D', 'Others']
    return subgenomes.every { sg ->
        file("${process_dir}/${sg}_test.plink2.pgen").exists()
    }
}

// Tuples matching merge_subgenome_test_pfile emit pfile shape: subgenome, chr_id, prefix, pgen, psam, pvar.
def listMergedSubgenomeTestPfileTuples(process_dir) {
    def subgenomes = ['A', 'B', 'D', 'Others']
    return subgenomes.collect { sg ->
        tuple(
            sg,
            "sub_${sg}",
            "${sg}_test.plink2",
            file("${process_dir}/${sg}_test.plink2.pgen"),
            file("${process_dir}/${sg}_test.plink2.psam"),
            file("${process_dir}/${sg}_test.plink2.pvar"),
        )
    }
}

// True when process_dir/variant already has PLINK2 --freq/--missing outputs from mk_plink_basic_info.
def hasPlinkBasicInfoForMergedTests(process_dir) {
    if (!process_dir) {
        return false
    }
    def subgenomes = ['A', 'B', 'D', 'Others']
    return subgenomes.every { sg ->
        file("${process_dir}/variant/${sg}.info.afreq").exists() &&
        file("${process_dir}/variant/${sg}.info.vmiss").exists()
    }
}

// Plot inputs for wheat_snp_qc: (afreq, vmiss, output_prefix) per merged subgenome.
// output_prefix is subgenome id only; job/mod live in publishDir.
def listMergedSubgenomeSnpQcPlotTuples(process_dir) {
    def subgenomes = ['A', 'B', 'D', 'Others']
    return subgenomes.collect { sg ->
        tuple(
            file("${process_dir}/variant/${sg}.info.afreq"),
            file("${process_dir}/variant/${sg}.info.vmiss"),
            "${sg}",
        )
    }
}

// Tuples matching merge_subgenome_test_pfile emit bfile shape.
def listMergedSubgenomeTestBfileTuples(process_dir) {
    def subgenomes = ['A', 'B', 'D', 'Others']
    return subgenomes.collect { sg ->
        tuple(
            sg,
            "sub_${sg}",
            "${sg}_test.plink",
            file("${process_dir}/${sg}_test.plink.bed"),
            file("${process_dir}/${sg}_test.plink.bim"),
            file("${process_dir}/${sg}_test.plink.fam"),
        )
    }
}
