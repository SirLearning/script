nextflow.enable.dsl=2

// Subgenomes that already have a merged *_test.plink2. A/B/D are required for
// reuse; Others is optional (exon-only extracts have no unplaced/organelle genes).
def mergedSubgenomeIds(process_dir) {
    if (!process_dir) {
        return []
    }
    return ['A', 'B', 'D', 'Others'].findAll { sg ->
        file("${process_dir}/${sg}_test.plink2.pgen").exists()
    }
}

// True when process_dir already holds merged A/B/D plink2 outputs.
def hasMergedSubgenomeTestPfiles(process_dir) {
    def found = mergedSubgenomeIds(process_dir)
    return found.containsAll(['A', 'B', 'D'])
}

// PLINK 1.9 bfiles for merged A/B/D test mods (Circos LD blocks, hom-SNP IBS).
def listMergedSubgenomeTestPlink1BedTuples(process_dir) {
    return ['A', 'B', 'D'].findAll { sg ->
        file("${process_dir}/${sg}_test.plink.bed").exists()
    }.collect { sg ->
        tuple(
            sg,
            "${sg}_test.plink",
            file("${process_dir}/${sg}_test.plink.bed"),
            file("${process_dir}/${sg}_test.plink.bim"),
            file("${process_dir}/${sg}_test.plink.fam"),
        )
    }
}

// Tuples matching merge_subgenome_test_pfile emit pfile shape: subgenome, chr_id, prefix, pgen, psam, pvar.
def listMergedSubgenomeTestPfileTuples(process_dir) {
    return mergedSubgenomeIds(process_dir).collect { sg ->
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
    def subgenomes = mergedSubgenomeIds(process_dir)
    if (!subgenomes) {
        return false
    }
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

// Tuples for sample_plink1_pca_plot: subgenome, chr_id, eigenvec, eigenval.
def listMergedSubgenomeTestPcaTuples(pca_dir) {
    def subgenomes = ['A', 'B', 'D', 'Others']
    return subgenomes.collect { sg ->
        tuple(
            sg,
            "sub_${sg}",
            file("${pca_dir}/${sg}.pca.eigenvec"),
            file("${pca_dir}/${sg}.pca.eigenval"),
        )
    }
}

// Tuples for sample_plink1_pca_ld_plot from process/{mod}/pca_ld/ (A/B/D only).
def listMergedSubgenomeTestPcaLdTuples(pca_ld_dir) {
    def subgenomes = ['A', 'B', 'D']
    return subgenomes.collect { sg ->
        tuple(
            sg,
            "sub_${sg}",
            file("${pca_ld_dir}/${sg}.pca.ld.eigenvec"),
            file("${pca_ld_dir}/${sg}.pca.ld.eigenval"),
        )
    }
}

// Tuples matching merge_subgenome_test_pfile emit bfile shape.
def listMergedSubgenomeTestBfileTuples(process_dir) {
    def subgenomes = mergedSubgenomeIds(process_dir)
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
