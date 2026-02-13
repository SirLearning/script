nextflow.enable.dsl=2

include { getTaxaBamMapFile_v1 } from './utils.nf'

// --- Genotype statistics workflow ---

workflow test_plink_stats {
    take:
    // Expect a channel: [ val(id), path(smiss) ]
    smiss
    scount
    vmiss
    gcount
    afreq
    hardy

    main:
    // prepare input
    def idxstats = file("${params.output_dir}/vmap4_v1_idxstat_summary.txt")
    def group_file = file("${params.output_dir}/sample_group.txt")
    def tbm_dir = file("${params.home_dir}/00data/05taxaBamMap")
    def dbone_dir = file("${params.user_dir}/git/DBone/Service/src/main/resources/raw/20251208")
    def tbm = smiss.map { id, chr, _smiss_path ->
        def subgenome_tbm = ''
        if (chr.startsWith('sub')) {
            if (id == "Others") {
                subgenome_tbm = "${params.home_dir}/00data/05taxaBamMap/all.ALL.taxaBamMap.txt"
            } else {
                subgenome_tbm = "${params.home_dir}/00data/05taxaBamMap/all.${id}.taxaBamMap.txt"
            }
        } else {
            subgenome_tbm = getTaxaBamMapFile_v1(chr, params.home_dir)
        }
        subgenome_tbm = file(subgenome_tbm)
        tuple(id, chr, subgenome_tbm)
    }
    // 1 sample stats
    def miss_out = sample_missing_stats(smiss)
    def cov_out = sample_coverage_stats(smiss, tbm)
    def het_out = sample_heterozygosity_stats(scount, smiss)
    def mr_out = sample_mapping_rate_stats(idxstats, group_file, smiss)
    // 1.1 kinship
    def ref_ibs_out = sample_ref_ibs_stats(group_file, idxstats, scount, smiss)
    // 1.2 germplasm dedup
    def dedup_out = sample_germplasm_dedup(smiss, tbm_dir, dbone_dir)

    // collect sample info and thresholds
    def sample_info = miss_out.info.combine(het_out.info).combine(mr_out.info).combine(ref_ibs_out.info).combine(cov_out.info)
    def sample_th = miss_out.th.combine(het_out.th).combine(dedup_out.th).combine(cov_out.th)
    
    // 2 variant stats
    // 2.1 Basic Variant Stats (Missing, MAF)
    def vmiss_out = variant_missing_stats(vmiss)
    def maf_out = variant_maf_stats(afreq)

    // emit:
}

workflow plink_stats {
    take:
    // Expect a channel: [ val(id), path(smiss) ]
    smiss
    scount
    vmiss
    gcount
    afreq
    hardy
    popdep

    main:
    // def stats_out = sample_missing_stats(smiss)
    def vmiss_out = variant_missing_stats(vmiss)
    def maf_out = variant_maf_stats(afreq)
    // def popdep_stats_out = variant_popdep_stats(popdep, vmiss)

    // emit:
}

// ... existing processes ...

process variant_missing_stats {
    tag "variant missing stats: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(vmiss)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.vmiss_ana.log", "w")
    sys.stderr = sys.stdout # f2

    from genetics.genomics.variant import ana_variant_missing
    print(f"Processing variant missing stats for ${chr}...")
    ana_variant_missing("${vmiss}", "${id}.variant.miss")
    """
}

process variant_maf_stats {
    tag "variant maf stats: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(afreq)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.maf_ana.log", "w")
    sys.stderr = sys.stdout # f2

    from genetics.genomics.variant import ana_minor_allele_frequency

    print(f"Processing variant MAF stats for ${chr}...")
    ana_minor_allele_frequency("${afreq}", "${id}.variant.maf")

    """
}

process variant_popdep_stats {
    tag "variant popdep stats: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(popdep)
    tuple val(id2), val(chr2), path(vmiss)

    output:
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.popdep_ana.log", "w")
    sys.stderr = sys.stdout 

    from genetics.genomics.variant.popdep import ana_popdep_curve_fit, ana_popdep_missing_reg

    print(f"Processing population depth Analysis 1 (Curve Fit) for ${chr}...")
    ana_popdep_curve_fit("${popdep}", "${id}.variant.popdep")

    print(f"Processing population depth Analysis 2 (vs Missing) for ${chr}...")
    ana_popdep_missing_reg("${popdep}", "${vmiss}", "${id}.variant.popdep")
    """
}

process sample_missing_stats {
    tag "compute missing rate threshold: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(smiss)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.smiss_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import ana_sample_missing

    print(f"Processing sample missing rate for ${chr}...")
    ana_sample_missing("${smiss}", "${id}.sample.miss")
    """
}

process sample_coverage_stats {
    tag "compute coverage stats: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(smiss)
    tuple val(id2), val(chr2), path(tbm)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.coverage_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import coverage_dist, reg_missing_coverage

    print(f"Processing coverage for ${chr}...")
    coverage_dist("${tbm}", "${id}.sample.cov")
    reg_missing_coverage("${tbm}", "${smiss}", "${id}.sample.cov")
    """
}

process sample_heterozygosity_stats {
    tag "compute heterozygosity rate: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(scount)
    tuple val(id2), val(chr2), path(smiss)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.het_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import ana_heterozygosity

    print(f"Processing heterozygosity rate for ${chr}...")
    ana_heterozygosity("${scount}", "${smiss}", "${id}.sample.het")
    """
}

process sample_mapping_rate_stats {
    tag "compute mapping rate stats: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path idxstats
    path group_file
    tuple val(id), val(chr), path(smiss)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.mr_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import reg_missing_idx_mapping

    print(f"Processing mapping rate for ${chr}...")
    reg_missing_idx_mapping("${idxstats}", "${smiss}", "${id}.sample.idx.mr", "${group_file}")
    """
}

process sample_ref_ibs_stats {
    tag "compute ref IBS stats"
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path group_file
    path idxstats
    tuple val(id), val(chr), path(scount)
    tuple val(id2), val(chr2), path(smiss)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.ref_ibs_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import ana_ref_ibs, ref_ibs_vs_mapping

    print(f"Processing reference IBS stats for ${chr}...")
    ana_ref_ibs("${scount}", "${smiss}", "${group_file}", "${id}.sample.ref_ibs")
    
    print(f"Processing reference IBS vs Mapping Rate for ${chr}...")
    # mapping_file, scount_file, group_file, output_prefix
    ref_ibs_vs_mapping("${idxstats}", "${scount}", "${group_file}", "${id}.sample.ref_ibs")
    """
}

process sample_king_stats {
    tag "compute king stats: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/thresholds", mode: 'copy', pattern: "*.stats.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(king_file)
    tuple val(id2), val(chr2), path(scount)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.stats.tsv"), emit: stats
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.king_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import ana_king_kinship, ana_derived_het, ana_het_vs_max_kinship

    print(f"Processing KING Kinship for ${chr}...")
    ana_king_kinship("${king_file}", "${id}.sample.king")
    
    print(f"Processing Derived Heterozygosity for ${chr}...")
    ana_derived_het("${king_file}", "${id}.sample.derived_het")
    
    print(f"Processing Het vs Max Kinship for ${chr}...")
    ana_het_vs_max_kinship("${scount}", "${king_file}", "${id}.sample.het_vs_kin")
    """
}

process sample_ibs_stats {
    tag "compute ibs stats: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path group_file
    tuple val(id), val(chr), path(ibs_matrix)
    tuple val(id2), val(chr2), path(ids_file)
    tuple val(id3), val(chr3), path(scount)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.ibs_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import ana_ibs_matrix, ana_het_vs_max_ibs, ana_ibs_trend

    print(f"Processing IBS Matrix Analysis for ${chr}...")
    # matrix_file, id_file, group_file, output_prefix
    ana_ibs_matrix("${ibs_matrix}", "${ids_file}", "${group_file}", "${id}.sample.ibs")
    
    print(f"Processing Het vs Max IBS for ${chr}...")
    # scount_file, matrix_file, id_file, output_prefix
    ana_het_vs_max_ibs("${scount}", "${ibs_matrix}", "${ids_file}", "${id}.sample.het_vs_ibs")

    print(f"Processing IBS Trend for ${chr}...")
    # matrix_file, id_file, group_file, output_prefix
    ana_ibs_trend("${ibs_matrix}", "${ids_file}", "${group_file}", "${id}.sample.ibs_trend")
    """
}

process sample_germplasm_dedup {
    tag "detect germplasm duplicates: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(smiss)
    path tbm_dir
    path dbone_dir

    output:
    tuple val(id), val(chr), path("${id}.sample.dedup.rm.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.germplasm.duplicates.log", "w")
    sys.stderr = sys.stdout

    from genetics.germplasm.sample import ana_duplication

    print(f"Detecting germplasm duplicates for ${chr}...")
    ana_duplication("${tbm_dir}", "${dbone_dir}", "${id}.sample.dedup")
    """
}

process plink_pca {
    tag "${id}"
    publishDir "${params.output_dir}/${params.job}/stats/pca", mode: 'copy'

    input:
    tuple val(id), path(bed), path(bim), path(fam)

    output:
    tuple val(id), path("${id}.eigenvec"), path("${id}.eigenval"), emit: pca

    script:
    """
    set -euo pipefail
    exec > ${id}.plink_pca.log 2>&1

    plink --bfile ${id} \\
        --pca ${params.pc_num} \\
        --out ${id} \\
        --allow-extra-chr \\
        --chr-set 42 \\
        --threads ${task.cpus}
    """
}

process PLOT_PLINK_PCA {
    tag "${id}"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy'

    input:
    tuple val(id), path(eigenvec), path(eigenval)
    path metadata

    output:
    tuple val(id), path("${id}.pca_plot.pdf"), emit: pca_plot

    script:
    def md_arg = metadata.name != 'NO_METADATA' ? "--metadata ${metadata}" : ""
    """
    set -euo pipefail
    Rscript ${params.src_dir}/r/genetics/plot_pca.r \\
        --eigenvec ${eigenvec} \\
        --eigenval ${eigenval} \\
        ${md_arg} \\
        --output ${id}.pca_plot.pdf
    """
}

process vcftools_vcf_qc_r {
    tag "${id}" ? "plot qc ${id}" : 'plot qc'
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy'

    input:
    tuple val(id), path(imiss), path(lmiss), path(het), path(frq), path(idepth), path(ldepth)

    output:
    path("${id}.qc_plots.pdf")

    script:
    """
    set -euo pipefail
    Rscript ${params.src_dir}/r/genetics/vcf_qc_plot.r \\
        --imiss ${imiss} \\
        --lmiss ${lmiss} \\
        --het ${het} \\
        --frq ${frq} \\
        --depth ${idepth} \\
        --site_depth ${ldepth} \\
        --output ${id}.qc_plots.png
    """
}

