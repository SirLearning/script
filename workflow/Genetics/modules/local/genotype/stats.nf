nextflow.enable.dsl=2

include { getTaxaBamMapFile_v1 } from './utils.nf'

// --- Genotype statistics workflow ---

workflow test_plink_stats {
    take:
    ld
    ld_cross
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
    variant_ld_decay_plot(ld)
    variant_ld_crosschr_plot(ld_cross)

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
    def popdep_maha_out = variant_popdep_mahalanobis(popdep)

    // emit:
}

// ... existing processes ...

process variant_missing_stats {
    tag "variant missing stats: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
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
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
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

process variant_popdep_mahalanobis {
    tag "variant popdep mahalanobis: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/thresholds", mode: "copy", pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: "copy", pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: "copy", pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(popdep)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.popdep_mahalanobis.log", "w")
    sys.stderr = sys.stdout 

    from genetics.genomics.variant.popdep import ana_popdep_mahalanobis

    print(f"Processing population depth Mahalanobis distance for ${chr}...")
    ana_popdep_mahalanobis("${popdep}", "${id}.variant.mahalanobis")
    """
}

process variant_ld_decay_plot {
    tag "variant ld decay: ${id}"
    label 'cpus_16'
    cache false
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(ld_vcor)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.ld_decay_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant import ana_ld_decay

    print("Running LD decay plotting for ${id} (${chr}) ...")
    ana_ld_decay(
        input_file="${ld_vcor}",
        output_prefix="${id}.variant.ld_decay",
        max_distance_kb=${params.ld_decay_max_kb},
        bin_size_kb=${params.ld_decay_bin_kb},
        scatter_max_points=${params.ld_decay_scatter_points},
        chunksize=${params.ld_decay_chunksize},
        median_sample_per_bin=${params.ld_decay_median_sample_per_bin},
        workers=${task.cpus},
        random_seed=${params.ld_decay_random_seed},
    )
    """
}

process variant_ld_crosschr_plot {
    tag "variant cross-chr ld baseline: ${id}"
    label 'cpus_16'
    cache false
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(ld_cross_vcor)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.ld_crosschr_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant import ana_ld_crosschr_baseline

    print("Running cross-chromosome LD baseline plotting for ${id} (${chr}) ...")
    ana_ld_crosschr_baseline(
        input_file="${ld_cross_vcor}",
        output_prefix="${id}.variant.ld_crosschr",
        chunksize=${params.ld_decay_chunksize},
    )
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
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
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
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
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
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
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
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
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
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
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

    from genetics.genomics.sample import safe_ref_ibs_stats

    print(f"Processing reference IBS stats for ${chr}...")
    safe_ref_ibs_stats(
        mapping_file="${idxstats}",
        scount_file="${scount}",
        missing_file="${smiss}",
        group_file="${group_file}",
        output_prefix="${id}.sample.ref_ibs",
    )
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
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
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

    from genetics.germplasm.sample import safe_ana_duplication

    print(f"Detecting germplasm duplicates for ${chr}...")
    safe_ana_duplication("${tbm_dir}", "${dbone_dir}", "${id}.sample.dedup")
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

// --- Integrated study plots (consume processor / awk outputs) ---

process plot_plink2_population_structure {
    tag "plot pca: ${id}"
    label 'cpus_2'
    cpus 8
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.wheat_plink_source_mod ?: params.mod}/${publish_mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.wheat_plink_source_mod ?: params.mod}/${publish_mod}/plots", mode: 'copy', pattern: "*.png"

    input:
    tuple val(id), val(chr), path(eigenvec), path(eigenval)
    val output_prefix
    val publish_mod

    output:
    path "*.tsv"
    path "*.png"

    script:
    """
    #!/usr/bin/env python
    from genetics.genomics.sample.pca_structure import plot_population_structure
    plot_population_structure(
        "${eigenvec}",
        "${eigenval}",
        "${output_prefix}",
        group_file="${params.output_dir}/sample_group.txt",
        tsne_n_input_pcs=${params.wheat_tsne_n_input_pcs},
        tsne_max_iter=${params.wheat_tsne_max_iter},
        tsne_random_state=${params.wheat_tsne_random_state},
        tsne_n_jobs=${task.cpus},
    )
    """
}

process report_plink2_tagsnp {
    tag "report tagsnp: ${id}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.wheat_plink_source_mod ?: params.mod}/${publish_mod}/info", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(id), val(chr), path(prune_in)
    val output_prefix
    val publish_mod

    output:
    path "*.tsv"

    script:
    """
    #!/usr/bin/env python
    from genetics.genomics.variant.tagsnp_report import report_tagsnp_selection
    report_tagsnp_selection("${prune_in}", "${output_prefix}", max_tags=${params.wheat_tagsnp_max_tags})
    """
}

process plot_plink2_snp_qc {
    tag "plot snp qc"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/plots", mode: 'copy', pattern: "*.png"

    input:
    path afreq
    path vmiss
    val output_prefix

    output:
    path "*.tsv"
    path "*.png"

    script:
    """
    #!/usr/bin/env python
    from genetics.genomics.variant.snp_qc_plot import summarize_and_plot_snp_qc
    summarize_and_plot_snp_qc(
        "${afreq}", "${vmiss}", "${output_prefix}",
        maf=${params.wheat_snp_qc_maf}, max_missing=${params.wheat_snp_qc_max_missing},
    )
    """
}

process plot_cnv_calls {
    tag "plot cnv"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/plots", mode: 'copy', pattern: "*.png"

    input:
    path cnv_tsv
    val output_prefix

    output:
    path "*.png"

    script:
    """
    #!/usr/bin/env python
    from genetics.genomics.variant.cnv_plot import plot_cnv_results
    plot_cnv_results("${cnv_tsv}", "${output_prefix}", del_z=${params.wheat_cnv_del_z}, dup_z=${params.wheat_cnv_dup_z})
    """
}

process plot_genetic_map {
    tag "plot genetic map"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/plots", mode: 'copy', pattern: "*.png"

    input:
    path map_tsv
    val output_prefix

    output:
    path "*.png"

    script:
    """
    #!/usr/bin/env python
    from genetics.genomics.variant.genetic_map_plot import plot_genetic_map
    plot_genetic_map("${map_tsv}", "${output_prefix}")
    """
}

process report_hapmap_table {
    tag "report hapmap"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/info", mode: 'copy', pattern: "*.tsv"

    input:
    path hapmap_tsv
    val output_prefix

    output:
    path "*.tsv"

    script:
    """
    #!/usr/bin/env python
    from genetics.genomics.variant.hapmap_report import report_hapmap_table
    report_hapmap_table("${hapmap_tsv}", "${output_prefix}")
    """
}

