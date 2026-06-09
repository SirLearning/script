nextflow.enable.dsl=2

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
