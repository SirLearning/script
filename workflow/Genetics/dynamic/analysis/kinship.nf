nextflow.enable.dsl=2

workflow kinship {
    take:
    // Expect a combined channel: [ val(meta), path(vcf), val(job_config) ]
    vcf_in

    main:
    // 1. Calculate IBS using PLINK
    ibs_results = calculate_ibs(vcf_in)

    // 2. Plot IBS distribution (Boxplot)
    // 需要一个分组文件。这里假设 job_config 中配置了 group_file，或者通过 params 传递
    // 为了演示，我们假设 params.group_file 存在
    if (params.group_file) {
        plot_ibs_boxplot(ibs_results.ibs_data, file(params.group_file))
    } else {
        log.warn "No group file provided (params.group_file), skipping IBS boxplot."
    }

    // 3. Calculate Kinship Matrix (King)
    kinship_results = calculate_kinship(vcf_in)

    emit:
    kinship = kinship_results.kinship_matrix
}

process calculate_ibs {
    tag { meta.id }
    publishDir "${params.output_dir}/${params.job}/genotype/kinship", mode: 'copy'

    input:
    tuple val(meta), path(vcf), val(job_config)

    output:
    tuple val(meta), path("${meta.id}.mibs"), path("${meta.id}.mibs.id"), emit: ibs_data

    script:
    def prefix = meta.id
    """
    # Calculate IBS distance matrix (square matrix)
    # --distance ibs: 1-IBS (distance)
    # --distance square: output square matrix
    plink --vcf ${vcf} \\
        --distance ibs square \\
        --allow-extra-chr \\
        --double-id \\
        --out ${prefix} \\
        --threads ${task.cpus}
    """
}

process plot_ibs_boxplot {
    tag { meta.id }
    publishDir "${params.output_dir}/${params.job}/genotype/kinship", mode: 'copy'

    input:
    tuple val(meta), path(mibs), path(ids)
    path group_file

    output:
    path("${meta.id}.ibs_boxplot.pdf"), emit: plot

    script:
    def prefix = meta.id
    // 默认参照样本为 CS，可以通过 params 修改
    def ref_sample = params.ref_sample ?: "CS"
    """
    Rscript ${params.src_dir}/r/genetics/static/genotype/ibs/plot_ibs_boxplot.r \\
        --ibs_files ${mibs} \\
        --id_files ${ids} \\
        --group_file ${group_file} \\
        --ref_sample "${ref_sample}" \\
        --output "${prefix}.ibs_boxplot.pdf"
    """
}

process calculate_kinship {
    tag { meta.id }
    publishDir "${params.output_dir}/${params.job}/genotype/kinship", mode: 'copy'
    
    input:
    tuple val(meta), path(vcf), val(job_config)

    output:
    tuple val(meta), path("*.kinship.txt"), emit: kinship_matrix

    script:
    def prefix = meta.id
    """
    # Use plink2 to calculate kinship (King-robust)
    plink2 \\
        --vcf ${vcf} \\
        --make-king-table \\
        --allow-extra-chr \\
        --out ${prefix} \\
        --threads ${task.cpus}

    mv ${prefix}.kin0 ${prefix}.kinship.txt
    """
}