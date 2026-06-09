nextflow.enable.dsl=2

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
