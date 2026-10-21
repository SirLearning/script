nextflow.enable.dsl=2

process sample_manifest_build {
    tag "sample manifest ${params.mod}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.keep"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.{manifest.tsv,qc.tsv,ref_v1.tsv}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'

    input:
    path group_file
    val taxa_bam_map_dir
    val home_dir

    output:
    path "assoc.keep", emit: assoc_keep
    path "wgs_core.keep", emit: wgs_keep
    path "*.manifest.tsv", emit: manifest
    path "*.qc.tsv", emit: qc
    path "*.th.tsv", emit: thresholds
    path "*.png", emit: plots, optional: true
    path "*.ref_v1.tsv", emit: ref_inventory
    path "*.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("sample_manifest.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.germplasm.sample.cohort_manifest import run_sample_manifest
    from genetics.genomics.ref.v1_assets import run_ref_v1_inventory

    run_sample_manifest(
        "${group_file}",
        "sample_cohort",
        taxa_bam_map_dir="${taxa_bam_map_dir}" if "${taxa_bam_map_dir}" != "null" else None,
        wgs_core_groups="${params.wgs_core_groups}",
        assoc_groups="${params.assoc_groups}",
    )
    extra = {}
    for key, path in (
        ("gff_lulab", "${params.ref_v1_gff}"),
        ("gff_hc", "${params.ref_v1_gff_hc}"),
        ("cds_hc", "${params.ref_v1_cds_hc}"),
        ("transcripts_hc", "${params.ref_v1_transcripts_hc}"),
    ):
        if path and path != "null":
            extra[key] = path
    run_ref_v1_inventory("${home_dir}", "sample_cohort", extra_assets=extra)
    """
}

process dual_af_from_afreq {
    tag "dual AF ${prefix}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'

    input:
    tuple val(prefix), path(wgs_afreq), path(wap_afreq)

    output:
    tuple val(prefix), path("*.dual_af.tsv"), emit: table
    path "*.th.tsv", emit: thresholds
    path "*.png", emit: plots, optional: true
    path "*.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("${prefix}.dual_af.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.dual_af import run_dual_af

    run_dual_af(
        "${wgs_afreq}",
        "${wap_afreq}",
        "${prefix}",
        rare_max=float("${params.dual_af_rare_max}"),
    )
    """
}

process coding_mask_from_vep {
    tag "coding mask ${prefix}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'

    input:
    tuple val(prefix), path(vep_table)

    output:
    tuple val(prefix), path("*.coding_mask.tsv"), emit: table
    path "*.th.tsv", emit: thresholds
    path "*.png", emit: plots, optional: true
    path "*.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("${prefix}.coding_mask.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.coding_mask import run_coding_mask

    run_coding_mask("${vep_table}", "${prefix}")
    """
}

process variant_sets_build {
    tag "variant sets ${prefix}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/variant", mode: 'copy', pattern: "*.extract"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.variant_sets.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'

    input:
    tuple val(prefix), path(coding_mask_table), path(dual_af_table), path(site_anno_table), path(gene_expr_table)

    output:
    path "*.variant_sets.tsv", emit: table
    path "*.extract", emit: extracts
    path "*.th.tsv", emit: thresholds
    path "*.png", emit: plots, optional: true
    path "*.log", emit: logs

    script:
    def set_names = params.variant_sets_names ? "\"${params.variant_sets_names}\"" : 'None'
    def site_arg = site_anno_table.name == 'NO_SITE_ANNO' ? 'None' : "\"${site_anno_table}\""
    def expr_arg = gene_expr_table.name == 'NO_GENE_EXPR' ? 'None' : "\"${gene_expr_table}\""
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("${prefix}.variant_sets.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.variant_sets import run_variant_sets

    run_variant_sets(
        "${coding_mask_table}",
        "${dual_af_table}",
        "${prefix}",
        set_names=${set_names},
        common_min=float("${params.variant_sets_common_min}"),
        site_anno_table=${site_arg},
        gene_expr_table=${expr_arg},
    )
    """
}
