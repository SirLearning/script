#!/usr/bin/env nextflow
/*
 * Same combine-group wiring as combine_ibs_ploidy_pca, but on the IBS-PCoA
 * coordinates instead of the variant PCA.
 *
 * One PCoA table per ploidy group, each built from the subgenome matrices that
 * group actually carries (AA -> A, AABB -> A+B, AABBDD -> A+B+D, DD -> D). The
 * tables are published alongside the variant PCA, then one figure is drawn per
 * cohort x group; cohorts that do not span a group (IPK/HZAU/Watkins is AABBDD
 * only) emit a single panel.
 *
 *   nextflow run .../subworkflows/tmp/combine_ibs_pcoa_plot.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_common_inter \
 *     --pcoa_dir /data/home/tusr1/01projects/vmap4/16stats.genome/21run_ibs_pcoa_groups
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mod = 'test_common_inter'
params.pcoa_dir = null
params.ploidies = 'AA,AABB,AABBDD,DD'

process publish_ibs_pcoa {
    tag "publish IBS-PCoA coordinates: ${group}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 1
    memory '4.GB'
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/pcoa_ibs", mode: 'copy', pattern: "*.{eigenvec,eigenval}"

    input:
    tuple val(group), path(eigenvec)

    output:
    tuple val(group), path("${group}.pcoa.eigenvec"), emit: coords

    script:
    """
    cp ${eigenvec} published_${group}.eigenvec
    mv published_${group}.eigenvec ${group}.pcoa.eigenvec
    """
}

process combine_ibs_pcoa_plot {
    tag "combine IBS PCoA: ${prefix} ${group}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '16.GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(prefix), val(threshold), val(group), path(pairs), path(eigenvec)

    output:
    path("*.png"), emit: plots, optional: true
    path("${prefix}.pcoa.${group}.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${prefix}.pcoa.${group}.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.sample.combine_ibs_validate import plot_combine_ibs_pca

    for categories in (None, ["pass"], ["fail"]):
        plot_combine_ibs_pca(
            "${pairs}",
            "${eigenvec}",
            "${prefix}",
            subgenome="${group}",
            threshold=${threshold},
            suffix="pcoa",
            allow_empty_edges=True,
            edge_categories=categories,
        )
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'combine_ibs_pcoa_plot: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir || !params.src_dir) {
        error 'combine_ibs_pcoa_plot: --user_dir and --src_dir are required'
    }
    if (!params.pcoa_dir) {
        error 'combine_ibs_pcoa_plot: --pcoa_dir is required'
    }
    def info = "${params.output_dir}/${params.job}/stats/${params.mod}/info"

    def requested = params.ploidies.toString().split(',').collect { it.trim() }
    def eigenvecs = Channel
        .fromPath("${params.pcoa_dir}/*.ibs_hom_snp.pcoa.eigenvec")
        .map { file -> tuple(file.name.split('\\.')[0], file) }
        .filter { group, file -> requested.contains(group) }

    def published = publish_ibs_pcoa(eigenvecs).coords

    // Cutoff per cohort, plus the populations that cohort actually spans.
    def cohorts = Channel.of(
        tuple('all.sample.combine_ibs', 0.995 as double, 'AA,AABB,AABBDD,DD'),
        tuple('all.sample.combine_ibs_ipk_hzau_watkins', 0.99 as double, 'AABBDD'),
    )
    def specs = published
        .combine(cohorts)
        .filter { item ->
            item[4].toString().split(',').collect { it.trim() }.contains(item[0])
        }
        .map { group, eigenvec, prefix, threshold, allowed ->
            tuple(
                prefix,
                threshold,
                group,
                file("${info}/${prefix}.pairs.tsv", checkIfExists: true),
                eigenvec,
            )
        }
    combine_ibs_pcoa_plot(specs)
}
