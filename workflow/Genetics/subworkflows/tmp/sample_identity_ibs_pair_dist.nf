#!/usr/bin/env nextflow
/*
 * Full-range PLINK and hom-SNP IBS histograms from the published
 * proportional pair sample (stats/{mod}/info/*.hom_snp.vs_plink.pairs.tsv).
 * Does not zoom to IBS>0.99.
 *
 *   cd /data/home/tusr1/01projects/vmap4/12stats.genome/<NNrun_ibs_pair_dist>
 *   nextflow run .../subworkflows/tmp/sample_identity_ibs_pair_dist.nf \
 *     -c .../nextflow.config \
 *     --user_dir /data/home/tusr1 \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mods test_thin,test_common_inter
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mods = 'test_thin,test_common_inter'

workflow {
    if (!params.output_dir || !params.job) {
        error 'sample_identity_ibs_pair_dist: --output_dir and --job are required'
    }
    if (!params.user_dir) {
        error 'sample_identity_ibs_pair_dist: --user_dir is required (Conda stats env)'
    }

    def mods = params.mods.toString().split(',').collect { it.trim() }.findAll { it }
    def rows = []
    mods.each { mod ->
        ['A', 'B', 'D'].each { sg ->
            def table = file("${params.output_dir}/${params.job}/stats/${mod}/info/${sg}.sample.identity.hom_snp.vs_plink.pairs.tsv")
            if (table.exists()) {
                rows << tuple(mod, sg, table)
            }
        }
    }
    if (!rows) {
        error "sample_identity_ibs_pair_dist: no vs_plink pair tables under mods ${mods}"
    }
    plot_ibs_sampled_pair_dist(channel.from(rows))
}

process plot_ibs_sampled_pair_dist {
    tag "ibs pair dist: ${mod} ${id}"
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(mod), val(id), path(pairs)

    output:
    tuple val(mod), val(id),
        path("${id}.sample.identity.ibs_pair.th.tsv"),
        path("${id}.sample.identity.ibs_pair.plink.dist.png"),
        path("${id}.sample.identity.ibs_pair.hom.dist.png"),
        path("${id}.sample.identity.ibs_pair.plink.dist_log.png"),
        path("${id}.sample.identity.ibs_pair.hom.dist_log.png")

    script:
    """
    #!/usr/bin/env python
    import sys
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")
    sys.stdout = open("${id}.identity_ibs_pair_dist.log", "w")
    sys.stderr = sys.stdout
    from genetics.genomics.sample.identity import ana_ibs_sampled_pair_dist
    ana_ibs_sampled_pair_dist(
        "${pairs}",
        "${id}.sample.identity.ibs_pair",
    )
    """
}
