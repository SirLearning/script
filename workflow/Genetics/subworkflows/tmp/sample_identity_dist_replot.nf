#!/usr/bin/env nextflow
/*
 * Republish identity pair-score histograms from existing
 * stats/{mod}/info/*.sample.identity.*.pairs.tsv
 *
 *   cd /data/home/tusr1/01projects/vmap4/12stats.genome/<NNrun_identity_dist_replot>
 *   nextflow run .../subworkflows/tmp/sample_identity_dist_replot.nf \
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
        error 'sample_identity_dist_replot: --output_dir and --job are required'
    }
    if (!params.user_dir) {
        error 'sample_identity_dist_replot: --user_dir is required (Conda stats env)'
    }

    def mods = params.mods.toString().split(',').collect { it.trim() }.findAll { it }
    def methods = ['king', 'jaccard', 'discord', 'ibd', 'fingerprint']
    def rows = []
    mods.each { mod ->
        ['A', 'B', 'D'].each { sg ->
            methods.each { method ->
                def table = file("${params.output_dir}/${params.job}/stats/${mod}/info/${sg}.sample.identity.${method}.pairs.tsv")
                if (table.exists()) {
                    rows << tuple(mod, sg, method, table)
                }
            }
        }
    }
    if (!rows) {
        error 'sample_identity_dist_replot: no identity pair tables found'
    }
    plot_identity_pair_dist(channel.from(rows))
}

process plot_identity_pair_dist {
    tag "identity dist: ${mod} ${id} ${method}"
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(mod), val(id), val(method), path(pairs)

    output:
    path "*.png", optional: true

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${id}.${method}.identity_dist.log", "w")
    sys.stderr = sys.stdout
    from genetics.genomics.sample.identity import plot_identity_pair_table
    plot_identity_pair_table(
        "${pairs}",
        "${id}.sample.identity.${method}",
        "${method}",
    )
    """
}
