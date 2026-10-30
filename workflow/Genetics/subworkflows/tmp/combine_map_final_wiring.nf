#!/usr/bin/env nextflow
/*
 * Draw the final WAP_PQ/VMap3 combine-map merges as edges on the published
 * variant PCA and IBS-PCoA coordinates.
 *
 * Edges come from combine_map.pair_decisions.tsv: only pairs with
 * ``merges = True`` whose sources are WAP/VMap3 (Watkins and IPK/HZAU are
 * excluded). PCA is drawn per subgenome (A/B/D/Others); PCoA per ploidy whose
 * coordinates hold at least one merged edge (AABBDD/AABB/DD).
 *
 * Publishes under {output_dir}/{job}/stats/{mod}/{plots,logs}:
 *   combine_map.final_merge.pca.{A,B,D,Others}.png
 *   combine_map.final_merge.pcoa.{AABBDD,AABB,DD}.png
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/combine_map_final_wiring.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_common_inter
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mod = 'test_common_inter'
params.decisions = '/data/home/tusr1/git/script/src/resource/sample_combine/combine_map.pair_decisions.tsv'
params.pca_views = 'A,B,D,Others'
params.pcoa_views = 'AABBDD,AABB,DD'

process combine_map_final_wiring {
    tag "final merge wiring: ${kind} ${axis_label}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '16.GB'
    errorStrategy 'terminate'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(kind), val(axis_label), path(coords), path(decisions)

    output:
    path("*.png"), emit: plots, optional: true
    path("*.log"), emit: logs

    script:
    def view = "${axis_label} ${kind == 'pca' ? 'PCA' : 'PCoA'}"
    def prefix = "combine_map.final_merge.${kind}.${axis_label}"
    """
    #!/usr/bin/env python
    import logging
    import sys
    sys.stdout = open("${prefix}.log", "w", buffering=1)
    sys.stderr = sys.stdout
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.sample.combine_ibs_validate import (
        plot_final_merge_wiring,
    )

    plot_final_merge_wiring(
        "${decisions}",
        "${coords}",
        "${prefix}",
        axis_label="${axis_label}",
        view="${view}",
    )
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'combine_map_final_wiring: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir || !params.src_dir) {
        error 'combine_map_final_wiring: --user_dir and --src_dir are required'
    }
    def process_dir = "${params.output_dir}/${params.job}/process/${params.mod}"
    def decisions = file(params.decisions, checkIfExists: true)

    def specs = []
    params.pca_views.toString().split(',').each { raw ->
        def sg = raw.trim()
        if (sg) {
            specs << tuple(
                'pca',
                sg,
                file("${process_dir}/pca/${sg}.pca.eigenvec", checkIfExists: true),
                decisions,
            )
        }
    }
    params.pcoa_views.toString().split(',').each { raw ->
        def group = raw.trim()
        if (group) {
            specs << tuple(
                'pcoa',
                group,
                file("${process_dir}/pcoa_ibs/${group}.pcoa.eigenvec", checkIfExists: true),
                decisions,
            )
        }
    }
    combine_map_final_wiring(Channel.fromList(specs))
}
