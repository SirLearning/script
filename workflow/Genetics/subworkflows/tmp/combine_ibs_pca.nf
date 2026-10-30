#!/usr/bin/env nextflow
/*
 * Draw one PCA projection per subgenome with every combine group wired up.
 * Points come from the published test_common_inter PCA eigenvector tables;
 * each combine-name pair becomes an edge, red above the cohort hom-SNP IBS
 * cutoff and blue below it, with the endpoints of the blue edges labelled.
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/combine_ibs_pca.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_common_inter
 *
 * Optional: --subgenomes A,B,D,Others selects the PCA views.
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mod = 'test_common_inter'
params.subgenomes = 'A,B,D,Others'

process combine_ibs_pca {
    tag "combine IBS PCA: ${prefix} ${subgenome}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '16.GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(prefix), val(threshold), val(subgenome), path(pairs), path(eigenvec)

    output:
    path("*.png"), emit: plots, optional: true
    path("${prefix}.pca.${subgenome}.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${prefix}.pca.${subgenome}.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.sample.combine_ibs_validate import plot_combine_ibs_pca

    plot_combine_ibs_pca(
        "${pairs}",
        "${eigenvec}",
        "${prefix}",
        subgenome="${subgenome}",
        threshold=${threshold},
    )
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'combine_ibs_pca: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir || !params.src_dir) {
        error 'combine_ibs_pca: --user_dir and --src_dir are required'
    }
    def info = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
    def pca = "${params.output_dir}/${params.job}/process/${params.mod}/pca"
    // One cutoff per cohort: WAP_PQ/VMap3 merges at 0.995, IPK/HZAU/Watkins at 0.99.
    def cohorts = [
        ['all.sample.combine_ibs', 0.995],
        ['all.sample.combine_ibs_ipk_hzau_watkins', 0.99],
    ]
    def specs = cohorts.collectMany { prefix, threshold ->
        params.subgenomes.toString().split(',').collect { subgenome ->
            tuple(
                prefix,
                threshold,
                subgenome.trim(),
                file("${info}/${prefix}.pairs.tsv", checkIfExists: true),
                file("${pca}/${subgenome.trim()}.pca.eigenvec", checkIfExists: true),
            )
        }
    }
    combine_ibs_pca(Channel.fromList(specs))
}
