#!/usr/bin/env nextflow
/*
 * LD-pruned PCA with every combine group wired up.
 *
 * The published test_common_inter PCA used the merged pfile directly, so
 * correlated blocks could dominate PC1/PC2. This entry reuses the pipeline's
 * LD-pruned PCA (process/<mod>/pca_ld/, --indep-pairwise 50 5 0.2, the same
 * recipe as plink1_pca_ld_prune) and draws one figure per cohort x subgenome.
 *
 * A/B/D pruned eigenvecs already exist and are reused. `Others` is off by
 * default: with --indep-pairwise 50 5 0.2 its 104k variants collapse to 24,
 * below the 22-per-PC minimum of --pca approx, so the pipeline also skips it
 * in its own stats path. Set --prune_others true (and loosen pca_ld_r2) to try.
 *
 *   nextflow run .../subworkflows/tmp/combine_ibs_pca_ld.nf \
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
params.pca_ld_window = 50
params.pca_ld_step = 5
params.pca_ld_r2 = 0.2
params.pc_num = 10
params.subgenomes = 'A,B,D'
params.prune_others = false

process others_pca_ld {
    tag "pca ld prune: ${prefix}"
    label 'cpus_8'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/pca_ld", mode: 'copy', pattern: "*.{eigenvec,eigenval}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/pca_ld/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(prefix), path("${prefix}.pca.ld.eigenvec"), path("${prefix}.pca.ld.eigenval"), emit: pca

    script:
    """
    set -euo pipefail
    exec > ${prefix}.plink2_pca_ld.log 2>&1

    plink2 --pfile ${prefix} \\
        --allow-extra-chr \\
        --indep-pairwise ${params.pca_ld_window} ${params.pca_ld_step} ${params.pca_ld_r2} \\
        --threads ${task.cpus} \\
        --out ${prefix}.pca.ld.prune

    plink2 --pfile ${prefix} \\
        --extract ${prefix}.pca.ld.prune.prune.in \\
        --mind 0.9999 \\
        --pca approx ${params.pc_num} \\
        --allow-extra-chr \\
        --threads ${task.cpus} \\
        --out ${prefix}.pca.ld
    """
}

process combine_ibs_pca_plot {
    tag "combine IBS PCA (LD): ${prefix} ${subgenome}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '16.GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(prefix), val(threshold), val(subgenome), path(pairs), path(eigenvec)

    output:
    path("*.png"), emit: plots, optional: true
    path("${prefix}.pca_ld.${subgenome}.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${prefix}.pca_ld.${subgenome}.log", "w", buffering=1)
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
        suffix="pca_ld",
    )
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'combine_ibs_pca_ld: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir || !params.src_dir) {
        error 'combine_ibs_pca_ld: --user_dir and --src_dir are required'
    }
    def info = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
    def process_dir = "${params.output_dir}/${params.job}/process/${params.mod}"
    def pca_ld = "${process_dir}/pca_ld"

    def others = params.prune_others
        ? others_pca_ld(
            Channel.of(
                tuple(
                    'Others_test.plink2',
                    file("${process_dir}/Others_test.plink2.pgen", checkIfExists: true),
                    file("${process_dir}/Others_test.plink2.psam", checkIfExists: true),
                    file("${process_dir}/Others_test.plink2.pvar", checkIfExists: true),
                )
            )
        ).pca.map { prefix, eigenvec, eigenval -> tuple('Others', eigenvec) }
        : Channel.empty()

    def on_disk = Channel.fromList(
        params.subgenomes.toString().split(',').collect { subgenome ->
            tuple(
                subgenome.trim(),
                file("${pca_ld}/${subgenome.trim()}.pca.ld.eigenvec", checkIfExists: true),
            )
        }
    )
    def eigenvectors = on_disk.mix(others)

    // WAP_PQ/VMap3 merges at 0.995, IPK/HZAU/Watkins at 0.99.
    def cohorts = Channel.of(
        tuple('all.sample.combine_ibs', 0.995 as double),
        tuple('all.sample.combine_ibs_ipk_hzau_watkins', 0.99 as double),
    )
    def specs = eigenvectors
        .combine(cohorts)
        .map { subgenome, eigenvec, prefix, threshold ->
            tuple(
                prefix,
                threshold,
                subgenome,
                file("${info}/${prefix}.pairs.tsv", checkIfExists: true),
                eigenvec,
            )
        }
    combine_ibs_pca_plot(specs)
}
