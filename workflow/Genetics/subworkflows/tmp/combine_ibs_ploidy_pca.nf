#!/usr/bin/env nextflow
/*
 * PCA per ploidy group (AA, AABB, AABBDD, DD) over the whole genome, not split
 * by subgenome. The per-chromosome common-inter pfiles are merged once, then
 * each group gets its own LD-pruned PCA (same recipe as plink1_pca_ld_prune:
 * --maf 0.05, --indep-pairwise 50 5 0.2, --pca approx 10) and the combine-group
 * wiring is drawn on top.
 *
 *   nextflow run .../subworkflows/tmp/combine_ibs_ploidy_pca.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_common_inter \
 *     --ploidies AA,AABB,AABBDD,DD
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mod = 'test_common_inter'
params.ploidies = 'AA,AABB,AABBDD,DD'
// Which column of the samples table defines a group: `ploidy` is the measured
// (BreachR) call that the IBS analysis also groups by; `nominal_ploidy` is the
// assignment the sample was collected as. They disagree for a few hundred
// samples (reason=breachr_ploidy_mismatch), e.g. 8 WAP `P####` accessions that
// measure as AABB and land in that panel below.
params.ploidy_column = 'ploidy'
params.pca_ld_window = 50
params.pca_ld_step = 5
params.pca_ld_r2 = 0.2
params.pca_maf = 0.05
params.pc_num = 10

process merge_inter_pfile {
    tag 'merge common-inter pfiles'
    label 'cpus_8'
    conda "${params.user_dir}/miniconda3/envs/stats"
    memory '32.GB'

    output:
    tuple val('all'), path('merged.pgen'), path('merged.psam'), path('merged.pvar'), emit: pfile

    script:
    """
    set -euo pipefail
    ls ${params.output_dir}/${params.job}/process/${params.mod}/chr*.common.inter.pgen \\
        | sed 's/\\.pgen\$//' > pfiles.txt
    wc -l pfiles.txt

    plink2 --pmerge-list pfiles.txt \\
        --make-pgen \\
        --allow-extra-chr \\
        --threads ${task.cpus} \\
        --out merged
    """
}

process ploidy_pca {
    tag "ploidy pca: ${ploidy}"
    label 'cpus_8'
    conda "${params.user_dir}/miniconda3/envs/stats"
    memory '48.GB'
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/pca_ploidy", mode: 'copy', pattern: "*.{eigenvec,eigenval}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/pca_ploidy/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), path(pgen), path(psam), path(pvar), val(ploidy), path(samples)

    output:
    tuple val(ploidy), path("${ploidy}.pca.eigenvec"), path("${ploidy}.pca.eigenval"), emit: pca
    path("${ploidy}.ploidy_pca.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > ${ploidy}.ploidy_pca.log 2>&1

    # The merged psam carries IID only, so --keep needs a single-column file.
    col=\$(head -1 ${samples} | tr '\\t' '\\n' \\
        | grep -n -x '${params.ploidy_column}' | cut -d: -f1)
    awk -v c="\$col" -v p="${ploidy}" 'NR > 1 && \$c == p {print \$1}' ${samples} > keep.txt
    wc -l keep.txt

    plink2 --pfile merged \\
        --keep keep.txt \\
        --maf ${params.pca_maf} \\
        --indep-pairwise ${params.pca_ld_window} ${params.pca_ld_step} ${params.pca_ld_r2} \\
        --allow-extra-chr \\
        --threads ${task.cpus} \\
        --out ${ploidy}.prune

    plink2 --pfile merged \\
        --keep keep.txt \\
        --extract ${ploidy}.prune.prune.in \\
        --pca approx ${params.pc_num} \\
        --allow-extra-chr \\
        --threads ${task.cpus} \\
        --out ${ploidy}.pca
    """
}

process combine_ibs_ploidy_plot {
    tag "combine IBS ploidy PCA: ${prefix} ${ploidy}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '16.GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(prefix), val(threshold), val(ploidy), path(pairs), path(eigenvec)

    output:
    path("*.png"), emit: plots, optional: true
    path("${prefix}.pca_ploidy.${ploidy}.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${prefix}.pca_ploidy.${ploidy}.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.sample.combine_ibs_validate import plot_combine_ibs_pca

    for categories in (None, ["pass"], ["fail"]):
        plot_combine_ibs_pca(
            "${pairs}",
            "${eigenvec}",
            "${prefix}",
            subgenome="${ploidy}",
            threshold=${threshold},
            suffix="pca_ploidy",
            allow_empty_edges=True,
            edge_categories=categories,
        )
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'combine_ibs_ploidy_pca: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir || !params.src_dir) {
        error 'combine_ibs_ploidy_pca: --user_dir and --src_dir are required'
    }
    def info = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
    def merged = merge_inter_pfile()
    def samples = file(
        "${info}/all.sample.identity.ploidy_group_ibs.samples.tsv",
        checkIfExists: true,
    )
    def ploidies = Channel.fromList(
        params.ploidies.toString().split(',').collect { it.trim() }
    ).map { ploidy -> tuple(ploidy, samples) }
    def pca = ploidy_pca(
        merged.pfile.combine(ploidies).map { id, pgen, psam, pvar, ploidy, sample_file ->
            tuple(id, pgen, psam, pvar, ploidy, sample_file)
        }
    ).pca

    // Cutoff per cohort, plus the populations that cohort actually spans:
    // IPK/HZAU/Watkins accessions are all AABBDD, so it gets one panel, while
    // WAP_PQ/VMap3 keeps the full set (AA included, even though the A_ samples
    // sit in singleton rows and therefore carry no edges).
    def cohorts = Channel.of(
        tuple('all.sample.combine_ibs', 0.995 as double, 'AA,AABB,AABBDD,DD'),
        tuple('all.sample.combine_ibs_ipk_hzau_watkins', 0.99 as double, 'AABBDD'),
    )
    def specs = pca
        .combine(cohorts)
        .filter { item ->
            item[5].toString().split(',').collect { it.trim() }.contains(item[0])
        }
        .map { ploidy, eigenvec, eigenval, prefix, threshold, allowed ->
            tuple(
                prefix,
                threshold,
                ploidy,
                file("${info}/${prefix}.pairs.tsv", checkIfExists: true),
                eigenvec,
            )
        }
    combine_ibs_ploidy_plot(specs)
}
