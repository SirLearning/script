#!/usr/bin/env nextflow
/*
 * Heterozygosity change of merged germplasms: before merge vs after merge.
 *
 * "Before" is the pre-merge pfile of the same chromosome (every source taxon as
 * its own sample, main_raw/chr006.plink2 published as chr006_all). "After" is
 * the scan3 pfile for the same library (one sample per final taxon, chr006_final).
 * PLINK2 --het runs on both, then each merged germplasm's heterozygosity is
 * compared with the median of its own sources, which removes the between-sample
 * coverage spread that dominates a merged-vs-unmerged comparison.
 *
 * Publishes under {output_dir}/{job}/stats/{mod}/{info,plots,logs}:
 *   <prefix>.merged.before_after.info.tsv           one row per germplasm
 *   <prefix>.merged.before_after.summary.info.tsv   overall + per source count
 *   <prefix>.merged.before_after.scatter.png        source median vs merged, y=x
 *   <prefix>.merged.before_after.delta.png          distribution of the change
 *   <prefix>.merged.before_after.log
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/scan3_merge_het_before_after.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_merge_chr6 \
 *     --pfile_all <chr006_all> --pfile_final <chr006_final> --merge_map <final.D.taxaMergeMap.txt>
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mod = 'test_merge_chr6'
params.pfile_all = null
params.pfile_final = null
params.merge_map = null
params.output_prefix = 'chr006'
params.title_prefix = ''

process scan3_merge_het_before_after {
    tag "merged heterozygosity before vs after"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 8
    memory '120.GB'
    errorStrategy 'terminate'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple path(all_pgen), path(all_pvar), path(all_psam), \
          path(final_pgen), path(final_pvar), path(final_psam), path(merge_map)

    output:
    path("*.info.tsv"), emit: info
    path("*.png"), emit: plots
    path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    plink2 --pfile chr006_all --het --threads ${task.cpus} --memory 110000 \
        --out het_all > het_all.log 2>&1
    plink2 --pfile chr006_final --het --threads ${task.cpus} --memory 110000 \
        --out het_final > het_final.log 2>&1

    cat > run_check.py <<'PY'
    import sys

    sys.path.insert(0, sys.argv[1])
    from genetics.genomics.sample.merge_het_check import run_merged_before_after_check

    run_merged_before_after_check(
        "het_all.het",
        "het_final.het",
        sys.argv[2],
        sys.argv[3],
        info_dir=".",
        plot_dir=".",
        title_prefix=sys.argv[4],
    )
    PY

    python run_check.py "${params.src_dir}/python" "${merge_map}" \
        "${params.output_prefix}" "${params.title_prefix}" \
        > "${params.output_prefix}.merged.before_after.log" 2>&1
    cat het_all.log >> "${params.output_prefix}.merged.before_after.log"
    cat het_final.log >> "${params.output_prefix}.merged.before_after.log"
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'scan3_merge_het_before_after: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir || !params.src_dir) {
        error 'scan3_merge_het_before_after: --user_dir and --src_dir are required'
    }
    if (!params.pfile_all || !params.pfile_final || !params.merge_map) {
        error 'scan3_merge_het_before_after: --pfile_all, --pfile_final, and --merge_map are required'
    }
    scan3_merge_het_before_after(
        Channel.value([
            file("${params.pfile_all}.pgen", checkIfExists: true),
            file("${params.pfile_all}.pvar", checkIfExists: true),
            file("${params.pfile_all}.psam", checkIfExists: true),
            file("${params.pfile_final}.pgen", checkIfExists: true),
            file("${params.pfile_final}.pvar", checkIfExists: true),
            file("${params.pfile_final}.psam", checkIfExists: true),
            file(params.merge_map, checkIfExists: true),
        ])
    )
}
