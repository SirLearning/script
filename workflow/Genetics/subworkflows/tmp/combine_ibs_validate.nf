#!/usr/bin/env nextflow
/*
 * Validate combine-name duplicate pairs with test_common_inter IBS.
 *
 * The combine-name table links WAP_PQ and VMap3 names for one germplasm.
 * Curation-excluded control slots (Jimai22 control plots, CK controls) are
 * removed before pairing. Genotype-available pairs are scored from the
 * published square PLINK and hom-SNP IBS matrices, assigned low/mid/high
 * IBS tiers, and compared with random AABBDD pairs.
 *
 * Sample depths come from the dedup snapshot. ``--combine_coverage_map``
 * (default final.ALL.taxaBamMap.txt) fills depths for genotyped taxa that the
 * dedup input list omitted, so ``coverage_min`` is defined for every pair.
 *
 * Publishes under {output_dir}/{job}/stats/{mod}/{info,plots,logs}:
 *   all.sample.combine_ibs.pairs.tsv
 *   all.sample.combine_ibs.summary.tsv
 *   all.sample.combine_ibs.vs_random.png
 *   all.sample.combine_ibs.by_dedup.png
 *   all.sample.combine_ibs.tiers.png
 *   all.sample.combine_ibs.log
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/combine_ibs_validate.nf \
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
params.combine_groups = '/data/home/tusr1/git/DBone/src/main/resources/curation_dict/curation_germplasm_combine_name_groups.tsv'
params.combine_excluded_untrusted = '/data/home/tusr1/git/DBone/src/main/resources/curation_dict/curation_germplasm_combine_sample_untrusted_excluded.tsv'
params.combine_excluded_conflict = '/data/home/tusr1/git/DBone/src/main/resources/curation_dict/curation_germplasm_combine_sample_conflict_ck_excluded.tsv'
params.combine_dedup = '/data/home/tusr1/git/script/src/resource/sample_combine/sample.dedup.checked.info.tsv'
params.combine_coverage_map = '/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap/final.ALL.taxaBamMap.txt'
params.combine_samples = null
params.combine_random_pairs = null
params.combine_ibs_dir = null
params.combine_hom_dir = null
params.combine_output_prefix = 'all.sample.combine_ibs'
params.combine_name_columns = 'wap_pq_combine_names,vmap3_combine_names'
params.combine_cohort_labels = 'wap,vmap3'
params.combine_ibs_metric = 'plink_ibs'
params.combine_tier_low_max = 0.975
params.combine_tier_high_min = 0.990

process combine_ibs_validate {
    tag "combine IBS validation: ${params.mod}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 4
    memory '16.GB'
    errorStrategy 'terminate'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(ibs_dir), val(hom_dir), path(combine_groups), path(excluded_untrusted), path(excluded_conflict), path(dedup_info), path(coverage_map), path(samples), path(random_pairs)

    output:
    path("${params.combine_output_prefix}.pairs.tsv"), emit: pairs
    path("${params.combine_output_prefix}.summary.tsv"), emit: summary
    path("${params.combine_output_prefix}.by_ploidy.tsv"), emit: by_ploidy
    path("${params.combine_output_prefix}.*.png"), emit: plots
    path("${params.combine_output_prefix}.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${params.combine_output_prefix}.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.sample.combine_ibs_validate import run_combine_ibs_validate

    run_combine_ibs_validate(
        combine_path="${combine_groups}",
        excluded_paths=["${excluded_untrusted}", "${excluded_conflict}"],
        dedup_path="${dedup_info}",
        samples_path="${samples}",
        ibs_dir="${ibs_dir}",
        hom_dir="${hom_dir}",
        random_pairs_path="${random_pairs}",
        output_prefix="${params.combine_output_prefix}",
        name_columns=[
            c.strip()
            for c in "${params.combine_name_columns}".split(',')
            if c.strip()
        ],
        cohort_labels=[
            c.strip()
            for c in "${params.combine_cohort_labels}".split(',')
            if c.strip()
        ],
        ibs_metric="${params.combine_ibs_metric}",
        tier_low_max=float("${params.combine_tier_low_max}"),
        tier_high_min=float("${params.combine_tier_high_min}"),
        coverage_map_path="${coverage_map}",
    )
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'combine_ibs_validate: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir || !params.src_dir) {
        error 'combine_ibs_validate: --user_dir and --src_dir are required'
    }

    def stats_dir = "${params.output_dir}/${params.job}/stats/${params.mod}"
    def samples_path = params.combine_samples ?: "${stats_dir}/info/all.sample.identity.ploidy_group_ibs.samples.tsv"
    def random_pairs_path = params.combine_random_pairs ?: "${stats_dir}/info/all.sample.identity.ploidy_group_ibs.pairs.tsv"
    def ibs_dir = params.combine_ibs_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}/ibs"
    def hom_dir = params.combine_hom_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}/ibs_hom_snp"

    def combine_groups = file(params.combine_groups, checkIfExists: true)
    def excluded_untrusted = file(params.combine_excluded_untrusted, checkIfExists: true)
    def excluded_conflict = file(params.combine_excluded_conflict, checkIfExists: true)
    def dedup_info = file(params.combine_dedup, checkIfExists: true)
    def coverage_map = file(params.combine_coverage_map, checkIfExists: true)
    def samples = file(samples_path, checkIfExists: true)
    def random_pairs = file(random_pairs_path, checkIfExists: true)
    if (!file(ibs_dir).exists()) {
        error "combine_ibs_validate: IBS directory not found: ${ibs_dir}"
    }
    if (!file(hom_dir).exists()) {
        error "combine_ibs_validate: hom-SNP IBS directory not found: ${hom_dir}"
    }

    ch_input = Channel.value([
        ibs_dir,
        hom_dir,
        combine_groups,
        excluded_untrusted,
        excluded_conflict,
        dedup_info,
        coverage_map,
        samples,
        random_pairs,
    ])
    combine_ibs_validate(ch_input)
}
