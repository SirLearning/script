#!/usr/bin/env nextflow
/*
 * Within-ploidy PLINK and hom-SNP IBS pair histograms.
 *
 * Previous ploidy is sample_group (AA/SS/DD/AABB/AABBDD). Breachr full-cohort
 * segments correct misassigned labels: drop a claimed subgenome when all 7
 * homoeologues are subgenome_absent or subgenome_reference_mismatch.
 * Expected nulls (AB D-null, diploid off-subgenomes) are not corrections.
 * Pairs are drawn inside each ploidy class on that subgenome matrix:
 * A → AA/AABB/AABBDD; B → SS/AABB/AABBDD; D → DD/AABBDD.
 *
 *   cd /data/home/tusr1/01projects/vmap4/12stats.genome/<NNrun_ibs_ploidy_dist>
 *   nextflow run .../subworkflows/tmp/sample_identity_ploidy_ibs_dist.nf \
 *     -c .../nextflow.config \
 *     --user_dir /data/home/tusr1 \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mods test_thin,test_common_inter
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mods = 'test_thin,test_common_inter'
params.sample_group = null
params.breachr_segments = null
params.n_pairs_per_class = 200000

workflow {
    if (!params.output_dir || !params.job) {
        error 'sample_identity_ploidy_ibs_dist: --output_dir and --job are required'
    }
    if (!params.user_dir) {
        error 'sample_identity_ploidy_ibs_dist: --user_dir is required (Conda stats env)'
    }

    def groupFile = params.sample_group
        ? file(params.sample_group)
        : file("${params.output_dir}/meta_data/sample_group.txt")
    def breachrFile = params.breachr_segments
        ? file(params.breachr_segments)
        : file("${params.output_dir}/benchmark/stats/aneuploidy_check/info/breachr_cohort/breachr_segments_all.tsv")
    if (!groupFile.exists()) {
        error "sample_identity_ploidy_ibs_dist: sample_group not found: ${groupFile}"
    }
    if (!breachrFile.exists()) {
        error "sample_identity_ploidy_ibs_dist: breachr segments not found: ${breachrFile}"
    }

    def mods = params.mods.toString().split(',').collect { it.trim() }.findAll { it }
    def rows = []
    mods.each { mod ->
        def proc = "${params.output_dir}/${params.job}/process/${mod}"
        ['A', 'B', 'D'].each { sg ->
            def hom_mibs = file("${proc}/ibs_hom_snp/${sg}.ibs_hom_snp.mibs")
            def hom_id = file("${proc}/ibs_hom_snp/${sg}.ibs_hom_snp.mibs.id")
            def plink_mibs = file("${proc}/ibs/${sg}.ibs.mibs")
            def plink_id = file("${proc}/ibs/${sg}.ibs.mibs.id")
            if (hom_mibs.exists() && hom_id.exists() && plink_mibs.exists() && plink_id.exists()) {
                rows << tuple(mod, sg, hom_mibs, hom_id, plink_mibs, plink_id, groupFile, breachrFile)
            }
        }
    }
    if (!rows) {
        error "sample_identity_ploidy_ibs_dist: no A/B/D hom+PLINK IBS matrices under mods ${mods}"
    }
    plot_ploidy_ibs_dist(channel.from(rows))
}

process plot_ploidy_ibs_dist {
    tag "ploidy IBS: ${mod} ${id}"
    cpus 2
    memory '16 GB'
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/info", mode: 'copy', pattern: "*.{samples,pairs,info}.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(mod), val(id), path(hom_mibs), path(hom_mibs_id), path(plink_mibs), path(plink_mibs_id), path(group_file), path(breachr_segments)

    output:
    tuple val(mod), val(id),
        path("${id}.sample.identity.ploidy_ibs.samples.tsv"),
        path("${id}.sample.identity.ploidy_ibs.pairs.tsv"),
        path("${id}.sample.identity.ploidy_ibs.info.tsv"),
        path("${id}.sample.identity.ploidy_ibs.th.tsv"),
        path("${id}.sample.identity.ploidy_ibs.plink.dist.png"),
        path("${id}.sample.identity.ploidy_ibs.hom.dist.png"),
        path("${id}.sample.identity.ploidy_ibs.plink.dist_log.png"),
        path("${id}.sample.identity.ploidy_ibs.hom.dist_log.png"),
        path("${id}.sample.identity.ploidy_ibs.plink.density.png"),
        path("${id}.sample.identity.ploidy_ibs.hom.density.png"),
        path("${id}.identity_ploidy_ibs.log")

    script:
    """
    #!/usr/bin/env python
    import sys
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")
    sys.stdout = open("${id}.identity_ploidy_ibs.log", "w")
    sys.stderr = sys.stdout
    from genetics.genomics.sample.identity import ana_ibs_ploidy_pair_dist
    ana_ibs_ploidy_pair_dist(
        "${plink_mibs}",
        "${plink_mibs_id}",
        "${hom_mibs}",
        "${hom_mibs_id}",
        "${group_file}",
        "${breachr_segments}",
        "${id}.sample.identity.ploidy_ibs",
        subgenome="${id}",
        n_pairs_per_class=${params.n_pairs_per_class},
    )
    """
}
