#!/usr/bin/env nextflow
/*
 * Within-ploidy PLINK and hom-SNP IBS using each class's own subgenomes.
 *
 * AA vs AA on A; SS vs SS on B; DD vs DD on D. AABB pairs take the mean of
 * A and B IBS; AABBDD pairs take the mean of A, B, and D. Two IBS types
 * (PLINK and hom-SNP) are scored that way, then plotted per class.
 * Ploidy labels: sample_group plus breachr full-cohort corrections.
 *
 *   cd /data/home/tusr1/01projects/vmap4/12stats.genome/<NNrun_ibs_ploidy_group>
 *   nextflow run .../subworkflows/tmp/sample_identity_ploidy_group_ibs_dist.nf \
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
        error 'sample_identity_ploidy_group_ibs_dist: --output_dir and --job are required'
    }
    if (!params.user_dir) {
        error 'sample_identity_ploidy_group_ibs_dist: --user_dir is required (Conda stats env)'
    }

    def groupFile = params.sample_group
        ? file(params.sample_group)
        : file("${params.output_dir}/meta_data/sample_group.txt")
    def breachrFile = params.breachr_segments
        ? file(params.breachr_segments)
        : file("${params.output_dir}/benchmark/stats/aneuploidy_check/info/breachr_cohort/breachr_segments_all.tsv")
    if (!groupFile.exists()) {
        error "sample_identity_ploidy_group_ibs_dist: sample_group not found: ${groupFile}"
    }
    if (!breachrFile.exists()) {
        error "sample_identity_ploidy_group_ibs_dist: breachr segments not found: ${breachrFile}"
    }

    def mods = params.mods.toString().split(',').collect { it.trim() }.findAll { it }
    def rows = []
    mods.each { mod ->
        def proc = "${params.output_dir}/${params.job}/process/${mod}"
        def paths = [:]
        def ok = true
        ['A', 'B', 'D'].each { sg ->
            def hom_mibs = file("${proc}/ibs_hom_snp/${sg}.ibs_hom_snp.mibs")
            def hom_id = file("${proc}/ibs_hom_snp/${sg}.ibs_hom_snp.mibs.id")
            def plink_mibs = file("${proc}/ibs/${sg}.ibs.mibs")
            def plink_id = file("${proc}/ibs/${sg}.ibs.mibs.id")
            if (!(hom_mibs.exists() && hom_id.exists() && plink_mibs.exists() && plink_id.exists())) {
                ok = false
            } else {
                paths[sg] = [plink_mibs, plink_id, hom_mibs, hom_id]
            }
        }
        if (ok) {
            rows << tuple(
                mod,
                paths['A'][0], paths['A'][1], paths['A'][2], paths['A'][3],
                paths['B'][0], paths['B'][1], paths['B'][2], paths['B'][3],
                paths['D'][0], paths['D'][1], paths['D'][2], paths['D'][3],
                groupFile,
                breachrFile,
            )
        }
    }
    if (!rows) {
        error "sample_identity_ploidy_group_ibs_dist: no complete A/B/D hom+PLINK IBS set under mods ${mods}"
    }
    plot_ploidy_group_ibs_dist(channel.from(rows))
}

process plot_ploidy_group_ibs_dist {
    tag "ploidy group IBS: ${mod}"
    cpus 2
    memory '32 GB'
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/info", mode: 'copy', pattern: "*.{samples,pairs,info}.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(mod),
        path(a_plink_mibs), path(a_plink_mibs_id), path(a_hom_mibs), path(a_hom_mibs_id),
        path(b_plink_mibs), path(b_plink_mibs_id), path(b_hom_mibs), path(b_hom_mibs_id),
        path(d_plink_mibs), path(d_plink_mibs_id), path(d_hom_mibs), path(d_hom_mibs_id),
        path(group_file), path(breachr_segments)

    output:
    tuple val(mod),
        path("all.sample.identity.ploidy_group_ibs.samples.tsv"),
        path("all.sample.identity.ploidy_group_ibs.pairs.tsv"),
        path("all.sample.identity.ploidy_group_ibs.info.tsv"),
        path("all.sample.identity.ploidy_group_ibs.th.tsv"),
        path("all.sample.identity.ploidy_group_ibs.plink.dist.png"),
        path("all.sample.identity.ploidy_group_ibs.hom.dist.png"),
        path("all.sample.identity.ploidy_group_ibs.plink.dist_log.png"),
        path("all.sample.identity.ploidy_group_ibs.hom.dist_log.png"),
        path("all.sample.identity.ploidy_group_ibs.plink.density.png"),
        path("all.sample.identity.ploidy_group_ibs.hom.density.png"),
        path("all.identity_ploidy_group_ibs.log")

    script:
    """
    #!/usr/bin/env python
    import sys
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")
    sys.stdout = open("all.identity_ploidy_group_ibs.log", "w")
    sys.stderr = sys.stdout
    from genetics.genomics.sample.identity import ana_ibs_ploidy_group_pair_dist
    ana_ibs_ploidy_group_pair_dist(
        "${a_plink_mibs}",
        "${a_plink_mibs_id}",
        "${a_hom_mibs}",
        "${a_hom_mibs_id}",
        "${b_plink_mibs}",
        "${b_plink_mibs_id}",
        "${b_hom_mibs}",
        "${b_hom_mibs_id}",
        "${d_plink_mibs}",
        "${d_plink_mibs_id}",
        "${d_hom_mibs}",
        "${d_hom_mibs_id}",
        "${group_file}",
        "${breachr_segments}",
        "all.sample.identity.ploidy_group_ibs",
        n_pairs_per_class=${params.n_pairs_per_class},
    )
    """
}
