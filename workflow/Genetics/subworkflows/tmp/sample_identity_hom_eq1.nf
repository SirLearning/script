#!/usr/bin/env nextflow
/*
 * Extract pairs whose hom-SNP IBS is 1 and plot their PLINK IBS.
 * Reads existing process/{mod}/ibs_hom_snp and process/{mod}/ibs matrices.
 *
 *   cd /data/home/tusr1/01projects/vmap4/12stats.genome/<NNrun_ibs_hom_eq1>
 *   nextflow run .../subworkflows/tmp/sample_identity_hom_eq1.nf \
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
        error 'sample_identity_hom_eq1: --output_dir and --job are required'
    }
    if (!params.user_dir) {
        error 'sample_identity_hom_eq1: --user_dir is required (Conda stats env)'
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
                rows << tuple(mod, sg, hom_mibs, hom_id, plink_mibs, plink_id)
            }
        }
    }
    if (!rows) {
        error "sample_identity_hom_eq1: no A/B/D hom+PLINK IBS matrices under mods ${mods}"
    }
    plot_hom_ibs_one(channel.from(rows))
}

process plot_hom_ibs_one {
    tag "hom IBS=1: ${mod} ${id}"
    label 'cpus_8'
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/info", mode: 'copy', pattern: "*.eq1.pairs.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/thresholds", mode: 'copy', pattern: "*.eq1.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(mod), val(id), path(hom_mibs), path(hom_mibs_id), path(plink_mibs), path(plink_mibs_id)

    output:
    tuple val(mod), val(id),
        path("${id}.sample.identity.hom_snp.eq1.pairs.tsv"),
        path("${id}.sample.identity.hom_snp.eq1.th.tsv"),
        path("${id}.sample.identity.hom_snp.eq1.plink.dist.png", optional: true),
        path("${id}.sample.identity.hom_snp.eq1.slope.png", optional: true)

    script:
    """
    #!/usr/bin/env python
    import sys
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")
    sys.stdout = open("${id}.identity_hom_eq1.log", "w")
    sys.stderr = sys.stdout
    from genetics.genomics.sample.identity import ana_hom_ibs_one_vs_plink
    ana_hom_ibs_one_vs_plink(
        "${hom_mibs}",
        "${hom_mibs_id}",
        "${plink_mibs}",
        "${plink_mibs_id}",
        "${id}.sample.identity.hom_snp",
    )
    """
}
