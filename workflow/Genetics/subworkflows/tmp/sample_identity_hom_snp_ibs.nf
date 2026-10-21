#!/usr/bin/env nextflow
/*
 * Extra IBS version: biallelic SNPs only, score a site only when both
 * samples are homozygous. Does not replace process/{mod}/ibs/*.ibs.mibs.
 *
 * Compute writes process/{mod}/ibs_hom_snp/*.ibs_hom_snp.mibs and is skipped
 * when that matrix already exists. Plots always republish from the matrix.
 *
 *   cd /data/home/tusr1/01projects/vmap4/12stats.genome/<NNrun_ibs_hom_snp>
 *   nextflow run .../subworkflows/tmp/sample_identity_hom_snp_ibs.nf \
 *     -c .../nextflow.config \
 *     --user_dir /data/home/tusr1 \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mods test_thin,test_common_inter
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mods = 'test_thin,test_common_inter'
params.hom_snp_min_ibs = 0.99
params.hom_snp_compare_ibs = 0.97
params.hom_snp_min_nobs = 0
params.hom_snp_chunk_sites = 2048
params.hom_snp_neighbor_k = 5
params.hom_snp_threads = 4
params.hom_snp_scatter_n = 200000

workflow {
    if (!params.output_dir || !params.job) {
        error 'sample_identity_hom_snp_ibs: --output_dir and --job are required'
    }
    if (!params.user_dir) {
        error 'sample_identity_hom_snp_ibs: --user_dir is required (Conda stats env)'
    }

    def mods = params.mods.toString().split(',').collect { it.trim() }.findAll { it }
    def need = []
    def ready = []
    mods.each { mod ->
        def proc = "${params.output_dir}/${params.job}/process/${mod}"
        ['A', 'B', 'D'].each { sg ->
            def bed = file("${proc}/${sg}_test.plink.bed")
            def bim = file("${proc}/${sg}_test.plink.bim")
            def fam = file("${proc}/${sg}_test.plink.fam")
            def plink_mibs = file("${proc}/ibs/${sg}.ibs.mibs")
            def plink_id = file("${proc}/ibs/${sg}.ibs.mibs.id")
            def hom_mibs = file("${proc}/ibs_hom_snp/${sg}.ibs_hom_snp.mibs")
            def hom_id = file("${proc}/ibs_hom_snp/${sg}.ibs_hom_snp.mibs.id")
            if (!(bed.exists() && bim.exists() && fam.exists() && plink_mibs.exists() && plink_id.exists())) {
                return
            }
            if (hom_mibs.exists() && hom_id.exists()) {
                ready << tuple(mod, sg, "sub_${sg}", hom_mibs, hom_id, plink_mibs, plink_id)
            } else {
                need << tuple(mod, sg, "sub_${sg}", "${sg}_test.plink", bed, bim, fam, plink_mibs, plink_id)
            }
        }
    }
    if (!need && !ready) {
        error "sample_identity_hom_snp_ibs: no A/B/D bfile+IBS inputs under mods ${mods}"
    }

    def ch_computed = need
        ? compute_hom_snp_ibs(channel.from(need))
        : channel.empty()
    def ch_ready = ready
        ? channel.from(ready)
        : channel.empty()
    plot_hom_snp_ibs(ch_ready.mix(ch_computed))
}

process compute_hom_snp_ibs {
    tag "hom-snp matrix: ${mod} ${id}"
    label 'cpus_8'
    maxForks 1
    errorStrategy 'retry'
    maxRetries 2
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/process/${mod}/ibs_hom_snp", mode: 'copy', pattern: "*.{mibs,mibs.id}"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(mod), val(id), val(chr), val(prefix), path(bed), path(bim), path(fam), path(plink_mibs), path(plink_mibs_id)

    output:
    tuple val(mod), val(id), val(chr), path("${id}.ibs_hom_snp.mibs"), path("${id}.ibs_hom_snp.mibs.id"), path(plink_mibs), path(plink_mibs_id)

    script:
    def nthread = [params.hom_snp_threads.toInteger(), task.attempt == 1 ? 4 : 1].min()
    """
    #!/usr/bin/env python
    import os
    import sys
    os.environ["OMP_NUM_THREADS"] = "${nthread}"
    os.environ["OPENBLAS_NUM_THREADS"] = "${nthread}"
    os.environ["MKL_NUM_THREADS"] = "${nthread}"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    sys.stdout = open("${chr}.identity_hom_snp_compute.log", "w")
    sys.stderr = sys.stdout
    from genetics.genomics.sample.identity import write_hom_snp_ibs_matrix
    write_hom_snp_ibs_matrix("${prefix}", "${id}.ibs_hom_snp", chunk_sites=${params.hom_snp_chunk_sites})
    """
}

process plot_hom_snp_ibs {
    tag "hom-snp plot: ${mod} ${id}"
    label 'cpus_8'
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/info", mode: 'copy', pattern: "*.{pairs.tsv,clusters.tsv,neighbors.tsv,neighbors.plink.tsv,neighbors.shift.tsv}"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(mod), val(id), val(chr), path(hom_mibs), path(hom_mibs_id), path(plink_mibs), path(plink_mibs_id)

    output:
    tuple val(mod), val(id),
        path("${id}.sample.identity.hom_snp.pairs.tsv"),
        path("${id}.sample.identity.hom_snp.clusters.tsv"),
        path("${id}.sample.identity.hom_snp.neighbors.tsv"),
        path("${id}.sample.identity.hom_snp.neighbors.plink.tsv"),
        path("${id}.sample.identity.hom_snp.neighbors.shift.tsv"),
        path("${id}.sample.identity.hom_snp.vs_plink.pairs.tsv"),
        path("${id}.sample.identity.hom_snp.th.tsv"),
        path("${id}.sample.identity.hom_snp.dist.png", optional: true),
        path("${id}.sample.identity.hom_snp.vs_plink.scatter.png", optional: true),
        path("${id}.sample.identity.hom_snp.vs_plink.delta.png", optional: true)

    script:
    """
    #!/usr/bin/env python
    import os
    import sys
    os.environ["OMP_NUM_THREADS"] = "4"
    os.environ["OPENBLAS_NUM_THREADS"] = "4"
    os.environ["MKL_NUM_THREADS"] = "4"
    sys.stdout = open("${chr}.identity_hom_snp_plot.log", "w")
    sys.stderr = sys.stdout
    from genetics.genomics.sample.identity import ana_hom_snp_ibs_from_mibs
    ana_hom_snp_ibs_from_mibs(
        "${hom_mibs}",
        "${hom_mibs_id}",
        "${id}.sample.identity.hom_snp",
        plink_mibs="${plink_mibs}",
        plink_mibs_id="${plink_mibs_id}",
        min_pair_ibs=${params.hom_snp_min_ibs},
        min_compare_ibs=${params.hom_snp_compare_ibs},
        min_nobs=${params.hom_snp_min_nobs},
        neighbor_k=${params.hom_snp_neighbor_k},
        scatter_n=${params.hom_snp_scatter_n},
    )
    """
}
