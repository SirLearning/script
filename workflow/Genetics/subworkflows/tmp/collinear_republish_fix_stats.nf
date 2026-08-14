#!/usr/bin/env nextflow
/*
 * Lean republish of test_collinear process artefacts from a mistaken source mod
 * (e.g. null) via publishDir copy (work on /data, publish on /data1), then re-run
 * the three ignored stats.
 *
 * Skips huge LD matrices (*.vcor). Stats tree can be republished separately or
 * skipped if already present under process|stats/{mod}.
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/collinear_republish_fix_stats.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink \
 *     --mod test_collinear --source_mod null
 */

nextflow.enable.dsl=2

include { sample_mapping_rate_stats } from '../../modules/local/genotype/stats/stats_sample.nf'
include {
    variant_maf_miss_group_stats
    variant_popdep_missing_reg
} from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = null
params.source_mod = 'null'
params.skip_republish = false
params.skip_stats_republish = false
params.skip_stats = false

process republish_one {
    tag "${kind}:${rel_dir}/${fname}"
    // work is on /data, publish on /data1 — hardlinks cannot cross devices.
    stageInMode 'symlink'
    publishDir "${params.output_dir}/${params.job}/${kind}/${params.mod}/${rel_dir}", mode: 'copy'
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'
    errorStrategy 'terminate'

    input:
    tuple val(kind), val(rel_dir), path(fname)

    output:
    path fname, emit: published

    script:
    """
    set -euo pipefail
    test -e "${fname}"
    """
}

process gate_fix_stats_inputs {
    tag "gate fix-stats inputs"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'
    errorStrategy 'terminate'

    input:
    val _ready
    val proc_dest

    output:
    path 'ready.txt', emit: ready
    path 'smiss.tsv', emit: smiss_manifest
    path 'maf_miss.tsv', emit: maf_miss_manifest
    path 'popdep.tsv', emit: popdep_manifest

    script:
    """
    set -euo pipefail
    PROC="${proc_dest}"
    test -d "\${PROC}"
    test -f "\${PROC}/A_test.plink2.pgen"
    test -f "\${PROC}/sample/A.info.smiss"
    test -f "\${PROC}/variant/A.info.afreq"
    test -f "\${PROC}/variant/A.popdep.info.tsv"

    : > ready.txt
    : > smiss.tsv
    : > maf_miss.tsv
    : > popdep.tsv
    for sg in A B D; do
      echo -e "\${sg}\\tsub_\${sg}\\t\${PROC}/sample/\${sg}.info.smiss" >> smiss.tsv
      echo -e "\${sg}\\tsub_\${sg}\\t\${PROC}/variant/\${sg}.info.afreq\\t\${PROC}/variant/\${sg}.info.vmiss" >> maf_miss.tsv
      echo -e "\${sg}\\tsub_\${sg}\\t\${PROC}/variant/\${sg}.popdep.info.tsv\\t\${PROC}/variant/\${sg}.info.vmiss" >> popdep.tsv
    done
    """
}

def list_republish_files(String kind, String srcRoot) {
    def rootStr = srcRoot.endsWith('/') ? srcRoot[0..-2] : srcRoot
    def root = file(rootStr)
    def out = []
    root.eachFileRecurse { f ->
        if (!f.isFile()) {
            return
        }
        // Skip huge LD pair matrices (tens of GB); keep plots/stats already published.
        if (f.name.endsWith('.vcor')) {
            return
        }
        def abs = f.toAbsolutePath().toString()
        if (!abs.startsWith(rootStr + '/')) {
            return
        }
        def rel = abs.substring(rootStr.length() + 1)
        def parent = rel.contains('/') ? rel.substring(0, rel.lastIndexOf('/')) : ''
        out << tuple(kind, parent, f)
    }
    return out
}

workflow {
    if (!params.output_dir || !params.job) {
        error 'collinear_republish_fix_stats: --output_dir and --job are required'
    }
    if (!params.mod) {
        error 'collinear_republish_fix_stats: --mod is required (e.g. test_collinear)'
    }
    if (!params.source_mod) {
        error 'collinear_republish_fix_stats: --source_mod is required (e.g. null)'
    }
    if (params.source_mod == params.mod) {
        error "collinear_republish_fix_stats: source_mod and mod must differ (got ${params.mod})"
    }

    def proc_dest = "${params.output_dir}/${params.job}/process/${params.mod}"
    def ch_ready

    if (!params.skip_republish) {
        def src_proc = "${params.output_dir}/${params.job}/process/${params.source_mod}"
        if (!file(src_proc).exists()) {
            error "Missing source process tree: ${src_proc}"
        }
        def items = list_republish_files('process', src_proc)
        if (!params.skip_stats_republish) {
            def src_stats = "${params.output_dir}/${params.job}/stats/${params.source_mod}"
            if (!file(src_stats).exists()) {
                error "Missing source stats tree: ${src_stats}"
            }
            items.addAll(list_republish_files('stats', src_stats))
        }
        log.info "Lean republish ${items.size()} files from ${params.source_mod} -> ${params.mod} (copy; skip *.vcor)"
        ch_ready = republish_one(Channel.fromList(items)).published.collect()
    } else {
        ch_ready = Channel.value(['skip_republish'])
    }

    if (!params.skip_stats) {
        def idxstats = file("${params.output_dir}/meta_data/vmap4_v1_idxstat_summary.txt")
        def group_file = file("${params.output_dir}/meta_data/sample_group.txt")
        if (!idxstats.exists()) {
            error "Missing idxstat summary: ${idxstats}"
        }
        if (!group_file.exists()) {
            error "Missing sample group file: ${group_file}"
        }

        def gated = gate_fix_stats_inputs(ch_ready, Channel.value(proc_dest))

        def ch_smiss = gated.smiss_manifest
            .splitCsv(sep: '\t')
            .map { row -> tuple(row[0], row[1], file(row[2])) }
        sample_mapping_rate_stats(idxstats, group_file, ch_smiss)

        def ch_maf_miss = gated.maf_miss_manifest
            .splitCsv(sep: '\t')
            .map { row -> tuple(row[0], row[1], file(row[2]), file(row[3])) }
        variant_maf_miss_group_stats(ch_maf_miss)

        def ch_popdep = gated.popdep_manifest
            .splitCsv(sep: '\t')
            .map { row -> tuple(row[0], row[1], file(row[2]), file(row[3])) }
        variant_popdep_missing_reg(ch_popdep)
    }
}
