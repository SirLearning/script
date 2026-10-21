#!/usr/bin/env nextflow
/*
 * Rebuild the WAP test pfiles at --mac 1 without touching the MAC≥2 libraries.
 *
 * Stages (one --wap_mac1_stage per nextflow invocation, because publishDir
 * follows params.mod):
 *   arm     process/test_wap_mac1 from main_raw: --keep WAP, --mac 1, --thin 0.005
 *   refchr  process/test_wap_mac1_refchr: recode the arm merge onto linear chroms,
 *           basic info, and plink2 --freq counts (assoc.acount). No LD.
 *   exon    process/test_wap_mac1_exon from main_raw: --keep WAP, exon BED, --mac 1,
 *           no thin. A/B/D only. Basic info and assoc.acount. No LD.
 *
 * Launch each stage from a vmap4 run folder after the previous stage has published.
 */

nextflow.enable.dsl=2

include {
    prepare_test_wap_keep
    subsampling_wap_pfile_for_test
    extract_wap_exon_pfile
    gff_feature_bed
    merge_subgenome_test_pfile
} from '../../modules/local/genotype/processor/processor_test.nf'
include {
    mk_plink_basic_info
    recode_subgenome_pfile_ref_chromosomes
} from '../../modules/local/genotype/processor/processor_plink2.nf'
include {
    hasMergedSubgenomeTestPfiles
    listMergedSubgenomeTestPfileTuples
} from '../../modules/local/infra/infra_plink_reuse.nf'

params.job = 'test_plink'
params.mod = null
params.process_dir = null
params.wap_mac1_stage = null
params.test_wap_min_mac = 1

process plink_assoc_acount {
    tag "assoc acount: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/variant", mode: 'copy', pattern: "*.assoc.acount"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.acount.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_16'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), path("${id}.assoc.acount"), emit: acount
    path "*.acount.log", emit: logs

    script:
    """
    set -euo pipefail
    plink2 --pfile ${prefix} \\
        --allow-extra-chr \\
        --freq counts \\
        --threads ${task.cpus} \\
        --out ${id}.assoc \\
        > ${id}.acount.log 2>&1
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod || !params.process_dir) {
        error 'wap_mac1_testsets: --output_dir, --job, --mod, and --process_dir are required'
    }
    if (!params.user_dir) {
        error 'wap_mac1_testsets: --user_dir is required (conda stats)'
    }
    def stage = params.wap_mac1_stage
    if (!(stage in ['arm', 'refchr', 'exon'])) {
        error "wap_mac1_testsets: --wap_mac1_stage must be arm, refchr, or exon (got ${stage})"
    }

    def map_subgenome = { id, chr, prefix, pgen, psam, pvar ->
        def subgenome = "Others"
        def c = chr.toString()
        if (["1","2","7","8","13","14","19","20","25","26","31","32","37","38"].contains(c)) subgenome = "A"
        else if (["3","4","9","10","15","16","21","22","27","28","33","34","39","40"].contains(c)) subgenome = "B"
        else if (["5","6","11","12","17","18","23","24","29","30","35","36","41","42"].contains(c)) subgenome = "D"
        return [ subgenome, id, chr, prefix, pgen, psam, pvar ]
    }

    if (stage == 'refchr') {
        if (!hasMergedSubgenomeTestPfiles(params.process_dir)) {
            error "wap_mac1_testsets refchr: missing *_test.plink2 under ${params.process_dir}"
        }
        def source = Channel.from(listMergedSubgenomeTestPfileTuples(params.process_dir))
        def recoded = recode_subgenome_pfile_ref_chromosomes(source)
        mk_plink_basic_info(recoded.pfile)
        plink_assoc_acount(recoded.pfile)
    }
    else {
        def group_file = file("${params.output_dir}/meta_data/sample_group.txt", checkIfExists: true)
        def source = Channel
            .fromPath("${params.process_dir}/chr*.plink2.pvar", checkIfExists: true)
            .map { pvar ->
                def id = pvar.getName().replace('.plink2.pvar', '')
                def chr = String.valueOf(Integer.parseInt(id.substring(3)))
                def prefix = pvar.toString().replaceAll(/\.pvar$/, '')
                tuple(id, chr, prefix, file("${prefix}.pgen"), file("${prefix}.psam"), pvar)
            }
        def keep = prepare_test_wap_keep(group_file)
        def extracted
        if (stage == 'arm') {
            extracted = subsampling_wap_pfile_for_test(keep.keep, source)
        }
        else {
            def gff = file(params.test_wap_exon_gff ?: params.ref_v1_gff, checkIfExists: true)
            def bed = gff_feature_bed(gff)
            def exon_source = source.filter { _id, chr, _prefix, _pgen, _psam, _pvar ->
                def n = chr.toString().toInteger()
                return n >= 1 && n <= 42
            }
            extracted = extract_wap_exon_pfile(keep.keep, bed.bed, exon_source)
        }
        def merged = merge_subgenome_test_pfile(extracted.pfile.map(map_subgenome).groupTuple(by: 0))
        mk_plink_basic_info(merged.pfile)
        if (stage == 'exon') {
            plink_assoc_acount(merged.pfile)
        }
    }
}
