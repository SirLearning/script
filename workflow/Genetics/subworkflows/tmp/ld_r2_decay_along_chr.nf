#!/usr/bin/env nextflow
/*
 * Along-chromosome LD: merge both PLINK arms of a vmap4 chromosome (e.g. 13+14
 * -> chr3A linear coordinates), sample evenly spaced focals, compute unphased r2,
 * and plot the distance at which mean r2 falls to --ld_along_chr_r2 (default 0.2).
 *
 * In each 5 Mb window sample 100 focals; r2 uses the full --ld-window-kb
 * (not the 5 Mb bin). Plot the mean decay-to-threshold distance per window.
 * Non-default r2 writes ld_r2_decay_along.r2_* so default 0.2 artefacts stay.
 *
 *   nextflow run .../subworkflows/tmp/ld_r2_decay_along_chr.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_common_inter_refchr \
 *     --process_dir /data1/.../test_plink/process/test_common_inter_refchr \
 *
 * Also used for --mod test_wap_common_refchr after recode.
 *     --ld_along_chr_chromosomes 1A,2A,3A,4A,5A,6A,7A,1B,2B,3B,4B,5B,6B,7B,1D,2D,3D,4D,5D,6D,7D
 *     --ld_along_chr_plot_only true --ld_along_chr_r2 0.5
 */

nextflow.enable.dsl=2

include { calc_plink_ld_focal_along_chr } from '../../modules/local/genotype/processor/processor_plink2.nf'
include {
    variant_ld_r2_decay_along_subgenome_plot
} from '../../modules/local/genotype/stats/stats_variant.nf'
include {
    hasMergedSubgenomeTestPfiles
    listMergedSubgenomeTestPfileTuples
} from '../../modules/local/infra/infra_plink_reuse.nf'

params.job = 'test_plink'
params.mod = 'test_common_inter'
params.process_dir = null
params.ld_along_chr_plot_only = false

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'ld_r2_decay_along_chr: --output_dir, --job, and --mod are required'
    }
    def process_dir = params.process_dir
        ? params.process_dir
        : "${params.output_dir}/${params.job}/process/${params.mod}"
    if (!hasMergedSubgenomeTestPfiles(process_dir)) {
        error "ld_r2_decay_along_chr: missing merged *_test.plink2 under ${process_dir}"
    }

    def chroms = "${params.ld_along_chr_chromosomes}"
        .split(',')
        .collect { it.trim() }
        .findAll { it }
    if (!chroms) {
        error 'ld_r2_decay_along_chr: --ld_along_chr_chromosomes is empty'
    }

    Channel
        .from(chroms)
        .map { chrom ->
            def token = chrom.replaceAll('^chr', '')
            tuple(token[-1], token)
        }
        .set { ch_chrom }

    if (params.ld_along_chr_plot_only) {
        Channel
            .from(chroms)
            .map { chrom ->
                def token = chrom.replaceAll('^chr', '')
                def id = token[-1]
                def vcor = file("${process_dir}/variant/${id}.${token}.focal.ld.vcor")
                def focals = file("${process_dir}/variant/${id}.${token}.focals.info.tsv")
                if (!vcor.exists() || !focals.exists()) {
                    error "ld_r2_decay_along_chr plot-only: missing ${vcor} or ${focals}"
                }
                tuple(id, token, vcor, focals)
            }
            .set { ch_ld }
    } else {
        Channel
            .from(listMergedSubgenomeTestPfileTuples(process_dir))
            .combine(ch_chrom, by: 0)
            .set { ch_pfile_chrom }

        calc_plink_ld_focal_along_chr(ch_pfile_chrom)
        ch_ld = calc_plink_ld_focal_along_chr.out.ld
    }

    ch_ld
        .filter { id, _chrom, _vcor, _focals -> id.toString() in ['A', 'B', 'D'] }
        .groupTuple(by: 0)
        .set { ch_sg }
    variant_ld_r2_decay_along_subgenome_plot(ch_sg)
}
