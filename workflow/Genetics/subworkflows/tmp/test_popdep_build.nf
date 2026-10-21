#!/usr/bin/env nextflow
/*
 * test_popdep: random 100M popdep sites + gap/repeat/GC/CNV/normal class plots.
 *
 * Does not extract variant pfiles or run variant stats. Publishes:
 *   {output_dir}/{job}/process/test_popdep/annotation/{gap,repeat,gc,cnv}.seg.bed
 *   {output_dir}/{job}/process/test_popdep/sites/chrNNN.popdep.sites.tsv
 *   {output_dir}/{job}/stats/test_popdep/{plots,info,logs}/
 *
 * Repeat: IWGSC RefSeq v1.0 CLARITE TE GFF.
 * GC: extreme-GC 1 kb windows from the v1.0 ABD FASTA (IWGSC does not ship a GC track).
 * CNV: 1 Mb mapq_bin_popdep bins with cohort cnv_frac >= threshold
 *      (IWGSC v1.0 annotations do not include a CS CNV track).
 *
 * Launch from a vmap4 run folder:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/14stats.genome/11run_test_popdep/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/14stats.genome/11run_test_popdep
 *   nextflow run .../subworkflows/tmp/test_popdep_build.nf -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink \
 *     --mod test_popdep \
 *     --popdep_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/_popdep_n500_gam_credible
 */

nextflow.enable.dsl=2

include {
    prepare_test_popdep_annotations
    prepare_test_popdep_chr_panel
} from '../../modules/local/genotype/processor/processor_test.nf'
include { test_popdep_class_plots } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = 'test_popdep'
params.popdep_dir = null

workflow {
    if (!params.mod) {
        error 'test_popdep_build: --mod test_popdep is required.'
    }
    if (!params.output_dir || !params.job) {
        error 'test_popdep_build: --output_dir and --job are required.'
    }
    if (!params.popdep_dir) {
        error 'test_popdep_build: --popdep_dir is required (frozen popdep BGZF root).'
    }

    ch_anno = prepare_test_popdep_annotations()
    ch_chr = ch_anno.alloc.splitCsv(header: true, sep: '\t').map { row ->
        tuple(row.Chr_id, row.n_sites.toInteger())
    }.filter { _chr_id, n_sites -> n_sites > 0 }

    ch_sites = prepare_test_popdep_chr_panel(ch_chr, ch_anno.annotation_dir)
    test_popdep_class_plots(ch_sites.sites.collect().map { true })
}
