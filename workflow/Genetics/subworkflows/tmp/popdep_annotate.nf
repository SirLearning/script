#!/usr/bin/env nextflow
/*
 * Tabix-annotate variant popdep from frozen BGZF grids (params.popdep_dir) onto merged test pfiles.
 *
 * Reads merged A/B/D/Others *_test.plink2 from params.process_dir; publishes
 *   {output_dir}/{job}/process/{mod}/variant/{sg}.popdep.info.tsv
 *
 * Launch from a vmap4 run folder (see doc/project_knowledge/workspace/vmap4_10stats_genome.yaml). Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/10stats.genome/33run_popdep_n500_annotate_thin/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/33run_popdep_n500_annotate_thin
 *   nextflow run .../subworkflows/tmp/popdep_annotate.nf -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink \
 *     --mod test_thin \
 *     --process_dir /data1/.../test_plink/process/test_thin \
 *     --popdep_dir /data1/.../test_plink/process/_popdep_n500_gam_credible
 */

nextflow.enable.dsl=2

include { annotate_subgenome_variant_popdep } from '../../modules/local/genotype/processor/processor_depth.nf'
include {
    hasMergedSubgenomeTestPfiles
    listMergedSubgenomeTestPfileTuples
} from '../../modules/local/infra/infra_plink_reuse.nf'

params.job = 'test_plink'
params.mod = null
params.process_dir = null
params.popdep_dir = null

workflow {
    if (!params.mod) {
        error 'popdep_annotate: --mod is required (e.g. test_thin, test_common_thin)'
    }
    if (!params.process_dir) {
        error 'popdep_annotate: --process_dir is required (merged *_test.plink2 under this dir)'
    }
    if (!params.popdep_dir) {
        error 'popdep_annotate: --popdep_dir is required (frozen chrNNN.popdep.txt.bgz + .tbi root)'
    }
    if (!params.output_dir || !params.job) {
        error 'popdep_annotate: --output_dir and --job are required'
    }
    if (!hasMergedSubgenomeTestPfiles(params.process_dir)) {
        error "popdep_annotate: process_dir missing merged subgenome pfiles: ${params.process_dir}"
    }

    merge_pfile = Channel.from(listMergedSubgenomeTestPfileTuples(params.process_dir))
    annotate_subgenome_variant_popdep(merge_pfile)
}
