#!/usr/bin/env nextflow
/*
 * Rewrite #CHROM of an already ref-chr merged mod from the first arm's numeric id
 * (chr2A labelled `7`) to the chromosome name (`chr2A`).
 *
 * Merging arms and naming them are separate concerns. The original recode did the
 * first and reused the first arm's id for the second, which is unambiguous but
 * reads as "chromosome 7" on every plot axis and in every derived table.
 *
 * This entry only renames. The pfiles it reads are already merged, so POS is
 * whole-chromosome linear and write_all_ref_merge_update_maps emits update rows
 * that change the chromosome and leave POS untouched. Variant IDs keep their
 * original {arm}-{local pos} form, so the VEP / coding_mask / dual_af / site_anno
 * / variant_sets tables built earlier still join.
 *
 * Rebuilding the mod from source would be wrong here: the common and thin mods are
 * random subsamples, so a rebuild would draw a different variant set and invalidate
 * their published LD results.
 *
 *   nextflow run .../subworkflows/tmp/refchr_rename.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_refchr
 *
 * Publishes back over the same mod, so the previous pfiles are replaced.
 */

nextflow.enable.dsl=2

include { recode_subgenome_pfile_ref_chromosomes } from '../../modules/local/genotype/processor/processor_plink2.nf'
include {
    hasMergedSubgenomeTestPfiles
    listMergedSubgenomeTestPfileTuples
} from '../../modules/local/infra/infra_plink_reuse.nf'

params.job = 'test_plink'
params.mod = null
params.process_dir = null

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'refchr_rename: --output_dir, --job, and --mod are required'
    }
    def process_dir = params.process_dir
        ? params.process_dir
        : "${params.output_dir}/${params.job}/process/${params.mod}"
    if (!hasMergedSubgenomeTestPfiles(process_dir)) {
        error "refchr_rename: missing merged *_test.plink2 under ${process_dir}"
    }
    log.info "${params.c_green}refchr rename:${params.c_reset} ${process_dir} (#CHROM -> chr1A..chr7D; POS and variant IDs unchanged)"

    Channel
        .from(listMergedSubgenomeTestPfileTuples(process_dir))
        .set { ch_pfile }

    recode_subgenome_pfile_ref_chromosomes(ch_pfile)
}
