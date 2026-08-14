nextflow.enable.dsl=2

/*
 * Build test_collinear: cached site lists + join + ~2M sample + full test stats pipeline.
 *
 * Cached under publish: {output_dir}/{job}/process/test_collinear/collinear/
 *   full_lineage_sites/*.tsv.gz
 *   common_variant_sites/chrNNN.tsv.gz
 *   extract/chrNNN.collinear.txt
 *
 * Full build:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/10stats.genome/70run_test_collinear_build/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/70run_test_collinear_build
 *   nextflow run .../subworkflows/tmp/collinear_test_build.nf -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink \
 *     --mod test_collinear
 *
 * Resume extract+stats from an existing prep tree (skip 5h prep):
 *   --collinear_skip_prep --collinear_root /path/to/collinear --mod test_collinear
 */

/*
 * params.mod: base.config defaults mod=null. Script-level params.mod=... cannot override
 * config params — always pass CLI --mod test_collinear (see launch examples above).
 */

include { check_input } from '../local/entry/genetics_helpers.nf'
include {
    collinear_sites_prep
    publish_collinear_root
} from '../../modules/local/genotype/processor/processor_test.nf'
include { RUN_TEST_COLLINEAR } from '../local/plink/plink_genotype_modes.nf'

workflow {
    if (!params.output_dir || !params.job) {
        error 'collinear_test_build: --output_dir and --job are required.'
    }
    if (!params.mod) {
        error 'collinear_test_build: --mod test_collinear is required (config mod=null; script assignment cannot override).'
    }
    params.main_raw_dir = params.main_raw_dir ?: params.popdep_dir
    params.process_dir = params.process_dir ?: "${params.output_dir}/${params.job}/process/test_collinear"
    // Preserve CLI overrides (e.g. reuse 67run prep work tree).
    params.collinear_root = params.collinear_root ?: "${params.process_dir}/collinear"
    params.collinear_extract_dir = params.collinear_extract_dir ?: "${params.collinear_root}/extract"
    params.syntenic_matrix_dir = params.syntenic_matrix_dir ?: '/data/home/tusr1/01projects/vmap4/00data/06maf/syntenic_sites_k2'
    params.full_lineage_summary_tsv = params.full_lineage_summary_tsv ?: '/data/home/tusr1/01projects/vmap4/10stats.genome/50run_maf_full_lineage_count/out/full_lineage_summary.tsv'

    def ch_collinear_root
    if (params.collinear_skip_prep) {
        def existing = file(params.collinear_root)
        if (!existing.exists() || !existing.isDirectory()) {
            error "collinear_test_build: --collinear_skip_prep requires existing --collinear_root (${params.collinear_root})"
        }
        log.info "Skipping collinear_sites_prep; publishing existing root: ${params.collinear_root}"
        ch_collinear_root = publish_collinear_root(Channel.value(existing))
    } else {
        ch_collinear_root = collinear_sites_prep()
    }

    def input_out = check_input()
    RUN_TEST_COLLINEAR(input_out.vcf, ch_collinear_root)
}
