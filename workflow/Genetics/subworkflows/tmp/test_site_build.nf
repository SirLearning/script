nextflow.enable.dsl=2

/*
 * Build test_site mod: random 10 Mb genome site panel + popdep extract + variant pfiles.
 *
 * Publish under {output_dir}/{job}/process/test_site/
 *   sites/panel.sites.tsv          — all sampled sites with popdep + Is_Variant
 *   sites/extract/chrNNN.variant.txt
 *   variant/chrNNN.site.pgen/psam/pvar  — variant sites only (from main_raw)
 *
 * Stats (popdep-only, variant highlighted): run test_site_popdep_stats.nf after build.
 *
 * Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/10stats.genome/76run_test_site_build/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/76run_test_site_build
 *   nextflow run .../subworkflows/tmp/test_site_build.nf -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink \
 *     --mod test_site \
 *     --popdep_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/_popdep_n500_gam_credible \
 *     --main_raw_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/main_raw
 */

include { check_input } from '../local/entry/genetics_helpers.nf'
include {
    prepare_test_site_panel
    extract_test_site_pfile_for_test
    publish_test_site_root
} from '../../modules/local/genotype/processor/processor_test.nf'

workflow {
    if (!params.output_dir || !params.job) {
        error 'test_site_build: --output_dir and --job are required.'
    }
    if (!params.mod) {
        error 'test_site_build: --mod test_site is required.'
    }
    if (!params.popdep_dir) {
        error 'test_site_build: --popdep_dir is required (frozen popdep BGZF root).'
    }

    params.main_raw_dir = params.main_raw_dir ?: "${params.output_dir}/${params.job}/process/main_raw"
    params.process_dir = params.process_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}"

    def ch_site_root
    if (params.test_site_skip_prep) {
        def existing = file("${params.process_dir}/sites/panel.sites.tsv")
        if (!existing.exists()) {
            error "test_site_build: --test_site_skip_prep requires ${params.process_dir}/sites/panel.sites.tsv"
        }
        ch_site_root = publish_test_site_root(Channel.value(file("${params.process_dir}")))
    } else {
        ch_site_root = prepare_test_site_panel()
    }

    def input_out = check_input()
    def main_raw = params.main_raw_dir
    source_pfiles = input_out.vcf.combine(ch_site_root).map { id, chr, _vcf, root ->
        def prefix = "${main_raw}/${id}.plink2"
        def extract = file("${root}/sites/extract/${id}.variant.txt")
        tuple(id, chr, prefix, file("${prefix}.pgen"), file("${prefix}.psam"), file("${prefix}.pvar"), extract)
    }.filter { id, chr, prefix, pgen, psam, pvar, extract ->
        extract.exists() && extract.size() > 0
    }

    extract_test_site_pfile_for_test(source_pfiles)
}
