nextflow.enable.dsl=2

/*
 * VEP consequence table for one VCF against the local IWGSC v1.1 GFF.
 * Feeds coding_mask_from_vep (partial_task variant_sets / coding_mask_chr_pilot.nf).
 *
 * Runs in Conda env `anno`. The GFF prepare step sorts and tabix-indexes a 200 MB+ GFF.
 *
 * Two input modes:
 *   --vep_vcf <file>      annotate one ready VCF
 *   --vep_from_mod <mod>  build one sites VCF per merged subgenome pfile from the
 *                         variant IDs, which carry <segment>-<segment_pos>
 *                         coordinates matching the GFF and FASTA
 *
 * Launch: partial_router.nf --partial_task annotate_vep --mod test_wap_refchr \
 *           --vep_from_mod test_wap_refchr
 */

include {
    hasMergedSubgenomeTestPfiles
    listMergedSubgenomeTestPfileTuples
} from '../../../modules/local/infra/infra_plink_reuse.nf'
include {
    pvar_sites_vcf
    vep_prepare_gff
    vep_annotate_vcf
} from '../../../modules/local/genotype/annotate/annotate_vep.nf'

workflow RUN_ANNOTATE_VEP {
    main:
    if (!params.output_dir || !params.job) {
        error 'RUN_ANNOTATE_VEP: params.output_dir and params.job are required.'
    }
    if (!params.user_dir) {
        error 'RUN_ANNOTATE_VEP: params.user_dir is required (Conda anno env).'
    }
    if (!params.mod) {
        error 'RUN_ANNOTATE_VEP: params.mod is required (publish folder).'
    }
    if (!params.vep_vcf && !params.vep_from_mod) {
        error 'RUN_ANNOTATE_VEP: give --vep_vcf (single VCF) or --vep_from_mod (build sites VCFs from a pfile mod).'
    }

    def fasta_path = params.vep_fasta
    if (!fasta_path) {
        if (!params.home_dir) {
            error 'RUN_ANNOTATE_VEP: set --vep_fasta, or --home_dir to use the default ABD RefSeq v1.0 FASTA.'
        }
        fasta_path = "${params.home_dir}/00data/03ref/03ABD/abd_iwgscV1.fa.gz"
    }
    def ch_fasta = channel.value(tuple(
        file(fasta_path, checkIfExists: true),
        file("${fasta_path}.fai", checkIfExists: true),
        file("${fasta_path}.gzi", checkIfExists: true),
    ))

    def ch_gff
    if (params.vep_prepared_gff) {
        def prepared = file(params.vep_prepared_gff, checkIfExists: true)
        def index = file("${params.vep_prepared_gff}.tbi", checkIfExists: true)
        log.info "${params.c_green}annotate_vep reusing prepared GFF ${prepared}${params.c_reset}"
        ch_gff = channel.value(tuple(prepared, index))
    }
    else {
        def gff = file(params.vep_gff ?: params.ref_v1_gff, checkIfExists: true)
        log.info "${params.c_green}annotate_vep preparing GFF ${gff}${params.c_reset}"
        ch_gff = vep_prepare_gff(gff).gff
    }

    def ch_vcf
    if (params.vep_from_mod) {
        def proc_dir = "${params.output_dir}/${params.job}/process/${params.vep_from_mod}"
        if (!hasMergedSubgenomeTestPfiles(proc_dir)) {
            error "RUN_ANNOTATE_VEP: missing merged *_test.plink2 under ${proc_dir}"
        }
        def pfiles = listMergedSubgenomeTestPfileTuples(proc_dir)
        if (params.vep_subgenomes) {
            def wanted = "${params.vep_subgenomes}".split(',').collect { it.trim() }.findAll { it }
            pfiles = pfiles.findAll { it[0] in wanted }
            if (!pfiles) {
                error "RUN_ANNOTATE_VEP: --vep_subgenomes ${params.vep_subgenomes} matched no pfile under ${proc_dir}"
            }
        }
        log.info "${params.c_green}annotate_vep building sites VCFs from ${proc_dir} for ${pfiles.collect { it[0] }.join(',')}${params.c_reset}"
        ch_vcf = pvar_sites_vcf(channel.from(pfiles), ch_gff).vcf
    }
    else {
        def vcf = file(params.vep_vcf, checkIfExists: true)
        def id = params.vep_id ?: vcf.name.replaceAll(/\.vcf(\.gz)?$/, '')
        log.info "${params.c_green}annotate_vep ${id} on ${vcf}${params.c_reset}"
        ch_vcf = channel.value(tuple(id, vcf))
    }

    vep_annotate_vcf(ch_vcf, ch_gff, ch_fasta)
}
