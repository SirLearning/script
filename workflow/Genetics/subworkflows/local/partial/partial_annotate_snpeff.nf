nextflow.enable.dsl=2

/*
 * SnpEff on the same sites VCFs VEP uses (variant IDs → Lulab GFF coordinates).
 *
 * Launch: partial_router.nf --partial_task annotate_snpeff --mod test_wap_exon \
 *           --vep_from_mod test_wap_exon
 *
 * Reuses process/{mod}/annotate/*.sites.vcf.gz when annotate_vep already wrote
 * them. Otherwise builds those VCFs the same way VEP does.
 */

include {
    hasMergedSubgenomeTestPfiles
    listMergedSubgenomeTestPfileTuples
} from '../../../modules/local/infra/infra_plink_reuse.nf'
include {
    pvar_sites_vcf
    vep_prepare_gff
} from '../../../modules/local/genotype/annotate/annotate_vep.nf'
include {
    snpeff_build_genome
    snpeff_annotate_vcf
    snpeff_vcf_to_tab
} from '../../../modules/local/genotype/annotate/annotate_snpeff.nf'

workflow RUN_ANNOTATE_SNPEFF {
    main:
    if (!params.output_dir || !params.job) {
        error 'RUN_ANNOTATE_SNPEFF: params.output_dir and params.job are required.'
    }
    if (!params.user_dir) {
        error 'RUN_ANNOTATE_SNPEFF: params.user_dir is required (Conda anno env).'
    }
    if (!params.mod) {
        error 'RUN_ANNOTATE_SNPEFF: params.mod is required (publish folder).'
    }
    def from_mod = params.vep_from_mod ?: params.mod
    def proc_dir = "${params.output_dir}/${params.job}/process/${from_mod}"
    if (!hasMergedSubgenomeTestPfiles(proc_dir)) {
        error "RUN_ANNOTATE_SNPEFF: missing merged *_test.plink2 under ${proc_dir}"
    }

    def fasta_path = params.vep_fasta
    if (!fasta_path) {
        if (!params.home_dir) {
            error 'RUN_ANNOTATE_SNPEFF: set --vep_fasta, or --home_dir for the ABD RefSeq v1.0 FASTA.'
        }
        fasta_path = "${params.home_dir}/00data/03ref/03ABD/abd_iwgscV1.fa.gz"
    }
    def ch_fasta = channel.value(tuple(
        file(fasta_path, checkIfExists: true),
        file("${fasta_path}.fai", checkIfExists: true),
        file("${fasta_path}.gzi", checkIfExists: true),
    ))

    def variant_dir = "${params.output_dir}/${params.job}/process/${params.mod}/variant"
    def published_done = ['A', 'B', 'D'].every { id ->
        file("${variant_dir}/${id}.snpeff.vcf.gz").exists() &&
            file("${variant_dir}/${id}.snpeff.tsv").exists()
    }
    def lock = params.home_dir
        ? file("${params.home_dir}/00data/03ref/07snpeff/snpeff.ann.lock")
        : null
    def lock_holder = (lock && lock.exists()) ? lock.text.trim() : ''
    def me = workflow.launchDir.toAbsolutePath().toString()
    if (published_done) {
        log.info "${params.c_green}annotate_snpeff skipping: published VCF+TSV already under ${variant_dir}${params.c_reset}"
        return
    }
    if (lock_holder && lock_holder != me) {
        log.info "${params.c_green}annotate_snpeff skipping: ${lock} held by ${lock_holder}${params.c_reset}"
        return
    }

    def prebuilt = params.home_dir ? "${params.home_dir}/00data/03ref/07snpeff" : null
    def genome_dir = params.snpeff_genome_dir
    if (!genome_dir && prebuilt && file("${prebuilt}/snpEff.config").exists() && file("${prebuilt}/data").exists()) {
        genome_dir = prebuilt
    }
    if (!genome_dir) {
        genome_dir = "${params.output_dir}/${params.job}/process/${params.mod}/annotate/snpeff_genome"
    }
    def cfg = file("${genome_dir}/snpEff.config")
    def data_dir = file("${genome_dir}/data")
    def ch_config
    def ch_data
    if (cfg.exists() && data_dir.exists()) {
        log.info "${params.c_green}annotate_snpeff reusing genome ${genome_dir}${params.c_reset}"
        ch_config = channel.value(cfg)
        ch_data = channel.value(data_dir)
    } else {
        def gff = file(params.vep_gff ?: params.ref_v1_gff, checkIfExists: true)
        log.info "${params.c_green}annotate_snpeff building ${params.snpeff_genome} from ${gff}${params.c_reset}"
        def built = snpeff_build_genome(gff, ch_fasta)
        ch_config = built.config
        ch_data = built.data
    }

    def existing = []
    def annotate_dir = file("${proc_dir}/annotate")
    if (annotate_dir.exists()) {
        existing = annotate_dir.listFiles()?.findAll { it.name.endsWith('.sites.vcf.gz') } ?: []
    }
    def ch_vcf
    if (existing) {
        log.info "${params.c_green}annotate_snpeff reusing sites VCFs under ${annotate_dir}${params.c_reset}"
        ch_vcf = channel.from(existing.collect { vcf ->
            tuple(vcf.name.replaceAll(/\.sites\.vcf\.gz$/, ''), vcf)
        })
    } else {
        def ch_gff
        if (params.vep_prepared_gff) {
            ch_gff = channel.value(tuple(
                file(params.vep_prepared_gff, checkIfExists: true),
                file("${params.vep_prepared_gff}.tbi", checkIfExists: true),
            ))
        } else {
            ch_gff = vep_prepare_gff(file(params.vep_gff ?: params.ref_v1_gff, checkIfExists: true)).gff
        }
        def pfiles = listMergedSubgenomeTestPfileTuples(proc_dir)
        log.info "${params.c_green}annotate_snpeff building sites VCFs from ${proc_dir}${params.c_reset}"
        ch_vcf = pvar_sites_vcf(channel.from(pfiles), ch_gff).vcf
    }

    def annotated = snpeff_annotate_vcf(ch_vcf, ch_config, ch_data)
    snpeff_vcf_to_tab(annotated.vcf)
}
