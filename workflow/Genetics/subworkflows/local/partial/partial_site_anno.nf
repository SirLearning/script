nextflow.enable.dsl=2

/*
 * Site annotation: join a mod's variants to the featHighRes conservation and
 * epigenome BEDs, then derive degeneracy / element / conserved / regulatory
 * classes that variant_sets turns into analysis sets.
 *
 * The BEDs cover chr1A-chr7D only, so the Others subgenome (chrUn, Mit, Chl)
 * has no annotation source and is skipped.
 *
 * Launch: partial_router.nf --partial_task site_anno --mod test_wap_refchr \
 *           --site_anno_from_mod test_wap_refchr
 */

include { getRefV1ChrName; getRefV1SubChr } from '../../../modules/local/infra/infra_ref_v1.nf'
include { site_anno_extract; site_anno_summary } from '../../../modules/local/genotype/annotate/annotate_site.nf'

// One entry per whole chromosome: take every other segment id (1, 7, 13, ... for
// subgenome A). The id doubles as the legacy #CHROM value of ref-chr pfiles made
// before the naming change; site_anno_extract matches either it or the chromosome
// name, so both pfile generations work.
def refchrLinearChrIds(subgenome) {
    def segments = getRefV1SubChr(subgenome)
    return (0..<segments.size()).findAll { it % 2 == 0 }.collect { segments[it] }
}

workflow RUN_SITE_ANNO {
    main:
    if (!params.output_dir || !params.job) {
        error 'RUN_SITE_ANNO: params.output_dir and params.job are required.'
    }
    if (!params.user_dir) {
        error 'RUN_SITE_ANNO: params.user_dir is required (Conda run and stats envs).'
    }
    if (!params.mod) {
        error 'RUN_SITE_ANNO: params.mod is required (publish folder).'
    }
    def source_mod = params.site_anno_from_mod ?: params.mod
    def process_dir = "${params.output_dir}/${params.job}/process/${source_mod}"
    def anno_dir = params.site_anno_dir ?: "${params.home_dir}/00data/11epi_anno"

    def wanted = params.site_anno_subgenomes
        ? params.site_anno_subgenomes.toString().split(',').collect { it.trim() }
        : ['A', 'B', 'D']
    def unsupported = wanted.findAll { !(it in ['A', 'B', 'D']) }
    if (unsupported) {
        log.warn "site_anno: no featHighRes BED for ${unsupported.join(',')}; only chr1A-chr7D are covered"
    }

    def jobs = []
    wanted.findAll { it in ['A', 'B', 'D'] }.each { sg ->
        def pvar = file("${process_dir}/${sg}_test.plink2.pvar")
        if (!pvar.exists()) {
            log.warn "site_anno: no pvar for ${sg} at ${pvar}, skipping"
            return
        }
        refchrLinearChrIds(sg).each { chr_id ->
            def chrom_name = getRefV1ChrName(chr_id)
            def bed = file("${anno_dir}/${chrom_name}_featHighRes.bed.gz")
            if (!bed.exists()) {
                log.warn "site_anno: missing ${bed}, skipping ${sg} ${chrom_name}"
                return
            }
            jobs << [sg, chr_id, chrom_name, pvar, bed]
        }
    }
    if (!jobs) {
        error "RUN_SITE_ANNO: found no {pvar, featHighRes BED} pair under ${process_dir} and ${anno_dir}"
    }

    log.info "${params.c_green}site_anno: ${jobs.size()} chromosome passes from ${anno_dir}${params.c_reset}"
    def chunks = site_anno_extract(channel.from(jobs)).chunk
    site_anno_summary(chunks.groupTuple())

    emit:
    table = site_anno_summary.out.table
}
