nextflow.enable.dsl=2

/*
 * Analysis variant sets: coding masks (VEP) intersected with cohort frequency
 * (WGS-core vs WAP). Emits PLINK2 --extract lists that burden and association
 * runs consume, so downstream tracks never re-derive set membership.
 *
 * Two input modes:
 *   --variant_sets_from_mod <mod>  derive everything per subgenome from a mod's
 *                                  publish tree: variant/{id}.vep.tsv plus the
 *                                  {id}.{cohort}.afreq pair written by cohort_freq
 *   --variant_sets_coding_mask + --variant_sets_dual_af  one ready pair of tables
 *
 * Optional extra inputs, both from partial_task site_anno / gene_expr:
 *   --variant_sets_site_anno_dir  {id}.site_anno.tsv per subgenome; adds the
 *                                 4-fold neutral, conserved and regulatory sets
 *   --variant_sets_gene_expr      one gene_expr.tsv; adds the expressed-gene sets
 *
 * Launch: partial_router.nf --partial_task variant_sets --mod test_wap_refchr \
 *           --variant_sets_from_mod test_wap_refchr
 */

include {
    coding_mask_from_vep
    dual_af_from_afreq
    variant_sets_build
} from '../../../modules/local/genotype/stats/stats_cohort.nf'

// site_anno tables are per subgenome; resolve by prefix (or the sentinel empty file).
def resolveSiteAnno(id, site_anno_dir, no_site_anno) {
    if (params.variant_sets_site_anno) {
        return file(params.variant_sets_site_anno, checkIfExists: true)
    }
    if (!site_anno_dir) {
        return no_site_anno
    }
    def candidate = file("${site_anno_dir}/${id}.site_anno.tsv")
    return candidate.exists() ? candidate : no_site_anno
}

workflow RUN_VARIANT_SETS {
    main:
    if (!params.output_dir || !params.job) {
        error 'RUN_VARIANT_SETS: params.output_dir and params.job are required.'
    }
    if (!params.user_dir) {
        error 'RUN_VARIANT_SETS: params.user_dir is required (Conda stats env).'
    }
    if (!params.mod) {
        error 'RUN_VARIANT_SETS: params.mod is required (publish folder).'
    }

    def no_site_anno = file("${moduleDir}/../../../assets/NO_SITE_ANNO", checkIfExists: true)
    def no_gene_expr = file("${moduleDir}/../../../assets/NO_GENE_EXPR", checkIfExists: true)
    def gene_expr = params.variant_sets_gene_expr
        ? file(params.variant_sets_gene_expr, checkIfExists: true)
        : no_gene_expr
    // site_anno tables are per subgenome, so resolve them by prefix inside the fan-out.
    def site_anno_dir = params.variant_sets_site_anno_dir
        ?: (params.variant_sets_from_mod
            ? "${params.output_dir}/${params.job}/stats/${params.variant_sets_from_mod}/info"
            : null)

    if (params.variant_sets_from_mod) {
        def variant_dir = "${params.output_dir}/${params.job}/process/${params.variant_sets_from_mod}/variant"
        def wgs_cohort = params.variant_sets_wgs_cohort
        def wap_cohort = params.variant_sets_wap_cohort

        def jobs = []
        file("${variant_dir}/*.vep.tsv").sort().each { vep ->
            def id = vep.name.replaceAll(/\.vep\.tsv$/, '')
            def wgs = file("${variant_dir}/${id}.${wgs_cohort}.afreq")
            def wap = file("${variant_dir}/${id}.${wap_cohort}.afreq")
            if (!wgs.exists() || !wap.exists()) {
                log.warn "variant_sets: no ${wgs_cohort}/${wap_cohort} afreq pair for ${id}, skipping"
                return
            }
            jobs << [id, vep, wgs, wap]
        }
        if (!jobs) {
            error "RUN_VARIANT_SETS: found no {vep.tsv + ${wgs_cohort}/${wap_cohort}.afreq} set under ${variant_dir}"
        }

        def annotated = jobs.findAll { row ->
            resolveSiteAnno(row[0], site_anno_dir, no_site_anno).name != 'NO_SITE_ANNO'
        }.collect { row -> row[0] }
        log.info "${params.c_green}variant_sets for ${jobs.collect { row -> row[0] }.join(',')} from ${variant_dir}" +
            " (site_anno: ${annotated ? annotated.join(',') : 'none'})${params.c_reset}"
        def ch_jobs = channel.from(jobs)
        def masks = coding_mask_from_vep(ch_jobs.map { id, vep, wgs, wap -> tuple(id, vep) }).table
        def afs = dual_af_from_afreq(ch_jobs.map { id, vep, wgs, wap -> tuple(id, wgs, wap) }).table
        variant_sets_build(
            masks.join(afs).map { id, mask, af ->
                tuple(id, mask, af, resolveSiteAnno(id, site_anno_dir, no_site_anno), gene_expr)
            }
        )
    }
    else {
        if (!params.variant_sets_coding_mask || !params.variant_sets_dual_af) {
            error 'RUN_VARIANT_SETS: give --variant_sets_from_mod, or both --variant_sets_coding_mask and --variant_sets_dual_af.'
        }
        def coding_mask = file(params.variant_sets_coding_mask, checkIfExists: true)
        def dual_af = file(params.variant_sets_dual_af, checkIfExists: true)
        def prefix = params.variant_sets_prefix ?: "${params.mod}"

        log.info "${params.c_green}variant_sets from ${coding_mask} x ${dual_af}${params.c_reset}"
        variant_sets_build(
            channel.value(tuple(prefix, coding_mask, dual_af, resolveSiteAnno(prefix, site_anno_dir, no_site_anno), gene_expr))
        )
    }
}
