nextflow.enable.dsl=2

/*
 * WAP association track: common-variant GLM now, rare-variant burden declared.
 *
 * Genotypes are the merged subgenome pfiles of the chosen mod, restricted to
 * assoc.keep (WAP only). Variant restriction comes from a variant_sets
 * *.extract list when given. Phenotypes come from partial_task phenotype_wap.
 *
 * Launch: partial_router.nf --partial_task assoc_wap --mod <pfile mod> \
 *           --assoc_wap_pheno ... --assoc_wap_trait ...
 */

include {
    hasMergedSubgenomeTestPfiles
    listMergedSubgenomeTestPfileTuples
} from '../../../modules/local/infra/infra_plink_reuse.nf'
include {
    plink2_wap_pca_covar
    plink2_wap_glm
    wap_burden_test
} from '../../../modules/local/static/gwas/gwas_wap.nf'
include { plot_wap_association; plot_wap_burden } from '../../../modules/local/static/gwas/gwas_plot.nf'

workflow RUN_ASSOC_WAP {
    main:
    if (!params.output_dir || !params.job) {
        error 'RUN_ASSOC_WAP: params.output_dir and params.job are required.'
    }
    if (!params.user_dir) {
        error 'RUN_ASSOC_WAP: params.user_dir is required (Conda stats env).'
    }
    if (!params.mod) {
        error 'RUN_ASSOC_WAP: params.mod is required (pfile + publish folder).'
    }
    if (!params.assoc_wap_pheno || !params.assoc_wap_trait) {
        error 'RUN_ASSOC_WAP: --assoc_wap_pheno and --assoc_wap_trait are required.'
    }

    def proc_dir = params.process_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}"
    if (!hasMergedSubgenomeTestPfiles(proc_dir)) {
        error "RUN_ASSOC_WAP: missing merged *_test.plink2 under ${proc_dir}"
    }

    def keep_dir = params.cohort_keep_dir ?: "${params.output_dir}/${params.job}/process/sample_manifest"
    def keep = file(params.assoc_wap_keep ?: "${keep_dir}/assoc.keep", checkIfExists: true)
    def pheno = file(params.assoc_wap_pheno, checkIfExists: true)
    // projectDir points at the entry script's folder, so resolve assets from this module instead.
    def assets_dir = "${moduleDir}/../../../assets"
    def covar = params.assoc_wap_covar
        ? file(params.assoc_wap_covar, checkIfExists: true)
        : file("${assets_dir}/NO_COVAR", checkIfExists: true)
    def extract = params.assoc_wap_extract
        ? file(params.assoc_wap_extract, checkIfExists: true)
        : file("${assets_dir}/NO_EXTRACT", checkIfExists: true)

    def ch_pfile = channel.from(listMergedSubgenomeTestPfileTuples(proc_dir))
    log.info "${params.c_green}assoc_wap ${params.assoc_wap_trait} on ${proc_dir}; keep ${keep}${params.c_reset}"

    // Structure PCs describe the individual, so compute them once from one pfile
    // and broadcast the same covariates to every subgenome. Burden needs them for
    // the same reason the single-variant track does.
    def ch_covar = channel.value(covar)
    if (params.assoc_wap_pca && !params.assoc_wap_covar) {
        def source_id = params.assoc_wap_pca_source
        def source = listMergedSubgenomeTestPfileTuples(proc_dir).find { it[0] == source_id }
        if (!source) {
            error "RUN_ASSOC_WAP: --assoc_wap_pca_source ${source_id} is not a merged subgenome pfile."
        }
        log.info "${params.c_green}assoc_wap structure PCs from ${source_id}${params.c_reset}"
        ch_covar = plink2_wap_pca_covar(channel.value(source), keep).covar.first()
    }

    if (params.assoc_wap_burden) {
        if (!params.assoc_wap_burden_set) {
            error 'RUN_ASSOC_WAP: --assoc_wap_burden needs --assoc_wap_burden_set (a variant_sets set name).'
        }
        def set_name = params.assoc_wap_burden_set
        def variant_dir = "${proc_dir}/variant"
        def sets_dir = params.assoc_wap_burden_sets_dir
            ?: "${params.output_dir}/${params.job}/stats/${params.mod}/info"

        // Both the extract list and the gene map are per subgenome, so pair them
        // up per pfile rather than broadcasting one file across all four.
        def jobs = []
        listMergedSubgenomeTestPfileTuples(proc_dir).each { pfile ->
            def id = pfile[0]
            def set_extract = params.assoc_wap_burden_extract
                ? file(params.assoc_wap_burden_extract)
                : file("${variant_dir}/${id}.${set_name}.extract")
            def sets_table = file("${sets_dir}/${id}.variant_sets.tsv")
            if (params.assoc_wap_burden_extract && !set_extract.name.startsWith("${id}.")) {
                return
            }
            if (!set_extract.exists()) {
                log.warn "assoc_wap burden: no ${set_extract} for ${id}, skipping"
                return
            }
            // plink2 --extract fails outright on an empty list; several sets are
            // legitimately empty for a given subgenome.
            if (set_extract.countLines() == 0) {
                log.warn "assoc_wap burden: ${set_extract.name} is empty for ${id}, skipping"
                return
            }
            if (!sets_table.exists()) {
                log.warn "assoc_wap burden: no ${sets_table} gene map for ${id}, skipping"
                return
            }
            // Nextflow tuples are immutable, so rebuild rather than append in place.
            jobs << (pfile.toList() + [set_name, set_extract, sets_table])
        }
        if (!jobs) {
            error "RUN_ASSOC_WAP: no subgenome has a non-empty ${set_name}.extract under ${variant_dir}"
        }
        log.info "${params.c_green}assoc_wap burden ${set_name} on " +
            "${jobs.collect { it[0] }.join(',')}${params.c_reset}"

        def burden = wap_burden_test(
            channel.from(jobs),
            keep,
            pheno,
            params.assoc_wap_trait,
            ch_covar,
        )
        def gff = file(params.assoc_wap_burden_gff ?: params.ref_v1_gff_hc, checkIfExists: true)
        plot_wap_burden(
            burden.burden.map { sg, trait_name, burden_set, table ->
                tuple(table, "${sg}.${trait_name}.${burden_set}")
            },
            gff,
        )
    }
    else {
        def glm = plink2_wap_glm(
            ch_pfile, keep, pheno, params.assoc_wap_trait, ch_covar, extract
        )
        // One plot per result file; the prefix keeps subgenomes apart in the publish tree.
        plot_wap_association(
            glm.glm
                .map { id, trait, results -> results }
                .flatten()
                .map { result -> tuple(result, result.name.replaceAll(/\.glm\..*$/, '')) }
        )
    }
}
