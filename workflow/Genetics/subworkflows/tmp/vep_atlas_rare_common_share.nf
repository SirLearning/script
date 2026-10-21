#!/usr/bin/env nextflow
/*
 * One 100% stacked bar per subgenome: common (red) under rare (blue).
 * The denominator is n_common + n_rare on that subgenome.
 *
 * Reads published atlas count tables. Does not rewrite them.
 * Writes stats/<mod>/plots/all.variant.vep_atlas.cv_reliable.keep{60,80}.rare_common_share.png
 *
 * Launch from a vmap4 run folder.
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mod = 'test_wap_mac1_refchr'

process variant_rare_common_share {
    tag "rare common share: ${params.mod}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    label 'cpus_1'

    input:
    path counts80
    path counts60

    output:
    path "all.variant.vep_atlas.cv_reliable.keep*.rare_common_share.png", emit: plots
    path "all.variant.vep_atlas.cv_reliable.keep*.rare_common_share.tsv", emit: tables
    path "rare_common_share.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("rare_common_share.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from infra.utils.io import load_df_from_tsv, save_df_to_tsv
    from genetics.genomics.variant.vep_atlas import plot_subgenome_rare_common_share

    for pct, path in (("80", "${counts80}"), ("60", "${counts60}")):
        counts = load_df_from_tsv(path)
        prefix = f"all.variant.vep_atlas.cv_reliable.keep{pct}.rare_common_share"
        prop = plot_subgenome_rare_common_share(
            counts,
            f"{prefix}.png",
            title=f"Rare and common share by subgenome (PopDep keep {pct}%)",
        )
        save_df_to_tsv(prop.reset_index(names="subgenome"), f"{prefix}.tsv")
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'vep_atlas_rare_common_share: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir || !params.src_dir) {
        error 'vep_atlas_rare_common_share: --user_dir and --src_dir are required'
    }
    def info = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
    def counts80 = file("${info}/all.variant.vep_atlas.cv_reliable.keep80.counts.tsv", checkIfExists: true)
    def counts60 = file("${info}/all.variant.vep_atlas.cv_reliable.keep60.counts.tsv", checkIfExists: true)
    variant_rare_common_share(counts80, counts60)
}
