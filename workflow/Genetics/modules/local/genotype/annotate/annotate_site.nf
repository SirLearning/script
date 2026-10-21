nextflow.enable.dsl=2

/*
 * Site-level annotation from the featHighRes conservation/epigenome BEDs
 * (00data/11epi_anno/chr*_featHighRes.bed.gz).
 *
 * Those files carry 30 per-base features on IWGSC RefSeq v1.0 whole-chromosome
 * coordinates, which is what the *_refchr pfiles use, so the join is a plain
 * (CHROM, POS) match with no lift-over. They are ~10 GB bgzip each and have no
 * tabix index; a sorted single pass costs less than indexing plus millions of
 * random seeks, so site_anno_extract streams one BED per chromosome and keeps
 * only the rows a pfile chromosome asks for.
 */

process site_anno_extract {
    tag "site anno ${subgenome} ${chrom_name}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/run"
    label 'site_anno_stream'

    input:
    tuple val(subgenome), val(chrom_id), val(chrom_name), path(pvar), path(bed)

    output:
    tuple val(subgenome), path("${subgenome}.${chrom_name}.site_anno.tsv.gz"), emit: chunk
    path "*.log", emit: logs

    script:
    """
    set -o pipefail
    AWK=\$(command -v mawk || command -v awk)

    bgzip -d -c -@ ${task.cpus} ${bed} \\
      | "\$AWK" -v want_chr='${chrom_id}' -v chrom_name='${chrom_name}' '
            FNR == NR {
                if (substr(\$0, 1, 1) == "#") next
                # A pfile labels a merged chromosome either by its name (chr4D) or,
                # before the naming change, by the first arm id (23). Accept both.
                if (\$1 != want_chr && \$1 != chrom_name) next
                matched = \$1
                id[\$2] = \$3
                n_want++
                next
            }
            {
                for (p = \$2 + 1; p <= \$3; p++) {
                    if (p in id) {
                        line = id[p] "\\t" chrom_name "\\t" p
                        for (i = 4; i <= NF; i++) line = line "\\t" \$i
                        print line
                        n_hit++
                    }
                }
            }
            END {
                printf "chrom=%s pfile_chr=%s variants=%d annotated=%d missing=%d\\n", \\
                    chrom_name, (matched == "" ? want_chr : matched), \\
                    n_want, n_hit, n_want - n_hit \\
                    > "${subgenome}.${chrom_name}.site_anno.log"
            }
        ' ${pvar} - \\
      | bgzip -c -@ ${task.cpus} > ${subgenome}.${chrom_name}.site_anno.tsv.gz
    """
}

process site_anno_summary {
    tag "site anno summary ${subgenome}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.site_anno.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'

    input:
    tuple val(subgenome), path(chunks)

    output:
    tuple val(subgenome), path("*.site_anno.tsv"), emit: table
    path "*.th.tsv", emit: thresholds
    path "*.png", emit: plots, optional: true
    path "*.log", emit: logs

    script:
    def chunk_list = chunks instanceof List ? chunks : [chunks]
    def quoted = chunk_list.collect { "\"${it}\"" }.join(', ')
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("${subgenome}.site_anno_summary.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.site_anno import run_site_anno

    run_site_anno(
        [${quoted}],
        "${subgenome}",
        phylop_min=float("${params.site_anno_phylop_min}"),
        phastcons_min=float("${params.site_anno_phastcons_min}"),
        gerp_min=float("${params.site_anno_gerp_min}"),
    )
    """
}
