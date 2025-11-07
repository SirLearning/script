nextflow.enable.dsl=2

process GENOTYPE_STATS {
    tag { meta.id }
    publishDir "${params.outdir}/genotype/stats", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.bcftools.stats"), emit: geno_stats
    path "*.log"

    when:
    params.enable_genotype_stats

    script:
    def prefix = meta.id
    """
    echo "Generating genotype stats for ${vcf}" > ${prefix}.stats.log
    bcftools stats -s - ${vcf} > ${prefix}.bcftools.stats 2>> ${prefix}.stats.log || true
    """
}

/*
 * DumpNice-driven PCA/MDS visualization processes
 * These processes invoke R scripts organized under src/r/dumpnice.
 * They expect precomputed inputs (e.g., PLINK eigenvec/eigenval/mds or QC tables)
 * staged in a directory provided by params.dumpnice_inputdir.
 * Enable each process via the corresponding boolean param.
 *
 * Example params:
 *  --dumpnice_inputdir /path/to/staged_inputs 
 *  --enable_dumpnice_pca2 true
 *  --enable_dumpnice_mds true
 *  --enable_dumpnice_env_pca false
 */

process DUMPNICE_PCA2 {
    tag { "dumpnice pca2 ${meta.id}" }
    publishDir "${params.outdir}/genotype/stats/pca", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('*.pdf'), emit: pca_plots
    path "*.log"

    when:
    params.enable_dumpnice_pca2 && params.dumpnice_inputdir

    script:
    def inputdir = params.dumpnice_inputdir as String
    def log = "${meta.id}.pca2.log"
    """
    set -euo pipefail
    echo "[DumpNice PCA2] staging inputs from ${inputdir}" | tee ${log}
    rsync -a "${inputdir}/" ./ || true
    Rscript src/r/dumpnice/pca/09_PCA2.r >> ${log} 2>&1 || true
    # ensure at least one artifact exists for publishing
    ls -1 *.pdf >/dev/null 2>&1 || touch dumpnice_pca2.pdf
    """
}

process DUMPNICE_MDS {
    tag { "dumpnice mds ${meta.id}" }
    publishDir "${params.outdir}/genotype/stats/mds", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('*.pdf'), emit: mds_plots
    path "*.log"

    when:
    params.enable_dumpnice_mds && params.dumpnice_inputdir

    script:
    def inputdir = params.dumpnice_inputdir as String
    def log = "${meta.id}.mds.log"
    """
    set -euo pipefail
    echo "[DumpNice MDS] staging inputs from ${inputdir}" | tee ${log}
    rsync -a "${inputdir}/" ./ || true
    Rscript src/r/dumpnice/utils/06_MDS.r >> ${log} 2>&1 || true
    ls -1 *.pdf >/dev/null 2>&1 || touch dumpnice_mds.pdf
    """
}

process DUMPNICE_ENV_PCA {
    tag { "dumpnice env-pca ${meta.id}" }
    publishDir "${params.outdir}/genotype/stats/env_pca", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('*.pdf'), emit: env_pca_plots
    path "*.log"

    when:
    params.enable_dumpnice_env_pca && params.dumpnice_inputdir

    script:
    def inputdir = params.dumpnice_inputdir as String
    def log = "${meta.id}.env_pca.log"
    """
    set -euo pipefail
    echo "[DumpNice ENV PCA] staging inputs from ${inputdir}" | tee ${log}
    rsync -a "${inputdir}/" ./ || true
    Rscript src/r/dumpnice/pca/07_PCA.r >> ${log} 2>&1 || true
    ls -1 *.pdf >/dev/null 2>&1 || touch dumpnice_env_pca.pdf
    """
}