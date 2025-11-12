nextflow.enable.dsl=2

workflow {
    take:
        ch_vcf // tuple val(meta), path(vcf)

    main:
        GENOTYPE_PLINK_STATS(ch_vcf)

    emit:
        bfiles: GENOTYPE_PLINK_STATS.out.bfiles
        missing: GENOTYPE_PLINK_STATS.out.missing
        freq: GENOTYPE_PLINK_STATS.out.freq
        het: GENOTYPE_PLINK_STATS.out.het
        pca: GENOTYPE_PLINK_STATS.out.pca
        plots: GENOTYPE_PLINK_STATS.out.plots
        pca_plot: GENOTYPE_PLINK_STATS.out.pca_plot
}

workflow GENOTYPE_PLINK_STATS {
    take:
        ch_vcf // tuple val(meta), path(vcf)

    main:
        ch_bfiles = PLINK_FROM_VCF(ch_vcf)
        ch_missing = PLINK_MISSING(ch_bfiles.bfiles)
        ch_freq = PLINK_FREQ(ch_bfiles.bfiles)
        ch_het = PLINK_HET(ch_bfiles.bfiles)
        ch_pca = PLINK_PCA(ch_bfiles.bfiles)

        if (params.enable_pca_plot) {
            // Try to locate optional metadata channel if provided; else use empty file
            def md = params.sample_metadata ? file(params.sample_metadata) : file('')
            ch_pca_plot = PLOT_PLINK_PCA(ch_pca.pca.map{ it -> tuple(it[0], it[1], it[2]) }, ch_pca.pca.map{ it -> it[3] }, ch_pca.pca.map{ it -> it[0] }, md)
        }

        if (params.enable_simple_plots) {
            ch_plots = PLOT_PLINK_QC(ch_missing.missing, ch_freq.freq, ch_het.het)
        }

    emit:
        bfiles: ch_bfiles.bfiles
        missing: ch_missing.missing
        freq: ch_freq.freq
        het: ch_het.het
        pca: ch_pca.pca
        plots: params.enable_simple_plots ? ch_plots.qc_plots : Channel.empty()
        pca_plot: params.enable_pca_plot ? ch_pca_plot.pca_plot : Channel.empty()
}

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
 * PLINK-based stats from VCF, aligned with note/00basic/03vcf_process.ipynb
 * Steps:
 *  - Convert VCF -> PLINK bed/bim/fam (chr-set default 42)
 *  - Missingness (imiss/lmiss)
 *  - MAF freq (.frq)
 *  - Heterozygosity (.het)
 *  - PCA (eigenvec/eigenval)
 *  - Optional quick plots (histograms) via inline R
 */

process PLINK_FROM_VCF {
    tag { "plink make-bed ${meta.id}" }
    publishDir "${params.outdir}/genotype/stats/plink", mode: 'copy'
    cpus 2

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.bed"), path("*.bim"), path("*.fam"), val(prefix), emit: bfiles
    path "*.log"

    script:
    def prefix = meta.id ?: 'plink'
    def chrset = params.plink_chr_set ?: 42
    def log = "${prefix}.plink_makebed.log"
    """
    set -euo pipefail
    echo "[PLINK] VCF->BED for ${vcf}" | tee ${log}
    plink \
      --vcf ${vcf} \
      --double-id \
      --make-bed \
      --chr-set ${chrset} \
      --out ${prefix} >> ${log} 2>&1
    echo ${prefix} > __prefix.txt
    """
}

process PLINK_MISSING {
    tag { "plink missing ${meta.id}" }
    publishDir "${params.outdir}/genotype/stats/plink", mode: 'copy'
    cpus 1

    input:
    tuple val(meta), path(bed), path(bim), path(fam), val(prefix)

    output:
    tuple val(meta), path("*.imiss"), path("*.lmiss"), emit: missing
    path "*.log"

    script:
    def log = "${prefix}.missing.log"
    """
    set -euo pipefail
    plink --bfile ${prefix} --missing --out ${prefix} > ${log} 2>&1
    """
}

process PLINK_FREQ {
    tag { "plink freq ${meta.id}" }
    publishDir "${params.outdir}/genotype/stats/plink", mode: 'copy'
    cpus 1

    input:
    tuple val(meta), path(bed), path(bim), path(fam), val(prefix)

    output:
    tuple val(meta), path("*.frq"), emit: freq
    path "*.log"

    script:
    def log = "${prefix}.freq.log"
    """
    set -euo pipefail
    plink --bfile ${prefix} --freq --out ${prefix} > ${log} 2>&1
    """
}

process PLINK_HET {
    tag { "plink het ${meta.id}" }
    publishDir "${params.outdir}/genotype/stats/plink", mode: 'copy'
    cpus 1

    input:
    tuple val(meta), path(bed), path(bim), path(fam), val(prefix)

    output:
    tuple val(meta), path("*.het"), emit: het
    path "*.log"

    script:
    def log = "${prefix}.het.log"
    """
    set -euo pipefail
    plink --bfile ${prefix} --het --out ${prefix} > ${log} 2>&1
    """
}

process PLINK_PCA {
    tag { "plink pca ${meta.id}" }
    publishDir "${params.outdir}/genotype/stats/plink", mode: 'copy'
    cpus 2

    input:
    tuple val(meta), path(bed), path(bim), path(fam), val(prefix)

    output:
    tuple val(meta), path("*.eigenvec"), path("*.eigenval"), emit: pca
    path "*.log"

    script:
    def ncomp = params.plink_pca_components ?: 10
    def log = "${prefix}.pca.log"
    """
    set -euo pipefail
    plink --bfile ${prefix} --pca ${ncomp} --out ${prefix} > ${log} 2>&1
    """
}

process PLOT_PLINK_PCA {
    tag { "plot plink pca ${meta.id}" }
    publishDir "${params.outdir}/genotype/stats/plots", mode: 'copy'
    cpus 1

    input:
    tuple val(meta), path(eigenvec), path(eigenval)
    val(meta2)
    path(metadata)

    output:
    tuple val(meta), path('*.pca_pc12.pdf'), emit: pca_plot
    path "*.log"

    when:
    params.enable_pca_plot

    script:
    def prefix = meta.id ?: 'pca'
    def log = "${prefix}.pca_plot.log"
    def mdflag = metadata ? "--metadata ${metadata}" : ""
    """
    set -euo pipefail
    echo "[PCA Plot] Generating PC1 vs PC2 plot" | tee ${log}
    Rscript src/r/dumpnice/wrappers/plot_pca_from_eigenvec.R \
        --eigenvec ${eigenvec} \
        --eigenval ${eigenval} \
        ${mdflag} \
        --out-prefix ${prefix} >> ${log} 2>&1 || true
    ls -1 *.pca_pc12.pdf >/dev/null 2>&1 || touch ${prefix}.pca_pc12.pdf
    """
}

process PLOT_PLINK_QC {
    tag { "plot qc ${meta.id}" }
    publishDir "${params.outdir}/genotype/stats/plots", mode: 'copy'
    cpus 1

    input:
    tuple val(meta), path(imiss), path(lmiss)
    tuple val(meta2), path(frq)
    tuple val(meta3), path(het)

    output:
    tuple val(meta), path('imiss.pdf'), path('lmiss.pdf'), path('maf_distribution.pdf'), path('heterozygosity.pdf'), emit: qc_plots

    when:
    params.enable_simple_plots

    script:
    """
    set -euo pipefail
    R --vanilla <<'RSCRIPT'
    indmiss <- read.table('${imiss}', header=TRUE)
    snpmiss <- read.table('${lmiss}', header=TRUE)
    maf_freq <- read.table('${frq}', header=TRUE)
    het <- read.table('${het}', header=TRUE)

    pdf('imiss.pdf')
    hist(indmiss[[6]], main='Histogram individual missingness', xlab='Individual missingness')
    dev.off()

    pdf('lmiss.pdf')
    hist(snpmiss[[5]], main='Histogram SNP missingness', xlab='SNP missingness')
    dev.off()

    pdf('maf_distribution.pdf')
    hist(maf_freq[[5]], main='MAF distribution', xlab='MAF')
    dev.off()

    pdf('heterozygosity.pdf')
    het$HET_RATE <- (het$N.NM. - het$O.HOM.) / het$N.NM.
    hist(het$HET_RATE, xlab='Heterozygosity Rate', ylab='Frequency', main='Heterozygosity Rate')
    dev.off()
RSCRIPT
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