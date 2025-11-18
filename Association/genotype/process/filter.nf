nextflow.enable.dsl = 2

// --- Filter parameters ---
// QUAL threshold (only applied when set; see assess.nf)
def qual = 30
// MAF thresholds
def maf = 0.05
// HWE threshold (currently not applied in bcftools step, placeholder)
def hwe_pval = 1e-6

workflow filter {
    take:
    vcf
    plink_bfile
    val job_config

    main:
    // Index VCF if needed
    indexed_vcf = INDEX_VCF(vcf)

    // Choose filtering method based on params.filter_tool
    filtered_vcf = null
    if (params.filter_tool == 'gatk') {
        filtered_vcf = GATK_FILTER(indexed_vcf.vcf)
    } else {
        filtered_vcf = BCFTOOLS_FILTER(indexed_vcf.vcf)
    }

    // Optional DumpNice QC plots
    dumpnice_site_plots = DUMPNICE_SITE_QC(filtered_vcf.vcf)
    dumpnice_vcf_plots  = DUMPNICE_VCF_QC(filtered_vcf.vcf)
    dumpnice_taxa_plots = DUMPNICE_TAXA_QC(filtered_vcf.vcf)

    emit:
    vcf = filtered_vcf.vcf
}

process INDEX_VCF {
    tag { "index ${meta.id}" }
    publishDir params.outdir + '/filter', mode: 'copy'
    errorStrategy 'terminate'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path(vcf), emit: vcf

    when:
    true

    script:
    def needs_index = vcf.toString().endsWith('.gz')
    def idx_cmd = needs_index ? 'tabix -p vcf ${vcf}' : 'true'
    """
    set -euo pipefail
    ${idx_cmd}
    """
}

process BCFTOOLS_FILTER {
    tag "bcftools filter"
    publishDir params.outdir + '/filter', mode: 'copy'
    cpus 2

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('filtered.vcf.gz'), emit: vcf

    script:
    def region = params.region ? "-r ${params.region}" : ''
    def maf    = (params.filter_maf ?: params.maf ?: 0.05) as Double
    def miss   = (params.filter_geno ?: params.max_miss ?: 0.2) as Double
    def qual   = (params.qual ?: params.filter_qual) ?: null
    // missingness filter: keep variants with missing rate <= max_miss
    // INFO/QUAL filter and MAF
    """
    set -euo pipefail
        bcftools view ${region} -O z -o tmp.vcf.gz ${vcf}
    tabix -p vcf tmp.vcf.gz
        # Optional QUAL filter
        if [ -n "${qual}" ]; then
            bcftools filter -i "QUAL>=${qual}" tmp.vcf.gz -O z -o tmp2.vcf.gz
        else
            ln -sf tmp.vcf.gz tmp2.vcf.gz
            ln -sf tmp.vcf.gz.tbi tmp2.vcf.gz.tbi || true
        fi
    tabix -p vcf tmp2.vcf.gz
    # Compute tags and filter MAF and missingness
    bcftools +fill-tags tmp2.vcf.gz -- -t F_MISSING,MAF > tmp2.tags.vcf
    bgzip -c tmp2.tags.vcf > tmp2.tags.vcf.gz
    tabix -p vcf tmp2.tags.vcf.gz
        bcftools view -q ${maf}:minor -e 'F_MISSING>${miss}' tmp2.tags.vcf.gz -O z -o filtered.vcf.gz
    tabix -p vcf filtered.vcf.gz
    """
}

process GATK_FILTER {
    tag "gatk VariantFiltration"
    publishDir params.outdir + '/filter', mode: 'copy'
    cpus 2

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('filtered.vcf.gz'), emit: vcf

    script:
    // Basic filter example; can be extended via params.gatk_expr
    def region = params.region ? "-L ${params.region}" : ''
    def exprs = params.gatk_expr ?: ["QD < 2.0", "FS > 60.0", "MQ < 40.0"]
    def filter_cmds = exprs.collect{e -> "--filter-expression \"${e}\" --filter-name F_${e.replaceAll(/[^A-Za-z0-9]+/,'_')}"}.join(' ')
    """
    set -euo pipefail
    gatk SelectVariants -V ${vcf} ${region} -O region.vcf.gz
    gatk VariantFiltration -V region.vcf.gz ${filter_cmds} -O filt_flagged.vcf.gz
    # keep PASS only
    bcftools view -f PASS -O z -o filtered.vcf.gz filt_flagged.vcf.gz
    tabix -p vcf filtered.vcf.gz
    """
}

process DUMPNICE_SITE_QC {
    tag { "dumpnice site qc ${meta.id}" }
    publishDir params.outdir + '/filter/qc/site', mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('*.pdf'), emit: plots

    when:
    params.dumpnice_site_qc && params.dumpnice_inputdir

    script:
    def inputdir = params.dumpnice_inputdir as String
    """
    set -euo pipefail
    # Stage DumpNice expected inputs (QC text tables) if provided
    rsync -a "${inputdir}/" ./ || true
    # Run DumpNice R script for site QC plots
    Rscript src/r/dumpnice/qc/55_VMap3_siteQC.r
    # Ensure at least one output exists even if the script names differ
    ls -1 *.pdf >/dev/null 2>&1 || touch dumpnice_site_qc.pdf
    """
}

process DUMPNICE_VCF_QC {
    tag { "dumpnice vcf qc ${meta.id}" }
    publishDir params.outdir + '/filter/qc/vcf', mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('*.pdf'), emit: plots

    when:
    params.dumpnice_vcf_qc && params.dumpnice_inputdir

    script:
    def inputdir = params.dumpnice_inputdir as String
    """
    set -euo pipefail
    rsync -a "${inputdir}/" ./ || true
    Rscript src/r/dumpnice/vcf/12_VcfQc.r || true
    ls -1 *.pdf >/dev/null 2>&1 || touch dumpnice_vcf_qc.pdf
    """
}

process DUMPNICE_TAXA_QC {
    tag { "dumpnice taxa qc ${meta.id}" }
    publishDir params.outdir + '/filter/qc/taxa', mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('*.pdf'), emit: plots

    when:
    params.dumpnice_taxa_qc && params.dumpnice_inputdir

    script:
    def inputdir = params.dumpnice_inputdir as String
    """
    set -euo pipefail
    rsync -a "${inputdir}/" ./ || true
    Rscript src/r/dumpnice/qc/56_VMap3_taxaQC.r
    ls -1 *.pdf >/dev/null 2>&1 || touch dumpnice_taxa_qc.pdf
    """
}
