nextflow.enable.dsl = 2

process INDEX_VCF {
    tag "index ${vcf.simpleName}"
    publishDir params.outdir + '/filter', mode: 'copy'
    errorStrategy 'terminate'

    input:
    path vcf

    output:
    path vcf, emit: vcf

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
    path vcf

    output:
    path 'filtered.vcf.gz', emit: vcf

    script:
    def region = params.region ? "-r ${params.region}" : ''
    def miss   = params.max_miss as Double
    // missingness filter: keep variants with missing rate <= max_miss
    // INFO/QUAL filter and MAF
    """
    set -euo pipefail
    bcftools view ${region} -O z -o tmp.vcf.gz ${vcf}
    tabix -p vcf tmp.vcf.gz
    # Filter QUAL
    bcftools filter -i 'QUAL>=${params.qual}' tmp.vcf.gz -O z -o tmp2.vcf.gz
    tabix -p vcf tmp2.vcf.gz
    # Compute tags and filter MAF and missingness
    bcftools +fill-tags tmp2.vcf.gz -- -t F_MISSING,MAF > tmp2.tags.vcf
    bgzip -c tmp2.tags.vcf > tmp2.tags.vcf.gz
    tabix -p vcf tmp2.tags.vcf.gz
    bcftools view -q ${params.maf}:minor -e 'F_MISSING>${miss}' tmp2.tags.vcf.gz -O z -o filtered.vcf.gz
    tabix -p vcf filtered.vcf.gz
    """
}

process GATK_FILTER {
    tag "gatk VariantFiltration"
    publishDir params.outdir + '/filter', mode: 'copy'
    cpus 2

    input:
    path vcf

    output:
    path 'filtered.vcf.gz', emit: vcf

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

workflow FILTER_VCF {
    take:
    ch_vcf

    main:
    ch_pre = INDEX_VCF(ch_vcf)
    if (params.filter_tool.toString().toLowerCase() == 'gatk') {
        ch_out = GATK_FILTER(ch_pre.vcf)
    } else {
        ch_out = BCFTOOLS_FILTER(ch_pre.vcf)
    }

    emit:
    ch_out.vcf
}
