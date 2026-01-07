nextflow.enable.dsl=2

include { getJavaSetupScript; getServerPopulations } from './utils.nf'

/*
 * Variant assessment metrics using bcftools/vcftools
 */

workflow assess {
    take:
        vcf_ch // [ id, vcf ]

    main:
        prepared_vcf_ch = prepare_vcf_gz(vcf_ch)
        counts_ch = quick_count(prepared_vcf_ch)
        vcftools_stats = vcftools_stats(prepared_vcf_ch)
        bcftools_tags = bcftools_tags(prepared_vcf_ch)
        
        // Optional DumpNice plots
        qc_plots = Channel.empty()
        if (params.enable_dumpnice_vcf_assess) {
            qc_plots = dumpnice_vcf_qc_assess(prepared_vcf_ch)
        }

    emit:
        counts = counts_ch
        maf_freq = vcftools_stats.maf_freq
        hwe = vcftools_stats.hwe
        miss_site = vcftools_stats.miss_site
        miss_indv = vcftools_stats.miss_indv
        depth_site = vcftools_stats.depth_site
        depth_indv = vcftools_stats.depth_indv
        site_qual = vcftools_stats.site_qual
        maf_missing = bcftools_tags.maf_missing
        gq_summary = bcftools_tags.gq_summary
        qc_plots = qc_plots
}

process prepare_vcf_gz {
    tag { id }
    publishDir "${params.output_dir}/${params.job}/assess", mode: 'copy'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}.prepared.vcf.gz"), emit: vcf

    script:
    """
    set -euo pipefail
    if [[ "${vcf}" == *.gz ]] || [[ "${vcf}" == *.bgz ]]; then
        ln -sf ${vcf} ${id}.prepared.vcf.gz
    else
        bgzip -c ${vcf} > ${id}.prepared.vcf.gz
    fi
    tabix -f -p vcf ${id}.prepared.vcf.gz
    """
}

process quick_count {
    tag { id }
    publishDir "${params.output_dir}/${params.job}/assess", mode: 'copy'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}.counts.tsv"), emit: counts

    script:
    """
    set -euo pipefail
    n_snp=\$(bcftools view -H ${vcf} | wc -l || echo 0)
    n_samples=\$(bcftools query -l ${vcf} | wc -l || echo 0)
    printf "id\\tn_variants\\tn_samples\\n%s\\t%s\\t%s\\n" "${id}" "\${n_snp}" "\${n_samples}" > ${id}.counts.tsv
    """
}

process vcftools_stats {
    tag { id }
    publishDir "${params.output_dir}/${params.job}/assess/vcftools", mode: 'copy'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}.frq"), emit: maf_freq
    tuple val(id), path("${id}.hwe"), emit: hwe
    tuple val(id), path("${id}.lmiss"), emit: miss_site
    tuple val(id), path("${id}.imiss"), emit: miss_indv
    tuple val(id), path("${id}.ldepth.mean"), emit: depth_site
    tuple val(id), path("${id}.idepth"), emit: depth_indv
    tuple val(id), path("${id}.lqual"), emit: site_qual
    path "*.log"

    script:
    """
    set -euo pipefail
    vcftools --gzvcf ${vcf} \\
        --freq --hardy --missing-site --missing-indv \\
        --site-mean-depth --depth --site-quality \\
        --out ${id} > ${id}.vcftools.log 2>&1 || true
    """
}

process bcftools_tags {
    tag { id }
    publishDir "${params.output_dir}/${params.job}/assess/bcftools", mode: 'copy'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}.maf_missing.tsv"), emit: maf_missing
    tuple val(id), path("${id}.gq_summary.tsv"), emit: gq_summary

    script:
    """
    set -euo pipefail
    bcftools +fill-tags ${vcf} -- -t MAF,F_MISSING,AC,AN | \\
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%MAF\\t%F_MISSING\\t%QUAL\\n' > ${id}.maf_missing.tsv

    bcftools query -f '[%GQ\\t]\\n' ${vcf} | tr '\\t' '\\n' | awk 'NF && \$1 != "."' > all.gq
    awk 'BEGIN{min=1e9;max=-1e9}{s+=\$1; n++; if(\$1<min)min=\$1; if(\$1>max)max=\$1} END{if(n==0){print "metric\\tvalue"; print "n\\t0"; exit} mean=s/n; print "metric\\tvalue"; print "n\\t"n; print "mean\\t"mean; print "min\\t"min; print "max\\t"max}' all.gq > ${id}.gq_summary.tsv
    """
}

process dumpnice_vcf_qc_assess {
    tag { id }
    publishDir "${params.output_dir}/${params.job}/assess/qc", mode: 'copy'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("*.pdf"), emit: qc_plots

    script:
    def inputdir = params.dumpnice_inputdir ?: "."
    """
    set -euo pipefail
    if [ -d "${inputdir}" ]; then
        rsync -a "${inputdir}/" ./ || true
    fi
    Rscript ${params.src_dir}/r/dumpnice/vcf/12_VcfQc.r || true
    ls -1 *.pdf >/dev/null 2>&1 || touch ${id}.vcf_qc_placeholder.pdf
    """
}
