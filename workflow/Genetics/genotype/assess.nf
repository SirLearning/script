nextflow.enable.dsl=2

include { getJavaSetupScript; getServerPopulations } from './utils.nf'

/*
 * Variant assessment metrics using bcftools/vcftools
 */

workflow plink_assess {
    take:
        smiss
        vmiss
        gcount
        afreq
        hardy

    main:
        // assess missing rate threshold
        // 1. sample assessment
        sample_missing = cpt_sample_missing(smiss)
    emit:
        missing_th = sample_missing.missing_th
}

process cpt_sample_missing {
    tag "compute missing rate threshold: ${id}"
    publishDir "${params.output_dir}/${params.job}/assess/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/assess/threshold", mode: 'copy', pattern: "*.smiss_th.tsv"
    publishDir "${params.output_dir}/${params.job}/assess/logs", mode: 'copy', pattern: "*.log"
    conda 'stats'

    input:
    tuple val(id), val(chr), path(smiss)

    output:
    tuple val(id), path("*.smiss_th.tsv"), emit: smiss_th
    tuple val(id), path("*.png"), emit: plots
    tuple val(id), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${id}.smiss_ana.log", "w")
    sys.stderr = sys.stdout

    from python_script.genomics.sample.smiss_ana import calculate_missing_threshold, plot_missing_dist

    print(f"Processing sample missing rate for ${id}...")
    plot_missing_dist("${smiss}", "${chr}.missing_dist")
    
    print(f"Calculating threshold...")
    calculate_missing_threshold("${smiss}", "${id}.smiss_th.tsv")
    """
}


// -- old code --

process quick_count {
    tag "quick count: ${id}"
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
    tag "vcftools stats: ${id}"
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
    path "${id}.vcftools.log"

    script:
    """
    set -euo pipefail
    vcftools --gzvcf ${vcf} \\
        --freq --hardy --missing-site --missing-indv \\
        --site-mean-depth --depth --site-quality \\
        --out ${id} > ${id}.vcftools.log 2>&1 || true
    """
}

process bcftools_qc_assess {
    tag "assess with bcftools: ${id}"
    publishDir "${params.output_dir}/${params.job}/assess", mode: 'copy'

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
    tag "dumpnice vcf qc assess: ${id}"
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
