nextflow.enable.dsl=2

include { getJavaSetupScript; getServerPopulations } from './utils.nf'

/*
 * Variant assessment: legacy bcftools/VCF paths plus PLINK2-native slice QC.
 * For merged test pfiles, prefer PLINK2 (--freq / --missing) over bcftools when possible.
 */

/**
 * Representative-chromosome slice on existing per-subgenome `*_test.plink2` under
 * `process/<mod>/`. Emits PLINK2 .afreq / .vmiss, site-count TSV, MAF-bin histogram,
 * and a GQ placeholder (FORMAT/GQ is not available on this path).
 */
process plink2_assess_debug_slice {
    tag "plink2 assess ${id} chr ${chr_list}"
    label 'cpus_8'
    cpus 8
    memory '32.GB'
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/assess/${params.mod ?: 'misc'}", mode: 'copy', pattern: "*.afreq"
    publishDir "${params.output_dir}/${params.job}/assess/${params.mod ?: 'misc'}", mode: 'copy', pattern: "*.vmiss"
    publishDir "${params.output_dir}/${params.job}/assess/${params.mod ?: 'misc'}", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(id), val(pfile_base), val(chr_list)

    output:
    tuple val(id),
        path("${id}.assess.afreq"),
        path("${id}.assess.vmiss"),
        path("${id}.counts.tsv"),
        path("${id}.mac_site_histogram.tsv"),
        path("${id}.gq_summary.tsv"),
        emit: bundle

    script:
    def procRoot = "${params.output_dir}/${params.job}/process/${params.mod}"
    """
    set -euo pipefail
    plink2 --threads ${task.cpus} \\
        --pfile ${procRoot}/${pfile_base} \\
        --chr ${chr_list} \\
        --freq \\
        --missing \\
        --out ${id}.assess

    n_var=\$(awk 'NR>1{c++} END{print c+0}' ${id}.assess.vmiss)
    n_samples=\$(awk 'NR>1{c++} END{print c+0}' ${procRoot}/${pfile_base}.psam)
    printf "id\\tn_variants\\tn_samples\\n%s\\t%s\\t%s\\n" "${id}" "\${n_var}" "\${n_samples}" > ${id}.counts.tsv

    awk -F'\\t' 'BEGIN {
        OFS="\\t"
        print "MAF_bin\\tn_sites"
        b["maf_eq_0"]=0; b["maf_0_1e3"]=0; b["maf_1e3_01"]=0; b["maf_ge_01"]=0; b["maf_missing"]=0
    }
    NR==1 { next }
    {
        af=\$5
        if (af == "." || af == "") { b["maf_missing"]++; next }
        af=af+0
        maf = (af > 0.5) ? (1 - af) : af
        if (maf <= 0) b["maf_eq_0"]++
        else if (maf < 0.001) b["maf_0_1e3"]++
        else if (maf < 0.01) b["maf_1e3_01"]++
        else b["maf_ge_01"]++
    }
    END {
        print "maf_eq_0\\t" b["maf_eq_0"]
        print "maf_0_1e3\\t" b["maf_0_1e3"]
        print "maf_1e3_01\\t" b["maf_1e3_01"]
        print "maf_ge_01\\t" b["maf_ge_01"]
        print "maf_missing\\t" b["maf_missing"]
    }' ${id}.assess.afreq > ${id}.mac_site_histogram.tsv

    printf 'metric\\tvalue\\nn\\t0\\nmean\\tNA\\nmin\\tNA\\nmax\\tNA\\n' > ${id}.gq_summary.tsv
    """
}

// -- legacy VCF/bcftools processes (kept for optional VCF-based tracks) --

process quick_count {
    tag "quick count: ${id}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/assess/${params.mod ?: 'misc'}", mode: 'copy'

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



process bcftools_qc_assess {
    tag "assess with bcftools: ${id}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    errorStrategy 'terminate'
    publishDir "${params.output_dir}/${params.job}/assess/${params.mod ?: 'misc'}", mode: 'copy'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}.maf_missing.tsv"), emit: maf_missing
    tuple val(id), path("${id}.gq_summary.tsv"), emit: gq_summary

    script:
    """
    set -euo pipefail
    # Do not request MAF in +fill-tags when INFO/MAF is already present (PLINK2 export).
    bcftools +fill-tags ${vcf} -- -t F_MISSING,AN,AC | \\
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/MAF\\t%INFO/F_MISSING\\t%QUAL\\n' > ${id}.maf_missing.tsv

    if bcftools view -h ${vcf} | grep -q '##FORMAT=<ID=GQ'; then
        bcftools query -f '[%GQ\\t]\\n' ${vcf} | tr '\\t' '\\n' | awk 'NF && \$1 != "."' > all.gq || true
        if [ ! -s all.gq ]; then
            printf 'metric\\tvalue\\nn\\t0\\nmean\\tNA\\nmin\\tNA\\nmax\\tNA\\n' > ${id}.gq_summary.tsv
        else
            awk 'BEGIN{min=1e9;max=-1e9}{s+=\$1; n++; if(\$1<min)min=\$1; if(\$1>max)max=\$1} END{if(n==0){print "metric\\tvalue"; print "n\\t0"; exit} mean=s/n; print "metric\\tvalue"; print "n\\t"n; print "mean\\t"mean; print "min\\t"min; print "max\\t"max}' all.gq > ${id}.gq_summary.tsv
        fi
    else
        printf 'metric\\tvalue\\nn\\t0\\nmean\\tNA\\nmin\\tNA\\nmax\\tNA\\n' > ${id}.gq_summary.tsv
    fi
    """
}

process dumpnice_vcf_qc_assess {
    tag "dumpnice vcf qc assess: ${id}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/assess/${params.mod ?: 'misc'}/qc", mode: 'copy'

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
