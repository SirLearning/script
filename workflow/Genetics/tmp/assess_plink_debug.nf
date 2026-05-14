nextflow.enable.dsl=2

/*
 * Debug assess track for test_thin / test_common_thin:
 * export a narrow PLINK2 chromosome slice from existing *_test.plink2, then run
 * quick_count + bcftools_qc_assess (assess.nf) plus MAF-bin site counts (§9 backlog).
 *
 * Per subgenome, a single representative PLINK chromosome is exported (small slice):
 *   A=1, B=3, D=5, Others=0 (matches *_test.plink2 CHROM layout).
 */
include { quick_count; bcftools_qc_assess } from '../genotype/assess.nf'

process plink2_export_subgenome_chr_vcf {
    tag "export vcf ${id} (${chr_list})"
    label 'cpus_8'
    cpus 8
    memory '32.GB'
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/assess/${params.mod}/export", mode: 'copy', pattern: "*.vcf.gz"

    input:
    tuple val(id), val(pfile_base), val(chr_list)

    output:
    tuple val(id), path("${id}.debug.vcf.gz"), emit: vcf

    script:
    def procRoot = "${params.output_dir}/${params.job}/process/${params.mod}"
    """
    set -euo pipefail
    plink2 --threads ${task.cpus} \\
        --pfile ${procRoot}/${pfile_base} \\
        --chr ${chr_list} \\
        --export vcf bgz \\
        --out ${id}.dbg
    test -s ${id}.dbg.vcf.gz
    mv ${id}.dbg.vcf.gz ${id}.debug.raw.gz
    bcftools view ${id}.debug.raw.gz | awk '!/^##chrSet=/' | bcftools view -Oz -o ${id}.debug.vcf.gz
    rm -f ${id}.debug.raw.gz
    test -s ${id}.debug.vcf.gz
    """
}

process bcftools_mac_site_histogram {
    tag "site MAF/AC ${id}"
    label 'cpus_4'
    cpus 4
    memory '16.GB'
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/assess/${params.mod}/info", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}.mac_site_histogram.tsv"), emit: tsv

    script:
    """
    set -euo pipefail
    bcftools query -f '%INFO/MAF\\t%INFO/AC\\n' ${vcf} | awk -F'\\t' '
    BEGIN {
        print "MAF_bin\\tn_sites"
        b["maf_eq_0"] = 0
        b["maf_0_1e3"] = 0
        b["maf_1e3_01"] = 0
        b["maf_ge_01"] = 0
        b["maf_missing"] = 0
    }
    {
        maf = \$1 + 0
        if (\$1 == "." || \$1 == "") { b["maf_missing"]++; next }
        if (maf <= 0) b["maf_eq_0"]++
        else if (maf < 0.001) b["maf_0_1e3"]++
        else if (maf < 0.01) b["maf_1e3_01"]++
        else b["maf_ge_01"]++
    }
    END {
        for (k in b) print k\"\\t\"b[k]
    }' | sort -k1,1 > ${id}.mac_site_histogram.tsv
    """
}

workflow {
    if (!(params.mod in ['test_thin', 'test_common_thin'])) {
        error "assess_plink_debug.nf: params.mod must be test_thin or test_common_thin (got: ${params.mod})."
    }
    if (!params.output_dir || !params.job) {
        error "params.output_dir and params.job are required."
    }

    def sub = ['A', 'B', 'D', 'Others']
    def chrPick = [
        'A'     : '1',
        'B'     : '3',
        'D'     : '5',
        'Others': '0',
    ]
    def ch_in = channel.from(sub.collect { sg -> tuple(sg, "${sg}_test.plink2", chrPick[sg]) })
    plink2_export_subgenome_chr_vcf(ch_in)
    def vcf_ch = plink2_export_subgenome_chr_vcf.out.vcf
    def id_ch = vcf_ch.map { id, vcf -> tuple("${params.mod}__${id}", vcf) }

    quick_count(id_ch)
    bcftools_qc_assess(id_ch)
    bcftools_mac_site_histogram(vcf_ch)
}
