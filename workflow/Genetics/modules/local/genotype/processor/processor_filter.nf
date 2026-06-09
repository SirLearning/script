nextflow.enable.dsl=2

process filter_sample_plink {
    tag "plink sample filter ${id}"
    publishDir "${params.output_dir}/${params.job}/process/logs", mode: 'copy', pattern: "*.log"
    publishDir "${params.output_dir}/${params.job}/process/filters", mode: 'copy', pattern: "*.{id}"

    input:
    tuple val(id), val(chr), path(prefix), path(pgen), path(psam), path(pvar)
    tuple val(id2), val(chr2), val(mind_th), path(dedup_id)

    output:
    tuple val(id), val(chr), val("${id}.filter.s"), path("${id}.filter.s.id"), emit: filter_s
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    set -euo pipefail
    exec > plink_sample_filter.${chr}.log 2>&1

    echo "Filtering PLINK samples..."
    echo "Threshold: duplicate germplasm samples, F_MISSING<=${mind_th}"

    prefix_out="${id}.filter.s"

    plink2 --pfile ${prefix} \\
        --mind ${mind_th} \\
        --remove ${dedup_id} \\
        --threads ${task.cpus} \\
        --out "\${prefix_out}"
    
    if [ -f "\${prefix_out}.mindrem.id" ]; then
        cat "\${prefix_out}.mindrem.id" > "\${prefix_out}.id"
    fi

    cat ${dedup_id} >> "\${prefix_out}.id"

    echo "PLINK sample filtering completed for ${id}."
    """
}

process filter_sample_after_variant_plink {
    tag "plink sample filter after variant ${id}"
    publishDir "${params.output_dir}/${params.job}/process/logs", mode: 'copy', pattern: "*.log"
    publishDir "${params.output_dir}/${params.job}/process/filters", mode: 'copy', pattern: "*.{id}"

    input:
    tuple val(id), val(chr), path(prefix), path(pgen), path(psam), path(pvar)
    tuple val(id2), val(chr2), val(filter_v_prefix), path(filter_in), path(filter_out)
    tuple val(mind_th), path(dedup_id)

    output:
    tuple val(id), val(chr), val("${id}.filter.v.s"), path("${id}.filter.v.s.id"), emit: filter_v_s
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    set -euo pipefail
    exec > plink_sample_filter_after_variant.${chr}.log 2>&1

    echo "Filtering PLINK samples after variant filtering..."
    echo "Removing samples in ${filter_v_prefix}"

    prefix_out="${id}.filter.v.s"

    plink2 --pfile ${prefix} \\
        --mind ${mind_th} \\
        --remove ${dedup_id} \\
        --extract ${filter_in} \\
        --threads ${task.cpus} \\
        --out "\${prefix_out}"

    if [ -f "\${prefix_out}.mindrem.id" ]; then
        cat "\${prefix_out}.mindrem.id" > "\${prefix_out}.id"
    fi

    cat ${dedup_id} >> "\${prefix_out}.id"

    echo "PLINK sample filtering after variant filtering completed for ${chr}."
    """
}

process filter_variant_plink {
    tag "plink filter ${id}"
    publishDir "${params.output_dir}/${params.job}/process/logs", mode: 'copy', pattern: "*.log"
    publishDir "${params.output_dir}/${params.job}/process/filters", mode: 'copy', pattern: "*.{in,out}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(maf_th), val(mac_th), val(geno_th), val(ld_th)
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), val("${id}.filter.v"), path("${id}.filter.v.in"), path("${id}.filter.v.out"), emit: filter_v
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    set -euo pipefail
    exec > ${chr}.plink_filter.log 2>&1

    echo "Filtering PLINK files..."
    echo "Thresholds: MAF>=${maf_th}, MAC>=${mac_th}, F_GENO>=${geno_th}, LD r2<=${ld_th}"

    prefix_out="${id}.filter.v"

    plink2 --pfile ${prefix} \\
        --maf ${maf_th} \\
        --mac ${mac_th} \\
        --geno ${geno_th} \\
        --indep-pairwise 50 5 ${ld_th} \\
        --threads ${task.cpus} \\
        --out "\${prefix_out}"

    cp "\${prefix_out}.prune.in" "\${prefix_out}.in"
    cp "\${prefix_out}.prune.out" "\${prefix_out}.out"

    echo "PLINK filtering completed for ${chr}."
    """
}

process filter_variant_after_sample_plink {
    tag "plink var filter after sample ${id}"
    publishDir "${params.output_dir}/${params.job}/process/logs", mode: 'copy', pattern: "*.log"
    publishDir "${params.output_dir}/${params.job}/process/filters", mode: 'copy', pattern: "*.{in,out}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(maf_th), val(mac_th), val(geno_th), val(ld_th)
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)
    tuple val(id2), val(chr2), val(filter_s_prefix), path(filter_id)

    output:
    tuple val(id), val(chr), val("${id}.filter.s.v"), path("${id}.filter.s.v.in"), path("${id}.filter.s.v.out"), emit: filter_s_v
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    set -euo pipefail
    exec > ${chr}.plink_var_filter_after_sample.log 2>&1

    echo "Filtering PLINK variants after sample filtering..."
    echo "Thresholds: remove samples in ${filter_s_prefix}, MAF>=${maf_th}, MAC>=${mac_th}, F_GENO>=${geno_th}, LD r2<=${ld_th}"

    prefix_out="${id}.filter.s.v"

    plink2 --pfile ${prefix} \\
        --remove ${filter_id} \\
        --maf ${maf_th} \\
        --mac ${mac_th} \\
        --geno ${geno_th} \\
        --indep-pairwise 50 5 ${ld_th} \\
        --threads ${task.cpus} \\
        --out "\${prefix_out}"

    cp "\${prefix_out}.prune.in" "\${prefix_out}.in"
    cp "\${prefix_out}.prune.out" "\${prefix_out}.out"

    echo "PLINK variant filtering after sample filtering completed for ${chr}."
    """
}

// v0 is rigid threshold filtering
process filter_vcf_v0_vcftools {
    tag "${id}" ? "vcftools filter ${id}" : 'vcftools filter'
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(id), val(chr), path(vcf)

    output:
    tuple val(id), val(chr), val("${id}.vcftools.flt.v0"), path("${id}.vcftools.flt.v0.vcf"), emit: vcf

    script:
    """
    set -euo pipefail
    exec > vcftools_filter.${chr}.log 2>&1

    in="${vcf}"
    out="${id}.vcftools.flt.v0"

    echo "Filtering VCF with vcftools... Thresholds are:"
    echo "removing indels"
    echo "MAF>=${params.maf}, MAC>=${params.mac}"
    echo "F_MISSING<=${params.max_missing}"
    echo "Min alleles: ${params.min_alleles}, Max alleles: ${params.max_alleles}"

    echo "run with vcftools..."
    vcftools --gzvcf \${in} \\
        --max-missing ${params.max_missing} \\
        --remove-indels \\
        --maf ${params.maf} \\
        --min-alleles ${params.min_alleles} \\
        --max-alleles ${params.max_alleles} \\
        --mac ${params.mac} \\
        --recode \\
        --recode-INFO-all \\
        --out "\${out}"

    echo "VCF filtering completed for ${chr}."
    mv "\${out}.recode.vcf" "\${out}.vcf"
    echo "Filtered and recoded VCF saved to \${out}.vcf"
    """
}

process filter_vcf_v0_bcftools {
    tag "${id}" ? "bcftools filter ${id}" : 'bcftools filter'
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(id), val(chr), path(vcf)

    output:
    tuple val(id), val(chr), val("${id}.bcftools.flt.v0"), path("${id}.bcftools.flt.v0.vcf"), emit: vcf

    script:
    """
    set -euo pipefail
    exec > bcftools_filter.${chr}.log 2>&1

    in="${vcf}"
    out="${id}.bcftools.flt.v0"

    echo "Filtering VCF with bcftools..."
    echo "thresholds: MAF>=${params.maf}, MAC>=${params.mac}, F_MISSING<=${params.max_missing}"

    echo "run with bcftools..."
    bcftools view \\
        -m2 -M2 -v snps \\
        -i 'MAF>=${params.maf} && AC>=${params.mac} && F_MISSING<=${params.max_missing}' \\
        -o ${id}.bcftools.flt.v0.vcf \\
        ${vcf}
    
    echo "VCF filtering completed for ${chr}."
    echo "Filtered VCF saved to ${id}.bcftools.flt.v0.vcf"
    """
}
