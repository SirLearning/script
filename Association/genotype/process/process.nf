nextflow.enable.dsl=2

// --- Filter Parameters ---
// standard filter parameters (std)
def maf = 0.05
def mac = 2
def min_alleles = 2
def max_alleles = 2
def max_missing = 0.05
// additional parameters can be added here
def qual = 30
def hwe_pval = 1e-6 // not sure if used

workflow process {
    take:
    // Expect a combined channel: [ val(meta), path(vcf), val(job_config) ]
    vcf_in

    main:
    // filter
    filtered_vcf = vcftools_std_filter(vcf_in)
    // Receive the full tuple and capture the single output channel directly
    gz_vcf = vcf_zip_idx(filtered_vcf.vcf)
    // Pass the resulting (meta, vcf.gz) tuples directly to PLINK
    plink_out = vcf_to_plink(gz_vcf)

    emit:
    vcf = gz_vcf.vcf
    plink_bfile = plink_out.plink_bfile
    plink_ped  = plink_out.plink_ped
}

process vcftools_std_filter {
    tag { meta && meta.id ? "vcftools filter ${meta.id}" : 'vcftools filter' }
    publishDir 'output/process', mode: 'copy'

    input:
    tuple val(meta), path(vcf), val(job_config)

    output:
    tuple val(meta) val("${meta.id}.std.filtered"), path("${meta.id}.std.filtered.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    in="${vcf}"
    out="${meta.id}.std.filtered"

    vcftools --gzvcf "\${in}" \\
        --max-missing ${max_missing} \\
        --remove-indels \\
        --maf ${maf} \\
        --min-alleles ${min_alleles} \\
        --max-alleles ${max_alleles} \\
        --mac ${mac} \\
        --recode \\
        --recode-INFO-all \\
        --out "\${out}"

    mv "\${out}.recode.vcf" "\${out}.vcf"
    """
}

process bcftools_std_filter {
    tag { meta && meta.id ? "bcftools filter ${meta.id}" : 'bcftools filter' }
    publishDir 'output/process', mode: 'copy'

    input:
    tuple val(meta), path(vcf), val(job_config)

    output:
    tuple val(meta), val("${meta.id}.std.filtered"), path("${meta.id}.std.filtered.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    bcftools view \\
        -m2 -M2 -v snps \\
        -i 'MAF>=${maf} && AC>=${mac} && && F_MISSING<=${max_missing}' \\
        -o ${meta.id}.std.filtered.vcf \\
        ${vcf}
    """
}

process vcf_zip_idx {
    tag { meta && meta.id ? "bgzip/tabix ${meta.id}" : 'bgzip/tabix' }
    // Use copy to keep files available for downstream processes (move would remove them before staging)
    publishDir 'output/process', mode: 'copy'

    input:
    tuple val(meta), val(prefix), path(vcf)

    output:
    tuple val(meta), val(prefix), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf

    script:
    """
    set -euo pipefail

    in=${vcf}
    out=${prefix}.vcf.gz
    echo "Processing VCF for work: ${meta.id}" >&2

    if [[ "\${in}" == *.vcf.gz ]]; then
        echo "Input VCF is already gzipped: \${in}" >&2
        in_base=\$(basename "\${in}")
        # Avoid self-linking when names are identical
        if [[ "\${in_base}" != "\${out}" ]]; then
            ln -sf "\${in}" "\${out}" || cp -f "\${in}" "\${out}"
        fi
    else
        echo "Compressing VCF with bgzip: \${in} -> \${out}" >&2
        bgzip -c "\${in}" > "\${out}"
    fi

    # Try to create index; if it fails (e.g., not BGZF), re-compress and retry
    if ! tabix -f -p vcf "\${out}"; then
        echo "tabix failed on \${out}; ensuring BGZF compression and retrying" >&2
        if [[ "\${in}" == *.vcf.gz ]]; then
            gunzip -c "\${in}" | bgzip -c > "\${out}.tmp"
        else
            bgzip -c "\${in}" > "\${out}.tmp"
        fi
        mv -f "\${out}.tmp" "\${out}"
        tabix -f -p vcf "\${out}"
    fi
    """
}

process vcf_to_plink {
    tag { meta && meta.id ? "plink ${meta.id}" : 'plink' }
    // Keep generated bed/bim/fam for further steps
    publishDir 'output/process', mode: 'copy'

    input:
    tuple val(meta), val(prefix), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.bed"), path("${prefix}.bim"), path("${prefix}.fam"), emit: plink_bfile
    path("${prefix}.ped"), emit: plink_ped
    script:
    """
    set -euo pipefail
    if [[ ! -s ${vcf} ]]; then
        echo "ERROR: VCF file missing or empty before PLINK: ${vcf}" >&2
        ls -l >&2 || true
        exit 201
    fi

    plink --vcf ${vcf} \\
        --autosome-num 42 \\
        --recode \\
        --double-id \\
        --allow-extra-chr \\
        --out ${prefix} \\
        --threads ${task.cpus}

    plink \\
        --file ${prefix}.ped \\
        --make-bed \\
        --chr-set 42 \\
        --out ${prefix} \\
        --threads ${task.cpus}
    """
}