nextflow.enable.dsl=2
// Using bcftools annotate to add functional annotations to VCF files

workflow annotate {
    take:
    // Expect a combined channel: [ val(meta), path(vcf), val(job_config) ]
    vcf_in

    main:
    // Annotate VCF using bcftools
    annotated_vcf = annotate_vcf_bcftools(vcf_in)

    emit:
    vcf = annotated_vcf.vcf
}

process annotate_vcf_bcftools {
    tag { id ? "bcftools annotate ${id}" : 'bcftools annotate' }
    publishDir 'output/annotate', mode: 'copy'

    input:
    tuple val(id), path(vcf), val(job_config)

    output:
    tuple val(id), val("${id}.annotated"), path("${id}.annotated.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    in="${vcf}"
    out="${id}.annotated"

    # Example annotation command, modify as needed
    bcftools annotate -a annotations.gff -c CHROM,FROM,TO,GENE -o \${out}.vcf \${in}

    """
}

process prepare_tab_file {
    tag { id ? "prepare tab file ${id}" : 'prepare tab file' }
    publishDir 'output/annotate', mode: 'copy'

    input:
    tuple val(id), path(vcf), val(job_config)

    output:
    tuple val(id), val("${id}.annotations.tab"), path("${id}.annotations.tab"), emit: tab

    script:
    """
    set -euo pipefail

    in="${vcf}"
    out="${id}.annotations.tab"

    # Extract annotations from VCF to a tab-delimited file
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/GENE\\n' \${in} > \${out}

    """
}

process generate_gene_bed {
    tag { id ? "generate gene bed ${id}" : 'generate gene bed' }
    publishDir 'output/annotate', mode: 'copy'

    input:
    tuple val(id), path(vcf), val(job_config)

    output:
    tuple val(id), val("${id}.gene.bed"), path("${id}.gene.bed"), emit: bed

    script:
    """
    set -euo pipefail

    in="${vcf}"
    out="${id}.gene.bed"

    # Extract gene regions from VCF and create BED file
    bcftools query -f '%CHROM\\t%POS\\t%END\\t%GENE\\n' \${in} | \
        awk '{print \$1, \$2-1, \$3, \$4}' OFS='\\t' > \${out}

    """
}
