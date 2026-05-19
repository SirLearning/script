nextflow.enable.dsl=2

workflow database {
    take:
    // Expect a channel: [ val(meta), path(vcf) ]
    vcf_in

    main:
    // Build Hail MatrixTable from VCF
    hail_mt = build_hail_from_vcf(vcf_in)

    emit:
    mt = hail_mt.mt
}

process build_hail_from_vcf {
    tag "build_hail_from_vcf"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.mt"), emit: mt

    script:
    """
    hail_script_build_from_vcf.py \\
        --input_vcf ${vcf} \\
        --output_mt ${meta}.mt
    """
}


process annotate_vcf_bcftools {
    tag "bcftools annotate ${id}"
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
    tag "prepare tab file ${id}"
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
    tag "generate gene bed ${id}"
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
