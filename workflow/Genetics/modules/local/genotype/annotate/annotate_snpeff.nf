nextflow.enable.dsl=2

/*
 * SnpEff against a custom IWGSC v1.1 Lulab genome (same GFF/FASTA as VEP).
 * Conda env `anno`. Chromosome 0 is renamed to params.vep_unplaced_name so the
 * VEP sites VCF and this database agree (Perl-style "0" is false in VEP only,
 * but the published sites VCF already uses the alias).
 */

process snpeff_build_genome {
    tag "snpeff build ${params.snpeff_genome}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/annotate/snpeff_genome", mode: 'copy'
    conda "${params.user_dir}/miniconda3/envs/anno"
    label 'snpeff_build'

    input:
    path gff
    tuple path(fasta), path(fasta_fai), path(fasta_gzi)

    output:
    path "snpEff.config", emit: config
    path "data", emit: data
    path "*.log", emit: logs

    script:
    def unplaced = params.vep_unplaced_name
    def genome = params.snpeff_genome
    """
    set -euo pipefail
    exec > snpeff_build.log 2>&1

    mkdir -p data/${genome}
    grep -v '^#' ${gff} \\
        | awk -F'\\t' -v OFS='\\t' -v unplaced='${unplaced}' \\
            '\$1 == "0" { \$1 = unplaced } { print }' \\
        > data/${genome}/genes.gff

    bgzip -dc ${fasta} \\
        | awk -v unplaced='${unplaced}' \\
            '/^>/ { if (\$0 == ">0") \$0 = ">" unplaced } { print }' \\
        > data/${genome}/sequences.fa

    cat > snpEff.config << EOF
data.dir = ./data/
${genome}.genome : IWGSC Chinese Spring RefSeq v1.1 Lulab
EOF

    snpEff -Xmx${params.snpeff_build_xmx} build -gff3 -noCheckCds -noCheckProtein -v ${genome}
    """
}

process snpeff_annotate_vcf {
    tag "snpeff: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/variant", mode: 'copy', pattern: "*.snpeff.vcf.gz"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/anno"
    label 'snpeff_ann'

    input:
    tuple val(id), path(vcf)
    path snpeff_config
    path snpeff_data

    output:
    tuple val(id), path("${id}.snpeff.vcf.gz"), emit: vcf
    path "*.log", emit: logs

    script:
    def jar = params.snpeff_jar
    if (!jar && params.home_dir) {
        def bundled = file("${params.home_dir}/00data/03ref/07snpeff/snpeff_4.3t_dist/snpEff/snpEff.jar")
        if (bundled.exists()) {
            jar = bundled.toString()
        }
    }
    def launcher = jar
        ? "java -Xmx${params.snpeff_ann_xmx} -jar ${jar}"
        : "snpEff -Xmx${params.snpeff_ann_xmx}"
    """
    set -euo pipefail
    exec > ${id}.snpeff.log 2>&1

    if [ ! -e data ]; then ln -s ${snpeff_data} data; fi

    # FastCall writes ALT I/D (unknown-sequence indel). SnpEff 4.3 treats a
    # single-letter ALT as an IUB SNP code and aborts on I. Recode, annotate,
    # then put the original allele back so the published VCF matches the sites VCF.
    bgzip -dc ${vcf} | awk -F'\\t' -v OFS='\\t' '
        /^#/ { print; next }
        \$5 == "I" {
            \$8 = (\$8 == "." ? "FC=I" : \$8 ";FC=I")
            \$5 = \$4 "N"
            print
            next
        }
        \$5 == "D" {
            \$8 = (\$8 == "." ? "FC=D" : \$8 ";FC=D")
            \$5 = "<DEL>"
            print
            next
        }
        { print }
    ' > recoded.vcf

    ${launcher} \\
        -c ${snpeff_config} \\
        -no-downstream -no-upstream -no-intergenic \\
        -noStats -noLog \\
        ${params.snpeff_genome} \\
        recoded.vcf \\
        | awk -F'\\t' -v OFS='\\t' '
            /^#/ { print; next }
            {
                if (\$8 ~ /(^|;)FC=I(;|\$)/) \$5 = "I"
                else if (\$8 ~ /(^|;)FC=D(;|\$)/) \$5 = "D"
                gsub(/;FC=[ID]/, "", \$8)
                gsub(/^FC=[ID];/, "", \$8)
                if (\$8 == "FC=I" || \$8 == "FC=D") \$8 = "."
                print
            }
        ' \\
        | bgzip -c -@ ${task.cpus} > ${id}.snpeff.vcf.gz
    """
}

process snpeff_vcf_to_tab {
    tag "snpeff tab: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/variant", mode: 'copy', pattern: "*.snpeff.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}.snpeff.tsv"), emit: table
    path "*.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("${id}.snpeff_tab.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.snpeff import run_snpeff_vcf_to_tab

    run_snpeff_vcf_to_tab("${vcf}", "${id}")
    """
}
