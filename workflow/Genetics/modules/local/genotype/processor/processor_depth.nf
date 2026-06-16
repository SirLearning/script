nextflow.enable.dsl=2

include { getPopDepTaxaBamFile_v1 } from '../../infra/infra_ref_v1.nf'
include { getRefV1ChrLength } from '../../infra/infra_ref_v1.nf'

def popdepVariantPublishDir() {
    def root = params.popdep_publish_dir ?: params.popdep_dir
    return root
        ? "${root}/variant"
        : "${params.output_dir}/${params.job}/process/variant"
}

def popdepLogPublishDir() {
    def root = params.popdep_publish_dir ?: params.popdep_dir
    return root
        ? "${root}/logs"
        : "${params.output_dir}/${params.job}/process"
}

process calc_population_depth {
    tag "prepare popdepth: ${chr}"
    publishDir "${popdepVariantPublishDir()}", mode: 'copy', pattern: "*.popdep.txt"
    publishDir "${popdepLogPublishDir()}", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/tiger"
    label 'popdep_tiger'

    input:
    tuple val(id), val(chr), path(vcf)
    tuple path(tiger_jar), val(app_name), val(java_version)

    output:
    tuple val(id), val(chr), path("${id}.popdep.txt"), emit: popdep

    script:
    def tb_file = getPopDepTaxaBamFile_v1(chr, params.home_dir)
    def chr_length = getRefV1ChrLength(chr)
    """
    set -euo pipefail
    exec > popdep_${chr}.log 2>&1

    echo "Starting population depth analysis for chromosome ${chr}..."
    echo "Using ${java_version} for ${app_name} process"
    echo "The environment java version is:"
    which java
    java -version
    
    echo "System resources before TIGER execution:"
    echo "Memory: \$(free -h | grep Mem)"
    echo "CPU cores (allocated): ${task.cpus}"
    echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g"
    echo "TIGER jar: ${tiger_jar}"

    java -Xmx${task.memory.toGiga()}G -jar ${tiger_jar} \\
        -app ${app_name} \\
        -a ${tb_file} \\
        -b ${chr} \\
        -c ${chr_length} \\
        -d samtools \\
        -e ${task.cpus} \\
        -f ${id}.popdep.txt.gz
    
    gzip -d ${id}.popdep.txt.gz

    echo "Population depth analysis completed for chromosome ${chr}."
    """
}

process annotate_subgenome_variant_popdep {
    tag "annotate variant popdep: ${chr}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/variant", mode: 'copy', pattern: "*.popdep.info.tsv"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.popdep.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), path("${id}.popdep.info.tsv"), emit: popdep
    path "${chr}.popdep.log", emit: log

    script:
    def popdep_workers = params.popdep_lookup_workers ?: 0
    def popdep_workers_py = (popdep_workers as int) > 0 ? "${popdep_workers}" : "None"
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.popdep.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.popdep import annotate_variants_popdep_from_pvar

    workers = ${popdep_workers_py}
    print(f"Annotating variant popdep for ${id} (${chr}) from ${params.popdep_dir} (workers={workers}) ...")
    annotate_variants_popdep_from_pvar(
        "${pvar}",
        "${params.popdep_dir}",
        "${id}.popdep.info.tsv",
        max_workers=workers,
    )
    """
}
