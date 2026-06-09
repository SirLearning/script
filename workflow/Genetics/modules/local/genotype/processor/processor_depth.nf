nextflow.enable.dsl=2

include { getPopDepTaxaBamFile_v1 } from '../../infra/infra_ref_v1.nf'
include { getRefV1ChrLength } from '../../infra/infra_ref_v1.nf'

process calc_population_depth {
    tag "prepare popdepth: ${chr}"
    publishDir "${params.output_dir}/${params.job}/process/variant", mode: 'copy', pattern: "*.popdep.txt"
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/tiger"
    label 'cpus_32'

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
