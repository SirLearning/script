nextflow.enable.dsl=2

process check_existing_files {
    tag "check_existing_files"
    memory '1g'
    cpus 1
    
    input:
    val output_dir
    val file_pattern
    val stage_name
    
    output:
    path "check_existing_files.log", emit: log
    
    script:
    """
    set -euo pipefail
    exec > check_existing_files.log 2>&1

    echo "Checking for existing ${stage_name} files..."
    find ${output_dir}/${stage_name} -name "${file_pattern}" 2>/dev/null | head -5
    count=\$(find ${output_dir}/${stage_name} -name "${file_pattern}" 2>/dev/null | wc -l)
    echo "Found \$count ${stage_name} files"
    """
}

process prepareTaxaBamMap {
    tag "prepare_taxa_bam_map_${pop}"
    memory '16g'
    cpus 4
    publishDir "${params.home_dir}/00data/05taxaBamMap", mode: 'copy'
    conda "${params.user_dir}/miniconda3/envs/dbone"
    
    input:
    tuple path(workshop_jar), val(java_version), val(java_lib_dir)
    tuple val(pop), val(depth_dir), val(bam_dir)
    
    output:
    path "${pop}.taxaBamMap.txt", emit: taxa_bam_map
    path "${pop}.taxaRunMap.txt", emit: taxa_run_map
    path "prepare_tbm_${pop}.log", emit: log
    
    script:
    def output_bam_file = "${pop}.taxaBamMap.txt"
    def output_run_file = "${pop}.taxaRunMap.txt"
    """
    set -euo pipefail
    exec > prepare_tbm_${pop}.log 2>&1

    echo "Preparing taxa-BAM mapping file for population ${pop}..." 
    echo "Using Java 21 for prepareTaxaBamMap process" 
    
    # Use TaxaBamMap.java only
    echo "Running TaxaBamMap from Workshop jar..." 
    
    java -Xmx16g -jar "${workshop_jar}" \\
        -d ${depth_dir} \\
        -b ${bam_dir} \\
        -o ${output_bam_file} \\
        -t ${output_run_file}
    
    echo "TaxaBamMap execution completed" 
    """
}

process concatTaxaBamMap {
    tag "concat_taxa_bam_map_${sub_genome}"
    memory '16g'
    cpus 4
    publishDir "${params.home_dir}/00data/05taxaBamMap", mode: 'copy'
    
    input:
    tuple path(taxa_bam_maps), val(job_name), val(sub_genome)
    
    output:
    tuple val(sub_genome), path("${job_name}.${sub_genome}.taxaBamMap.txt"), emit: tbm
    path "${sub_genome}.concat.log", emit: log
    
    script:
    """
    set -euo pipefail
    exec > ${sub_genome}.concat.log 2>&1

    echo "Concatenating taxa-BAM mapping files (${taxa_bam_maps.size()})..."

    cat ${taxa_bam_maps.join(' ')} > ${job_name}.${sub_genome}.taxaBamMap.txt
    echo -e "Taxa\tCoverage-Of-All-Bams\tBams(A list of bams of the taxon, seperated by the delimiter of Tab)" > ${job_name}.${sub_genome}.taxaBamMap.txt.header
    cat ${job_name}.${sub_genome}.taxaBamMap.txt >> ${job_name}.${sub_genome}.taxaBamMap.txt.header
    mv ${job_name}.${sub_genome}.taxaBamMap.txt.header ${job_name}.${sub_genome}.taxaBamMap.txt

    echo "Concatenation completed. Output: ${job_name}.${sub_genome}.taxaBamMap.txt" 
    """
}
