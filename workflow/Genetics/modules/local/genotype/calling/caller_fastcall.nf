nextflow.enable.dsl=2

process fastcall3_disc {
    tag "disc_${chromosome}"
    publishDir("${params.output_dir}/${params.job}/ing", mode: 'copy', overwrite: true)
    conda "${params.user_dir}/miniconda3/envs/tiger"
    label 'fastcall3_disc'
    
    input:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map)
    tuple val(tiger_jar), val(app_name), val(java_version), val(java_lib_dir)
    val samtools_path
    val ing_path
    
    output:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map), emit: chr_config
    path "disc_${chromosome}.log", emit: log
    
    script:
    """
    set -euo pipefail
    exec > disc_${chromosome}.log 2>&1
    
    echo "Starting disc analysis for chromosome ${chromosome}..." 
    echo "Using ${java_version} for ${app_name} disc process"
    
    # Monitor system resources
    echo "System resources before TIGER execution:" 
    echo "Memory: \$(free -h | grep Mem)" 
    echo "CPU cores (allocated): ${task.cpus}" 
    echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g" 
    echo "TIGER jar: ${tiger_jar}" 
    echo "FastCall version: ${app_name}" 
    
    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod disc \\
        -a ${reference} \\
        -b ${taxa_bam_map} \\
        -c 0 \\
        -d ${params.disc_min_MQ} \\
        -e ${params.disc_min_BQ} \\
        -f ${params.disc_min_DC} \\
        -g ${params.disc_min_DR} \\
        -h ${params.disc_max_DR} \\
        -i ${params.disc_HoR} \\
        -j ${params.disc_HeR} \\
        -k ${params.disc_TDR} \\
        -l ${chromosome} \\
        -m ${task.cpus} \\
        -n ${ing_path} \\
        -o ${samtools_path}

    echo "Disc analysis for chromosome ${chromosome} completed"
    """
}

process fastcall3_blib {
    tag "blib_${chromosome}"
    publishDir("${params.output_dir}/${params.job}/vLib", mode: 'copy', overwrite: true)
    conda "${params.user_dir}/miniconda3/envs/tiger"
    label 'fastcall3_blib'

    input:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map)
    tuple val(tiger_jar), val(app_name), val(java_version), val(java_lib_dir)
    val ing_path
    val vLib_path
    
    output:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map), path("${chromosome}_*.lib.gz"), emit: chr_config_lib
    path "blib_${chromosome}.log", emit: log
    
    script:
    """
    set -euo pipefail
    exec > blib_${chromosome}.log 2>&1
    
    echo "Starting blib generation for chromosome ${chromosome}..." 
    echo "Using ${java_version} for ${app_name} blib process" 
    
    # Check if disc files exist
    # Check if disc files exist under ing_path for this chromosome (search recursively)
    if ! find ${ing_path} -type f -name "${chromosome}_*.ing.gz" | grep -q . ; then
        echo "Error: No .ing.gz files found for chromosome ${chromosome} under ${ing_path}"
        exit 1
    fi
    echo "Counting input disc files under ${ing_path}..."
    cnt_ing=\$(find ${ing_path} -type f -name "${chromosome}_*.ing.gz" | wc -l)
    echo "Found \$cnt_ing .ing.gz chunks for chromosome ${chromosome}"
    
    echo "Input disc files (showing up to 5 examples):" 
    find ${ing_path} -type f -name "${chromosome}_*.ing.gz" | head -n 5 || true
    
    # Monitor system resources
    echo "System resources before TIGER execution:" 
    echo "Memory: \$(free -h | grep Mem)" 
    echo "CPU cores (allocated): ${task.cpus}" 
    echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g" 
    echo "TIGER jar: ${tiger_jar}" 
    echo "FastCall version: ${app_name}" 
    
    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod blib \\
        -a ${reference} \\
        -b ${chromosome} \\
        -c ${params.blib_MAO} \\
        -d ${task.cpus} \\
        -e ${ing_path} \\
        -f ${vLib_path}
    
    echo "Blib generation for chromosome ${chromosome} completed"

    # Symlink produced lib for this chromosome back into workdir so Nextflow can capture it
    if ls ${vLib_path}/${chromosome}_*.lib.gz >/dev/null 2>&1; then
        ln -sf ${vLib_path}/${chromosome}_*.lib.gz . || true
    else
        # If naming differs, just link the latest .lib.gz
        latest_lib=\$(ls -1t ${vLib_path}/*.lib.gz 2>/dev/null | head -n1)
        [ -n "\$latest_lib" ] && ln -sf "\$latest_lib" . || true
    fi
    """
}

process fastcall3_scan {
    tag "scan_${chromosome}"
    publishDir("${params.output_dir}/${params.job}/gen", mode: 'move', overwrite: false)
    conda "${params.user_dir}/miniconda3/envs/tiger"
    label 'fastcall3_scan'

    input:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map), path(blib_file)
    tuple val(tiger_jar), val(app_name), val(java_version), val(java_lib_dir)
    val samtools_path
    val gen_path
    
    output:
    // tuple val(chromosome), path("*.vcf*"), emit: vcf_files
    path "scan_${chromosome}.log", emit: log
    
    script:
    // FIX: Typo 'param.home_dir' -> 'params.home_dir'
    // def chr_config = getChrConfig(chromosome, params.home_dir)
    def chr_int = chromosome.toString().toInteger()
    def vcf_name = (chr_int < 10) ? "chr00${chr_int}.vcf" : "chr0${chr_int}.vcf"
    """
    set -euo pipefail
    exec > scan_${chromosome}.log 2>&1
    
    echo "Starting scan analysis for chromosome ${chromosome}..."
    echo "Using ${java_version} for ${app_name} scan process"
    
    # Resolve the lib file (staged by Nextflow) for this chromosome
    lib_file_path="${blib_file}"
    if [ ! -f "\$lib_file_path" ]; then
        echo "Error: Library file \$lib_file_path not found"
        exit 1
    fi
    echo "Input library file: \$lib_file_path"
    echo "Library file size: \$(stat -c%s \"\$lib_file_path\") bytes"
    
    # Monitor system resources
    echo "System resources before TIGER execution:"
    echo "Memory: \$(free -h | grep Mem)"
    echo "CPU cores (allocated): ${task.cpus}"
    echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g"
    echo "TIGER jar: ${tiger_jar}"
    echo "FastCall version: ${app_name}"

    java -Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod scan \\
        -a ${reference} \\
        -b ${taxa_bam_map} \\
        -c "\$lib_file_path" \\
        -d ${chromosome} \\
        -e 0 \\
        -f ${params.scan_min_MQ} \\
        -g ${params.scan_min_BQ} \\
        -h ${params.scan_error_rate} \\
        -i ${samtools_path} \\
        -j ${task.cpus} \\
        -k ${gen_path}
    
    bgzip -f -@ ${task.cpus} ${gen_path}/VCF/${vcf_name}

    echo "Scan analysis for chromosome ${chromosome} completed"
    """
}

process fastcall3_scan2 {
    tag "scan2_${chromosome}"
    publishDir("${params.output_dir}/${params.job}/gen", mode: 'move', overwrite: false)
    conda "${params.user_dir}/miniconda3/envs/tiger"
    label 'fastcall3_scan2'
    
    input:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map), path(blib_file)
    tuple val(tiger_jar), val(app_name), val(java_version), val(java_lib_dir)
    val samtools_path
    val gen_path
    
    output:
    // tuple val(chromosome), path("*.vcf*"), emit: vcf_files
    path "scan_${chromosome}.log", emit: log
    
    script:
    // FIX: Typo 'param.home_dir' -> 'params.home_dir'
    // def chr_config = getChrConfig(chromosome, params.home_dir)
    def chr_int = chromosome.toString().toInteger()
    def vcf_name = (chr_int < 10) ? "chr00${chr_int}.vcf" : "chr0${chr_int}.vcf"
    """
    set -euo pipefail
    exec > scan_${chromosome}.log 2>&1

    echo "Starting scan analysis for chromosome ${chromosome}..."
    echo "Using ${java_version} for ${app_name} scan process"
    
    # Resolve the lib file (staged by Nextflow) for this chromosome
    lib_file_path="${blib_file}"
    if [ ! -f "\$lib_file_path" ]; then
        echo "Error: Library file \$lib_file_path not found"
        exit 1
    fi
    echo "Input library file: \$lib_file_path"
    echo "Library file size: \$(stat -c%s \"\$lib_file_path\") bytes"
    
    # Monitor system resources
    echo "System resources before TIGER execution:"
    echo "Memory: \$(free -h | grep Mem)"
    echo "CPU cores (allocated): ${task.cpus}"
    echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g"
    echo "TIGER jar: ${tiger_jar}"
    echo "FastCall version: ${app_name}"

    java -Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod scan2 \\
        -a ${reference} \\
        -b ${taxa_bam_map} \\
        -c "\$lib_file_path" \\
        -d ${chromosome} \\
        -e 0 \\
        -f ${params.scan_min_MQ} \\
        -g ${params.scan_min_BQ} \\
        -h ${params.scan_error_rate} \\
        -i ${samtools_path} \\
        -j ${task.cpus} \\
        -k ${gen_path}
    
    bgzip -f -@ ${task.cpus} ${gen_path}/VCF/${vcf_name}

    echo "Scan analysis for chromosome ${chromosome} completed"
    """
}
