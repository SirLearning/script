nextflow.enable.dsl=2

def chrVcfFilename(chromosome) {
    def chr_int = chromosome.toString().toInteger()
    return (chr_int < 10) ? "chr00${chr_int}.vcf" : "chr0${chr_int}.vcf"
}

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

    {
        echo "Starting disc analysis for chromosome ${chromosome}..."
        echo "Using ${java_version} for ${app_name} disc process"
        echo "System resources before TIGER execution:"
        echo "Memory: \$(free -h | grep Mem)"
        echo "CPU cores (allocated): ${task.cpus}"
        echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g"
        echo "TIGER jar: ${tiger_jar}"
        echo "FastCall version: ${app_name}"
    } > disc_${chromosome}.log

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
        -o ${samtools_path} \\
        2>> disc_${chromosome}.log

    echo "Disc analysis for chromosome ${chromosome} completed" >> disc_${chromosome}.log
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

    {
        echo "Starting blib generation for chromosome ${chromosome}..."
        echo "Using ${java_version} for ${app_name} blib process"
        if ! find ${ing_path} -type f -name "${chromosome}_*.ing.gz" | grep -q . ; then
            echo "Error: No .ing.gz files found for chromosome ${chromosome} under ${ing_path}"
            exit 1
        fi
        echo "Counting input disc files under ${ing_path}..."
        cnt_ing=\$(find ${ing_path} -type f -name "${chromosome}_*.ing.gz" | wc -l)
        echo "Found \$cnt_ing .ing.gz chunks for chromosome ${chromosome}"
        echo "Input disc files (showing up to 5 examples):"
        find ${ing_path} -type f -name "${chromosome}_*.ing.gz" | head -n 5 || true
        echo "System resources before TIGER execution:"
        echo "Memory: \$(free -h | grep Mem)"
        echo "CPU cores (allocated): ${task.cpus}"
        echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g"
        echo "TIGER jar: ${tiger_jar}"
        echo "FastCall version: ${app_name}"
    } > blib_${chromosome}.log

    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod blib \\
        -a ${reference} \\
        -b ${chromosome} \\
        -c ${params.blib_MAO} \\
        -d ${task.cpus} \\
        -e ${ing_path} \\
        -f ${vLib_path} \\
        2>> blib_${chromosome}.log

    echo "Blib generation for chromosome ${chromosome} completed" >> blib_${chromosome}.log

    if ls ${vLib_path}/${chromosome}_*.lib.gz >/dev/null 2>&1; then
        ln -sf ${vLib_path}/${chromosome}_*.lib.gz . || true
    else
        latest_lib=\$(ls -1t ${vLib_path}/*.lib.gz 2>/dev/null | head -n1)
        [ -n "\$latest_lib" ] && ln -sf "\$latest_lib" . || true
    fi
    """
}

process fastcall3_vlib {
    tag "vlib_${chromosome}"
    publishDir("${params.output_dir}/${params.job}/vLib", mode: 'copy', overwrite: true)
    conda "${params.user_dir}/miniconda3/envs/tiger"
    label 'fastcall3_vlib'

    input:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map), path(blib_file)
    tuple val(tiger_jar), val(app_name), val(java_version), val(java_lib_dir)
    val vLib_path

    output:
    path "vlib_${chromosome}.log", emit: log

    script:
    """
    set -euo pipefail
    lib_file_path="${blib_file}"
    txt_out="${vLib_path}/\$(basename "\$lib_file_path" .lib.gz).txt"

    {
        echo "Starting vlib (library to text) for chromosome ${chromosome}..."
        echo "Using ${java_version} for ${app_name} vlib process"
        echo "Input library file: \$lib_file_path"
        echo "Text output: \$txt_out"
        echo "TIGER jar: ${tiger_jar}"
    } > vlib_${chromosome}.log

    if [ ! -f "\$lib_file_path" ]; then
        echo "Error: Library file \$lib_file_path not found" >> vlib_${chromosome}.log
        exit 1
    fi

    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod vlib \\
        -a "\$lib_file_path" \\
        -b "\$txt_out" \\
        2>> vlib_${chromosome}.log

    echo "Vlib for chromosome ${chromosome} completed" >> vlib_${chromosome}.log
    """
}

process fastcall3_clib {
    tag "clib_${chromosome}"
    publishDir("${params.output_dir}/${params.job}/vLib_custom", mode: 'copy', overwrite: true)
    conda "${params.user_dir}/miniconda3/envs/tiger"
    label 'fastcall3_clib'

    input:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map), path(blib_file)
    tuple val(tiger_jar), val(app_name), val(java_version), val(java_lib_dir)
    path positions
    val clib_path

    output:
    path "clib_${chromosome}.log", emit: log

    script:
    """
    set -euo pipefail
    lib_file_path="${blib_file}"
    mkdir -p ${clib_path}
    custom_lib="${clib_path}/${chromosome}_custom.lib.gz"

    {
        echo "Starting clib (position subset library) for chromosome ${chromosome}..."
        echo "Using ${java_version} for ${app_name} clib process"
        echo "Input library file: \$lib_file_path"
        echo "Positions: ${positions}"
        echo "Output library: \$custom_lib"
        echo "TIGER jar: ${tiger_jar}"
    } > clib_${chromosome}.log

    if [ ! -f "\$lib_file_path" ]; then
        echo "Error: Library file \$lib_file_path not found" >> clib_${chromosome}.log
        exit 1
    fi

    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod clib \\
        -a "\$lib_file_path" \\
        -b ${positions} \\
        -c "\$custom_lib" \\
        2>> clib_${chromosome}.log

    echo "Clib for chromosome ${chromosome} completed" >> clib_${chromosome}.log
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
    path "scan_${chromosome}.log", emit: log
    
    script:
    def vcf_name = chrVcfFilename(chromosome)
    def bgzip_threads = Math.min(8, task.cpus as int)
    """
    set -euo pipefail
    lib_file_path="${blib_file}"
    mkdir -p ${gen_path}/VCF
    vcf_gz="${gen_path}/VCF/${vcf_name}.gz"

    {
        echo "Starting scan analysis for chromosome ${chromosome}..."
        echo "Using ${java_version} for ${app_name} scan process"
        echo "Input library file: \$lib_file_path"
        if [ -f "\$lib_file_path" ]; then
            echo "Library file size: \$(stat -c%s "\$lib_file_path") bytes"
        fi
        echo "Scan output dir (-k, indiCounts): ${gen_path}"
        echo "VCF stdout -> bgzip -c -@ ${bgzip_threads} -> \$vcf_gz"
        echo "MiDR (-m): ${params.scan_midr}"
        echo "System resources before TIGER execution:"
        echo "Memory: \$(free -h | grep Mem)"
        echo "CPU cores (allocated): ${task.cpus}"
        echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g"
        echo "TIGER jar: ${tiger_jar}"
        echo "FastCall version: ${app_name}"
    } > scan_${chromosome}.log

    if [ ! -f "\$lib_file_path" ]; then
        echo "Error: Library file \$lib_file_path not found" >> scan_${chromosome}.log
        exit 1
    fi

    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod scan \\
        -a ${reference} \\
        -b ${taxa_bam_map} \\
        -c "\$lib_file_path" \\
        -d ${chromosome} \\
        -e ${params.scan_baq} \\
        -f ${params.scan_min_MQ} \\
        -g ${params.scan_min_BQ} \\
        -h ${params.scan_error_rate} \\
        -i ${samtools_path} \\
        -j ${task.cpus} \\
        -k ${gen_path} \\
        -m ${params.scan_midr} \\
        2>> scan_${chromosome}.log \\
        | bgzip -c -@ ${bgzip_threads} > "\${vcf_gz}.tmp"
    mv -f "\${vcf_gz}.tmp" "\${vcf_gz}"
    echo "Wrote \${vcf_gz} (VCF from stdout, no uncompressed VCF on disk)" >> scan_${chromosome}.log
    echo "Scan analysis for chromosome ${chromosome} completed" >> scan_${chromosome}.log
    """
}

process fastcall3_scan_iac {
    tag "scan_iac_${chromosome}"
    publishDir("${params.output_dir}/${params.job}/gen", mode: 'move', overwrite: false)
    conda "${params.user_dir}/miniconda3/envs/tiger"
    label 'fastcall3_scan_iac'

    input:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map), path(blib_file)
    tuple val(tiger_jar), val(app_name), val(java_version), val(java_lib_dir)
    val samtools_path
    val gen_path

    output:
    path "scan_iac_${chromosome}.log", emit: log

    script:
    """
    set -euo pipefail
    lib_file_path="${blib_file}"
    mkdir -p ${gen_path}/indiCounts

    {
        echo "Starting scan_iac (indiCounts only) for chromosome ${chromosome}..."
        echo "Using ${java_version} for ${app_name} scan_iac process"
        echo "Input library file: \$lib_file_path"
        if [ -f "\$lib_file_path" ]; then
            echo "Library file size: \$(stat -c%s "\$lib_file_path") bytes"
        fi
        echo "Scan output dir (-k): ${gen_path}"
        echo "IAC files: ${gen_path}/indiCounts/<taxon>/${chromosome}_<binStart>_<binEnd>.iac.gz"
        echo "No genotype call and no VCF."
        echo "MiDR (-m): ${params.scan_midr}"
        echo "System resources before TIGER execution:"
        echo "Memory: \$(free -h | grep Mem)"
        echo "CPU cores (allocated): ${task.cpus}"
        echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g"
        echo "TIGER jar: ${tiger_jar}"
        echo "FastCall version: ${app_name}"
    } > scan_iac_${chromosome}.log

    if [ ! -f "\$lib_file_path" ]; then
        echo "Error: Library file \$lib_file_path not found" >> scan_iac_${chromosome}.log
        exit 1
    fi

    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod scan_iac \\
        -a ${reference} \\
        -b ${taxa_bam_map} \\
        -c "\$lib_file_path" \\
        -d ${chromosome} \\
        -e ${params.scan_baq} \\
        -f ${params.scan_min_MQ} \\
        -g ${params.scan_min_BQ} \\
        -h ${params.scan_error_rate} \\
        -i ${samtools_path} \\
        -j ${task.cpus} \\
        -k ${gen_path} \\
        -m ${params.scan_midr} \\
        >> scan_iac_${chromosome}.log 2>&1

    # grep -q closes the pipe early; with pipefail that becomes a false failure after IAC is written.
    iac_hit=\$(find ${gen_path}/indiCounts -type f -name "${chromosome}_*.iac.gz" -print -quit || true)
    if [ -z "\$iac_hit" ]; then
        echo "Error: scan_iac wrote no ${chromosome}_*.iac.gz under ${gen_path}/indiCounts" >> scan_iac_${chromosome}.log
        exit 1
    fi
    echo "Found IAC example: \$iac_hit" >> scan_iac_${chromosome}.log
    echo "Scan_iac for chromosome ${chromosome} completed" >> scan_iac_${chromosome}.log
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
    val gen_path
    
    output:
    path "scan2_${chromosome}.log", emit: log
    
    script:
    def vcf_name = chrVcfFilename(chromosome)
    def bgzip_threads = Math.min(8, task.cpus as int)
    """
    set -euo pipefail
    lib_file_path="${blib_file}"
    mkdir -p ${gen_path}/VCF
    vcf_gz="${gen_path}/VCF/${vcf_name}.gz"

    {
        echo "Starting scan2 (VCF from existing .iac) for chromosome ${chromosome}..."
        echo "Using ${java_version} for ${app_name} scan2 process"
        echo "Input library file: \$lib_file_path"
        echo "Scan output dir (-k, must contain indiCounts/): ${gen_path}"
        echo "VCF stdout -> bgzip -c -@ ${bgzip_threads} -> \$vcf_gz"
        echo "MiDR (-m): ${params.scan_midr}"
        echo "System resources before TIGER execution:"
        echo "Memory: \$(free -h | grep Mem)"
        echo "CPU cores (allocated): ${task.cpus}"
        echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g"
        echo "TIGER jar: ${tiger_jar}"
        echo "FastCall version: ${app_name}"
    } > scan2_${chromosome}.log

    if [ ! -f "\$lib_file_path" ]; then
        echo "Error: Library file \$lib_file_path not found" >> scan2_${chromosome}.log
        exit 1
    fi
    if [ ! -d "${gen_path}/indiCounts" ]; then
        echo "Error: scan2 requires ${gen_path}/indiCounts from a prior scan" >> scan2_${chromosome}.log
        exit 1
    fi

    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod scan2 \\
        -a ${reference} \\
        -b ${taxa_bam_map} \\
        -c "\$lib_file_path" \\
        -d ${chromosome} \\
        -h ${params.scan_error_rate} \\
        -j ${task.cpus} \\
        -k ${gen_path} \\
        -m ${params.scan_midr} \\
        2>> scan2_${chromosome}.log \\
        | bgzip -c -@ ${bgzip_threads} > "\${vcf_gz}.tmp"
    mv -f "\${vcf_gz}.tmp" "\${vcf_gz}"
    echo "Wrote \${vcf_gz} (VCF from stdout, no uncompressed VCF on disk)" >> scan2_${chromosome}.log
    echo "Scan2 for chromosome ${chromosome} completed" >> scan2_${chromosome}.log
    """
}

process fastcall3_scan3 {
    tag "scan3_${chromosome}"
    publishDir("${params.output_dir}/${params.job}/gen3", mode: 'move', overwrite: false)
    conda "${params.user_dir}/miniconda3/envs/tiger"
    label 'fastcall3_scan3'

    input:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map), path(blib_file), path(merge_map)
    tuple val(tiger_jar), val(app_name), val(java_version), val(java_lib_dir)
    val source_dir
    val out_dir

    output:
    path "scan3_${chromosome}.log", emit: log

    script:
    def vcf_name = chrVcfFilename(chromosome)
    def bgzip_threads = Math.min(8, task.cpus as int)
    """
    set -euo pipefail
    lib_file_path="${blib_file}"
    mkdir -p ${out_dir}/VCF
    vcf_gz="${out_dir}/VCF/${vcf_name}.gz"

    {
        echo "Starting scan3 (merge taxa .iac by AD) for chromosome ${chromosome}..."
        echo "Using ${java_version} for ${app_name} scan3 process"
        echo "Input library file: \$lib_file_path"
        echo "Source scan dir (-o): ${source_dir}"
        echo "Merge map (-p): ${merge_map}"
        echo "Output dir (-k): ${out_dir}"
        echo "VCF stdout -> bgzip -c -@ ${bgzip_threads} -> \$vcf_gz"
        echo "MiDR (-m): ${params.scan_midr}"
        echo "System resources before TIGER execution:"
        echo "Memory: \$(free -h | grep Mem)"
        echo "CPU cores (allocated): ${task.cpus}"
        echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g"
        echo "TIGER jar: ${tiger_jar}"
        echo "FastCall version: ${app_name}"
    } > scan3_${chromosome}.log

    if [ ! -f "\$lib_file_path" ]; then
        echo "Error: Library file \$lib_file_path not found" >> scan3_${chromosome}.log
        exit 1
    fi
    if [ ! -d "${source_dir}/indiCounts" ]; then
        echo "Error: scan3 requires ${source_dir}/indiCounts from a prior scan" >> scan3_${chromosome}.log
        exit 1
    fi

    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod scan3 \\
        -a ${reference} \\
        -b ${taxa_bam_map} \\
        -c "\$lib_file_path" \\
        -d ${chromosome} \\
        -h ${params.scan_error_rate} \\
        -j ${task.cpus} \\
        -k ${out_dir} \\
        -m ${params.scan_midr} \\
        -o ${source_dir} \\
        -p ${merge_map} \\
        2>> scan3_${chromosome}.log \\
        | bgzip -c -@ ${bgzip_threads} > "\${vcf_gz}.tmp"
    mv -f "\${vcf_gz}.tmp" "\${vcf_gz}"
    echo "Wrote \${vcf_gz} (VCF from stdout, no uncompressed VCF on disk)" >> scan3_${chromosome}.log
    echo "Scan3 for chromosome ${chromosome} completed" >> scan3_${chromosome}.log
    """
}

process fastcall3_scan4 {
    tag "scan4_${chromosome}"
    publishDir("${params.output_dir}/${params.job}/gen", mode: 'move', overwrite: false)
    conda "${params.user_dir}/miniconda3/envs/tiger"
    label 'fastcall3_scan4'

    input:
    tuple val(chromosome), path(reference), path(fai), path(gzi), path(taxa_bam_map), path(blib_file)
    tuple val(tiger_jar), val(app_name), val(java_version), val(java_lib_dir)
    val gen_path

    output:
    path "scan4_${chromosome}.log", emit: log

    script:
    def vcf_name = chrVcfFilename(chromosome)
    def bgzip_threads = Math.min(8, task.cpus as int)
    """
    set -euo pipefail
    lib_file_path="${blib_file}"
    mkdir -p ${gen_path}/VCF
    vcf_gz="${gen_path}/VCF/${vcf_name}.gz"

    {
        echo "Starting scan4 (windowed VCF from existing .iac) for chromosome ${chromosome}..."
        echo "Using ${java_version} for ${app_name} scan4 process"
        echo "Input library file: \$lib_file_path"
        echo "Scan output dir (-k, must contain indiCounts/): ${gen_path}"
        echo "VCF stdout -> bgzip -c -@ ${bgzip_threads} -> \$vcf_gz"
        echo "MiDR (-m): ${params.scan_midr}"
        echo "System resources before TIGER execution:"
        echo "Memory: \$(free -h | grep Mem)"
        echo "CPU cores (allocated): ${task.cpus}"
        echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g"
        echo "TIGER jar: ${tiger_jar}"
        echo "FastCall version: ${app_name}"
    } > scan4_${chromosome}.log

    if [ ! -f "\$lib_file_path" ]; then
        echo "Error: Library file \$lib_file_path not found" >> scan4_${chromosome}.log
        exit 1
    fi
    if [ ! -d "${gen_path}/indiCounts" ]; then
        echo "Error: scan4 requires ${gen_path}/indiCounts from a prior scan" >> scan4_${chromosome}.log
        exit 1
    fi

    java -Xmx${task.memory.toGiga()}g -jar ${tiger_jar} \\
        -app ${app_name} \\
        -mod scan4 \\
        -a ${reference} \\
        -b ${taxa_bam_map} \\
        -c "\$lib_file_path" \\
        -d ${chromosome} \\
        -h ${params.scan_error_rate} \\
        -j ${task.cpus} \\
        -k ${gen_path} \\
        -m ${params.scan_midr} \\
        2>> scan4_${chromosome}.log \\
        | bgzip -c -@ ${bgzip_threads} > "\${vcf_gz}.tmp"
    mv -f "\${vcf_gz}.tmp" "\${vcf_gz}"
    echo "Wrote \${vcf_gz} (VCF from stdout, no uncompressed VCF on disk)" >> scan4_${chromosome}.log
    echo "Scan4 for chromosome ${chromosome} completed" >> scan4_${chromosome}.log
    """
}
