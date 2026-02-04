nextflow.enable.dsl=2

include { getTigerJarConfig } from './utils.nf'
include { getPopDepTaxaBamFile_v1 } from './utils.nf'
include { getRefV1ChrLength } from './utils.nf'
include { getRefV1ChrName } from './utils.nf'
include { getRefV1ChrOffset } from './utils.nf'

/*
    Module: genotype/processor.nf
    Description: Processes genotype VCF files, generating information, filtering, formatting, and QC statistics.
*/

workflow test_plink {
    take:
    // Expect a channel: [ val(id), val(chr), path(vcf) ]
    vcf_in

    main:
    preprocess_out = plink_preprocess(vcf_in)

    // subsample pfile for testing
    def subsampling_out = subsampling_pfile_for_test(preprocess_out.pfile, params.thin_rate)
    // Map pfiles to include subgenome info
    def grouped_pfile = subsampling_out.pfile.map { id, chr, prefix, pgen, psam, pvar ->
        def subgenome = "Others"
        def c = chr.toString()
        if (["1","2","7","8","13","14","19","20","25","26","31","32","37","38"].contains(c)) subgenome = "A"
        else if (["3","4","9","10","15","16","21","22","27","28","33","34","39","40"].contains(c)) subgenome = "B"
        else if (["5","6","11","12","17","18","23","24","29","30","35","36","41","42"].contains(c)) subgenome = "D"
        return [ subgenome, id, chr, prefix, pgen, psam, pvar ]
    }.groupTuple(by: 0)

    // Generate merge lists
    merge_out = merge_subgenome_test_pfile(grouped_pfile)

    emit:
    vcf = preprocess_out.gz_vcf
    merged_bfile = merge_out.bfile
    merged_pfile = merge_out.pfile
}

workflow plink_processor {
    take:
    // Expect a channel: [ val(id), val(chr), path(vcf) ]
    vcf_in

    main:
    // 1. preprocess VCF to PLINK format
    preprocess_out = plink_preprocess(vcf_in)

    // 2. make sample and variant basic information
    basic_info_out = mk_plink_basic_info(preprocess_out.pfile)

    // 3. calculate population depth using TIGER
    def pd_config = getTigerJarConfig("TIGER_PD_20260130.jar", params.home_dir)
    ch_tiger_config = channel.of([ pd_config.path, pd_config.app_name, pd_config.java_version ])
    popdep_out = calc_population_depth(preprocess_out.gz_vcf, ch_tiger_config)

    emit:
    vcf = preprocess_out.gz_vcf
    plink_bfile = preprocess_out.bfile
    plink_pfile  = preprocess_out.pfile
    smiss = basic_info_out.smiss
    vmiss = basic_info_out.vmiss
    scount = basic_info_out.scount
    gcount = basic_info_out.gcount
    afreq = basic_info_out.afreq
    hardy = basic_info_out.hardy
    popdep = popdep_out.popdep
}

workflow plink_preprocess {
    take:
    // Expect a channel: [ val(id), val(chr), path(vcf) ]
    vcf_in

    main:
    gz_vcf = format_vcf_bgzip(vcf_in)
    plink_out = format_vcf_plink(gz_vcf.vcf)

    emit:
    gz_vcf = gz_vcf.vcf
    bfile = plink_out.bfile
    pfile = plink_out.pfile
}

process format_vcf_bgzip {
    tag "bgzip ${chr}"
    publishDir "${params.output_dir}/${params.job}/process", mode: 'symlink', pattern: "*.vcf.gz"
    publishDir "${params.output_dir}/${params.job}/process/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), path(vcf)

    output:
    tuple val(id), val(chr), path("${id}.vcf.gz"), emit: vcf

    script:
    """
    set -euo pipefail
    exec > bgzip.${chr}.log 2>&1

    in=${vcf}
    out=${id}.vcf.gz

    if [[ "\${in}" == *.vcf.gz ]]; then
        echo "Input VCF is already gzipped, linking or copying to output..."
        in_base=\$(basename "\${in}")
        if [[ "\${in_base}" != "\${out}" ]]; then
            ln -sf "\${in}" "\${out}" || cp -f "\${in}" "\${out}"
        fi
    else
        echo "Compressing VCF with bgzip..."
        bgzip -c -@ ${task.cpus} "\${in}" > "\${out}"
    fi

    echo "BGZF compression completed for ${id}."
    """
}

workflow arrange_merge_vcf {
    take:
    // Expect a channel: [ val(id), val(chr), path(vcf) ]
    vcf_in

    main:
    arrange_out = arrange_vcf_wheat_chr_by_awk(vcf_in)
    arrange_ch = arrange_out.vcf.groupTuple(by: 2)
    merge_out = merge_arranged_vcf(arrange_ch)

    emit:
    vcf = merge_out.vcf
}

process format_vcf_bgzip_idx {
    tag "bgzip/tabix ${chr}"
    publishDir "${params.output_dir}/${params.job}/process", mode: 'symlink', pattern: "*.vcf.gz"
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy', pattern: "*.vcf.gz.tbi"
    publishDir "${params.output_dir}/${params.job}/process/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), path(vcf)

    output:
    tuple val(id), val(chr), path("${id}.vcf.gz"), path("${id}.vcf.gz.tbi"), emit: vcf

    script:
    def tmp_dir = "${params.output_dir}/${params.job}/tmp"
    """
    set -euo pipefail
    exec > bgzip_tabix.${chr}.log 2>&1

    in=${vcf}
    out=${id}.vcf.gz

    if [[ "\${in}" == *.vcf.gz ]]; then
        echo "Input VCF is already gzipped, linking or copying to output..."
        in_base=\$(basename "\${in}")
        if [[ "\${in_base}" != "\${out}" ]]; then
            ln -sf "\${in}" "\${out}" || cp -f "\${in}" "\${out}"
        fi
    else
        echo "Compressing VCF with bgzip..."
        bgzip -c -@ ${task.cpus} "\${in}" > "\${out}"
    fi

    if [ -f "${tmp_dir}/\${out}.tbi" ]; then
        echo "Index file already exists for \${out}, skipping tabix indexing."
        ln -sf "${tmp_dir}/\${out}.tbi" ./
    else
        echo "Indexing VCF with tabix..."
        if ! tabix -@ ${task.cpus} -f -p vcf "\${out}"; then
            echo "tabix failed on \${out}; ensuring BGZF compression and retrying" >&2
            if [[ "\${in}" == *.vcf.gz ]]; then
                gunzip -c "\${in}" | bgzip -c -@ ${task.cpus} > "\${out}.tmp"
            else
                bgzip -c -@ ${task.cpus} "\${in}" > "\${out}.tmp"
            fi
            mv -f "\${out}.tmp" "\${out}"
            tabix -@ ${task.cpus} -f -p vcf "\${out}"
        fi
    fi

    echo "BGZF compression and indexing completed for ${id}."
    """
}

process arrange_vcf_wheat_chr_by_awk {
    tag "arrange chrom pos awk ${chr}"
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy', pattern: "*.vcf.gz"
    publishDir "${params.output_dir}/${params.job}/process/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), path(vcf)

    output:
    tuple val(id), val(chr), val(getRefV1ChrName(chr)), path("${chr}.arr_chr.vcf.gz"), emit: vcf
    path "arrange_chr_pos_${chr}.log", emit: log

    script:
    def chr_name = getRefV1ChrName(chr)
    def chr_offset = getRefV1ChrOffset(chr)
    """
    set -euo pipefail
    exec > arrange_chr_pos_${chr}.log 2>&1

    echo "Arranging chromosome names and positions for ${chr}..."
    bgzip -dc ${vcf} | awk -v new_chr="${chr_name}" -v offset=${chr_offset} '
    BEGIN { OFS="\\t" }
    /^##contig/ {
        print \$0; next 
    }
    /^#/ { print \$0; next }
    {
        \$1 = new_chr;
        \$2 = \$2 + offset;
        print \$0;
    }
    ' | bgzip -c -@ ${task.cpus} > ${id}.arr_chr.vcf.gz

    echo "Chromosome arrangement completed for ${chr}."
    """
}

process merge_arranged_vcf {
    tag "merge arranged vcf: ${chr_name}"
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy', pattern: "*.vcf.gz"
    publishDir "${params.output_dir}/${params.job}/process/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(chr_name), val(ids), val(chrs), path(vcfs)

    output:
    tuple val(chr_name), path("${chr_name}.vcf.gz"), emit: vcf
    path "merge_arranged_vcf_${chr_name}.log", emit: log

    script:
    """
    set -euo pipefail
    exec > merge_arranged_vcf_${chr_name}.log 2>&1

    echo "Merging arranged VCFs for ${chr_name}..."
    ls ${vcfs.join(' ')} > ${chr_name}_vcf_list.txt

    bcftools merge -m all -O z -o ${chr_name}.vcf.gz -l ${chr_name}_vcf_list.txt --threads ${task.cpus}

    echo "VCF merging completed for ${chr_name}."
    """
}

process format_vcf_plink {
    tag "format plink ${chr}"
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy', pattern: "*.{bed,bim,fam,pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), path(vcf)

    output:
    tuple val(id), val(chr), val("${id}.plink"), path("${id}.plink.bed"), path("${id}.plink.bim"), path("${id}.plink.fam"), emit: bfile
    tuple val(id), val(chr), val("${id}.plink2"), path("${id}.plink2.pgen"), path("${id}.plink2.psam"), path("${id}.plink2.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    def tmp_dir = "${params.output_dir}/${params.job}/tmp"
    """
    set -euo pipefail
    exec > format_plink.${chr}.log 2>&1
    
    plink2_out="${id}.plink2"
    plink_out="${id}.plink"
    
    EXIST_PLINK1=false
    EXIST_PLINK2=false
    if [ -f "${tmp_dir}/${id}.bed" ] && [ -f "${tmp_dir}/${id}.bim" ] && [ -f "${tmp_dir}/${id}.fam" ]; then
        EXIST_PLINK1=true
    fi
    if [ -f "${tmp_dir}/${id}.plink2.pgen" ] && [ -f "${tmp_dir}/${id}.plink2.psam" ] && [ -f "${tmp_dir}/${id}.plink2.pvar" ]; then
        EXIST_PLINK2=true
    fi

    if [ "\$EXIST_PLINK2" = true ]; then
        echo "PLINK2 files for ${id} already exist in ${tmp_dir}, linking them..."
        ln -sf "${tmp_dir}/${id}.plink2.pgen" ./
        ln -sf "${tmp_dir}/${id}.plink2.psam" ./
        ln -sf "${tmp_dir}/${id}.plink2.pvar" ./
    else
        if [[ ! -s ${vcf} ]]; then
            echo "ERROR: VCF file missing or empty before PLINK: ${vcf}" >&2
            exit 201
        fi
        plink2 --vcf ${vcf} \\
            --chr-set 44 \\
            --allow-extra-chr \\
            --make-pgen \\
            --max-alleles 2 \\
            --threads ${task.cpus} \\
            --out \${plink2_out}
    fi
    
    if [ "\$EXIST_PLINK1" = true ]; then
        echo "PLINK1 files for ${id} already exist in ${tmp_dir}, linking them..."
        ln -sf "${tmp_dir}/${id}.plink.bed" ./
        ln -sf "${tmp_dir}/${id}.plink.bim" ./
        ln -sf "${tmp_dir}/${id}.plink.fam" ./
    else
        plink2 --pfile \${plink2_out} \\
            --make-bed \\
            --threads ${task.cpus} \\
            --out \${plink_out}
    fi
    """
}

process subsampling_pfile_for_test {
    tag "subsample pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)
    val thin_rate

    output:
    tuple val(id), val(chr), val("${id}.thin"), path("${id}.thin.pgen"), path("${id}.thin.psam"), path("${id}.thin.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > subsample_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --thin ${thin_rate} \\
        --seed \$(date +%s) \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.thin
    """
}

process merge_subgenome_test_pfile {
    tag "merge plink2 test subgenome pfile: ${subgenome}"
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy', pattern: "*.{bed,bim,fam,pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    // Receives lists of files/values grouped by subgenome
    tuple val(subgenome), val(ids), val(chrs), val(prefixes), path(pgens), path(psams), path(pvars)

    output:
    tuple val(subgenome), val(chrs), val("${subgenome}_test.plink"), path("${subgenome}_test.plink.bed"), path("${subgenome}_test.plink.bim"), path("${subgenome}_test.plink.fam"), emit: bfile
    tuple val(subgenome), val(chrs), val("${subgenome}_test.plink2"), path("${subgenome}_test.plink2.pgen"), path("${subgenome}_test.plink2.psam"), path("${subgenome}_test.plink2.pvar"), emit: pfile
    tuple val(subgenome), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > merge_subgenome_test.${subgenome}.log 2>&1
    
    echo "${prefixes[1..-1].join('\n')}" > ${subgenome}.merge_list.txt
    
    plink2 --pfile ${prefixes[0]} \\
        --pmerge-list ${subgenome}.merge_list.txt \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${subgenome}_test.plink2

    plink2 --pfile ${subgenome}_test.plink2 \\
        --make-bed \\
        --threads ${task.cpus} \\
        --out ${subgenome}_test.plink
    """
}

process mk_plink_basic_info {
    tag "make plink basic info: ${chr}"
    publishDir "${params.output_dir}/${params.job}/process/sample", mode: 'copy', pattern: "*.{smiss,scount}"
    publishDir "${params.output_dir}/${params.job}/process/variant", mode: 'copy', pattern: "*.{vmiss,gcount,afreq,hardy}"
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'
    
    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), path("${id}.info.smiss"), emit: smiss
    tuple val(id), val(chr), path("${id}.info.vmiss"), emit: vmiss
    tuple val(id), val(chr), path("${id}.info.scount"), emit: scount
    tuple val(id), val(chr), path("${id}.info.gcount"), emit: gcount
    tuple val(id), val(chr), path("${id}.info.afreq"), emit: afreq
    tuple val(id), val(chr), path("${id}.info.hardy"), emit: hardy
    tuple val(id), val(chr), path("${id}.info.log"), emit: log

    script:
    """
    set -euo pipefail
    plink2 --pfile ${prefix} \\
        --allow-extra-chr \\
        --missing \\
        --sample-counts \\
        --geno-counts \\
        --freq \\
        --hardy \\
        --threads ${task.cpus} \\
        --out ${id}.info > ${chr}.info.log 2>&1
    """
}

process mk_vcftools_basic_info {
    tag "vcftools stats: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/sample", mode: 'copy', pattern: "*.{imiss, idepth}"
    publishDir "${params.output_dir}/${params.job}/process/variant", mode: 'copy', pattern: "*.{frq,hwe,lmiss,ldepth.mean,lqual}"
    publishDir "${params.output_dir}/${params.job}/assess/vcftools/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), path(vcf)

    output:
    tuple val(id), val(chr), path("${id}.frq"), emit: frq
    tuple val(id), val(chr), path("${id}.hwe"), emit: hwe
    tuple val(id), val(chr), path("${id}.lmiss"), emit: lmiss
    tuple val(id), val(chr), path("${id}.imiss"), emit: imiss
    tuple val(id), val(chr), path("${id}.ldepth.mean"), emit: ldepth_mean
    tuple val(id), val(chr), path("${id}.idepth"), emit: idepth
    tuple val(id), val(chr), path("${id}.lqual"), emit: lqual
    path "${chr}.vcftools.log", emit: log

    script:
    """
    set -euo pipefail
    exec > ${chr}.vcftools.log 2>&1

    vcftools --gzvcf ${vcf} \\
        --freq \\
        --hardy \\
        --missing-site \\
        --missing-indv \\
        --site-mean-depth \\
        --depth \\
        --site-quality \\
        --out ${id}.vcftools
    """
}

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

    echo "VCF filtering completed for ${id}."
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
    
    echo "VCF filtering completed for ${id}."
    echo "Filtered VCF saved to ${id}.bcftools.flt.v0.vcf"
    """
}
