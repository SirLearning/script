nextflow.enable.dsl=2

/**
 * Genetics Genotype Processor Workflow
 * 
 * 该 workflow 主要负责 VCF 文件的标准化处理、过滤、质控 (QC) 以及格式转换。
 * 主要步骤包括：
 * 1. 使用 vcftools 进行标准参数过滤 (MAF, MAC, Alleles, Missingness)。
 * 2. 对过滤后的 VCF 进行 bgzip 压缩并建立 tabix 索引，确保后续工具的兼容性和高效访问。
 * 3. 生成 QC 统计数据并利用 R 脚本绘图。
 * 4. 将 VCF 转换为 PLINK (bed/bim/fam, ped) 格式。
 * 
 * 注意：所有的 gzip 压缩均采用 bgzip 格式，以支持 tabix 随机访问。
 */

workflow processor {
    take:
    // Expect a channel: [ val(id), path(vcf) ]
    vcf_in

    main:
    // --- filtering steps ---
    log.info """\
    Starting VCF filtering using standard parameters:
        MAF >= ${params.maf}
        MAC >= ${params.mac}
        Min Alleles = ${params.min_alleles}
        Max Alleles = ${params.max_alleles}
        Max Missing = ${params.max_missing}
    """
    
    // 1. 首先进行过滤
    // filter_vcftools_std expects [ id, vcf ]
    filtered_vcf = filter_vcftools_std(vcf_in)

    // 2. 对过滤后的 VCF 进行 bgzip 压缩并建立 tabix 索引
    // filtered_vcf.vcf emits [ id, prefix, vcf_path ]
    gz_vcf = format_vcf_bgzip_idx(filtered_vcf.vcf)
    
    // 3. 将压缩索引后的 VCF 传递给 PLINK 进行格式转换
    // gz_vcf.vcf emits [ id, prefix, vcf_gz, tbi ]
    plink_out = format_vcf_plink(gz_vcf.vcf)
    
    // 4. QC Stats
    // 使用过滤后的 VCF 进行统计
    stats_out = vcf_stats(filtered_vcf.vcf.map{ id, prefix, vcf -> tuple(id, vcf) })
    plot_vcf_qc(stats_out.stats)

    emit:
    vcf = gz_vcf.vcf
    plink_bfile = plink_out.plink_bfile
    plink_ped  = plink_out.plink_ped
}

workflow plink_processor {
    take:
    // Expect a channel: [ val(id), path(vcf) ]
    vcf_in

    main:
    // 1. 对输入的 VCF 进行 bgzip 压缩并建立 tabix 索引
    gz_vcf = format_vcf_bgzip_idx(vcf_in)

    // 2. 将压缩索引后的 VCF 传递给 PLINK 进行格式转换
    plink_out = format_vcf_plink(gz_vcf.vcf)

    emit:
    vcf = gz_vcf.vcf
    plink_bfile = plink_out.plink_bfile
    plink_ped  = plink_out.plink_ped
}

process vcf_stats {
    tag "${id}" ? "vcf stats ${id}" : 'vcf stats' 
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}.imiss"), path("${id}.lmiss"), path("${id}.het"), path("${id}.frq"), path("${id}.idepth"), path("${id}.ldepth.mean"), emit: stats

    script:
    """
    set -euo pipefail
    vcftools --vcf ${vcf} --missing-indv --out ${id}
    vcftools --vcf ${vcf} --missing-site --out ${id}
    vcftools --vcf ${vcf} --het --out ${id}
    vcftools --vcf ${vcf} --freq2 --out ${id}
    vcftools --vcf ${vcf} --depth --out ${id}
    vcftools --vcf ${vcf} --site-mean-depth --out ${id}
    """
}

process plot_vcf_qc {
    tag "${id}" ? "plot qc ${id}" : 'plot qc'
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy'

    input:
    tuple val(id), path(imiss), path(lmiss), path(het), path(frq), path(idepth), path(ldepth)

    output:
    path("${id}.qc_plots.pdf")

    script:
    """
    set -euo pipefail
    Rscript ${params.src_dir}/r/genetics/vcf_qc_plot.r \\
        --imiss ${imiss} \\
        --lmiss ${lmiss} \\
        --het ${het} \\
        --frq ${frq} \\
        --depth ${idepth} \\
        --site_depth ${ldepth} \\
        --output ${id}.qc_plots.pdf
    """
}

process filter_vcftools_std {
    tag "${id}" ? "vcftools filter ${id}" : 'vcftools filter'
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), val("${id}.std.filtered"), path("${id}.std.filtered.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    in="${vcf}"
    out="${id}.std.filtered"

    # 如果输入是压缩的，使用 --gzvcf
    if [[ "\${in}" == *.gz ]] || [[ "\${in}" == *.bgz ]]; then
        VCF_OPT="--gzvcf"
    else
        VCF_OPT="--vcf"
    fi

    vcftools \${VCF_OPT} "\${in}" \\
        --max-missing ${params.max_missing} \\
        --remove-indels \\
        --maf ${params.maf} \\
        --min-alleles ${params.min_alleles} \\
        --max-alleles ${params.max_alleles} \\
        --mac ${params.mac} \\
        --recode \\
        --recode-INFO-all \\
        --out "\${out}"

    mv "\${out}.recode.vcf" "\${out}.vcf"
    """
}

process filter_bcftools_std {
    tag "${id}" ? "bcftools filter ${id}" : 'bcftools filter'
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), val("${id}.std.filtered"), path("${id}.std.filtered.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    bcftools view \\
        -m2 -M2 -v snps \\
        -i 'MAF>=${params.maf} && AC>=${params.mac} && F_MISSING<=${params.max_missing}' \\
        -o ${id}.std.filtered.vcf \\
        ${vcf}
    """
}

process format_vcf_bgzip_idx {
    tag "${prefix}" ? "bgzip/tabix ${prefix}" : 'bgzip/tabix'
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(id), val(prefix), path(vcf)

    output:
    tuple val(id), val(prefix), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf

    script:
    """
    set -euo pipefail

    in=${vcf}
    out=${prefix}.vcf.gz

    # 检查输入是否已经是 gzipped。
    if [[ "\${in}" == *.vcf.gz ]]; then
        in_base=\$(basename "\${in}")
        if [[ "\${in_base}" != "\${out}" ]]; then
            ln -sf "\${in}" "\${out}" || cp -f "\${in}" "\${out}"
        fi
    else
        bgzip -c "\${in}" > "\${out}"
    fi

    # 尝试创建 tabix 索引。
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

process format_vcf_plink {
    tag "${prefix}" ? "plink ${prefix}" : 'plink'
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(id), val(prefix), path(vcf), path(tbi)

    output:
    tuple val(id), path("${prefix}.bed"), path("${prefix}.bim"), path("${prefix}.fam"), emit: plink_bfile
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

// only for text files
process arrange_wheat_chr_by_awk {
    tag "${id}" ? "arrange chrom pos awk ${id}" : 'arrange chrom pos awk'
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(id), path(input_file)
    path map_tsv
    tuple val(chrom_col), val(pos_col)

    output:
    tuple val(id), path("${input_file.baseName}.arr_chr.txt"), emit: vcf

    script:
    """
    # Process file using awk
    awk -v c_arg="${chrom_col}" -v p_arg="${pos_col}" '
    BEGIN { FS=OFS="\t" }
    
    # Load mapping (File 1: map_tsv)
    NR==FNR {
        map_chrom[\$1] = \$2
        map_offset[\$1] = \$3
        next
    }

    # Process Header (File 2: input_file, Line 1)
    FNR==1 {
        print \$0
        
        c_idx = -1
        p_idx = -1
        
        # Check if arguments are numeric (1-based index)
        if (c_arg ~ /^[0-9]+\$/) c_idx = c_arg
        if (p_arg ~ /^[0-9]+\$/) p_idx = p_arg
        
        # Check header names (overrides numeric if match found)
        for (i=1; i<=NF; i++) {
            if (\$i == c_arg) c_idx = i
            if (\$i == p_arg) p_idx = i
        }
        
        if (c_idx == -1 || p_idx == -1) {
            print "Error: Columns " c_arg " or " p_arg " not found." > "/dev/stderr"
            exit 1
        }
        next
    }

    # Process Data
    {
        # Skip empty lines
        if (NF < 1) next

        # Check if chromosome needs mapping
        if (\$c_idx in map_chrom) {
            old_chr = \$c_idx
            \$c_idx = map_chrom[old_chr]
            \$p_idx = \$p_idx + map_offset[old_chr]
        }
        print \$0
    }
    ' "${map_tsv}" "${input_file}" > "${input_file.baseName}.arr_chr.txt"
    """
}

// only for text files
process arrange_wheat_chr_by_python {
    tag "${id}" ? "arrange chrom pos ${id}" : 'arrange chrom pos'
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(id), path(input_file)
    path map_json
    tuple val(chrom_col), val(pos_col)

    output:
    tuple val(id), path("${input_file.baseName}.arr_chr.txt"), emit: vcf

    script:
    """
    python ${params.src_dir}/python/infra/utils/wheat_ref_v1.py \\
        -m "${map_json}" \\
        -i "${input_file}" \\
        -o "${input_file.baseName}.arr_chr.txt" \\
        -c "${chrom_col}" \\
        -p "${pos_col}"
    """
}
