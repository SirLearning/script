nextflow.enable.dsl=2

include { plink_assess as PLINK_ASSESS } from './assess.nf'
include { getTigerJarConfig } from './utils.nf'

/*
    Module: genotype/processor.nf
    Description: Processes genotype VCF files, generating information, filtering, formatting, and QC statistics.
*/

workflow plink_processor {
    take:
    // Expect a channel: [ val(id), val(chr), path(vcf) ]
    vcf_in

    main:
    // 1. bgzip + tabix (check and ensure)
    gz_vcf = format_vcf_bgzip_idx(vcf_in)

    // 2. tranfer to PLINK format
    plink_out = format_vcf_plink(gz_vcf.vcf)

    // 3. make sample and variant basic information
    basic_info_out = mk_plink_basic_info(plink_out.pfile)

    // 4. assess with plink
    assess_out = PLINK_ASSESS(
        basic_info_out.smiss, 
        basic_info_out.vmiss,
        basic_info_out.scount,
        basic_info_out.gcount,
        basic_info_out.afreq,
        basic_info_out.hardy)

    emit:
    vcf = gz_vcf.vcf
    plink_bfile = plink_out.bfile
    plink_pfile  = plink_out.pfile
}


process format_vcf_bgzip_idx {
    tag "${chr}" ? "bgzip/tabix ${chr}" : 'bgzip/tabix'
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'
    conda 'stats'

    input:
    tuple val(id), val(chr), path(vcf)

    output:
    tuple val(id), val(chr), path("${id}.vcf.gz"), path("${id}.vcf.gz.tbi"), emit: vcf

    script:
    """
    set -euo pipefail

    in=${vcf}
    out=${id}.vcf.gz

    if [[ "\${in}" == *.vcf.gz ]]; then
        in_base=\$(basename "\${in}")
        if [[ "\${in_base}" != "\${out}" ]]; then
            ln -sf "\${in}" "\${out}" || cp -f "\${in}" "\${out}"
        } else
        bgzip -c -@ ${task.cpus} "\${in}" > "\${out}"
    fi

    if ! tabix -f -p vcf "\${out}"; then
        echo "tabix failed on \${out}; ensuring BGZF compression and retrying" >&2
        if [[ "\${in}" == *.vcf.gz ]]; then
            gunzip -c "\${in}" | bgzip -c -@ ${task.cpus} > "\${out}.tmp"
        else
            bgzip -c -@ ${task.cpus} "\${in}" > "\${out}.tmp"
        fi
        mv -f "\${out}.tmp" "\${out}"
        tabix -@ ${task.cpus} -f -p vcf "\${out}"
    fi
    """
}


process format_vcf_plink {
    tag "${chr}" ? "plink ${chr}" : 'plink'
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'
    conda 'stats'

    input:
    tuple val(id), val(chr), path(vcf), path(tbi)

    output:
    tuple val(id), val(chr), val("${id}.plink"), path("${id}.bed"), path("${id}.bim"), path("${id}.fam"), emit: bfile
    tuple val(id), val(chr), val("${id}.plink2"), path("${id}.pgen"), path("${id}.psam"), path("${id}.pvar"), emit: pfile

    script:
    """
    set -euo pipefail
    if [[ ! -s ${vcf} ]]; then
        echo "ERROR: VCF file missing or empty before PLINK: ${vcf}" >&2
        ls -l >&2 || true
        exit 201
    fi

    plink2_out=${id}.plink2
    plink_out=${id}.plink

    plink2 --vcf ${vcf} \\
        --allow-extra-chr \\
        --make-pgen \\
        --max-alleles 2 \\
        --thread ${task.cpus} \\
        --out \${plink2_out} > ${chr}.plk2.log 2>&1
    
    plink2 --pfile \${plink2_out} \\
        --make-bed \\
        --thread ${task.cpus} \\
        --out \${plink_out} > ${chr}.plk.log 2>&1

    """
}

process mk_plink_basic_info {
    tag "make plink basic info: ${chr}"
    publishDir "${params.output_dir}/${params.job}/process/sample", mode: 'copy', pattern: "*.{smiss,scount}"
    publishDir "${params.output_dir}/${params.job}/process/variant", mode: 'copy', pattern: "*.{vmiss,gcount,afreq,hardy}"
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy', pattern: "*.log"
    conda 'stats'

    input:
    tuple val(id), val(chr), path(pgen), path(psam), path(pvar)

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
    plink2 --pfile ${pgen.baseName} \\
        --allow-extra-chr \\
        --missing \\
        --sample-counts \\
        --geno-counts \\
        --freq \\   MAF
        --hardy \\  heterozygosity
        --threads ${task.cpus} \\
        --out ${id}.info > ${chr}.info.log 2>&1
    """
}

process prepare_popdepth {
    tag "prepare popdepth: ${chr}"
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'
    conda 'run'

    input:
    tuple val(id), val(chr), path(vcf)

    output:
    tuple val(id), val(chr), path("${id}.popdepth.txt"), emit: popdepth

    script:
    """
    set -euo pipefail
    bcftools query -f '%CHROM\\t%POS\\t%INFO/DP\\n' ${vcf} > ${id}.popdepth.txt
    """
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
