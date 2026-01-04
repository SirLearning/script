nextflow.enable.dsl=2

// --- Filter Parameters ---
// standard filter parameters (std)
def maf = 0.05
def mac = 2
def min_alleles = 2
def max_alleles = 2
def max_missing = 0.05
// additional parameters can be added here
def qual = 30
def hwe_pval = 1e-6 // not sure if used, maybe not now

workflow processor {
    take:
    // Expect a combined channel: [ val(), path(vcf), val(job_config) ]
    vcf_in

    main:
    // --- format VCF ---
    // Receive the full tuple and capture the single output channel directly
    gz_vcf = format_vcf_bgzip_idx(filtered_vcf.vcf)
    // Pass the resulting (, vcf.gz) tuples directly to PLINK
    plink_out = format_vcf_plink(gz_vcf)
    
    // --- filtering steps ---
    log.info """\
    Starting VCF filtering using standard parameters:
        MAF >= ${maf}
        MAC >= ${mac}
        Min Alleles = ${min_alleles}
        Max Alleles = ${max_alleles}
        Max Missing = ${max_missing}
    """
    // filter
    filtered_vcf = filter_vcftools_std(vcf_in)

    // --- QC Stats ---
    stats_out = vcf_stats(filtered_vcf.vcf)
    plot_vcf_qc(stats_out.stats)

    // --- Convert to MT ---
    mt_out = vcf_to_mt_hail(filtered_vcf.vcf)

    emit:
    vcf = gz_vcf.vcf
    plink_bfile = plink_out.plink_bfile
    plink_ped  = plink_out.plink_ped
}

process vcf_stats {
    tag { id ? "vcf stats ${id}" : 'vcf stats' }
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
    tag { id ? "plot qc ${id}" : 'plot qc' }
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
    tag { id ? "vcftools filter ${id}" : 'vcftools filter' }
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(), path(vcf), val(job_config)

    output:
    tuple val(), val("${id}.std.filtered"), path("${id}.std.filtered.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    in="${vcf}"
    out="${id}.std.filtered"

    vcftools --gzvcf "\${in}" \\
        --max-missing ${max_missing} \\
        --remove-indels \\
        --maf ${maf} \\
        --min-alleles ${min_alleles} \\
        --max-alleles ${max_alleles} \\
        --mac ${mac} \\
        --recode \\
        --recode-INFO-all \\
        --out "\${out}"

    mv "\${out}.recode.vcf" "\${out}.vcf"
    """
}

process filter_bcftools_std {
    tag { id ? "bcftools filter ${id}" : 'bcftools filter' }
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(), path(vcf), val(job_config)

    output:
    tuple val(), val("${id}.std.filtered"), path("${id}.std.filtered.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    bcftools view \\
        -m2 -M2 -v snps \\
        -i 'MAF>=${maf} && AC>=${mac} && && F_MISSING<=${max_missing}' \\
        -o ${id}.std.filtered.vcf \\
        ${vcf}
    """
}

process format_vcf_bgzip_idx {
    tag { id ? "bgzip/tabix ${id}" : 'bgzip/tabix' }
    // Use copy to keep files available for downstream processes (move would remove them before staging)
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(), val(prefix), path(vcf)

    output:
    tuple val(), val(prefix), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf

    script:
    """
    set -euo pipefail

    in=${vcf}
    out=${prefix}.vcf.bgz
    echo "Processing VCF for work: ${.id}" >&2

    if [[ "\${in}" == *.vcf.gz ]]; then
        echo "Input VCF is already gzipped: \${in}" >&2
        in_base=\$(basename "\${in}")
        # Avoid self-linking when names are identical
        if [[ "\${in_base}" != "\${out}" ]]; then
            ln -sf "\${in}" "\${out}" || cp -f "\${in}" "\${out}"
        fi
    else
        echo "Compressing VCF with bgzip: \${in} -> \${out}" >&2
        bgzip -c "\${in}" > "\${out}"
    fi

    # Try to create index; if it fails (e.g., not BGZF), re-compress and retry
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
    tag { id ? "plink ${id}" : 'plink' }
    // Keep generated bed/bim/fam for further steps
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(), val(prefix), path(vcf), path(tbi)

    output:
    tuple val(), path("${prefix}.bed"), path("${prefix}.bim"), path("${prefix}.fam"), emit: plink_bfile
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
    tag { id ? "arrange chrom pos awk ${id}" : 'arrange chrom pos awk' }
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(), path(input_file)
    path map_tsv
    tuple val(chrom_col), val(pos_col)

    output:
    tuple val(), path("${input_file.baseName}.arr_chr.txt"), emit: vcf

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
    tag { id ? "arrange chrom pos ${id}" : 'arrange chrom pos' }
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

