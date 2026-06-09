nextflow.enable.dsl=2

include { getRefV1ChrName } from '../../infra/infra_ref_v1.nf'
include { getRefV1ChrOffset } from '../../infra/infra_ref_v1.nf'

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

process arrange_sample_info {
    tag "arrange sample info: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/sample_info", mode: 'copy', pattern: "*.sample_info.tsv"
    publishDir "${params.output_dir}/${params.job}/process/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), path(vcf)

    output:
    tuple val(id), val(chr), path("${id}.sample_info.tsv"), emit: sample_info
    path "${id}.sample_info.log", emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${id}.sample_info.log", "w")
    sys.stderr = sys.stdout

    from python_script.genomics.germplasm.sample.ana import ana_sample_info

    print(f"Processing sample information for ${id}...")
    ana_sample_info("${vcf}", "${id}.sample_info.tsv")
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
    if [ -f "${tmp_dir}/${id}.plink.bed" ] && [ -f "${tmp_dir}/${id}.plink.bim" ] && [ -f "${tmp_dir}/${id}.plink.fam" ]; then
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

