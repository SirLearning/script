nextflow.enable.dsl=2

process mk_plink_basic_info_camp_pop_with_filter {
    tag "make plink basic info for camp pop: ${chr}"
    publishDir "${params.output_dir}/${params.job}/process/sample", mode: 'copy', pattern: "*.{smiss,scount}"
    publishDir "${params.output_dir}/${params.job}/process/variant", mode: 'copy', pattern: "*.{vmiss,gcount,afreq,hardy}"
    publishDir "${params.output_dir}/${params.job}/process/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'
    
    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)
    path camp_vmap4_map_tsv

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
    exec > ${chr}.info.log 2>&1

    # Extract sample IDs from column 2 (Bams) and column 3 (Bams-2), skipping the header
    awk -F'\\t' 'NR>1 {
        if(\$2 != "") print \$2;
        if(NF>=3 && \$3 != "") print \$3;
    }' ${camp_vmap4_map_tsv} > camp_pop.id

    plink2 --pfile ${prefix} \\
        --keep camp_pop.id \\
        --allow-extra-chr \\
        --geno 0.1 \\
        --maf 0.01 \\
        --missing \\
        --sample-counts \\
        --geno-counts \\
        --freq \\
        --hardy \\
        --threads ${task.cpus} \\
        --out ${id}.info
    """
}

process mk_plink_basic_info_camp_pop {
    tag "make plink basic info for camp pop: ${chr}"
    publishDir "${params.output_dir}/${params.job}/process/sample", mode: 'copy', pattern: "*.{smiss,scount}"
    publishDir "${params.output_dir}/${params.job}/process/variant", mode: 'copy', pattern: "*.{vmiss,gcount,afreq,hardy}"
    publishDir "${params.output_dir}/${params.job}/process/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'
    
    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)
    path camp_vmap4_map_tsv

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
    exec > ${chr}.info.log 2>&1

    # Extract sample IDs from column 2 (Bams) and column 3 (Bams-2), skipping the header
    awk -F'\\t' 'NR>1 {
        if(\$2 != "") print \$2;
        if(NF>=3 && \$3 != "") print \$3;
    }' ${camp_vmap4_map_tsv} > camp_pop.id

    plink2 --pfile ${prefix} \\
        --keep camp_pop.id \\
        --allow-extra-chr \\
        --missing \\
        --sample-counts \\
        --geno-counts \\
        --freq \\
        --hardy \\
        --threads ${task.cpus} \\
        --out ${id}.info
    """
}

process plink2_pca {
    tag "plink2 pca: ${id}"
    label 'cpus_8'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/process/${params.wheat_plink_source_mod ?: params.mod}/integrated/pca", mode: 'copy', pattern: "*.{eigenvec,eigenval,log}"
    publishDir "${params.output_dir}/${params.job}/process/${params.wheat_plink_source_mod ?: params.mod}/integrated/pca/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), path("${id}.pca.eigenvec"), path("${id}.pca.eigenval"), emit: pca

    script:
    """
    set -euo pipefail
    exec > ${chr}.plink2_pca.log 2>&1
    plink2 --pfile ${prefix} \\
        --allow-extra-chr \\
        --pca approx ${params.wheat_pca_n_pcs} \\
        --threads ${task.cpus} \\
        --out ${id}.pca
    """
}

process plink2_tagsnp_prune {
    tag "plink2 tagsnp: ${id}"
    label 'cpus_8'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/process/${params.wheat_plink_source_mod ?: params.mod}/integrated/tagsnp", mode: 'copy', pattern: "*.{prune.in,prune.out,log}"
    publishDir "${params.output_dir}/${params.job}/process/${params.wheat_plink_source_mod ?: params.mod}/integrated/tagsnp/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), path("${id}.tagsnp.prune.in"), emit: prune

    script:
    """
    set -euo pipefail
    exec > ${chr}.plink2_tagsnp.log 2>&1
    plink2 --pfile ${prefix} \\
        --allow-extra-chr \\
        --indep-pairwise 50 5 ${params.wheat_tagsnp_ld_threshold} \\
        --threads ${task.cpus} \\
        --out ${id}.tagsnp
    """
}

process awk_depth_cnv_call {
    tag "awk cnv"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/compute", mode: 'copy', pattern: "*.cnv.tsv"

    input:
    path depth_matrix
    val output_prefix

    output:
    path "${output_prefix}.cnv.tsv", emit: cnv

    script:
    """
    set -euo pipefail
    exec > ${output_prefix}.cnv.log 2>&1
    awk -F'\\t' -v del_z=${params.wheat_cnv_del_z} -v dup_z=${params.wheat_cnv_dup_z} '
    NR==FNR && FNR==1 {
        for (i=4; i<=NF; i++) { sn[i]=\$i; sum[i]=0; sumsq[i]=0; n[i]=0 }
        next
    }
    NR==FNR && FNR>1 {
        for (i=4; i<=NF; i++) {
            v=\$i+0
            if (v==v) { sum[i]+=v; sumsq[i]+=v*v; n[i]++ }
        }
        next
    }
    FNR==1 && NR!=FNR {
        for (i=4; i<=NF; i++) {
            if (n[i]>1) { mu[i]=sum[i]/n[i]; sd[i]=sqrt((sumsq[i]-n[i]*mu[i]*mu[i])/(n[i]-1)) }
            else if (n[i]==1) { mu[i]=sum[i]; sd[i]=0 }
            else { mu[i]=0; sd[i]=0 }
        }
        print "Sample\\tCHR\\tSTART\\tEND\\tZScore\\tCNVType"
        next
    }
    FNR>1 && NR!=FNR {
        chr=\$1; st=\$2; en=\$3
        for (i=4; i<=NF; i++) {
            v=\$i+0
            if (v!=v) continue
            z=(sd[i]>0)?(v-mu[i])/sd[i]:0
            state="NORMAL"
            if (z<=del_z) state="DELETION"
            else if (z>=dup_z) state="DUPLICATION"
            if (state!="NORMAL")
                print sn[i]"\\t"chr"\\t"st"\\t"en"\\t"z"\\t"state
        }
    }' "${depth_matrix}" "${depth_matrix}" > ${output_prefix}.cnv.tsv
    """
}

/*
 * Single PLINK2 pass for per-dataset summary statistics on pfiles.
 * Flags: --missing, --sample-counts, --geno-counts, --freq, --hardy.
 * Downstream modules (stats, wheat_snp_qc, assess) must reuse
 * process/<mod>/sample/*.info.{smiss,scount} and variant/*.info.{vmiss,gcount,afreq,hardy}
 * instead of re-running --freq/--missing. Mode-specific flags (--pca, --indep-pairwise,
 * --r2-unphased, --glm) stay in dedicated processes.
 */
process mk_plink_basic_info {
    tag "make plink basic info: ${chr}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/sample", mode: 'copy', pattern: "*.{smiss,scount}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/variant", mode: 'copy', pattern: "*.{vmiss,gcount,afreq,hardy}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.log"
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

process calc_plink_ld_unphased {
    tag "make plink LD info: ${chr}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/variant", mode: 'copy', pattern: "*.info.ld.vcor"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.ld.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), path("${id}.info.ld.vcor"), emit: ld
    tuple val(id), val(chr), path("${id}.info.ld.log"), emit: log

    script:
    """
    set -euo pipefail
    exec > ${chr}.ld.log 2>&1

    plink2 --pfile ${prefix} \\
        --allow-extra-chr \\
        --r2-unphased \\
        --ld-window-kb ${params.ld_window_kb} \\
        --ld-window ${params.ld_window} \\
        --ld-window-r2 ${params.ld_window_r2} \\
        --threads ${task.cpus} \\
        --out ${id}.info.ld
    """
}

process calc_plink_ld_crosschr_random {
    tag "make cross-chr plink LD info: ${chr}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/variant", mode: 'copy', pattern: "*.info.ld.crosschr.vcor"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.ld.crosschr.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), path("${id}.info.ld.crosschr.vcor"), emit: ld
    tuple val(id), val(chr), path("${id}.info.ld.crosschr.log"), emit: log

    script:
    """
    set -euo pipefail
    exec > ${chr}.ld.crosschr.log 2>&1

    awk -v seed=${params.ld_decay_random_seed} -v rate=${params.ld_crosschr_sample_rate} 'BEGIN{srand(seed)} NR>1 { if (rand() <= rate) print \$3 }' ${pvar} > ${id}.crosschr.sample.id
    sampled_n=\$(wc -l < ${id}.crosschr.sample.id || true)
    if [ "\$sampled_n" -lt 2 ]; then
        awk 'NR>1 {print \$3}' ${pvar} | head -n ${params.ld_crosschr_min_variants} > ${id}.crosschr.sample.id
    fi

    plink2 --pfile ${prefix} \\
        --allow-extra-chr \\
        --extract ${id}.crosschr.sample.id \\
        --r2-unphased inter-chr \\
        --ld-window-r2 ${params.ld_window_r2} \\
        --threads ${task.cpus} \\
        --out ${id}.info.ld.crosschr
    """
}
