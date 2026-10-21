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

process plink2_het_inbreeding {
    tag "plink2 het inbreeding: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.het_source_mod}/sample", mode: 'copy', pattern: "*.het"
    publishDir "${params.output_dir}/${params.job}/process/${params.het_source_mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 8
    memory '16.GB'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), path("${id}.info.het"), emit: het
    tuple val(id), val(chr), path("${id}.info.het.log"), emit: log

    script:
    """
    set -euo pipefail
    exec > ${id}.info.het.log 2>&1

    plink2 --pfile ${prefix} \\
        --allow-extra-chr \\
        --het \\
        --threads ${task.cpus} \\
        --out ${id}.info
    """
}

process plink2_het_extract {
    tag "plink2 het extract: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.het"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 8
    memory '16.GB'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar), path(extract)

    output:
    tuple val(id), val(chr), path("${id}.all_het.info.het"), emit: het
    tuple val(id), val(chr), path("${id}.all_het.info.het.log"), emit: log

    when:
    extract.size() > 0

    script:
    """
    set -euo pipefail
    exec > ${id}.all_het.info.het.log 2>&1

    plink2 --pfile ${prefix} \\
        --allow-extra-chr \\
        --extract ${extract} \\
        --het \\
        --threads ${task.cpus} \\
        --out ${id}.all_het.info
    """
}

process calc_plink_ld_focal_along_chr {
    tag "focal LD along chr${chrom}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/variant", mode: 'copy', pattern: "*.{vcor,focals.info.tsv,focals.id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar), val(chrom)

    output:
    tuple val(id), val(chrom), path("${id}.${chrom}.focal.ld.vcor"), path("${id}.${chrom}.focals.info.tsv"), emit: ld
    tuple val(id), val(chrom), path("${id}.${chrom}.focal.ld.log"), emit: log

    script:
    """
    set -euo pipefail
    exec > ${id}.${chrom}.focal.ld.log 2>&1

    python - <<'PY'
from genetics.genomics.variant.ld_along_chr import (
    sample_focals_per_spatial_bin_from_pvar,
    write_ref_merge_update_maps,
)

write_ref_merge_update_maps(
    "${pvar}",
    "${chrom}",
    "update.chr",
    "update.map",
)
sample_focals_per_spatial_bin_from_pvar(
    "${pvar}",
    ["${chrom}"],
    n_per_bin=${params.ld_along_chr_n_focal_per_bin},
    bin_mb=${params.ld_along_chr_spatial_bin_mb},
    output_prefix="${id}.${chrom}",
)
PY

    CHR_IDS=\$(cat ${id}.${chrom}.plink_chr_ids.txt)
    plink2 --pfile ${prefix} \\
        --allow-extra-chr \\
        --chr \${CHR_IDS} \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.${chrom}.extract

    plink2 --pfile ${id}.${chrom}.extract \\
        --allow-extra-chr \\
        --update-chr update.chr \\
        --update-map update.map \\
        --sort-vars \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.${chrom}.refmerge

    plink2 --pfile ${id}.${chrom}.refmerge \\
        --allow-extra-chr \\
        --ld-snp-list ${id}.${chrom}.focals.id \\
        --r2-unphased \\
        --ld-window-kb ${params.ld_along_chr_max_kb} \\
        --ld-window ${params.ld_window} \\
        --ld-window-r2 ${params.ld_window_r2} \\
        --threads ${task.cpus} \\
        --out ${id}.${chrom}.focal.ld
    """
}

process recode_subgenome_pfile_ref_chromosomes {
    tag "ref-chr recode: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{bed,bim,fam,pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), val("${id}_test.plink2"), path("${id}_test.plink2.pgen"), path("${id}_test.plink2.psam"), path("${id}_test.plink2.pvar"), emit: pfile
    tuple val(id), val(chr), val("${id}_test.plink"), path("${id}_test.plink.bed"), path("${id}_test.plink.bim"), path("${id}_test.plink.fam"), emit: bfile
    tuple val(id), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > recode_refchr.${id}.log 2>&1

    python - <<'PY'
from genetics.genomics.variant.ld_along_chr import write_all_ref_merge_update_maps

write_all_ref_merge_update_maps("${pvar}", "update.chr", "update.map")
PY

    if [ -s update.chr ]; then
        plink2 --pfile ${prefix} \\
            --allow-extra-chr \\
            --update-chr update.chr \\
            --update-map update.map \\
            --sort-vars \\
            --make-pgen \\
            --threads ${task.cpus} \\
            --out recoded
    else
        plink2 --pfile ${prefix} \\
            --allow-extra-chr \\
            --make-pgen \\
            --threads ${task.cpus} \\
            --out recoded
    fi

    plink2 --pfile recoded \\
        --allow-extra-chr \\
        --make-bed \\
        --threads ${task.cpus} \\
        --out recoded_bed

    mv recoded.pgen ${id}_test.plink2.pgen
    mv recoded.psam ${id}_test.plink2.psam
    mv recoded.pvar ${id}_test.plink2.pvar
    mv recoded_bed.bed ${id}_test.plink.bed
    mv recoded_bed.bim ${id}_test.plink.bim
    mv recoded_bed.fam ${id}_test.plink.fam
    """
}

process plink2_vcftools_tajima_d {
    tag "tajima D: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.Tajima.D"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 16
    memory '128.GB'
    maxForks 3
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    path "${id}.variant.tajima_d.${window_tag}.Tajima.D", emit: tajima
    path "${id}.variant.tajima_d.${window_tag}.run.log", emit: log

    script:
    def window_bp = (params.tajima_d_window_bp ?: 100000) as int
    window_tag = "${window_bp.intdiv(1000)}kb"
    """
    set -euo pipefail
    exec > ${id}.variant.tajima_d.${window_tag}.run.log 2>&1

    plink2 --pfile ${prefix} \\
        --allow-extra-chr \\
        --snps-only just-acgt \\
        --export vcf-4.2 bgz \\
        --threads ${task.cpus} \\
        --out ${id}.export

    vcftools --gzvcf ${id}.export.vcf.gz \\
        --TajimaD ${window_bp} \\
        --out ${id}.variant.tajima_d.${window_tag}

    rm -f ${id}.export.vcf.gz
    """
}

process plink1_ld_dprime_window {
    tag "plink1 r2 dprime window: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.circos_ld_mod ?: params.mod}/variant", mode: 'copy', pattern: "*.ld"
    publishDir "${params.output_dir}/${params.job}/process/${params.circos_ld_mod ?: params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'
    memory '64.GB'

    input:
    tuple val(id), val(prefix), path(bed), path(bim), path(fam)

    output:
    tuple val(id), path("${id}.info.ld_block.ld"), emit: ld
    tuple val(id), path("${id}.info.ld_block.log"), emit: log

    script:
    def snp_window = params.ld_block_snp_window ?: 20
    def window_kb = params.ld_block_window_kb ?: 2000
    """
    set -euo pipefail

    plink --bfile ${prefix} \\
        --allow-extra-chr \\
        --r2 dprime \\
        --ld-window ${snp_window} \\
        --ld-window-kb ${window_kb} \\
        --ld-window-r2 0 \\
        --out ${id}.info.ld_block
    """
}
