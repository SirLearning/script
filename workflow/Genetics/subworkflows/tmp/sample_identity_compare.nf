#!/usr/bin/env nextflow
/*
 * Five sample-identity methods on merged A/B/D test pfiles/bfiles:
 *   1. KING kinship          (plink2 --make-king-table)
 *   2. Rare-allele Jaccard   (--max-maf, then ALT-carrier Jaccard — not rare IBS)
 *   3. Genotype discordance  (hardcall discord on a common SNP panel)
 *   4. PLINK IBD             (--indep-pairwise, then --genome PI_HAT)
 *   5. Somalier-style panel  (same common SNP panel; somalier is not installed)
 *
 * Restricted to test_thin / test_common_inter so the pairwise work stays tractable.
 *
 *   cd /data/home/tusr1/01projects/vmap4/12stats.genome/<NNrun_sample_identity_*>
 *   nextflow run .../subworkflows/tmp/sample_identity_compare.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_thin
 */

nextflow.enable.dsl=2

include {
    listMergedSubgenomeTestPfileTuples
    listMergedSubgenomeTestBfileTuples
} from '../../modules/local/infra/infra_plink_reuse.nf'

params.job = 'test_plink'
params.mod = null
params.identity_king_min = 0.177
params.identity_rare_max_maf = 0.01
params.identity_jaccard_min = 0.5
params.identity_jaccard_min_shared = 5
params.identity_jaccard_max_carriers = 400
params.identity_fp_maf = 0.1
params.identity_fp_snps = 8000
params.identity_discord_max = 0.05
params.identity_somalier_min = 0.4
params.identity_fp_min_nobs = 100
params.identity_pi_hat_min = 0.25
params.identity_ibd_max_snps = 25000
params.identity_min_methods = 2

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'sample_identity_compare: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir) {
        error 'sample_identity_compare: --user_dir is required (Conda stats env)'
    }
    if (!(params.mod in ['test_thin', 'test_common_inter'])) {
        error "sample_identity_compare: unsupported mod ${params.mod} (use test_thin or test_common_inter)"
    }

    def proc = params.process_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}"
    def ch_pfile = channel.from(listMergedSubgenomeTestPfileTuples(proc))
        .filter { id, chr, prefix, pgen, psam, pvar -> id in ['A', 'B', 'D'] }
    def ch_bfile = channel.from(listMergedSubgenomeTestBfileTuples(proc))
        .filter { id, chr, prefix, bed, bim, fam -> id in ['A', 'B', 'D'] }

    def king_raw = plink2_identity_king(ch_pfile)
    def rare_vcf = plink2_identity_rare_vcf(ch_pfile)
    def fp_raw = plink2_identity_fingerprint_raw(ch_pfile)
    def genome_raw = plink1_identity_genome(ch_bfile)

    def king_pairs = ana_identity_king(king_raw)
    def jac_pairs = ana_identity_jaccard(rare_vcf)
    def fp_pairs = ana_identity_fingerprint(fp_raw)
    def ibd_pairs = ana_identity_genome(genome_raw)

    ana_identity_consensus_sg(
        king_pairs
            .join(jac_pairs, by: [0, 1])
            .join(fp_pairs, by: [0, 1])
            .join(ibd_pairs, by: [0, 1])
    )
}

process plink2_identity_king {
    tag "identity king: ${id}"
    label 'cpus_8'
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/identity", mode: 'copy', pattern: "*.kin0"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/identity/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), path("${id}.identity.king.kin0")

    script:
    """
    set -euo pipefail
    exec > ${chr}.identity_king.log 2>&1
    plink2 --pfile ${prefix} \\
        --chr-set 44 --allow-extra-chr \\
        --make-king-table \\
        --king-table-filter ${params.identity_king_min} \\
        --threads ${task.cpus} \\
        --out ${id}.identity.king
    """
}

process plink2_identity_rare_vcf {
    tag "identity rare vcf: ${id}"
    label 'cpus_8'
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/identity", mode: 'copy', pattern: "*.vcf"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/identity/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), path("${id}.identity.rare.vcf")

    script:
    """
    set -euo pipefail
    exec > ${chr}.identity_rare_vcf.log 2>&1
    if plink2 --pfile ${prefix} \\
        --chr-set 44 --allow-extra-chr \\
        --max-maf ${params.identity_rare_max_maf} \\
        --export vcf \\
        --threads ${task.cpus} \\
        --out ${id}.identity.rare
    then
        :
    else
        echo "No rare sites at --max-maf ${params.identity_rare_max_maf}; writing empty VCF"
        printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\n' \\
            > ${id}.identity.rare.vcf
    fi
    """
}

process plink2_identity_fingerprint_raw {
    tag "identity fingerprint: ${id}"
    label 'cpus_8'
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/identity", mode: 'copy', pattern: "*.raw"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/identity/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), path("${id}.identity.fp.raw")

    script:
    """
    set -euo pipefail
    exec > ${chr}.identity_fp.log 2>&1
    if ! plink2 --pfile ${prefix} \\
        --chr-set 44 --allow-extra-chr \\
        --maf ${params.identity_fp_maf} \\
        --thin-count ${params.identity_fp_snps} \\
        --export A \\
        --threads ${task.cpus} \\
        --out ${id}.identity.fp
    then
        echo "thin-count ${params.identity_fp_snps} exceeded remaining SNPs; exporting all MAF>=${params.identity_fp_maf}"
        plink2 --pfile ${prefix} \\
            --chr-set 44 --allow-extra-chr \\
            --maf ${params.identity_fp_maf} \\
            --export A \\
            --threads ${task.cpus} \\
            --out ${id}.identity.fp
    fi
    """
}

process plink1_identity_genome {
    tag "identity ibd: ${id}"
    label 'cpus_16'
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/identity", mode: 'copy', pattern: "*.genome"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/identity/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), val(prefix), path(bed), path(bim), path(fam)

    output:
    tuple val(id), val(chr), path("${id}.identity.ibd.genome")

    script:
    """
    set -euo pipefail
    exec > ${chr}.identity_ibd.log 2>&1
    plink --bfile ${prefix} \\
        --chr-set 44 --allow-extra-chr --allow-no-sex \\
        --indep-pairwise ${params.pca_ld_window} ${params.pca_ld_step} ${params.pca_ld_r2} \\
        --threads ${task.cpus} \\
        --out ${id}.identity.prune
    n_prune=\$(wc -l < ${id}.identity.prune.prune.in)
    if [ "\${n_prune}" -gt ${params.identity_ibd_max_snps} ]; then
        plink2 --bfile ${prefix} \\
            --chr-set 44 --allow-extra-chr \\
            --extract ${id}.identity.prune.prune.in \\
            --thin-count ${params.identity_ibd_max_snps} \\
            --make-bed \\
            --threads ${task.cpus} \\
            --out ${id}.identity.ibd.in
    else
        plink2 --bfile ${prefix} \\
            --chr-set 44 --allow-extra-chr \\
            --extract ${id}.identity.prune.prune.in \\
            --make-bed \\
            --threads ${task.cpus} \\
            --out ${id}.identity.ibd.in
    fi
    plink --bfile ${id}.identity.ibd.in \\
        --chr-set 44 --allow-extra-chr --allow-no-sex \\
        --genome --min ${params.identity_pi_hat_min} \\
        --threads ${task.cpus} \\
        --out ${id}.identity.ibd
    """
}

process ana_identity_king {
    tag "ana king: ${id}"
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.pairs.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), path(king_file)

    output:
    tuple val(id), val(chr), path("${id}.sample.identity.king.pairs.tsv")

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.identity_king_ana.log", "w")
    sys.stderr = sys.stdout
    from genetics.genomics.sample.identity import ana_king_identity
    ana_king_identity("${king_file}", "${id}.sample.identity.king", min_kinship=${params.identity_king_min})
    """
}

process ana_identity_jaccard {
    tag "ana jaccard: ${id}"
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.pairs.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), path(vcf_file)

    output:
    tuple val(id), val(chr), path("${id}.sample.identity.jaccard.pairs.tsv")

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.identity_jaccard_ana.log", "w")
    sys.stderr = sys.stdout
    from genetics.genomics.sample.identity import ana_rare_jaccard_vcf
    ana_rare_jaccard_vcf(
        "${vcf_file}",
        "${id}.sample.identity.jaccard",
        min_shared=${params.identity_jaccard_min_shared},
        min_jaccard=${params.identity_jaccard_min},
        max_carriers=${params.identity_jaccard_max_carriers},
    )
    """
}

process ana_identity_fingerprint {
    tag "ana fingerprint: ${id}"
    label 'cpus_8'
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.pairs.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), path(raw_file)

    output:
    tuple val(id), val(chr), path("${id}.sample.identity.discord.pairs.tsv"), path("${id}.sample.identity.fingerprint.pairs.tsv")

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.identity_fp_ana.log", "w")
    sys.stderr = sys.stdout
    from genetics.genomics.sample.identity import ana_fingerprint_additive
    ana_fingerprint_additive(
        "${raw_file}",
        "${id}.sample.identity",
        max_discord=${params.identity_discord_max},
        min_relatedness=${params.identity_somalier_min},
        min_nobs=${params.identity_fp_min_nobs},
    )
    """
}

process ana_identity_genome {
    tag "ana ibd: ${id}"
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.pairs.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), path(genome_file)

    output:
    tuple val(id), val(chr), path("${id}.sample.identity.ibd.pairs.tsv")

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.identity_ibd_ana.log", "w")
    sys.stderr = sys.stdout
    from genetics.genomics.sample.identity import ana_genome_ibd
    ana_genome_ibd("${genome_file}", "${id}.sample.identity.ibd", min_pi_hat=${params.identity_pi_hat_min})
    """
}

process ana_identity_consensus_sg {
    tag "ana consensus: ${id}"
    errorStrategy 'terminate'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), path(king_pairs), path(jaccard_pairs), path(discord_pairs), path(fingerprint_pairs), path(ibd_pairs)

    output:
    tuple val(id), val(chr), path("${id}.sample.identity.consensus.pairs.tsv"), path("${id}.sample.identity.consensus.clusters.tsv")

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.identity_consensus_ana.log", "w")
    sys.stderr = sys.stdout
    from infra.utils.io import load_df_from_tsv
    from genetics.genomics.sample.identity import ana_identity_consensus
    tables = {
        "king": load_df_from_tsv("${king_pairs}"),
        "jaccard": load_df_from_tsv("${jaccard_pairs}"),
        "discord": load_df_from_tsv("${discord_pairs}"),
        "fingerprint": load_df_from_tsv("${fingerprint_pairs}"),
        "ibd": load_df_from_tsv("${ibd_pairs}"),
    }
    ana_identity_consensus(tables, "${id}.sample.identity", min_methods=${params.identity_min_methods})
    """
}
