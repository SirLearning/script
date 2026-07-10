#!/usr/bin/env nextflow
/*
 * Subset per-chromosome MAF alignments to a chosen species list (UCSC mafSpeciesSubset).
 *
 * Input MAFs: {maf_dir}/{chr}_MZ_{n}_triticeae.maf.gz (default 06maf triticeae 29-way).
 * Species list: one prefix per line (text before the first '.' in MAF source names), comments with #.
 *
 * Launch from a vmap4 run folder. Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/10stats.genome/01run_maf_species_subset/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/01run_maf_species_subset
 *   nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/maf_species_subset.nf \
 *     -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
 *     --maf_dir /data/home/tusr1/01projects/vmap4/00data/06maf \
 *     --species_list /data/home/tusr1/01projects/vmap4/00data/06maf/species_wheat_aegilops_synteny.txt \
 *     --chromosomes 1A,1B,1D \
 *     --publish_dir /data/home/tusr1/01projects/vmap4/00data/06maf/subset_wheat_aegilops
 */

nextflow.enable.dsl=2

params.maf_dir = "${params.home_dir ?: '/data/home/tusr1/01projects/vmap4'}/00data/06maf"
params.species_list = null
params.chromosomes = '1A'
params.maf_tag = '29_triticeae'
params.publish_dir = null

process maf_species_subset {
    tag "${chr}"

    conda "${params.user_dir}/miniconda3/envs/stats"

    publishDir "${params.publish_dir}", mode: 'copy', pattern: '*.maf.gz'
    publishDir "${params.publish_dir}", mode: 'copy', pattern: '*.log'

    input:
    tuple val(chr), path(maf_gz), path(species_txt)

    output:
    tuple val(chr), path("${chr}_subset.maf.gz"), emit: maf

    script:
    def stem = maf_gz.name.replaceAll(/\\.maf\\.gz\$/, '')
    """
    set -euo pipefail
    exec > maf_species_subset_${chr}.log 2>&1

    echo "Input MAF: ${maf_gz}"
    echo "Species list: ${species_txt}"

    zcat ${maf_gz} > ${chr}_full.maf
    mafSpeciesSubset ${chr}_full.maf ${species_txt} ${chr}_subset.maf

    echo "Species retained:"
    grep '^s ' ${chr}_subset.maf | cut -f2 | cut -d '.' -f1 | sort | uniq -c | sort -rn

    gzip -f ${chr}_subset.maf
    rm -f ${chr}_full.maf
    """
}

workflow {
    if (!params.species_list) {
        error 'maf_species_subset: --species_list is required (one species prefix per line)'
    }
    if (!params.publish_dir) {
        error 'maf_species_subset: --publish_dir is required (output directory for subset MAFs)'
    }

    def species_file = file(params.species_list)
    def chr_list = (params.chromosomes as String)
        .split(',')
        .collect { it.trim() }
        .findAll { it }

    ch = channel.from(chr_list).map { chr ->
        def maf_path = "${params.maf_dir}/${chr}_MZ_${params.maf_tag}.maf.gz"
        tuple(chr, file(maf_path), species_file)
    }

    maf_species_subset(ch)
}
