#!/usr/bin/env nextflow
/*
 * Extract k-of-n syntenic alignment columns from per-chromosome subset MAFs.
 *
 * Prerequisite: maf_species_subset outputs under {maf_dir}/subset_wheat_aegilops/.
 *
 * Launch from a vmap4 run folder. Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/10stats.genome/03run_maf_syntenic_sites/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/03run_maf_syntenic_sites
 *   nextflow run .../subworkflows/tmp/maf_syntenic_sites.nf -c .../nextflow.config \
 *     --user_dir /data/home/tusr1 --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --maf_dir /data/home/tusr1/01projects/vmap4/00data/06maf/subset_wheat_aegilops \
 *     --species_list /data/home/tusr1/01projects/vmap4/00data/06maf/species_wheat_aegilops_synteny.txt \
 *     --chromosomes 1A,1B,1D \
 *     --min_species 2 \
 *     --publish_dir /data/home/tusr1/01projects/vmap4/00data/06maf/syntenic_sites_k2
 */

nextflow.enable.dsl=2

params.maf_dir = "${params.home_dir ?: '/data/home/tusr1/01projects/vmap4'}/00data/06maf/subset_wheat_aegilops"
params.species_list = "${params.home_dir ?: '/data/home/tusr1/01projects/vmap4'}/00data/06maf/species_wheat_aegilops_synteny.txt"
params.chromosomes = '1A'
params.min_species = 2
params.min_score = 0.0
params.filter_cross_chr = true
params.publish_dir = null

process maf_syntenic_sites {
    tag "${chr}"

    conda "${params.user_dir}/miniconda3/envs/stats"

    publishDir "${params.publish_dir}", mode: 'copy', pattern: '*.matrix.tsv.gz'
    publishDir "${params.publish_dir}", mode: 'copy', pattern: '*.summary.tsv'
    publishDir "${params.publish_dir}", mode: 'copy', pattern: '*.th.tsv'
    publishDir "${params.publish_dir}", mode: 'copy', pattern: '*.log'

    input:
    tuple val(chr), path(maf_gz), path(species_txt)

    output:
    tuple val(chr), path("${chr}.syntenic.matrix.tsv.gz"), emit: matrix

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.syntenic.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.alignment.maf_syntenic import extract_syntenic_sites_from_maf

    print("Extracting syntenic sites for ${chr} (min_species=${params.min_species})...")
    summary = extract_syntenic_sites_from_maf(
        "${maf_gz}",
        "${chr}.syntenic",
        "${chr}",
        species_list_path="${species_txt}",
        min_species=${params.min_species},
        min_score=${params.min_score},
        filter_cross_chr=${params.filter_cross_chr ? 'True' : 'False'},
    )
    for key, value in summary.items():
        print(f"{key}={value}")
    """
}

workflow {
    if (!params.publish_dir) {
        error 'maf_syntenic_sites: --publish_dir is required'
    }

    def species_file = file(params.species_list)
    def chr_list = (params.chromosomes as String)
        .split(',')
        .collect { it.trim() }
        .findAll { it }

    ch = channel.from(chr_list).map { chr ->
        tuple(chr, file("${params.maf_dir}/${chr}_subset.maf.gz"), species_file)
    }

    maf_syntenic_sites(ch)
}
