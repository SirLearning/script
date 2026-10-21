#!/usr/bin/env nextflow
/*
 * Plot species LD decay distance vs genome size from 00data/10meta_data/species_LD.csv.
 *
 * Publishes to params.output_dir/params.job/stats/params.mod/{info,plots,logs}.
 * Default job=meta_data, mod=species_ld under vmap4.VCF.v1 publish tree.
 *
 * Launch from a vmap4 run folder; log in doc/progress/YYYY-MM-DD.md (#### NF replay). Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/12stats.genome/92run_species_ld_vs_genome
 *   cd /data/home/tusr1/01projects/vmap4/12stats.genome/92run_species_ld_vs_genome
 *   nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/species_ld_plot.nf \\
 *     -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \\
 *     --home_dir /data/home/tusr1/01projects/vmap4 \\
 *     --user_dir /data/home/tusr1 \\
 *     --src_dir /data/home/tusr1/git/script/src \\
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \\
 *     --job meta_data \\
 *     --mod species_ld \\
 *     --species_ld_csv /data/home/tusr1/01projects/vmap4/00data/10meta_data/species_LD.csv
 */

nextflow.enable.dsl=2

include { species_ld_plot } from '../../modules/local/germplasm/species_ld_plot.nf'

params.job = 'meta_data'
params.mod = 'species_ld'
params.species_ld_csv = null

workflow {
    if (!params.output_dir) {
        error 'species_ld_plot: --output_dir is required'
    }
    if (!params.job || !params.mod) {
        error 'species_ld_plot: --job and --mod are required'
    }
    if (!params.species_ld_csv) {
        error 'species_ld_plot: --species_ld_csv is required'
    }

    species_ld_csv = file(params.species_ld_csv)
    species_ld_plot(species_ld_csv)
}
