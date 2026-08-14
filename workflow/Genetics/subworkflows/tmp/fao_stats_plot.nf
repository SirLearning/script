#!/usr/bin/env nextflow
/*
 * Plot global wheat FAOSTAT production and food-supply share trends.
 *
 * Publishes to params.output_dir/params.job/stats/params.mod/{info,plots,logs}.
 * Default job=meta_data, mod=wheat_fao under vmap4.VCF.v1 publish tree.
 *
 * Launch from a vmap4 run folder; log in doc/progress/YYYY-MM-DD.md (#### NF replay). Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/00data/07fao_stats/02run_meta_data_publish
 *   cd /data/home/tusr1/01projects/vmap4/00data/07fao_stats/02run_meta_data_publish
 *   nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/fao_stats_plot.nf \\
 *     -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \\
 *     --home_dir /data/home/tusr1/01projects/vmap4 \\
 *     --user_dir /data/home/tusr1 \\
 *     --src_dir /data/home/tusr1/git/script/src \\
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \\
 *     --job meta_data \\
 *     --mod wheat_fao \\
 *     --fao_production_csv /data/home/tusr1/01projects/vmap4/00data/07fao_stats/FAOSTAT_data_en_8-10-2026_production_full.csv \\
 *     --fao_supply_csv /data/home/tusr1/01projects/vmap4/00data/07fao_stats/FAOSTAT_data_en_8-10-2026_supply.csv
 */

nextflow.enable.dsl=2

include { fao_stats_plot } from '../../modules/local/germplasm/fao_stats_plot.nf'

params.job = 'meta_data'
params.mod = 'wheat_fao'
params.fao_production_csv = null
params.fao_supply_csv = null

workflow {
    if (!params.output_dir) {
        error 'fao_stats_plot: --output_dir is required'
    }
    if (!params.job || !params.mod) {
        error 'fao_stats_plot: --job and --mod are required'
    }
    if (!params.fao_production_csv) {
        error 'fao_stats_plot: --fao_production_csv is required'
    }
    if (!params.fao_supply_csv) {
        error 'fao_stats_plot: --fao_supply_csv is required'
    }

    production_csv = file(params.fao_production_csv)
    supply_csv = file(params.fao_supply_csv)
    fao_stats_plot(production_csv, supply_csv)
}
