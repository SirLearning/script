#!/usr/bin/env nextflow
/*
 * Plot Wu et al. 2025 LD block size distribution from Table S26 CSV.
 *
 * Publishes to params.output_dir/params.job/stats/params.mod/{info,plots,logs}.
 * Default job=meta_data, mod=ld_block_wu2025 under vmap4.VCF.v1 publish tree.
 *
 * Launch from a vmap4 run folder; log in doc/progress/YYYY-MM-DD.md (#### NF replay). Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/00data/10meta_data/01run_ld_block_plot
 *   cd /data/home/tusr1/01projects/vmap4/00data/10meta_data/01run_ld_block_plot
 *   nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/ld_block_size_plot.nf \\
 *     -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \\
 *     --home_dir /data/home/tusr1/01projects/vmap4 \\
 *     --user_dir /data/home/tusr1 \\
 *     --src_dir /data/home/tusr1/git/script/src \\
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \\
 *     --job meta_data \\
 *     --mod ld_block_wu2025 \\
 *     --ld_block_csv /data/home/tusr1/01projects/vmap4/00data/10meta_data/Wu_2025_LD_block_size.csv
 */

nextflow.enable.dsl=2

include { ld_block_size_plot } from '../../modules/local/germplasm/ld_block_size_plot.nf'

params.job = 'meta_data'
params.mod = 'ld_block_wu2025'
params.ld_block_csv = null

workflow {
    if (!params.output_dir) {
        error 'ld_block_size_plot: --output_dir is required'
    }
    if (!params.job || !params.mod) {
        error 'ld_block_size_plot: --job and --mod are required'
    }
    if (!params.ld_block_csv) {
        error 'ld_block_size_plot: --ld_block_csv is required'
    }

    ld_block_csv = file(params.ld_block_csv)
    ld_block_size_plot(ld_block_csv)
}
