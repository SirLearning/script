nextflow.enable.dsl=2

/*
 * Replot all-het sample heterozygous burden from published PLINK2 --het tables.
 *
 *     nextflow run .../mac_het_rh_sample_replot.nf \
 *       -c .../nextflow.config \
 *       --user_dir /data/home/tusr1 --src_dir /data/home/tusr1/git/script/src \
 *       --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink --mod test_rare_only
 */

include { variant_mac_het_rh_sample_burden } from '../../modules/local/genotype/stats/stats_variant.nf'

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'mac_het_rh_sample_replot: --output_dir, --job, and --mod are required.'
    }
    def info_dir = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
    ch = channel.fromList(['A', 'B', 'D', 'Others']).map { sg ->
        tuple(sg, "sub_${sg}", file("${info_dir}/${sg}.all_het.info.het", checkIfExists: true))
    }
    variant_mac_het_rh_sample_burden(ch)
}
