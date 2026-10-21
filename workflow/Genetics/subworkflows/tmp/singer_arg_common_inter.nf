#!/usr/bin/env nextflow
/*
 * WAP-China subset (~1300) of test_common_inter_refchr: 64 spread 2 Mb windows,
 * selfing-phase GT, Beagle-impute, SINGER ARG per window.
 *
 *   nextflow run .../subworkflows/tmp/singer_arg_common_inter.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_common_inter_2mb_wap_china_refchr
 *
 * Failed-window retry (reuse published Beagle VCFs; our random -seed, not
 * singer_master -debug):
 *     --singer_from_beagle true
 */

nextflow.enable.dsl=2

include {
    prepare_wap_china_singer_inputs
    extract_common_inter_window_pfile
    vcf_selfing_phase_slash_to_pipe
    beagle_impute_vcf
    singer_sample_arg
} from '../../modules/local/genotype/processor/processor_singer.nf'

params.job = 'test_plink'
params.mod = 'test_common_inter_2mb_wap_china_refchr'
params.process_dir = null

workflow {
    if (!params.output_dir || !params.job || !params.user_dir || !params.home_dir) {
        error 'singer_arg_common_inter: --output_dir, --job, --user_dir, and --home_dir are required'
    }
    def beagle_jar = file(
        params.beagle_jar
            ?: "${params.home_dir}/00data/software/beagle/beagle.22Jul22.46e.jar"
    )
    def singer_bin_dir = file(
        params.singer_bin_dir
            ?: "${params.home_dir}/00data/software/singer-0.1.9-linux"
    )
    def singer_engine = file("${singer_bin_dir}/singer")
    def have_singer = singer_engine.exists() && singer_engine.size() > 1_000_000
    if (!have_singer) {
        log.warn "SINGER engine missing under ${singer_bin_dir}; will extract, selfing-phase, and Beagle-impute only."
    }

    if (params.singer_from_beagle) {
        if (!have_singer) {
            error "singer_arg_common_inter: --singer_from_beagle needs the singer engine under ${singer_bin_dir}"
        }
        def pub = "${params.output_dir}/${params.job}/process/${params.mod}"
        def windows_tsv = file("${pub}/singer_windows.tsv")
        def export_dir = file("${pub}/export")
        def singer_out = file("${pub}/singer")
        if (!windows_tsv.exists()) {
            error "singer_arg_common_inter: missing ${windows_tsv}"
        }
        log.info "${params.c_green}SINGER ARG retry${params.c_reset}: published Beagle VCFs under ${export_dir}; skip windows with singer_*_nodes_0.txt"
        Channel
            .fromPath(windows_tsv.toString())
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                def start_bp = row.start_bp as Integer
                def end_bp = row.end_bp as Integer
                def id = "${row.subgenome}.${row.ref_chrom}.${start_bp}_${end_bp}"
                tuple(
                    id,
                    row.ref_chrom,
                    start_bp,
                    end_bp,
                    file("${export_dir}/${id}.beagle.vcf.gz"),
                    file("${singer_out}/singer_${id}_nodes_0.txt"),
                )
            }
            .filter { id, ref_chrom, start_bp, end_bp, vcf, nodes ->
                if (!vcf.exists()) {
                    log.warn "singer_arg_common_inter: skip ${id}; missing Beagle VCF ${vcf}"
                    return false
                }
                if (nodes.exists()) {
                    return false
                }
                return true
            }
            .map { id, ref_chrom, start_bp, end_bp, vcf, nodes ->
                tuple(id, ref_chrom, start_bp, end_bp, vcf)
            }
            .set { ch_beagle }
        singer_sample_arg(ch_beagle, singer_bin_dir)
    } else {
    if (!beagle_jar.exists()) {
        error "singer_arg_common_inter: Beagle jar not found: ${beagle_jar}"
    }

    def source_dir = params.process_dir
        ? params.process_dir
        : "${params.output_dir}/${params.job}/process/${params.singer_source_mod}"
    def germplasm_tsv = file(
        params.singer_germplasm_tsv
            ?: "${params.user_dir}/git/DBone/src/main/resources/db/tmp/wap_r1_germplasm_table.tsv"
    )
    if (!germplasm_tsv.exists()) {
        error "singer_arg_common_inter: germplasm table not found: ${germplasm_tsv}"
    }

    def pfiles = [:]
    ['A', 'B', 'D'].each { sg ->
        def prefix = "${source_dir}/${sg}_test.plink2"
        def pgen = file("${prefix}.pgen")
        def psam = file("${prefix}.psam")
        def pvar = file("${prefix}.pvar")
        if (!pgen.exists() || !psam.exists() || !pvar.exists()) {
            error "singer_arg_common_inter: missing ${sg}_test.plink2 under ${source_dir}"
        }
        pfiles[sg] = [pgen: pgen, psam: psam, pvar: pvar]
    }

    log.info "${params.c_green}SINGER ARG${params.c_reset}: WAP ${params.singer_country} x ${params.singer_n_windows} windows of ${params.singer_window_bp} bp from ${source_dir}"

    prepare_wap_china_singer_inputs(
        germplasm_tsv,
        pfiles['A'].psam,
        pfiles['A'].pvar,
        pfiles['B'].pvar,
        pfiles['D'].pvar,
        params.singer_n_windows,
        params.singer_window_bp,
        params.singer_min_snps,
        params.singer_country,
    )

    def ch_jobs = prepare_wap_china_singer_inputs.out.windows
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sg = row.subgenome
            def pf = pfiles[sg]
            if (pf == null) {
                error "singer_arg_common_inter: unknown subgenome in windows table: ${sg}"
            }
            def start_bp = row.start_bp as Integer
            def end_bp = row.end_bp as Integer
            tuple(
                "${sg}.${row.ref_chrom}.${start_bp}_${end_bp}",
                row.plink_chr,
                row.ref_chrom,
                start_bp,
                end_bp,
                pf.pgen,
                pf.psam,
                pf.pvar,
            )
        }

    extract_common_inter_window_pfile(ch_jobs, prepare_wap_china_singer_inputs.out.keep)
    vcf_selfing_phase_slash_to_pipe(extract_common_inter_window_pfile.out.vcf)
    beagle_impute_vcf(vcf_selfing_phase_slash_to_pipe.out.vcf, beagle_jar)
    if (have_singer) {
        singer_sample_arg(beagle_impute_vcf.out.vcf, singer_bin_dir)
    }
    }
}
