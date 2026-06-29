#!/usr/bin/env nextflow
/*
 * Parallel mosdepth aneuploidy scan from *.mosdepth.summary.txt under depth cohort folders.
 *
 * One Nextflow task per summary file (parallel). Publishes under the benchmark job tree:
 *   {output_dir}/{job}/stats/{mod}/{info,logs,plots}/…
 *
 * Run cwd (level 1 module under vmap4, not under 00data):
 *   /data/home/tusr1/01projects/vmap4/11aneuploidy/01run_benchmark_full
 *
 * Publish tree:
 *   {output_dir}/{job}/stats/{mod}/info/{sample_aneuploidy_summary,flagged_samples}.tsv
 *   {output_dir}/{job}/stats/{mod}/info/chr_depth/{chr1A,...}.depth.tsv
 *   {output_dir}/{job}/stats/{mod}/plots/*.png   (mosdepth_aneuploidy_replot.nf only)
 *   {output_dir}/{job}/stats/{mod}/logs/
 *
 * Example (full cohort):
 *   mkdir -p /data/home/tusr1/01projects/vmap4/11aneuploidy/01run_benchmark_full/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/11aneuploidy/01run_benchmark_full
 *   nextflow run .../subworkflows/tmp/mosdepth_aneuploidy_scan.nf \\
 *     -c .../nextflow.config -c .../conf/mosdepth_aneuploidy.config \\
 *     --home_dir /data/home/tusr1/01projects/vmap4 \\
 *     --user_dir /data/home/tusr1 \\
 *     --src_dir /data/home/tusr1/git/script/src \\
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \\
 *     --job benchmark \\
 *     --mod aneuploidy_check \\
 *     --depth_root /data/home/tusr1/01projects/vmap4/00data/04depth \\
 *     --mosdepth_cohorts 01A \\
 *     --mosdepth_aneuploidy_max_forks 16
 */

nextflow.enable.dsl=2

include {
    mosdepth_aneuploidy_sample
    mosdepth_aneuploidy_collect
} from '../../modules/local/germplasm/mosdepth_aneuploidy.nf'

params.job = 'benchmark'
params.mod = 'aneuploidy_check'
params.depth_root = "${params.home_dir}/00data/04depth"
params.mosdepth_cohorts = '01A,02AB,03ABD,04D,05HZNU,06Nature,07S,08WAP,09Watkins'
params.mosdepth_aneuploidy_max_forks = 16

def stripMosdepthSummarySuffixes(String s) {
    // Longest suffixes first; repeated tail match mirrors Python normalize_sample_id().
    def stripPattern = /(?:\.mosdepth\.summary\.txt|\.rmdup\.bam|_deduped\.bam|\.bam|_deduped|\.rmdup)+$/
    return s.replaceAll(stripPattern, '')
}

def normalizeSampleFromSummaryName(String name) {
    return stripMosdepthSummarySuffixes(name)
}

workflow {
    if (!params.depth_root) {
        error 'mosdepth_aneuploidy_scan: --depth_root is required'
    }
    if (!params.output_dir || !params.job || !params.mod) {
        error 'mosdepth_aneuploidy_scan: --output_dir, --job, and --mod are required'
    }

    cohort_list = params.mosdepth_cohorts.split(',').collect { it.trim() }.findAll { it }

    summary_in = Channel
        .from(cohort_list)
        .flatMap { cohort ->
            def dir = file("${params.depth_root}/${cohort}", type: 'dir')
            if (!dir.exists()) {
                log.warn "mosdepth aneuploidy: cohort dir missing: ${dir}"
                return []
            }
            dir.listFiles({ f -> f.name.endsWith('.mosdepth.summary.txt') })
                .collect { summary ->
                    tuple(cohort, normalizeSampleFromSummaryName(summary.name), summary)
                }
        }

    per_sample = mosdepth_aneuploidy_sample(summary_in)
    mosdepth_aneuploidy_collect(
        per_sample.chr.collect(),
        per_sample.sample.collect(),
    )
}
