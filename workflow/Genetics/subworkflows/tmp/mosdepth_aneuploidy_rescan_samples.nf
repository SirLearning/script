#!/usr/bin/env nextflow
/*
 * Re-scan mosdepth aneuploidy for a sample subset and patch the published benchmark tables.
 *
 * Run cwd:
 *   /data/home/tusr1/01projects/vmap4/11aneuploidy/12run_abd_rescan_patch
 *
 * Sample list TSV (no header): sample<TAB>cohort
 *
 * Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/11aneuploidy/12run_abd_rescan_patch/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/11aneuploidy/12run_abd_rescan_patch
 *   nextflow run .../subworkflows/tmp/mosdepth_aneuploidy_rescan_samples.nf \\
 *     -c .../nextflow.config \\
 *     --home_dir /data/home/tusr1/01projects/vmap4 \\
 *     --user_dir /data/home/tusr1 \\
 *     --src_dir /data/home/tusr1/git/script/src \\
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \\
 *     --job benchmark \\
 *     --mod aneuploidy_check \\
 *     --depth_root /data/home/tusr1/01projects/vmap4/00data/04depth \\
 *     --rescan_samples_file /path/to/rescan_samples.tsv
 */

nextflow.enable.dsl=2

include {
    mosdepth_aneuploidy_sample
    mosdepth_aneuploidy_patch_publish
} from '../../modules/local/germplasm/mosdepth_aneuploidy.nf'

params.job = 'benchmark'
params.mod = 'aneuploidy_check'
params.depth_root = "${params.home_dir}/00data/04depth"
params.rescan_samples_file = null

def stripMosdepthSummarySuffixes(String s) {
    def stripPattern = /(?:\.mosdepth\.summary\.txt|\.rmdup\.bam|_deduped\.bam|\.bam|_deduped|\.rmdup)+$/
    return s.replaceAll(stripPattern, '')
}

def normalizeSampleFromSummaryName(String name) {
    return stripMosdepthSummarySuffixes(name)
}

def resolveSummaryFile(String cohort, String sample) {
    def cohortDir = file("${params.depth_root}/${cohort}", type: 'dir')
    if (!cohortDir.exists()) {
        error "mosdepth aneuploidy rescan: cohort dir missing: ${cohortDir}"
    }
    def exact = file("${cohortDir}/${sample}.bam.mosdepth.summary.txt")
    if (exact.exists()) {
        return exact
    }
    def matches = cohortDir.listFiles({ f ->
        f.name.endsWith('.mosdepth.summary.txt') &&
            normalizeSampleFromSummaryName(f.name) == sample
    })
    if (matches.size() == 1) {
        return matches[0]
    }
    error "mosdepth aneuploidy rescan: summary not found for ${cohort}:${sample}"
}

workflow {
    if (!params.depth_root) {
        error 'mosdepth_aneuploidy_rescan_samples: --depth_root is required'
    }
    if (!params.rescan_samples_file) {
        error 'mosdepth_aneuploidy_rescan_samples: --rescan_samples_file is required'
    }
    if (!params.output_dir || !params.job || !params.mod) {
        error 'mosdepth_aneuploidy_rescan_samples: --output_dir, --job, and --mod are required'
    }

    def stats_base = "${params.output_dir}/${params.job}/stats/${params.mod}"

    summary_in = Channel
        .fromPath(params.rescan_samples_file)
        .splitCsv(sep: '\t', header: false)
        .map { row ->
            def sample = row[0].toString().trim()
            def cohort = row[1].toString().trim()
            tuple(cohort, sample, resolveSummaryFile(cohort, sample))
        }

    per_sample = mosdepth_aneuploidy_sample(summary_in)
    mosdepth_aneuploidy_patch_publish(
        per_sample.chr.collect(),
        per_sample.sample.collect(),
        stats_base,
    )
}
