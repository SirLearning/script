#!/usr/bin/env nextflow
/*
 * Per-segment sample missing rate (.smiss) from test_thin chrNNN.thin pfiles.
 *
 * Input pfiles (read-only):
 *   {output_dir}/{thin_job}/process/{thin_mod}/chrNNN.thin.{pgen,psam,pvar}
 *
 * Publish:
 *   {output_dir}/{job}/process/{mod}/sample/chrNNN.smiss
 *   {output_dir}/{job}/process/{mod}/info/sample_chr_miss.{segments,long}.tsv
 *   {output_dir}/{job}/process/{mod}/logs/
 *
 * Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/13sample_chr_miss/01run_benchmark
 *   cd /data/home/tusr1/01projects/vmap4/13sample_chr_miss/01run_benchmark
 *   nextflow run .../subworkflows/tmp/sample_chr_miss.nf \\
 *     -c .../nextflow.config -c .../conf/sample_chr_miss.config \\
 *     --user_dir /data/home/tusr1 \\
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \\
 *     --job benchmark \\
 *     --mod sample_chr_miss \\
 *     --thin_job test_plink \\
 *     --thin_mod test_thin
 */

nextflow.enable.dsl=2

include {
    plink2_sample_missing_per_chr
    collect_sample_chr_miss
} from '../../modules/local/genotype/processor/processor_sample_chr_miss.nf'

params.job = 'benchmark'
params.mod = 'sample_chr_miss'
params.thin_job = 'test_plink'
params.thin_mod = 'test_thin'

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'sample_chr_miss: --output_dir, --job, and --mod are required'
    }

    def thin_root = "${params.output_dir}/${params.thin_job}/process/${params.thin_mod}"

    chr_pfiles = Channel
        .fromPath("${thin_root}/chr*.thin.pgen", checkIfExists: true)
        .map { pgen ->
            def base = pgen.getName().replaceAll(/\.pgen$/, '')
            def m = (base =~ /^chr(\d+)\.thin$/)
            if (!m.matches()) {
                error "Unexpected thin pfile name: ${pgen.name}"
            }
            def chr_num = m[0][1]
            tuple(
                chr_num,
                base,
                pgen,
                file("${pgen.parent}/${base}.psam"),
                file("${pgen.parent}/${base}.pvar"),
            )
        }

    smiss_out = plink2_sample_missing_per_chr(chr_pfiles)
    collect_sample_chr_miss(
        smiss_out.smiss.map { _chr, smiss -> smiss }.collect(),
        params.thin_job,
        params.thin_mod,
    )
}
