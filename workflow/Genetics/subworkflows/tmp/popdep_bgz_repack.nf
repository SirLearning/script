#!/usr/bin/env nextflow
/*
 * Re-pack existing TIGER *.popdep.txt.gz into BGZF + tabix (no TIGER recompute).
 *
 * Preserves the full TIGER header (including RelativeDepth_Mean / RelativeDepth_SD)
 * and prepends Chrom. Publishes to params.popdep_publish_dir/variant/.
 *
 * Example (N500 panel — CrossChr gz + PopDepFull chr32 gz):
 *   mkdir -p /data/home/tusr1/01projects/vmap4/10stats.genome/37run_popdep_n500_bgz_repack/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/37run_popdep_n500_bgz_repack
 *   nextflow run .../subworkflows/tmp/popdep_bgz_repack.nf -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --popdep_publish_dir /data1/.../test_plink/process/_popdep_n500_gam_credible \
 *     --popdep_gz_dirs '/path/crosschr_out,/path/chr032_gz_dir'
 */

nextflow.enable.dsl=2

include { popdep_tiger_gz_to_bgzip_tabix } from '../../modules/local/genotype/processor/processor_depth.nf'

params.popdep_gz_dirs = null
params.popdep_publish_dir = null

def popdepTupleFromGzFile(gz) {
    def name = gz.getName()
    if (!name.endsWith('.popdep.txt.gz')) {
        throw new IllegalArgumentException("popdep_bgz_repack: unexpected gz name: ${name}")
    }
    def base = name.substring(0, name.length() - '.popdep.txt.gz'.length())
    def chr = base.replaceFirst('(?i)^chr', '')
    if (!chr.isInteger()) {
        throw new IllegalArgumentException("popdep_bgz_repack: cannot parse chromosome from gz name: ${name}")
    }
    def id = String.format('chr%03d', chr.toInteger())
    return tuple(id, chr, gz)
}

workflow {
    if (!params.popdep_gz_dirs) {
        error 'popdep_bgz_repack: --popdep_gz_dirs is required (comma-separated directories of *.popdep.txt.gz)'
    }
    if (!params.popdep_publish_dir) {
        error 'popdep_bgz_repack: --popdep_publish_dir is required'
    }

    def dir_list = params.popdep_gz_dirs.split(',').collect { it.trim() }.findAll { it }

    ch_gz = channel
        .from(dir_list)
        .flatMap { dir_path ->
            def dir = file(dir_path, type: 'dir')
            if (!dir.exists()) {
                log.warn "popdep_bgz_repack: missing directory ${dir_path}"
                return []
            }
            dir.listFiles({ f -> f.name.endsWith('.popdep.txt.gz') })
                .collect { gz -> popdepTupleFromGzFile(gz) }
        }

    popdep_tiger_gz_to_bgzip_tabix(ch_gz)
}
