#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    ct_md5_wtk_fq
    ct_md5_wtk_bam
    ck_md5_wtk
    collect_md5_results
} from './align_md5.nf'

include {
    prepare_reference_genome
    align_with_bwa_mem2
    filter_bam_with_samtools
    filter_bam_with_samtools_no_rmdup
    prepare_bam_index_with_samtools
    ct_depth_with_mosdepth
} from './align_bwa.nf'

include { cp_bams_based_on_usb_size } from './align_transfer.nf'

def helpMessage() {
    log.info """
    Usage: nextflow run align.nf [options]

    Options:

    Example:
    screen -dmS cp_115 bash -c "\
        cd /data/home/tuser1/run && \
        source ~/.bashrc && conda activate run && \
        nextflow run /data/home/tuser1/git/script/workflow/Genetics/modules/local/genotype/align/align.nf \
            --user_dir /data/home/tuser1 \
            --output_dir /data/home/tuser1 \
            --bams_dir /data/home/tuser1/01bam \
            --usb_mnt_dir /mnt \
            --server s115 \
            --usb_dirs usb,usb-2,usb-3 \
            -resume "
    screen -dmS cp_204 bash -c "\
        cd /data/home/tusr1/01projects/vmap4/06transfer.watkins/01run && \
        source ~/.bashrc && conda activate run && \
        nextflow run /data/home/tusr1/git/script/workflow/Genetics/modules/local/genotype/align/align.nf \
            --user_dir /data/home/tusr1 \
            --output_dir /data/home/tusr1/01projects/vmap4/06transfer.watkins \
            --bams_dir /data1/dazheng_tusr1/204 \
            --usb_mnt_dir /mnt \
            --server s243 \
            --usb_dirs usb,usb-2,usb-3 \
            -resume "
    """
}

workflow RUN_ALIGN_USB_TRANSFER {
    main:
    files_ch = channel.fromPath(params.bams_dir + "/CRR*.rmdup.bam", type: 'file').map { bam ->
        def id = bam.baseName.replace(".rmdup.bam", "")
        def bai = file("${bam}.bai")
        def md5 = file("${bam}.md5")
        [params.server, id, bam, bai, md5]
    }

    files = files_ch.groupTuple(by: 0)

    usb_dirs_list = params.usb_dirs.split(',')
    dirs = channel.from(usb_dirs_list).map { dir ->
        def dir_path = file("${params.usb_mnt_dir}/${dir}")
        [params.server, dir, dir_path]
    }.groupTuple(by: 0)

    _cp_out = cp_bams_based_on_usb_size(files, dirs)

    md5_check_ch = _cp_out.copied_manifest.flatMap { server, manifest ->
        manifest.readLines().collect { line ->
            def cols = line.split('\t')
            if (cols.size() >= 3) {
                def bam_file = cols[0].trim()
                def bai_file = cols[1].trim()
                def md5_file = cols[2].trim()
                def bam_id = new File(bam_file).getName()
                return [server, bam_id, bam_file, bai_file, md5_file]
            } else if (cols.size() == 2) {
                def bam_file = cols[0].trim()
                def md5_file = cols[1].trim()
                def bam_id = new File(bam_file).getName()
                return [server, bam_id, bam_file, "MISSING_BAI", md5_file]
            }
            return null
        }.findAll { it != null }
    }

    individual_checks_ch = ck_md5_wtk(md5_check_ch)
    collect_md5_results(individual_checks_ch.groupTuple(by: 0))
}

workflow {
    RUN_ALIGN_USB_TRANSFER()
}
