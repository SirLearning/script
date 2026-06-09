nextflow.enable.dsl=2

process cp_bams_based_on_usb_size {
    tag "cp_bams"
    publishDir "${params.output_dir}/logs", mode: 'copy'
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(server), val(bam_ids), path(bams), path(bai_files), path(md5_files)
    tuple val(server1), val(dir_names), val(dir_paths)

    output:
    path "usb_transfer.txt", emit: transfer_log
    tuple val(server1), path("copied_manifest.tsv"), emit: copied_manifest


    script:
    """
    #!/usr/bin/env python3
    import sys
    sys.stdout = open("cp_bams.log", "w")
    sys.stderr = sys.stdout # f2

    from infra.server.cp import run_copy_process

    bams_str = "${bams}"
    bams = bams_str.split()

    bais_str = "${bai_files}"
    bais = bais_str.split()

    md5s_str = "${md5_files}"
    md5s = md5s_str.split()
    
    # Group files (bam, bai, md5) together
    files = []
    # groupTuple preserves alignment of list elements
    for i in range(len(bams)):
        group = [bams[i], bais[i], md5s[i]]
        files.append(group)

    usb_dirs_str = "${dir_paths}"
    # Remove brackets if present
    usb_dirs_str = usb_dirs_str.replace('[', '').replace(']', '').replace(',', ' ')
    usb_dirs_list = usb_dirs_str.split()
    
    # Append server subdir to paths
    import os
    server_name = "${server1}"
    usb_dirs = []
    for d in usb_dirs_list:
        path = os.path.join(d, server_name)
        # Create subdir if not exists (handled in cp.py sort of, but good to be explicit or leave to cp.py logic)
        usb_dirs.append(path)

    transfer_file = "usb_transfer.txt"
    copied_manifest_file = "copied_manifest.tsv"
    threads = int("${task.cpus}")

    print(f"Running copy process for bams with {threads} threads...")
    run_copy_process(files, usb_dirs, transfer_file, threads, copied_manifest_file)
    """
}
