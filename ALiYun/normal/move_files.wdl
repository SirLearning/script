version 1.0

workflow move_files {
    input {
        Array[File] input_files
        String output_dir = "oss://your-bucket/reference"
        String docker = "registry.cn-hangzhou.aliyuncs.com/damo-gene/dna-variant-calling@sha256:e94797b919ac4d3d13652d7c1b472552d33348316bb118a7f26de6473eecf5e7"
    }

    call move_to_directory {
        input:
            files = input_files,
            target_directory = output_dir,
            docker = docker
    }

    output {
        Array[File] moved_files = move_to_directory.moved_files
    }
}

task move_to_directory {
    input {
        Array[File] files
        String target_directory

        Int cpu = 32
        String memory = "128G"
        String disk = "100G"
        String docker
    }

    command <<<
        # Check if output directory is OSS path
        if [[ ${target_directory} == oss://* ]]; then
            # For OSS output, copy files to local temp directory first
            mkdir -p temp_output
            for file in ~{sep=' ' files}; do
                cp $file temp_output/
            done
            # Then upload to OSS using ossutil
            ossutil cp -r temp_output/ ${target_directory}/
        else
            # For local output
            mkdir -p ${target_directory}
            for file in ~{sep=' ' files}; do
                cp $file ${target_directory}/
            done
        fi
    >>>

    output {
        Array[File] moved_files = if (sub(target_directory, "oss://.*", "oss") == "oss")
            then glob("temp_output/*")
            else glob("${target_directory}/*")
    }

    runtime {
        cpu: cpu
        memory: memory
        disk: disk
        docker: docker
        continueOnReturnCode: 0
    }
}

