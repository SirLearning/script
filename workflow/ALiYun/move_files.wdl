version 1.0

workflow move_files {
    input {
        # 使用OSS bucket挂载数据
        String bucket = "填入您当前账号下的需要挂载的bucket名字"
        Array[String] input_files = ["如/easygene_data/reference/hg38.fa", "/easygene_data/reference/hg38.fa.fai", "/easygene_data/reference/other_file.txt"]
        String target_directory = "reference/"

        String docker = "registry.cn-hangzhou.aliyuncs.com/damo-gene/dna-variant-calling-cpu@sha256:64c131e1c08fed5c7399ddef7b9aeb87eea0eaacaa8f5dde0ab93019b5a7230e"
    }

    scatter (single_file in input_files) {
        call CopyFiles {
            input:
                bucket = bucket,
                input_file = single_file,  # 传入单个文件
                target_directory = target_directory,
                docker = docker
        }
    }

    output {
    }
}

task CopyFiles {
    input {
        String bucket
        String input_file
        String target_directory
        
        Int cpu = 32
        String memory = "128G"
        String disk = "100G"
        String docker
    }

    
    command <<<
        set -euxo pipefail
        
        echo "Processing file: ~{input_file}"
        
        if [[ -f "~{input_file}" ]]; then
            echo "Copying file: ~{input_file} to ~{target_directory}"
            cp "~{input_file}" "~{target_directory}"
            echo "✓ Successfully copied: ~{input_file}"
        else
            echo "✗ File not found: ~{input_file}"
            exit 1
        fi
        
        echo "File copied to: ~{target_directory}/$(basename ~{input_file})"
        ls -la "~{target_directory}"
    >>>

    runtime {
        cpu: cpu
        memory: memory
        disk: disk
        docker: docker
        env: {
            "BUCKET": bucket
        }
    }
}
