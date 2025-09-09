version 1.0
# WDL is meant to be a human readable and writable way to express tasks and workflows. 
# https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#introduction


# 在此定义您的分析流程
# https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#workflow-definition

workflow decompress_tar_gz {
    
    # 流程输入参数
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#types
    
    input {
        File compressed_file = "oss://gstor-default-workspace-cn-hangzhou-507d42e2/abd_bwa2Lib.tar.gz"
        String output_dir = "abd_bwa2Lib"
    }
    
    # 调用任务或者子流程, 提供对应的输入参数
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#call-statement
    call extract_tar_gz {
        input:
            input_file = compressed_file,
            output_directory = output_dir
    }
    
    # 流程分析最终结果
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#outputs
    output {
        Array[File] extracted_files = extract_tar_gz.result_files
        String extraction_info = extract_tar_gz.info
    }
}

# 在此定义分析流程中的任务
# https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#task-definition

task extract_tar_gz {

    # 必选，任务输入参数，支持多种类型变量声明。
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#types

    input {
        File input_file
        String output_directory
        # Int i = 0
        # Float f = 27.3
        # Boolean b = true
        # File f = "oss://bucket/to/file"
        # Array[File] = ["oss://bucket/to/file1","oss://bucket/to/file2"]
    }
    
    # 可选，无需输入的变量声明
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#non-input-declarations
    String info_file = "extraction_info.txt"


    # 任务命令行，为工具的实际执行脚本
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#command-section

    command <<<
        # 使用声明的输入文件和参数，组成实际执行的命令行
        echo "Starting extraction of ~{input_file}..." > ~{info_file}
        echo "Output directory: ~{output_directory}" >> ~{info_file}
        echo "Start time: $(date)" >> ~{info_file}
        
        # 创建输出目录
        mkdir -p ~{output_directory}
        
        # 解压tar.gz文件
        tar -xzf ~{input_file} -C ~{output_directory}
        
        # 检查解压是否成功
        if [ -d "~{output_directory}" ]; then
            echo "Extraction completed successfully!" >> ~{info_file}
            echo "Extracted files:" >> ~{info_file}
            ls -la ~{output_directory} >> ~{info_file}
            echo "Total files: $(find ~{output_directory} -type f | wc -l)" >> ~{info_file}
        else
            echo "Extraction failed!" >> ~{info_file}
            exit 1
        fi
        
        echo "End time: $(date)" >> ~{info_file}
    >>>
    
    # 可选，默认提供1核2G，40G硬盘的执行环境
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#runtime-section
    # https://help.aliyun.com/document_detail/257724.html
    runtime {
        cpu: 2
        memory: "4G"
        disks: "local-disk 100 SSD"
        # docker: "registry-vpc.cn-beijing.aliyuncs.com/easygene/genomes-in-the-cloud:1.0"
        # software: "sentieon:202010.02"
        continueOnReturnCode: [0, 1]
    }

    # 任务运行输出结果
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#outputs-section
    output {
        # 将工作目录下生产的结果文件上传到OSS中
        Array[File] result_files = glob("~{output_directory}/**")
        String info = read_string(info_file)
    }
    
}
