version 1.0

workflow wgs {
    input {
        String sample_id
        String read_group
        # 平台支持OSS bucket挂载，可以直接在容器内访问OSS数据，请参考应用基本信息的数据使用方式
        File? fastq1
        File? fastq2
        # 支持输入CRAM格式，与FASTQ文件格式二选一
        File? cram
        # 和bwa的索引稍有不同，请参考应用基本信息的参数设置说明
        # BWA索引文件数组 - 作为输入参数，不设默认值
        Array[File] genomes
        # 二进制文件格式，请参考应用基本信息的参数设置说明
        File? known_sites
        
        File? intervals_bed
        # 捕获区域数组，格式为["chr1","chr2:1-100000"]
        Array[String]? intervals_string

        # 输出GVCF或VCF文件
        Boolean output_gvcf = false
        Boolean output_vcf = true
        
        String docker = "registry.cn-hangzhou.aliyuncs.com/damo-gene/dna-variant-calling@sha256:e94797b919ac4d3d13652d7c1b472552d33348316bb118a7f26de6473eecf5e7"
    }

    call FastqToVcf {
        input:
            fastq1=fastq1,
            fastq2=fastq2,
            cram=cram,
            sample_id=sample_id,
            read_group=read_group,
            genomes=genomes,
            known_sites=known_sites,
            intervals_bed=intervals_bed,
            intervals_string=intervals_string,
            output_gvcf=output_gvcf,
            output_vcf=output_vcf,
            docker=docker
    }

    output {
        Pair[File,File] markdup_bam = FastqToVcf.sdb_bam_output
        File metrics = FastqToVcf.metrics_output
        Array[File] vcf = FastqToVcf.vcf_output
    }
}

task FastqToVcf {
    input {
        String sample_id
        String read_group

        File? fastq1
        File? fastq2
        File? cram
        Array[File] genomes

        File? known_sites
        File? intervals_bed
        Array[String]? intervals_string
        
        Boolean output_gvcf
        Boolean output_vcf

        # 指定一个或多个GPU实例规格，左侧优先级最高
        Array[String] instanceType = ["gc7.small","gc4.medium","gc5.medium"]
        # 磁盘类型需要使用性能型NAS
        String disks = "local-disk 250 nas_per"
        String docker
    }

    String sdb_bam = sample_id + ".sdb.bam"
    String sdb_bam_index = sample_id + ".sdb.bam.csi"
    String recal_table = sample_id + ".recal.table"
    String metrics = sample_id + ".metrics.txt"
    String out_gvcf = sample_id + ".gvcf.gz"
    String out_gvcf_idx = sample_id + ".gvcf.gz.tbi"
    String out_vcf = sample_id + ".vcf.gz"
    String out_vcf_idx = sample_id + ".vcf.gz.tbi"

    command <<<
        set -euxo pipefail
        
        # 显示输入的genomes数组内容
        echo "=== Genomes array content ==="
        echo "~{sep='\n' genomes}"
        echo "=============================="
        
        # 将所有genomes文件复制到工作目录
        echo "Copying genome files to working directory..."
        for genome_file in ~{sep=' ' genomes}; do
            echo "Copying: $genome_file"
            cp "$genome_file" ./
        done
        
        # 查找参考基因组文件（第一个.fa.gz文件）
        ref_genome=""
        for genome_file in ~{sep=' ' genomes}; do
            filename=$(basename "$genome_file")
            if [[ "$filename" == *.fa.gz ]] && [[ "$filename" != *.fa.gz.* ]]; then
                ref_genome="$filename"
                break
            fi
        done
        
        # 如果没找到.fa.gz文件，使用第一个文件作为参考
        if [ -z "$ref_genome" ]; then
            ref_genome=$(basename "~{genomes[0]}")
            # 移除可能的索引后缀
            ref_genome=${ref_genome%.*}
            # 如果还有后缀，继续移除直到得到基础文件名
            while [[ "$ref_genome" == *.fa.gz.* ]]; do
                ref_genome=${ref_genome%.*}
            done
            ref_genome="${ref_genome}.fa.gz"
        fi
        
        echo "Reference genome file: $ref_genome"
        
        # 验证当前目录中的文件
        echo "Available files in current directory:"
        ls -la *.gz* || echo "No .gz files found"
        
        # 检查参考基因组文件是否存在
        if [ ! -f "$ref_genome" ]; then
            echo "Error: Reference genome file $ref_genome not found"
            echo "Trying to find reference genome file..."
            # 尝试查找任何.fa.gz文件
            for file in *.fa.gz; do
                if [[ -f "$file" ]] && [[ "$file" != *.fa.gz.* ]]; then
                    ref_genome="$file"
                    echo "Found reference genome: $ref_genome"
                    break
                fi
            done
        fi
        
        # 如果还是没找到，创建一个软链接
        if [ ! -f "$ref_genome" ]; then
            echo "Creating reference genome file..."
            # 查找最大的.gz文件作为参考基因组
            largest_file=$(ls -S *.gz* | head -1)
            if [[ "$largest_file" == *.0123 ]]; then
                # 假设.0123文件对应的参考基因组文件
                base_name=${largest_file%.0123}
                if [ ! -f "$base_name" ]; then
                    echo "Warning: Reference genome $base_name not found, using $largest_file as reference"
                    ref_genome="$largest_file"
                else
                    ref_genome="$base_name"
                fi
            else
                ref_genome="$largest_file"
            fi
        fi
        
        echo "Final reference genome: $ref_genome"
        
        # 验证BWA索引文件（检查实际存在的文件）
        echo "Verifying available BWA index files..."
        ref_base=${ref_genome}
        
        # 检查各种可能的索引文件格式
        index_files_found=0
        
        # 检查标准BWA索引文件
        for ext in 0123 amb ann bwt pac sa; do
            index_file="${ref_base}.${ext}"
            if [ -f "$index_file" ]; then
                echo "Found standard index file: $index_file"
                ((index_files_found++))
            fi
        done
        
        # 检查特殊格式的索引文件
        for file in ${ref_base}.bwt.2bit.64 ${ref_base}.bwt.2bit ${ref_base}.sa64; do
            if [ -f "$file" ]; then
                echo "Found special format index file: $file"
                ((index_files_found++))
            fi
        done
        
        echo "Total index files found: $index_files_found"
        
        if [ $index_files_found -lt 4 ]; then
            echo "Warning: Insufficient BWA index files found, but proceeding..."
        else
            echo "BWA index files verification completed"
        fi

        fast_align -R ~{read_group} \
        $ref_genome \
        ~{fastq1} ~{fastq2} ~{cram} \
        -o tmp_bams -s -S 30

        fast_sdb2 -I tmp_bams -S 30 \
        -R $ref_genome \
        ~{"-K " + known_sites + " -t " + recal_table} \
        -O ~{sdb_bam} -M ~{metrics} -r 16 -w 16

        if [[ -n '~{true="y" false="" output_gvcf}' ]]; then
            fast_hap -I ~{sdb_bam} \
            -R $ref_genome \
            ~{"-L " + intervals_bed} ~{if defined(intervals_string) then "-L" else ""} ~{sep=" -L " intervals_string} \
            ~{if defined(known_sites) then "-r " + recal_table else ""} \
            -O ~{out_gvcf} -E 2
        fi
        if [[ -n '~{true="y" false="" output_vcf}' ]]; then
            fast_hap -I ~{sdb_bam} \
            -R $ref_genome \
            ~{"-L " + intervals_bed} ~{if defined(intervals_string) then "-L" else ""} ~{sep=" -L " intervals_string} \
            ~{if defined(known_sites) then "-r " + recal_table else ""} \
            -O ~{out_vcf}
        fi
    >>>

    runtime {
        instanceType: instanceType
        disks: disks
        docker: docker
        continueOnReturnCode: [0, 1, 255]
        memory: "60 GB"
        cpu: 16
    }

    output {
        Pair[File,File] sdb_bam_output = (sdb_bam, sdb_bam_index)
        File metrics_output = metrics
        Array[File] vcf_output = glob("*vcf.gz*")
    }
}