version 1.0

workflow wgs {

    input {
        String sample_id
        String reads_group

        File fastq1
        File fastq2

        Array[File] genome_indexes
        File reference_fasta

        String docker_img = "registry-vpc.cn-hangzhou.aliyuncs.com/easy-gene/genomes-in-the-cloud:1.0"
    }

    call create_sequence_dict {
        input:
            reference_fasta = reference_fasta,
            docker_img = docker_img
    }

    call bwa_align {
        input: 
            sample_id=sample_id,
            reads_group=reads_group,
            fastq1=fastq1,
            fastq2=fastq2,
            genome_indexes=genome_indexes,
            docker_img=docker_img
    }

    call dedup {
        input:
            sample_id=sample_id,
            sorted_bam = bwa_align.sorted_bam_output,
            docker_img = docker_img
    }

    # call gen_bqsr_scatter {
    #     input:
    #         genome_dict=create_sequence_dict.genome_dict
    # }

    # scatter (sequence_group_intervals in gen_bqsr_scatter.sequence_grouping) {
    #     call base_recalibrator {
    #         input:
    #             grouping=sequence_group_intervals,
    #             sample_id=sample_id,
    #             dedup_bam=dedup.dedup_bam_output,
    #             genome_indexes=genome_indexes,
    #             docker_img=docker_img
    #     }
    # }

    # call gather_bqsr_report {
    #     input:
    #         recalibration_reports=base_recalibrator.recalibration_report_output,
    #         sample_id=sample_id,
    #         docker_img=docker_img
    # }

    # scatter (subgroup in gen_bqsr_scatter.sequence_grouping_with_unmapped) {
    #     call apply_bqsr {
    #         input:
    #             sample_id=sample_id,
    #             calling_region=subgroup,
    #             dedup_bam=dedup.dedup_bam_output,
    #             genome_indexes=genome_indexes,
    #             bqsr_grp=gather_bqsr_report.bqsr_grp_output,
    #             docker_img=docker_img
    #     }
    # }

    # call gather_bam_files {
    #     input:
    #         sample_id=sample_id,
    #         separate_recalibrated_bams = apply_bqsr.separate_recalibrated_bam,
    #         docker_img=docker_img
    # }

    call gen_hc_scatter {
        input:
            genome_dict=create_sequence_dict.genome_dict
    }

    scatter (intervals in gen_hc_scatter.scatter_interval_list) {
        call haplotype_caller {
            input:
                genome_dict=create_sequence_dict.genome_dict,
                sample_id=sample_id,
                # input_bam=gather_bam_files.recalibrated_bam_output,
                input_bam=dedup.dedup_bam_output,
                calling_region=intervals,
                reference_indices=bwa_align.reference_indices,
                docker_img=docker_img
        }
    }

    call merge_vcfs {
        input:
            sample_id=sample_id,
            docker_img=docker_img,
            vcfs=haplotype_caller.hc_block_gvcf
    }

    call genotype_gvcfs {
        input: 
            genome_dict=create_sequence_dict.genome_dict,
            sample_id=sample_id,
            reference_indices=bwa_align.reference_indices,
            gvcf=merge_vcfs.gvcf_output,
            docker_img=docker_img
    }

    output {
        # Pair[File,File] output_bam = gather_bam_files.recalibrated_bam_output
        Array[File] reference_indices = bwa_align.reference_indices  # 引用task的输出
        Pair[File,File] output_bam = dedup.dedup_bam_output
        Pair[File,File] output_gvcf = merge_vcfs.gvcf_output
        File output_vcf = genotype_gvcfs.vcf_output
    }

}


task create_sequence_dict {
    input {
        File reference_fasta
        String gatk_Launcher = "/usr/local/bin/gatk-4.1.4.1/gatk"
        String docker_img
    }

    String dict_filename = basename(reference_fasta, ".fa.gz") + ".dict"
    
    command <<<
        ~{gatk_Launcher} CreateSequenceDictionary \
        -R ~{reference_fasta} \
        -O ~{dict_filename}
    >>>

    runtime {
        docker: docker_img
        cpu: 16
        memory: "64G"
        disk: "250G"
    }

    output {
        File genome_dict = dict_filename
    }
}


task bwa_align {

    input {
        # reads pair
        File fastq1
        File fastq2

        # pre-built genome index files
        Array[File] genome_indexes

        # experiments info
        String reads_group
        String sample_id

        # bwa parameters
        String bwa_release = "/usr/local/bin/bwakit-0.7.15/"
        String samtools_release = "/usr/local/bin/samtools-1.9/"
        String bwa_opts = "-K 10000000 -M -Y"
    
        # Resource
        Int cpu = 32
        String memory = "128G"
        String disk = "1000G"

        # docker image
        String docker_img
    }

    String sorted_bam = sample_id + ".sorted.bam"
    String sorted_bam_index = sorted_bam + ".bai"

    # 预定义索引文件名
    String ref_basename = basename(genome_indexes[0])
    String amb_file = ref_basename + ".amb"
    String ann_file = ref_basename + ".ann"
    String bwt_file = ref_basename + ".bwt"
    String pac_file = ref_basename + ".pac"
    String sa_file = ref_basename + ".sa"
    String fai_file = ref_basename + ".fai"

    command <<<
        cpu_cores=$(nproc)

        # 找到参考基因组文件
        REF_GENOME=""
        for file in ~{sep=" " genome_indexes}; do
            if [[ "$file" == *.fa || "$file" == *.fasta || "$file" == *.fa.gz || "$file" == *.fasta.gz ]]; then
                REF_GENOME="$file"
                break
            fi
        done
        
        if [[ -z "$REF_GENOME" ]]; then
            REF_GENOME="~{genome_indexes[0]}"
        fi
        
        echo "Reference genome: $REF_GENOME"
        
        # 检查是否需要构建索引
        BASE_NAME=$(basename "$REF_GENOME")

        # 将参考基因组文件复制到当前工作目录
        echo "Copying reference genome to local directory..."
        cp "$REF_GENOME" "./$BASE_NAME"    # 修复：使用正确的变量名


        # # 将所有索引文件复制到当前工作目录
        # echo "Copying all index files to local directory..."
        # for index_file in ~{sep=" " genome_indexes}; do
        #     index_basename=$(basename "$index_file")
        #     cp "$index_file" "./$index_basename"
        #     echo "Copied: $index_basename"
        # done

        if [[ ! -f "${BASE_NAME}.bwt" ]]; then
            echo "Building BWA index..."
            ~{bwa_release}/bwa index "./$BASE_NAME"
        else
            echo "BWA index already exists"
        fi

        # 检查复制后的文件
        echo "Files in current directory:"
        ls -la

        # 验证所有BWA索引文件是否存在
        echo "Checking BWA index files:"
        for ext in .amb .ann .bwt .pac .sa; do
            index_file="${BASE_NAME}${ext}"
            if [[ -f "$index_file" ]]; then
                echo "✓ Found: $index_file"
            else
                echo "✗ Missing: $index_file"
                exit 1
            fi
        done

        # 构建正确的 read group 格式
        if [[ "~{reads_group}" == @RG* ]]; then
            # 如果已经是正确格式，直接使用
            RG_LINE="~{reads_group}"
        else
            # 如果不是，构建标准格式
            RG_LINE="@RG\tID:~{sample_id}\tSM:~{sample_id}\tPL:ILLUMINA\tLB:~{reads_group}"
        fi
        
        echo "Using read group: $RG_LINE"
        echo "Using reference: ./$BASE_NAME"

        ~{bwa_release}/bwa mem \
        ~{bwa_opts} \
        -t $cpu_cores \
        -R "$RG_LINE" \
        "./$BASE_NAME" \
        ~{fastq1} ~{fastq2} \
        | ~{samtools_release}/samtools sort -@ $cpu_cores -o ~{sorted_bam}

        ~{samtools_release}/samtools index ~{sorted_bam}
    >>>

    runtime {
        cpu: cpu
        memory: memory
        disk: disk
        docker: docker_img
    }

    output {
        Pair[File, File] sorted_bam_output = (sorted_bam, sorted_bam_index)
        Array[File] reference_indices = [ref_basename, amb_file, ann_file, bwt_file, pac_file, sa_file, fai_file]
    }
}

task dedup {
    input {
        String sample_id
        Pair[File, File] sorted_bam
        String gatk_Launcher = "/usr/local/bin/gatk-4.1.4.1/gatk"
        String java_opts = "'-Xms4000m -Djava.io.tmpdir=/mnt -Dsamjdk.compression_level=2'"

        Boolean remove_duplicates = true
        Boolean create_index = true

        # Resource
        Int cpu = 16
        String memory = "64G"
        String disk = "1000G"
        
        # docker 
        String docker_img
    }

    String dedup_bam = sample_id + ".deduplicated.bam"
    String dedup_bam_index = sample_id + ".deduplicated.bai"
    String dedup_metrics = sample_id + ".deduplicated.metrics"


    command <<<
        ~{gatk_Launcher} --java-options ~{java_opts} \
        MarkDuplicates \
        -I ~{sorted_bam.left} \
        -O ~{dedup_bam} \
        -M ~{dedup_metrics} \
        ~{true="--REMOVE_DUPLICATES true" false="" remove_duplicates} \
        ~{true="--CREATE_INDEX true" false="" create_index}
    >>>

    runtime {
        cpu: cpu
        memory: memory
        disk: disk
        docker: docker_img
    }

    output {
        Pair[File,File] dedup_bam_output = (dedup_bam, dedup_bam_index)
        File dedup_metrics_output = dedup_metrics
    }
}

# Split the whole genome for BQSR scatter-gather
task gen_bqsr_scatter {

    input {
        File genome_dict
        Int scatter_count = 20
    }

    command <<<
        # 从字典文件提取染色体信息
        grep "^@SQ" ~{genome_dict} | \
        awk '{split($2,a,":"); print "-L " a[2]}' > sequence_grouping.txt
        
        # 复制并添加unmapped
        cp sequence_grouping.txt sequence_grouping_with_unmapped.txt
        echo "-L unmapped" >> sequence_grouping_with_unmapped.txt
    >>>

    runtime {
        docker: "registry-vpc.cn-hangzhou.aliyuncs.com/easy-gene/ubuntu:24_04-python"
    }

    output {
        Array[String] sequence_grouping = read_lines("sequence_grouping.txt")
        Array[String] sequence_grouping_with_unmapped = read_lines("sequence_grouping_with_unmapped.txt")
    }
}



task base_recalibrator {
    input {

        String grouping
        
        Pair[File,File] dedup_bam

        Array[File] genome_indexes
        
        String sample_id
        String gatk_Launcher = "/usr/local/bin/gatk-4.1.4.1/gatk"
        String java_opts = "'-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m -Djava.io.tmpdir=/mnt'"

        # Resource
        Int cpu = 16
        String memory = "64G"
        String disk = "1000G"
        
        # docker 
        String docker_img
    }

    String recalibration_report = sample_id + "_bqsr.table"

    command <<<

        ~{gatk_Launcher} --java-options ~{java_opts} \
        BaseRecalibrator \
        -I ~{dedup_bam.left} \
        -R ~{genome_indexes[0]} \
        ~{grouping} \
        -O ~{recalibration_report}
    >>>

    runtime {
        cpu:cpu
        memory:memory
        disk:disk
        docker: docker_img
    }

    output {
        File recalibration_report_output = recalibration_report
    }

}

task gather_bqsr_report {
    input {
        Array[File] recalibration_reports
        String sample_id
        String gatk_Launcher = "/usr/local/bin/gatk-4.1.4.1/gatk"
        String java_opts = "'-Xms3000m -Djava.io.tmpdir=/mnt'"

        # Resource
        Int cpu = 16
        String memory = "64G"
        String disk = "1000G"
        
        # docker 
        String docker_img

    }

    String bqsr_grp = sample_id + "-bqsr.grp"

    command <<<
        ~{gatk_Launcher} --java-options ~{java_opts} \
        GatherBQSRReports \
        -I ~{sep=" -I " recalibration_reports} \
        -O ~{bqsr_grp}
    >>>

    runtime {
        cpu:cpu
        memory:memory
        disk:disk
        docker: docker_img
    }

    output {
        File bqsr_grp_output = bqsr_grp
    }
}

task apply_bqsr {
    input {
        Pair[File,File] dedup_bam
        Array[File] genome_indexes
        String calling_region
        File bqsr_grp
        String sample_id

        String gatk_Launcher = "/usr/local/bin/gatk-4.1.4.1/gatk"
        String java_opts = "'-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m -Djava.io.tmpdir=/mnt'"

        # Resource
        Int cpu = 32
        String memory = "128G"
        String disk = "1000G"
        
        # docker 
        String docker_img

    }

    String tmp_recalibrated_bam = sample_id + "_tmp_recaled.bam"

    command <<<
        ~{gatk_Launcher} --java-options ~{java_opts} \
        ApplyBQSR \
        -R ~{genome_indexes[0]} \
        -I ~{dedup_bam.left} \
        --bqsr-recal-file ~{bqsr_grp} \
        ~{calling_region} \
        -O ~{tmp_recalibrated_bam}
    >>>

    runtime {
        cpu:cpu
        memory:memory
        disk:disk
        docker: docker_img
    }

    output {
        File separate_recalibrated_bam = tmp_recalibrated_bam
    }
}

task gather_bam_files {

    input {
        Array[File] separate_recalibrated_bams
        String sample_id
        String gatk_Launcher="/usr/local/bin/gatk-4.1.4.1/gatk"
        String java_opts = "'-Xms3000m -Djava.io.tmpdir=/mnt -Dsamjdk.compression_level=2'"
        Boolean create_index = true

        # Resource
        Int cpu = 16
        String memory = "64G"
        String disk = "1000G"
        
        # docker 
        String docker_img
    }

    String recalibrated_bam = sample_id + ".recalibrated.bam"
    String recalibrated_bam_index = sample_id + ".recalibrated.bai"

    command <<<
        ~{gatk_Launcher} --java-options ~{java_opts} \
        GatherBamFiles \
        -I ~{sep=" -I " separate_recalibrated_bams} \
        ~{true="--CREATE_INDEX true " false = "" create_index} \
        -O ~{recalibrated_bam}
    >>>

    runtime {
        cpu:cpu
        memory:memory
        disk:disk
        docker: docker_img
    }

    output {
        Pair[File,File] recalibrated_bam_output = (recalibrated_bam, recalibrated_bam_index)
    }
}


# Break the calling interval_list into sub-intervals
# Perform variant calling on the sub-intervals, and then gather the results

task gen_hc_scatter_old {

    input {
        File wgs_calling_regions
    }

    command <<<

        awk 'BEGIN {
        prev_total = 0
        frag = 1
        container = ""
        } 
        { if ( $1 !~ /^@/ ) 
        {
            len = ($3 - $2 + 1)
            if ( prev_total + len < 324860607 ) {
            prev_total += len
            container = container sprintf("-L %s:%d-%d ", $1, $2, $3)
            }
            else {
            a1 = prev_total + len - 324860607
            a2 = 324860607 - prev_total
            if ( a1 > a2 ) { print container; container = sprintf("-L %s:%d-%d ", $1, $2, $3); prev_total = len}
            else { container = container sprintf("-L %s:%d-%d ", $1, $2, $3); print container; container = ""; prev_total = 0}
            frag += 1
            }
        }
        }
        END {
        if ( container ) { print container }
        }' ~{wgs_calling_regions} > "ScatterIntervalList.txt"

    >>>

    output {
        Array[String] scatter_interval_list = read_lines("ScatterIntervalList.txt")
    }
}

task gen_hc_scatter {

    input {
        File genome_dict
        Int scatter_count = 20
    }

    command <<<
        python - ~{genome_dict} ~{scatter_count} << EOF
import sys
import math

with open(sys.argv[1], "r") as ref_dict_file:
    sequence_list = []
    total_length = 0
    
    for line in ref_dict_file:
        if line.startswith("@SQ"):
            line_split = line.split("\t")
            seq_name = line_split[1].split("SN:")[1]
            seq_length = int(line_split[2].split("LN:")[1])
            sequence_list.append((seq_name, seq_length))
            total_length += seq_length

scatter_count = int(sys.argv[2])
target_size = total_length // scatter_count

with open("ScatterIntervalList.txt", "w") as output_file:
    current_size = 0
    current_intervals = []
    
    for seq_name, seq_length in sequence_list:
        if current_size + seq_length <= target_size or not current_intervals:
            current_intervals.append(f"-L {seq_name}")
            current_size += seq_length
        else:
            # 输出当前组合
            output_file.write(" ".join(current_intervals) + "\n")
            # 开始新的组合
            current_intervals = [f"-L {seq_name}"]
            current_size = seq_length
    
    # 输出最后一组
    if current_intervals:
        output_file.write(" ".join(current_intervals) + "\n")
EOF
    >>>

    runtime {
        docker: "registry-vpc.cn-hangzhou.aliyuncs.com/easy-gene/ubuntu:24_04-python"
    }

    output {
        Array[String] scatter_interval_list = read_lines("ScatterIntervalList.txt")
    }
}


# Call variants in parallel over WGS calling intervals

task haplotype_caller {
    input {
        File genome_dict

        String sample_id
        Pair[File,File] input_bam

        String calling_region 
        Array[File] reference_indices

        String gatk_Launcher = "/usr/local/bin/gatk-4.1.4.1/gatk"
        String java_opts = "'-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/mnt'"
        
        # Resource
        Int cpu = 32
        String memory = "128G"
        String disk = "1000G"
        
        # docker 
        String docker_img    


    }

    String hc_gvcf = sample_id + "_hc.g.vcf.gz"    

    command <<<
        # 复制参考基因组文件到本地
        REF_BASENAME=$(basename ~{reference_indices[0]})
        cp ~{reference_indices[0]} ./$REF_BASENAME
        # 复制参考基因组索引文件到本地
        for ext in .bwt .pac .sa; do
            index_file="${REF_BASENAME}${ext}"
            cp ~{reference_indices[0]}$ext ./$index_file
        done

        cp ~{genome_dict} ./$(basename ~{reference_indices[0]} .fa.gz).dict

        ~{gatk_Launcher} --java-options ~{java_opts} \
        HaplotypeCaller \
        -R ./$REF_BASENAME \
        -I ~{input_bam.left} \
        ~{calling_region} \
        --ERC GVCF \
        -O ~{hc_gvcf}
    >>>

    runtime {
        cpu:cpu
        memory:memory
        disk:disk
        docker: docker_img
    }

    output {
        File hc_block_gvcf = hc_gvcf
    }
}


task merge_vcfs {

    input {
        String sample_id
        Array[File] vcfs
        String gatk_Launcher = "/usr/local/bin/gatk-4.1.4.1/gatk"
        String java_opts = "'-Xms2000m -Djava.io.tmpdir=/mnt'"
        
        # Resource
        Int cpu = 16
        String memory = "64G"
        String disk = "1000G"
        
        # docker 
        String docker_img    
    }


    #test_hc.g.vcf.gz.tbi

    String gvcf = sample_id + ".g.vcf.gz"
    String gvcf_idx = sample_id + ".g.vcf.gz.tbi"

    command <<<
        ~{gatk_Launcher} --java-options ~{java_opts} \
        MergeVcfs \
        -I ~{sep=" -I " vcfs} \
        -O ~{gvcf}
    >>>

    runtime {
        cpu:cpu
        memory:memory
        disk:disk
        docker: docker_img
    }

    output {
        Pair[File, File] gvcf_output = (gvcf, gvcf_idx)
    }
}


task genotype_gvcfs {
    input {
        File genome_dict

        String sample_id
        Pair[File, File] gvcf
        
        Array[File] reference_indices

        String gatk_Launcher = "/usr/local/bin/gatk-4.1.4.1/gatk"
        String java_opts = "'-Xms4000m -Djava.io.tmpdir=/mnt'"
       
        # Resource
        Int cpu = 16
        String memory = "64G"
        String disk = "1000G"
        
        # docker 
        String docker_img
    }

    String vcf = sample_id + ".vcf.gz"

    command <<<
        # 复制参考基因组文件到本地
        REF_BASENAME=$(basename ~{reference_indices[0]})
        cp ~{reference_indices[0]} ./$REF_BASENAME
        # 复制参考基因组索引文件到本地
        for ext in .bwt .pac .sa; do
            index_file="${REF_BASENAME}${ext}"
            cp ~{reference_indices[0]}$ext ./$index_file
        done

        cp ~{genome_dict} ./$(basename ~{reference_indices[0]} .fa.gz).dict

        ~{gatk_Launcher} --java-options ~{java_opts} \
        GenotypeGVCFs \
        -R ~{reference_indices[0]} \
        -V ~{gvcf.left} \
        -O ~{vcf}
    >>>

    runtime {
        cpu:cpu
        memory:memory
        disk:disk
        docker: docker_img
    }

    output {
        File vcf_output = vcf
    }

}
