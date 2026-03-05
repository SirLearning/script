nextflow.enable.dsl=2

/*
    Module: static/llm.nf
    Description: Convert VCF to sequences using a JAR, then calculate Evo score using a model.
*/

workflow eval_evo_score {
    take:
    vcf_in // expect channel: [ val(id), val(chr), path(vcf) ]

    main:
    // 1. Convert VCF variant to sequence
    seq_out = convert_vcf_to_seq(vcf_in)

    // 2. Call evo model to calculate score
    score_out = calc_evo_score(seq_out.seq)

    emit:
    fasta = seq_out.seq
    score = score_out.score
}

// Ensure your Nextflow config (`nextflow.config`) has docker/singularity enabled.

process convert_vcf_to_seq {
    tag "vcf2seq ${id}"
    publishDir "${params.output_dir}/${params.job}/llm/fasta", mode: "copy", pattern: "*.fasta"
    publishDir "${params.output_dir}/${params.job}/llm/logs", mode: "copy", pattern: "*.log"
    label "cpus_4"
    // Use an environment or container config if needed

    input:
    tuple val(id), val(chr), path(vcf)

    output:
    tuple val(id), val(chr), path("${id}_${chr}.fasta"), emit: seq
    path "vcf_to_seq_${id}_${chr}.log", emit: log

    script:
    def vcf2seq_jar = params.vcf2seq_jar ?: "default.jar"
    """
    set -euo pipefail
    exec > vcf_to_seq_${id}_${chr}.log 2>&1

    echo "Converting VCF to sequences for ${id} ${chr}..."
    java -Xmx${task.memory.toGiga()}G -jar ${vcf2seq_jar} \\
        -i ${vcf} \\
        -c ${chr} \\
        -o ${id}_${chr}.fasta

    echo "Conversion completed."
    """
}

process calc_evo_score {
    tag "evo ${id}"
    // Container configuration: specify the Docker image containing the Evo model
    container { params.evo_docker ?: "evo-model:latest" }
    publishDir "${params.output_dir}/${params.job}/llm/score", mode: "copy", pattern: "*.score.txt"
    publishDir "${params.output_dir}/${params.job}/llm/logs", mode: "copy", pattern: "*.log"
    label "cpus_8"
    // Consider adding label "gpu" if computation requires GPU

    input:
    tuple val(id), val(chr), path(fasta)

    output:
    tuple val(id), val(chr), path("${id}_${chr}.score.txt"), emit: score
    path "calc_evo_score_${id}_${chr}.log", emit: log

    script:
    """
    set -euo pipefail
    exec > calc_evo_score_${id}_${chr}.log 2>&1

    echo "Calculating Evo scores for ${id} ${chr}..."
    
    # Run the Evo model inside the container. Modify the command based on your model's execution method.
    python /app/run_evo.py --input ${fasta} --output ${id}_${chr}.score.txt

    echo "Score calculation completed."
    """
}

