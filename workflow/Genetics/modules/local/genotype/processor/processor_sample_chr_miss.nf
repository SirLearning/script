nextflow.enable.dsl=2

/*
 * Per-segment PLINK2 sample missing rate (.smiss) from chrNNN.thin pfiles.
 * Segment id N maps to ref v1 names via get_ref_v1_chr_name (see ref_v1.py).
 */

process plink2_sample_missing_per_chr {
    tag "smiss_chr${chr_num}"

    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/sample", mode: 'copy', pattern: "*.smiss"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(chr_num), val(pfile_base), path(pgen), path(psam), path(pvar)

    output:
    tuple val(chr_num), path("chr${chr_num}.smiss"), emit: smiss
    path "chr${chr_num}.smiss.log", emit: log

    script:
    """
    set -euo pipefail
    exec > chr${chr_num}.smiss.log 2>&1

    echo "Computing sample missing rate for PLINK segment chr${chr_num} (${pfile_base})..."
    plink2 --pfile ${pfile_base} \\
        --allow-extra-chr \\
        --missing \\
        --threads ${task.cpus} \\
        --out chr${chr_num}
    echo "Wrote chr${chr_num}.smiss"
    """
}

process collect_sample_chr_miss {
    tag "collect sample_chr_miss tables"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path(smiss_files, stageAs: 'smiss/*')
    val(thin_job)
    val(thin_mod)

    output:
    path("sample_chr_miss.segments.tsv"), emit: segments
    path("sample_chr_miss.long.tsv"), emit: long_table
    path("collect_sample_chr_miss.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    from pathlib import Path

    sys.stdout = open("collect_sample_chr_miss.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample.sample_chr_miss import collect_segment_smiss_tables

    paths = collect_segment_smiss_tables(
        smiss_dir="smiss",
        thin_job="${thin_job}",
        thin_mod="${thin_mod}",
    )
    print("Wrote {} segment rows; long table samples={}".format(
        paths["segments"].stat().st_size,
        paths["long_table"].stat().st_size,
    ))
    """
}
