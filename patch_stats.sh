sed -i '/process variant_popdep_stats {/i \
process variant_popdep_mahalanobis {\n\
    tag "variant popdep mahalanobis: ${chr}"\n\
    publishDir "${params.output_dir}/${params.job}/stats/thresholds", mode: "copy", pattern: "*.info.tsv"\n\
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: "copy", pattern: "*.png"\n\
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: "copy", pattern: "*.log"\n\
    conda "${params.user_dir}/miniconda3/envs/stats"\n\
\n\
    input:\n\
    tuple val(id), val(chr), path(popdep)\n\
\n\
    output:\n\
    tuple val(id), val(chr), path("*.info.tsv"), emit: info\n\
    tuple val(id), val(chr), path("*.png"), emit: plots\n\
    tuple val(id), val(chr), path("*.log"), emit: logs\n\
\n\
    script:\n\
    """\n\
    #!/usr/bin/env python\n\
    import sys\n\
    sys.stdout = open("${chr}.popdep_mahalanobis.log", "w")\n\
    sys.stderr = sys.stdout \n\
\n\
    from genetics.genomics.variant.popdep import ana_popdep_mahalanobis\n\
\n\
    print(f"Processing population depth Mahalanobis distance for ${chr}...")\n\
    ana_popdep_mahalanobis("${popdep}", "${id}.variant.mahalanobis")\n\
    """\n\
}\n\
' /data/home/tusr1/git/script/workflow/Genetics/genotype/stats.nf
