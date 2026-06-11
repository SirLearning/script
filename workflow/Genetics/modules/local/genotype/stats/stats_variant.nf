nextflow.enable.dsl=2

process variant_missing_stats {
    tag "variant missing stats: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(vmiss)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.vmiss_ana.log", "w")
    sys.stderr = sys.stdout # f2

    from genetics.genomics.variant import ana_variant_missing
    print(f"Processing variant missing stats for ${chr}...")
    ana_variant_missing("${vmiss}", "${id}.variant.miss")
    """
}

process variant_maf_stats {
    tag "variant maf stats: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(afreq)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.maf_ana.log", "w")
    sys.stderr = sys.stdout # f2

    from genetics.genomics.variant import ana_minor_allele_frequency

    print(f"Processing variant MAF stats for ${chr}...")
    ana_minor_allele_frequency("${afreq}", "${id}.variant.maf")

    """
}

process variant_mac_stats {
    tag "variant mac stats: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(gcount)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.mac_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.mac import ana_mac_stats

    print(f"Processing variant MAC stats for ${chr}...")
    ana_mac_stats("${gcount}", "${id}.variant.mac")

    """
}

process variant_mac_maf_reg {
    tag "variant mac vs maf reg: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(gcount)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.mac_maf_reg.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.mac import ana_mac_maf_reg

    print(f"Processing MAC vs MAF regression for ${chr}...")
    ana_mac_maf_reg("${gcount}", "${id}.variant.mac_maf")
    """
}

process variant_mac_missing_reg {
    tag "variant mac vs missing reg: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(gcount), path(vmiss)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.mac_miss_reg.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.mac import ana_mac_missing_reg

    print(f"Processing MAC vs missing regression for ${chr}...")
    ana_mac_missing_reg("${gcount}", "${vmiss}", "${id}.variant.mac_miss")
    """
}

process variant_mac_missing_reg_bin50_sample {
    tag "variant mac_miss bin50 sample: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.bin50sample.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.reg.bin50s.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.reg.mac_an.bin50s.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.R_mac.bin50s.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(gcount), path(vmiss)

    output:
    tuple val(id), val(chr), path("*.reg.bin50s.png"), emit: plots_reg
    tuple val(id), val(chr), path("*.reg.mac_an.bin50s.png"), emit: plots_mac_an_bin50s
    tuple val(id), val(chr), path("*.R_mac.bin50s.png"), emit: plots_R_mac
    tuple val(id), val(chr), path("*.bin50sample.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.mac_miss_bin50s.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.mac import ana_mac_missing_reg_bin50_sample

    print(f"MAC vs missing bin50-sample plots for ${chr}...")
    ana_mac_missing_reg_bin50_sample("${gcount}", "${vmiss}", "${id}.variant.mac_miss")
    """
}

process variant_mac_dist_log_redraw {
    tag "variant mac dist log: ${id} (${plink_mod})"
    publishDir "${params.output_dir}/${params.job}/stats/${plink_mod}/plots", mode: 'copy', pattern: "*.log.png"
    publishDir "${params.output_dir}/${params.job}/stats/${plink_mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(plink_mod), val(id), val(chr), path(mac_info)

    output:
    tuple val(plink_mod), val(id), val(chr), path("*.log.png"), emit: plots
    tuple val(plink_mod), val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.mac_dist_log.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.mac import redraw_mac_dist_log_from_info

    print(f"Redraw MAC log distribution for ${id} (${plink_mod})...")
    redraw_mac_dist_log_from_info("${mac_info}", "${id}.variant.mac")
    """
}

process variant_popdep_mahalanobis {
    tag "variant popdep mahalanobis: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/thresholds", mode: "copy", pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: "copy", pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: "copy", pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(popdep)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.popdep_mahalanobis.log", "w")
    sys.stderr = sys.stdout 

    from genetics.genomics.variant.popdep import ana_popdep_mahalanobis

    print(f"Processing population depth Mahalanobis distance for ${chr}...")
    ana_popdep_mahalanobis("${popdep}", "${id}.variant.mahalanobis")
    """
}

process variant_ld_decay_plot {
    tag "variant ld decay: ${id}"
    label 'cpus_16'
    cache false
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(ld_vcor)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.ld_decay_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant import ana_ld_decay

    print("Running LD decay plotting for ${id} (${chr}) ...")
    ana_ld_decay(
        input_file="${ld_vcor}",
        output_prefix="${id}.variant.ld_decay",
        max_distance_kb=${params.ld_decay_max_kb},
        bin_size_kb=${params.ld_decay_bin_kb},
        scatter_max_points=${params.ld_decay_scatter_points},
        chunksize=${params.ld_decay_chunksize},
        median_sample_per_bin=${params.ld_decay_median_sample_per_bin},
        workers=${task.cpus},
        random_seed=${params.ld_decay_random_seed},
    )
    """
}

process variant_ld_crosschr_plot {
    tag "variant cross-chr ld baseline: ${id}"
    label 'cpus_16'
    cache false
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(ld_cross_vcor)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.ld_crosschr_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant import ana_ld_crosschr_baseline

    print("Running cross-chromosome LD baseline plotting for ${id} (${chr}) ...")
    ana_ld_crosschr_baseline(
        input_file="${ld_cross_vcor}",
        output_prefix="${id}.variant.ld_crosschr",
        chunksize=${params.ld_decay_chunksize},
    )
    """
}


process variant_popdep_stats {
    tag "variant popdep stats: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(popdep)
    tuple val(id2), val(chr2), path(vmiss)

    output:
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.popdep_ana.log", "w")
    sys.stderr = sys.stdout 

    from genetics.genomics.variant.popdep import ana_popdep_curve_fit, ana_popdep_missing_reg

    print(f"Processing population depth Analysis 1 (Curve Fit) for ${chr}...")
    ana_popdep_curve_fit("${popdep}", "${id}.variant.popdep")

    print(f"Processing population depth Analysis 2 (vs Missing) for ${chr}...")
    ana_popdep_missing_reg("${popdep}", "${vmiss}", "${id}.variant.popdep")
    """
}
