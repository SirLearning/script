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

process variant_mq_missing_reg {
    tag "variant mq vs missing reg: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(mq), path(vmiss)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.mq_miss_reg.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.mq import ana_mq_missing_reg

    print(f"Processing MQ vs missing regression for ${chr}...")
    ana_mq_missing_reg("${mq}", "${vmiss}", "${id}.variant.mq_miss")
    """
}

process variant_popdep_qc_plots {
    tag "variant popdep qc plots: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(popdep), path(vmiss)

    output:
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.popdep_qc.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.popdep import ana_popdep_qc_plots

    print(f"Processing extended popdep QC plots for ${chr}...")
    ana_popdep_qc_plots("${popdep}", "${vmiss}", "${id}.variant.popdep_qc")
    """
}

process variant_popdep_missing_reg {
    tag "variant popdep vs missing reg: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(popdep), path(vmiss)

    output:
    // Popdep sidecars may have all-null depth (e.g. annotate miss); allow no PNGs.
    tuple val(id), val(chr), path("*.png"), optional: true, emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.popdep_miss_reg.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.popdep import ana_popdep_missing_reg

    print(f"Processing popdep vs missing regression for ${chr}...")
    ana_popdep_missing_reg("${popdep}", "${vmiss}", "${id}.variant.popdep_miss")
    """
}

process variant_popdep_mahalanobis {
    tag "variant popdep mahalanobis: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
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
    ana_popdep_mahalanobis("${popdep}", "${id}.variant.mahalanobis", relative_depth=False)
    ana_popdep_mahalanobis("${popdep}", "${id}.variant.mahalanobis_reldepth", relative_depth=True)
    """
}

process test_site_popdep_plots {
    tag "test_site popdep plots"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    val(process_dir)

    output:
    path "*.png", emit: plots
    path "*.tsv", emit: info
    path "*.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("test_site_popdep_plots.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.popdep_site_panel import ana_test_site_popdep_plots

    print("test_site popdep-only plots (variant sites highlighted)...")
    ana_test_site_popdep_plots(
        "${process_dir}",
        "test_site.popdep",
    )
    """
}

process popdep_loo_influence_summary {
    tag "popdep loo influence summary"
    publishDir "${publish_root}/stats", mode: 'copy', pattern: "*"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(n_cases), val(publish_root), val(baseline_grid), val(grid_id), path(manifest), val(summary_prefix)

    output:
    path "*.info.tsv", emit: info
    path "*.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${summary_prefix}.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.popdep import summarize_popdep_leave_one_out

    print(f"PopDep LOO influence summary: n_cases=${n_cases} publish_root=${publish_root}")
    print(f"Baseline grid=${baseline_grid} chr_id=${grid_id}")
    summarize_popdep_leave_one_out(
        "${baseline_grid}",
        "${publish_root}",
        "${summary_prefix}",
        chr_id="${grid_id}",
        min_baseline_mean=0.0,
    )
    """
}

process popdep_loo_aneuploidy_samechr_compare {
    tag "popdep loo aneuploidy samechr compare"
    publishDir "${params.popdep_loo_aneuploidy_compare_publish_root ?: "${params.output_dir}/${params.job}/process/_popdep_n500_loo_aneuploidy_samechr"}/stats", mode: 'copy', pattern: "*"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path manifest_tsv
    path chr1_info_tsv

    output:
    path "*.info.tsv", emit: info
    path "*.log", emit: logs

    script:
    def publish_base = params.popdep_loo_aneuploidy_compare_publish_root
        ?: "${params.output_dir}/${params.job}/process/_popdep_n500_loo_aneuploidy_samechr"
    def out_prefix = params.popdep_loo_aneuploidy_compare_prefix ?: 'n500_aneuploidy_samechr_vs_chr1'
    def bench_bp = params.popdep_loo_bench_max_bp ?: 10000000
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${out_prefix}.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.popdep import compare_aneuploidy_samechr_loo_influence

    compare_aneuploidy_samechr_loo_influence(
        "${manifest_tsv}",
        "${chr1_info_tsv}",
        "${publish_base}",
        "${out_prefix}",
        bench_max_bp=${bench_bp},
    )
    """
}

process popdep_loo_aneuploidy_samechr_plots {
    tag "popdep loo aneuploidy samechr plots"
    publishDir "${params.popdep_loo_aneuploidy_compare_publish_root ?: "${params.output_dir}/${params.job}/process/_popdep_n500_loo_aneuploidy_samechr"}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.popdep_loo_aneuploidy_compare_publish_root ?: "${params.output_dir}/${params.job}/process/_popdep_n500_loo_aneuploidy_samechr"}/stats/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.popdep_loo_aneuploidy_compare_publish_root ?: "${params.output_dir}/${params.job}/process/_popdep_n500_loo_aneuploidy_samechr"}/stats/logs", mode: 'copy', pattern: "*.log"
    publishDir "${params.popdep_loo_aneuploidy_compare_publish_root ?: "${params.output_dir}/${params.job}/process/_popdep_n500_loo_aneuploidy_samechr"}", mode: 'copy', pattern: "chr*/plots/*.png"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    val compare_info_tsv
    val publish_root
    path chr1_info_tsv
    path manifest_tsv

    output:
    path "*.png", emit: plots
    path "chr*/plots/*.png", emit: chr_plots
    path "*.info.tsv", emit: info
    path "*.log", emit: logs

    script:
    def out_prefix = params.popdep_loo_aneuploidy_compare_prefix ?: 'n500_aneuploidy_samechr_vs_chr1'
    def bench_bp = params.popdep_loo_bench_max_bp ?: 10000000
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${out_prefix}.aneuploidy_plots.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.popdep import plot_aneuploidy_samechr_loo_batch

    print("PopDep aneuploidy same-chr LOO plots (split absolute / pattern influence)")
    print(f"  compare=${compare_info_tsv}")
    print(f"  manifest=${manifest_tsv}")
    print(f"  publish_root=${publish_root}")
    print(f"  chr1_info=${chr1_info_tsv}")
    plot_aneuploidy_samechr_loo_batch(
        "${publish_root}",
        "${compare_info_tsv}",
        "${out_prefix}",
        bench_max_bp=${bench_bp},
        chr1_info_tsv="${chr1_info_tsv}",
        per_chr_plot_root=".",
        manifest_tsv="${manifest_tsv}",
    )
    """
}

process popdep_loo_influence_plots {
    tag "popdep loo influence plots"
    publishDir "${params.popdep_loo_publish_root ?: "${params.output_dir}/${params.job}/process/${params.popdep_loo_publish_slug}"}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.popdep_loo_publish_root ?: "${params.output_dir}/${params.job}/process/${params.popdep_loo_publish_slug}"}/stats/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.popdep_loo_publish_root ?: "${params.output_dir}/${params.job}/process/${params.popdep_loo_publish_slug}"}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path info_tsv

    output:
    path "*.png", emit: plots
    path "*.info.tsv", emit: info
    path "*.log", emit: logs

    script:
    def out_prefix = params.popdep_loo_summary_prefix
        ?: "n500_chr${params.popdep_loo_chr}_${params.popdep_loo_bench_max_bp}_loo"
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${out_prefix}.influence_plots.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.popdep import plot_popdep_loo_influence_summary

    print(f"PopDep LOO influence plots from ${info_tsv}")
    plot_popdep_loo_influence_summary(
        "${info_tsv}",
        "${out_prefix}",
        top_n=${params.popdep_loo_plot_top_n ?: 15},
    )
    """
}

process popdep_loo_outlier_review {
    tag "popdep loo outlier review"
    publishDir "${params.popdep_loo_publish_root ?: "${params.output_dir}/${params.job}/process/${params.popdep_loo_publish_slug}"}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.popdep_loo_publish_root ?: "${params.output_dir}/${params.job}/process/${params.popdep_loo_publish_slug}"}/stats/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.popdep_loo_publish_root ?: "${params.output_dir}/${params.job}/process/${params.popdep_loo_publish_slug}"}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple path(info_tsv), path(sample_qc_tsv)

    output:
    path "*.png", emit: plots
    path "*.info.tsv", emit: info
    path "*.log", emit: logs

    script:
    def out_prefix = params.popdep_loo_summary_prefix
        ?: "n500_chr${params.popdep_loo_chr}_${params.popdep_loo_bench_max_bp}_loo"
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${out_prefix}.outlier_review.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.popdep import plot_popdep_loo_outlier_review

    print(f"PopDep LOO outlier review: info=${info_tsv} qc=${sample_qc_tsv}")
    plot_popdep_loo_outlier_review(
        "${info_tsv}",
        "${sample_qc_tsv}",
        "${out_prefix}",
        top_n_highlight=${params.popdep_loo_review_top_n ?: 15},
        top_pct=${params.popdep_loo_review_top_pct ?: 5.0},
    )
    """
}

process popdep_loo_quadrant_review {
    tag "popdep loo quadrant review"
    publishDir "${params.popdep_loo_publish_root ?: "${params.output_dir}/${params.job}/process/${params.popdep_loo_publish_slug}"}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.popdep_loo_publish_root ?: "${params.output_dir}/${params.job}/process/${params.popdep_loo_publish_slug}"}/stats/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.popdep_loo_publish_root ?: "${params.output_dir}/${params.job}/process/${params.popdep_loo_publish_slug}"}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple path(info_tsv), path(sample_qc_tsv)

    output:
    path "*.png", emit: plots
    path "*.info.tsv", emit: info
    path "*.log", emit: logs

    script:
    def out_prefix = params.popdep_loo_summary_prefix
        ?: "n500_chr${params.popdep_loo_chr}_${params.popdep_loo_bench_max_bp}_loo"
    def inf_thr = params.popdep_loo_quadrant_influence_threshold
    def qual_thr = params.popdep_loo_quadrant_quality_threshold
    def inf_col = params.popdep_loo_quadrant_influence_col ?: 'influence_score'
    def manifest = params.popdep_loo_aneuploidy_manifest
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${out_prefix}.quadrant_review.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.popdep import plot_popdep_loo_influence_quality_quadrant

    print(f"PopDep LOO quadrant review: info=${info_tsv} qc=${sample_qc_tsv} col=${inf_col}")
    plot_popdep_loo_influence_quality_quadrant(
        "${info_tsv}",
        "${sample_qc_tsv}",
        "${out_prefix}",
        influence_col="${inf_col}",
        influence_threshold=${inf_thr != null ? inf_thr : 'None'},
        quality_threshold=${qual_thr != null ? qual_thr : 'None'},
        highlight_manifest_tsv=${manifest ? "\"${manifest}\"" : 'None'},
    )
    """
}

process popdep_loo_aneuploidy_samechr_pattern_quadrant {
    tag "popdep loo aneuploidy samechr pattern quadrant"
    publishDir "${params.popdep_loo_aneuploidy_compare_publish_root ?: "${params.output_dir}/${params.job}/process/_popdep_n500_loo_aneuploidy_samechr"}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.popdep_loo_aneuploidy_compare_publish_root ?: "${params.output_dir}/${params.job}/process/_popdep_n500_loo_aneuploidy_samechr"}/stats/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.popdep_loo_aneuploidy_compare_publish_root ?: "${params.output_dir}/${params.job}/process/_popdep_n500_loo_aneuploidy_samechr"}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path manifest_tsv
    val samechr_publish_root
    path sample_qc_tsv
    val pattern_threshold
    val quality_threshold

    output:
    path "*.png", emit: plots
    path "*.info.tsv", emit: info
    path "*.log", emit: logs

    script:
    def out_prefix = params.popdep_loo_aneuploidy_compare_prefix ?: 'n500_aneuploidy_samechr_vs_chr1'
    def bench_bp = params.popdep_loo_bench_max_bp ?: 10000000
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${out_prefix}.aneuploidy_samechr_pattern_quadrant.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.popdep import plot_aneuploidy_samechr_pattern_quality_quadrant

    plot_aneuploidy_samechr_pattern_quality_quadrant(
        "${manifest_tsv}",
        "${samechr_publish_root}",
        "${sample_qc_tsv}",
        "${out_prefix}",
        bench_max_bp=${bench_bp},
        influence_threshold=${pattern_threshold != 'null' ? pattern_threshold : 'None'},
        quality_threshold=${quality_threshold != 'null' ? quality_threshold : 'None'},
    )
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

process variant_maf_miss_group_stats {
    tag "variant maf/miss group plots: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_8'

    input:
    tuple val(id), val(chr), path(afreq), path(vmiss)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.variant_maf_miss_group.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.maf_miss_group import ana_variant_maf_miss_by_group

    pfile_prefix = "${params.output_dir}/${params.job}/process/${params.mod}/${id}_test.plink2"
    group_file = "${params.output_dir}/meta_data/sample_group.txt"
    print(f"Processing variant MAF/miss group plots for ${chr}...")
    ana_variant_maf_miss_by_group(
        "${afreq}",
        "${vmiss}",
        pfile_prefix,
        group_file,
        "${id}.variant.maf_miss",
    )
    """
}
