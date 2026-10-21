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

process test_popdep_class_plots {
    tag "test_popdep class plots"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_16'

    input:
    val sites_ready

    output:
    path "*.png", emit: plots
    path "*.tsv", emit: info
    path "*.log", emit: logs

    script:
    def sites_dir = params.test_popdep_sites_dir ?: "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_popdep'}/sites"
    """
    #!/usr/bin/env python
    import sys
    from pathlib import Path

    sys.stdout = open("test_popdep_class_plots.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.popdep_site_panel import ana_test_popdep_class_plots

    scatter_tsv = "${params.test_popdep_scatter_tsv ?: ''}"
    if scatter_tsv and scatter_tsv not in ("null", "None"):
        print(f"Redrawing JointGrids from scatter table {scatter_tsv} (joints_only, no 100M stream)...")
        ana_test_popdep_class_plots(scatter_tsv, "test_popdep.popdep", joints_only=True)
    else:
        tables = sorted(Path("${sites_dir}").glob("*.popdep.sites.tsv"))
        print(f"Streaming {len(tables)} chromosome popdep panels from ${sites_dir} (pre-binned hist, no 100M concat)...")
        ana_test_popdep_class_plots(tables, "test_popdep.popdep")
    """
}

process test_popdep_density_keep_plots {
    tag "test_popdep density keep vs drop"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_16'

    input:
    tuple path(scatter), path(masks)

    output:
    path "*.png", emit: plots
    path "*.tsv", emit: info
    path "*.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("test_popdep_density_keep_plots.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.popdep_density_filter import ana_test_popdep_density_keep_plots

    print("Colouring test_popdep scatter by genome 60%/80% trusted cells ...", flush=True)
    ana_test_popdep_density_keep_plots(
        "${scatter}",
        "${masks}",
        "test_popdep.popdep.density",
    )
    """
}

process test_popdep_log_density {
    tag "test_popdep log10 RelDepth density"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/info", mode: 'copy', pattern: "*.{tsv,npz}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_16'

    input:
    val sites_ready

    output:
    path "*.png", emit: plots
    path "*.tsv", emit: info
    path "*.npz", emit: npz
    path "*.log", emit: logs

    script:
    def sites_dir = params.test_popdep_sites_dir ?: "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_popdep'}/sites"
    def scatter = params.test_popdep_scatter_tsv ?: "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/info/test_popdep.popdep.scatter.sites.tsv"
    def n_cells = params.popdep_log_n_cells ?: 90000
    def counts_npz = params.popdep_log_counts_npz ?: ''
    def plot_keep_drop = params.popdep_log_plot_keep_drop
    if (plot_keep_drop == null) {
        plot_keep_drop = true
    }
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("test_popdep_log_density.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.popdep_density_filter import ana_test_popdep_log_density

    counts_npz = "${counts_npz}"
    kwargs = dict(n_cells=int("${n_cells}"), scatter="${scatter}")
    if counts_npz and counts_npz not in ("null", "None", ""):
        kwargs["counts_npz"] = counts_npz
        print(f"Reuse log histogram {counts_npz}; site-coverage 60%/80% ...", flush=True)
    else:
        print("log10 RelDepth density on test_popdep sites, n_cells=${n_cells} ...", flush=True)
    kwargs["plot_keep_drop"] = str("${plot_keep_drop}").strip().lower() in ("true", "1", "yes")
    if not kwargs["plot_keep_drop"]:
        print("Skip keep/drop JointGrids; write site-coverage curve.", flush=True)
    ana_test_popdep_log_density(
        "${sites_dir}",
        "test_popdep.popdep.log_density",
        **kwargs,
    )
    """
}

process test_popdep_linear_density {
    tag "test_popdep linear RelDepth density"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/info", mode: 'copy', pattern: "*.{tsv,npz}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_16'

    input:
    val sites_ready

    output:
    path "*.png", emit: plots
    path "*.tsv", emit: info
    path "*.npz", emit: npz
    path "*.log", emit: logs

    script:
    def sites_dir = params.test_popdep_sites_dir ?: "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_popdep'}/sites"
    def scatter = params.test_popdep_scatter_tsv ?: "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/info/test_popdep.popdep.scatter.sites.tsv"
    def n_cells = params.popdep_lin_n_cells ?: 90000
    def counts_npz = params.popdep_lin_counts_npz ?: ''
    def plot_keep_drop = params.popdep_lin_plot_keep_drop
    if (plot_keep_drop == null) {
        plot_keep_drop = true
    }
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("test_popdep_linear_density.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.popdep_density_filter import ana_test_popdep_linear_density

    counts_npz = "${counts_npz}"
    kwargs = dict(n_cells=int("${n_cells}"), scatter="${scatter}")
    if counts_npz and counts_npz not in ("null", "None", ""):
        kwargs["counts_npz"] = counts_npz
        print(f"Reuse linear histogram {counts_npz}; site-coverage 60%/80%/90% ...", flush=True)
    else:
        print("linear RelDepth density on test_popdep sites, n_cells=${n_cells} range=[0.02, 6] ...", flush=True)
    kwargs["plot_keep_drop"] = str("${plot_keep_drop}").strip().lower() in ("true", "1", "yes")
    if not kwargs["plot_keep_drop"]:
        print("Skip keep/drop JointGrids; write site-coverage curve.", flush=True)
    ana_test_popdep_linear_density(
        "${sites_dir}",
        "test_popdep.popdep.lin_density",
        **kwargs,
    )
    """
}

process test_popdep_cv_density {
    tag "test_popdep RelMean vs CV density"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/info", mode: 'copy', pattern: "*.{tsv,npz}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_16'

    input:
    val sites_ready

    output:
    path "*.png", emit: plots
    path "*.tsv", emit: info
    path "*.npz", emit: npz
    path "*.log", emit: logs

    script:
    def sites_dir = params.test_popdep_sites_dir ?: "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_popdep'}/sites"
    def scatter = params.test_popdep_scatter_tsv ?: "${params.output_dir}/${params.job}/stats/${params.mod ?: 'test_popdep'}/info/test_popdep.popdep.scatter.sites.tsv"
    def n_cells = params.popdep_cv_n_cells ?: 90000
    def counts_npz = params.popdep_cv_counts_npz ?: ''
    def plot_keep_drop = params.popdep_cv_plot_keep_drop
    if (plot_keep_drop == null) {
        plot_keep_drop = true
    }
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("test_popdep_cv_density.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.popdep_density_filter import ana_test_popdep_cv_density

    counts_npz = "${counts_npz}"
    kwargs = dict(n_cells=int("${n_cells}"), scatter="${scatter}")
    if counts_npz and counts_npz not in ("null", "None", ""):
        kwargs["counts_npz"] = counts_npz
        print(f"Reuse mean-CV histogram {counts_npz}; site-coverage 60-99% ...", flush=True)
    else:
        print("RelMean vs CV density on test_popdep sites, n_cells=${n_cells} mean=[0, 6] CV=observed ...", flush=True)
    kwargs["plot_keep_drop"] = str("${plot_keep_drop}").strip().lower() in ("true", "1", "yes")
    if not kwargs["plot_keep_drop"]:
        print("Skip keep/drop JointGrids; write site-coverage curve and class plot.", flush=True)
    ana_test_popdep_cv_density(
        "${sites_dir}",
        "test_popdep.popdep.cv_density",
        **kwargs,
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

process variant_ld_decay_r2_threshold_plot {
    tag "variant ld decay r2 threshold: ${id}"
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
    tuple val(id), val(chr), path("${id}.variant.ld_decay.r2_*.info.tsv"), emit: info
    tuple val(id), val(chr), path("${id}.variant.ld_decay.r2_*.th.tsv"), emit: th
    tuple val(id), val(chr), path("${id}.variant.ld_decay.r2_*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    def r2txt = "${params.ld_decay_r2_threshold}".replace('.', 'p')
    def bintxt = "${params.ld_decay_fine_bin_kb}".replace('.', 'p')
    """
    #!/usr/bin/env python
    import sys
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")
    sys.stdout = open("${chr}.ld_decay_r2_threshold.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.ld_decay import ana_ld_decay_mean_r2_threshold

    print("Fine-bin all-pairs mean r2 to ${params.ld_decay_r2_threshold} for ${id} ...")
    ana_ld_decay_mean_r2_threshold(
        input_file="${ld_vcor}",
        output_prefix="${id}.variant.ld_decay.r2_${r2txt}.${bintxt}kb",
        r2_threshold=${params.ld_decay_r2_threshold},
        max_distance_kb=${params.ld_decay_fine_max_kb},
        bin_size_kb=${params.ld_decay_fine_bin_kb},
        scatter_max_points=${params.ld_decay_scatter_points},
        chunksize=${params.ld_decay_chunksize},
        median_sample_per_bin=${params.ld_decay_median_sample_per_bin},
        workers=${task.cpus},
        random_seed=${params.ld_decay_random_seed},
    )
    """
}

process variant_ld_decay_rare_common_plot {
    tag "variant ld decay rare/common: ${id}"
    label 'cpus_16'
    cache false
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(ld_vcor), path(afreq)

    output:
    tuple val(id), val(chr), path("${id}.variant.ld_decay.rare_common.info.tsv"), emit: info
    tuple val(id), val(chr), path("${id}.variant.ld_decay.rare_common.th.tsv"), emit: th
    tuple val(id), val(chr), path("${id}.variant.ld_decay.rare_common.maf_scatter.info.tsv"), emit: maf_scatter
    tuple val(id), val(chr), path("${id}.variant.ld_decay.rare_common.dprime.info.tsv"), emit: dprime
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.ld_decay_rare_common_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant import ana_ld_decay_rare_common

    print("Running rare-rare / rare-common / common-common LD decay plotting for ${id} (${chr}) ...")
    ana_ld_decay_rare_common(
        input_file="${ld_vcor}",
        afreq_file="${afreq}",
        output_prefix="${id}.variant.ld_decay.rare_common",
        maf_threshold=${params.hf_maf},
        max_distance_kb=${params.ld_decay_max_kb},
        bin_size_kb=${params.ld_decay_bin_kb},
        scatter_max_points=${params.ld_decay_scatter_points},
        maf_scatter_points=${params.ld_decay_maf_scatter_points},
        chunksize=${params.ld_decay_chunksize},
        median_sample_per_bin=${params.ld_decay_median_sample_per_bin},
        workers=${task.cpus},
        random_seed=${params.ld_decay_random_seed},
        pfile_prefix="${params.output_dir}/${params.job}/process/${params.mod}/${id}_test.plink2",
    )
    """
}

process variant_ld_decay_maf_bin_plot {
    tag "variant ld decay maf bins: ${id}"
    label 'cpus_16'
    cache false
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(ld_vcor), path(afreq)

    output:
    tuple val(id), val(chr), path("${id}.variant.ld_decay.maf_bin.info.tsv"), emit: info
    tuple val(id), val(chr), path("${id}.variant.ld_decay.maf_bin.th.tsv"), emit: th
    tuple val(id), val(chr), path("${id}.variant.ld_decay.maf_bin.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.ld_decay_maf_bin_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.ld_decay import ana_ld_decay_maf_bins

    print("Running within-bin MAF LD decay for ${id} (${chr}) ...")
    ana_ld_decay_maf_bins(
        input_file="${ld_vcor}",
        afreq_file="${afreq}",
        output_prefix="${id}.variant.ld_decay.maf_bin",
        max_distance_kb=${params.ld_decay_max_kb},
        bin_size_kb=${params.ld_decay_bin_kb},
        chunksize=${params.ld_decay_chunksize},
    )
    """
}

process variant_ld_maf_bin_decay_replot {
    tag "ld maf-bin decay replot: ${id}"
    cache false
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), path(info)

    output:
    tuple val(id), path("*.png"), emit: plots
    tuple val(id), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${id}.ld_maf_bin_decay_replot.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.ld_decay import plot_ld_maf_bin_decay_from_table

    print("Replotting within-bin MAF LD decay from ${info} for ${id} ...")
    plot_ld_maf_bin_decay_from_table(
        "${info}",
        "${id}.variant.ld_decay.maf_bin",
        max_distance_kb=${params.ld_decay_max_kb},
    )
    """
}

process variant_ld_maf_bin_decay_all {
    tag "ld maf-bin decay pooled ABD: ${params.mod}"
    cache false
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path info_files

    output:
    path "all.variant.ld_decay.maf_bin.info.tsv", emit: info
    path "all.variant.ld_decay.maf_bin.th.tsv", emit: th
    path "all.variant.ld_decay.maf_bin.png", emit: plots
    path "*.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    from pathlib import Path
    sys.stdout = open("all.ld_maf_bin_decay_pooled.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.ld_decay import combine_ld_maf_bin_decay_tables

    paths = sorted(Path(".").glob("*.variant.ld_decay.maf_bin.info.tsv"))
    print(f"Pooling {len(paths)} subgenome MAF-bin LD decay tables...")
    combine_ld_maf_bin_decay_tables(
        paths,
        "all.variant.ld_decay.maf_bin",
        max_distance_kb=${params.ld_decay_max_kb},
    )
    """
}

process variant_ld_rare_common_decay_replot {
    tag "ld rare/common decay replot: ${id}"
    cache false
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), path(info)

    output:
    tuple val(id), path("*.png"), emit: plots
    tuple val(id), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${id}.ld_rare_common_decay_replot.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.ld_decay import plot_ld_rare_common_decay_from_table

    print("Replotting class LD decay from ${info} for ${id} ...")
    plot_ld_rare_common_decay_from_table(
        "${info}",
        "${id}.variant.ld_decay.rare_common",
        max_distance_kb=${params.ld_decay_max_kb},
    )
    """
}

process variant_ld_pplus_vs_q_overall_replot {
    tag "ld pplus vs q overall: ${id} ≤${params.pplus_vs_q_max_kb} kb"
    cache false
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), path(info)

    output:
    tuple val(id), path("*.png"), emit: plots
    tuple val(id), path("*.info.tsv"), emit: info
    tuple val(id), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${id}.ld_pplus_vs_q_overall_replot.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.ld_decay import plot_signed_dprime_pplus_vs_q_overall

    print(
        "Plotting P+ vs q overall (distance ≤ ${params.pplus_vs_q_max_kb} kb) "
        "from ${info} for ${id} ..."
    )
    plot_signed_dprime_pplus_vs_q_overall(
        "${info}",
        "${id}.variant.ld_decay.rare_common",
        maf_threshold=${params.hf_maf},
        max_distance_kb=${params.pplus_vs_q_max_kb},
        filename_suffix="${params.pplus_vs_q_suffix}",
    )
    """
}

process variant_ld_maf_scatter_r2max_replot {
    tag "ld maf scatter r_max contours: ${id}"
    label 'cpus_16'
    cache false
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), path(info)

    output:
    tuple val(id), path("*.png"), emit: plots
    tuple val(id), path("*.info.tsv"), emit: info
    tuple val(id), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${id}.ld_maf_scatter_r2max_replot.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.ld_decay import plot_ld_maf_r2_scatter_from_table

    print("Replotting MAF scatter and dosage-based |D'| for ${id} ...")
    plot_ld_maf_r2_scatter_from_table(
        "${info}",
        "${id}.variant.ld_decay.rare_common",
        maf_max=${params.hf_maf},
        maf_threshold=${params.hf_maf},
        max_distance_kb=${params.ld_decay_max_kb},
        bin_size_kb=${params.ld_decay_bin_kb},
        pfile_prefix="${params.output_dir}/${params.job}/process/${params.mod}/${id}_test.plink2",
        workers=${task.cpus},
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


process variant_ld_homoeolog_compare_plot {
    tag "ld homoeolog ABD compare"
    cache false
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path a_by_chr
    path b_by_chr
    path d_by_chr

    output:
    path "*.info.tsv", emit: info
    path "*.png", emit: plots
    path "*.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("ABD.ld_homoeolog_compare.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.ld_decay import plot_ld_homoeolog_subgenome_compare

    print("Plotting homoeolog LD compare across A/B/D ...")
    plot_ld_homoeolog_subgenome_compare(
        "${a_by_chr}",
        "${b_by_chr}",
        "${d_by_chr}",
        "ABD.variant.ld_decay.homoeolog_compare",
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

process variant_mac_het_excess_stats {
    tag "variant mac het excess: ${id}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '16.GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.sweep.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.sweep_by_mac.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*max_allowed_het*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), path(gcount), path(het)

    output:
    tuple val(id), val(chr), path("*.sweep.tsv"), emit: sweep
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    def f_grid = params.mac_het_f_grid ?: '0.90,0.95,0.98,0.99'
    def p_thresholds = params.mac_het_p_thresholds ?: '0.05,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-12'
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.mac_het_excess.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.mac_het_excess import ana_mac_het_excess_sweep

    print(f"Excess-het probability sweep for ${id} (no filtered pfile)...")
    ana_mac_het_excess_sweep(
        "${gcount}",
        "${id}.variant.mac_het",
        het_path="${het}",
        f_grid="${f_grid}",
        thresholds="${p_thresholds}",
    )
    """
}

process variant_mac_het_excess_combine {
    tag "variant mac het excess combine"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 1
    memory '4.GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.sweep.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    path sweep_files

    output:
    path("all.variant.mac_het.sweep.tsv"), emit: sweep
    path("*.png"), emit: plots
    path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("all.mac_het_excess.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.mac_het_excess import combine_mac_het_excess_sweeps
    import glob

    files = sorted(glob.glob("*.sweep.tsv"))
    print(f"Combining {len(files)} subgenome excess-het sweeps...")
    combine_mac_het_excess_sweeps(files, "all.variant.mac_het")
    """
}

process variant_mac_hetfrac_heatmap_retained {
    tag "mac hetfrac heatmap retained: ${id}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '32.GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.het_retain.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), path(gcount), path(het)

    output:
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.het_retain.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    def alpha = params.mac_het_retain_alpha ?: '1e-6'
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.mac_hetfrac_heatmap_retained.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.mac_het_excess import ana_mac_hetfrac_heatmap_retained

    print(f"MAC het-fraction heatmap for ${id} retained at P_excess>=${alpha} (median F)...")
    ana_mac_hetfrac_heatmap_retained(
        "${gcount}",
        "${id}.variant.mac",
        het_path="${het}",
        alpha=float("${alpha}"),
    )
    """
}

process variant_popdep_het_filter_scatter {
    tag "popdep het-filter scatter: ${id}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '32.GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), path(gcount), path(het), path(popdep)

    output:
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    def alpha = params.mac_het_retain_alpha ?: '1e-6'
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.popdep_het_filter_scatter.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.mac_het_excess import ana_popdep_scatter_by_het_filter

    print(f"Popdep het-filter scatters for ${id} at P_excess>=${alpha} (median F)...")
    ana_popdep_scatter_by_het_filter(
        "${gcount}",
        "${popdep}",
        "${id}.variant.popdep.het_filter",
        het_path="${het}",
        alpha=float("${alpha}"),
    )
    """
}

process variant_mac_het_rh_diagnostics {
    tag "mac het R_H diagnostics: ${id}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '32.GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.extract"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), path(gcount), path(het), path(popdep), path(popdep_mapq30, stageAs: 'popdep_mapq30.info.tsv'), val(has_mapq30), path(mq), path(vmiss)

    output:
    tuple val(id), val(chr), path("${id}.variant.mac_het.rh.th.tsv"), emit: th
    tuple val(id), val(chr), path("${id}.variant.mac_het.rh.all_het.extract"), emit: extract
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.png"), emit: plots, optional: true
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    def mq30_py = has_mapq30 ? "\"${popdep_mapq30}\"" : 'None'
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.mac_het_rh_diagnostics.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.mac_het_rh import ana_mac_het_rh_diagnostics

    print(f"MAC–H R_H band diagnostics for ${id} ...")
    ana_mac_het_rh_diagnostics(
        "${gcount}",
        "${id}.variant.mac_het.rh",
        het_path="${het}",
        popdep_file="${popdep}",
        popdep_mapq30_file=${mq30_py},
        mq_file="${mq}",
        vmiss_file="${vmiss}",
    )
    """
}

process variant_mac_het_rh_combine {
    tag "mac het R_H combine"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 1
    memory '4.GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    path th_files

    output:
    path("all.variant.mac_het.rh.th.info.tsv"), emit: info
    path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("all.mac_het_rh.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.mac_het_rh import combine_mac_het_rh_thresholds
    import glob

    files = sorted(glob.glob("*.th.tsv"))
    print(f"Combining {len(files)} subgenome R_H threshold tables...")
    combine_mac_het_rh_thresholds(files, "all.variant.mac_het.rh")
    """
}

process variant_mac_het_rh_sample_burden {
    tag "mac het all_het sample burden: ${id}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 1
    memory '8.GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), val(chr), path(het)

    output:
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots, optional: true
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.mac_het_rh_sample_burden.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.mac_het_rh import ana_allhet_sample_burden

    print(f"All-het sample heterozygous burden for ${id} ...")
    ana_allhet_sample_burden("${het}", "${id}.variant.mac_het.rh_sample")
    """
}

process variant_ld_r2_decay_along_chr_plot {
    tag "along-chr r2 decay: ${id} chr${chrom}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 8
    memory '64.GB'

    input:
    tuple val(id), val(chrom), path(vcor), path(focals)

    output:
    tuple val(id), val(chrom), path("${id}.variant.ld_r2_decay_along*.${chrom}.info.tsv"), emit: info
    tuple val(id), val(chrom), path("${id}.variant.ld_r2_decay_along*.${chrom}.th.tsv"), emit: th
    tuple val(id), val(chrom), path("${id}.variant.ld_r2_decay_along*.${chrom}.png"), emit: plots
    tuple val(id), val(chrom), path("${chrom}.ld_r2_decay_along*.log"), emit: logs
    path "${id}.variant.ld_r2_decay_along*.${chrom}.focal_decay.info.tsv", optional: true, emit: focal_decay

    script:
    """
    #!/usr/bin/env python
    import sys
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")
    from genetics.genomics.variant.ld_along_chr import (
        ana_ld_r2_decay_along_chr,
        ld_along_chr_output_stem,
    )

    stem = ld_along_chr_output_stem(${params.ld_along_chr_r2})
    sys.stdout = open(f"${chrom}.{stem}.log", "w")
    sys.stderr = sys.stdout

    print("Along-chr mean r2 decay for ${id} chr${chrom} ...")
    ana_ld_r2_decay_along_chr(
        input_file="${vcor}",
        focals_file="${focals}",
        output_prefix=f"${id}.variant.{stem}.${chrom}",
        r2_threshold=${params.ld_along_chr_r2},
        max_distance_kb=${params.ld_along_chr_max_kb},
        bin_size_kb=${params.ld_along_chr_bin_kb},
        min_pairs_per_bin=${params.ld_along_chr_min_pairs_per_bin},
    )
    """
}

process variant_ld_r2_decay_along_subgenome_plot {
    tag "along-chr r2 decay: ${id} chr1-chr7 stacked"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 8
    memory '64.GB'

    input:
    tuple val(id), val(chroms), path(vcors), path(focals)

    output:
    tuple val(id), path("${id}.variant.ld_r2_decay_along*.png"), emit: plots
    path "${id}.variant.ld_r2_decay_along*.info.tsv", emit: info
    path "${id}.variant.ld_r2_decay_along*.th.tsv", emit: th
    tuple val(id), path("${id}.ld_r2_decay_along*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")
    from genetics.genomics.variant.ld_along_chr import (
        ana_ld_r2_decay_along_subgenome,
        collect_chrom_ld_files,
        ld_along_chr_output_stem,
    )

    stem = ld_along_chr_output_stem(${params.ld_along_chr_r2})
    sys.stdout = open(f"${id}.{stem}.log", "w")
    sys.stderr = sys.stdout

    print("Along-chr mean r2 decay stacked chr1-chr7 for ${id} ...")
    pairs = collect_chrom_ld_files(".", "${id}")
    if not pairs:
        raise SystemExit("no chromosome vcor/focals pairs staged in the work dir")
    print(f"[Info] Chromosomes: {', '.join(chrom for chrom, _, _ in pairs)}")
    ana_ld_r2_decay_along_subgenome(
        pairs,
        f"${id}.variant.{stem}",
        subgenome="${id}",
        r2_threshold=${params.ld_along_chr_r2},
        max_distance_kb=${params.ld_along_chr_max_kb},
        bin_size_kb=${params.ld_along_chr_bin_kb},
        min_pairs_per_bin=${params.ld_along_chr_min_pairs_per_bin},
    )
    """
}

process variant_filter_cascade_maf_plot {
    tag "filter cascade maf plot: ${params.mod}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*filter_steps*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.filtered_sites.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*subgenome_maf_rare_common_prop.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*stepwise_maf_class*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*genome_pipeline_maf_class*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 4
    memory '32.GB'

    input:
    path wap_keep

    output:
    path("all.variant.filter_cascade.maf_rare_common.png"), emit: plot
    path("all.variant.filter_cascade.maf_rare_common.th.tsv"), emit: maf_th
    path("all.variant.filter_cascade.filter.th.tsv"), emit: filter_th
    path("all.variant.filter_cascade.filter_steps_total.tsv"), emit: steps_total
    path("all.variant.filter_cascade.filter_steps_by_subgenome.tsv"), emit: steps_by_sg
    path("all.variant.filter_cascade.filtered_sites.tsv"), emit: filtered_sites
    path("all.variant.filter_cascade.subgenome_maf_rare_common_prop.png"), emit: sg_prop_plot
    path("all.variant.filter_cascade.subgenome_maf_rare_common_prop.tsv"), emit: sg_prop_tsv
    path("all.variant.filter_cascade.stepwise_maf_class_counts.tsv"), emit: maf_class_counts
    path("all.variant.filter_cascade.stepwise_maf_class_retention.tsv"), emit: maf_class_retention
    path("all.variant.filter_cascade.genome_pipeline_maf_class_projection.tsv"), emit: genome_proj
    path("all.variant.filter_cascade.genome_pipeline_maf_class_compare.tsv"), emit: genome_compare
    path("filter_cascade_maf_plot.log"), emit: logs

    script:
    def het_mod = params.filter_cascade_het_mod ?: 'test_common_only'
    def het_base = "${params.output_dir}/${params.job}/process/${het_mod}/sample"
    def alpha = params.filter_cascade_het_alpha ?: '1e-6'
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("filter_cascade_maf_plot.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.filter_cascade_maf_plot import run_filter_cascade_maf_plot

    process_dir = "${params.output_dir}/${params.job}/process/${params.mod}"
    het_paths = {
        sg: f"${het_base}/{sg}.info.het"
        for sg in ("A", "B", "D", "Others")
    }
    print(f"Filter cascade MAF plot for mod=${params.mod} (het_mod=${het_mod})...")
    run_filter_cascade_maf_plot(
        process_dir,
        "all.variant.filter_cascade",
        het_paths=het_paths,
        wap_keep_path="${wap_keep}",
        alpha=float("${alpha}"),
        plink2="plink2",
    )
    """
}

process variant_subgenome_maf_plot {
    tag "subgenome MAF stacked: ${params.mod}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 4
    memory '16.GB'

    input:
    path afreq_files

    output:
    path "all.variant.maf.subgenome.dist.png", emit: dist
    path "all.variant.maf.subgenome.dist.log.png", emit: dist_log
    path "all.variant.maf.subgenome.th.tsv", emit: summary
    path "subgenome_maf_plot.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    from pathlib import Path
    sys.stdout = open("subgenome_maf_plot.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.subgenome_maf_tajima import run_subgenome_maf_plots

    paths = sorted(Path(".").glob("*.info.afreq"))
    print(f"Stacked MAF plots from {len(paths)} afreq tables...")
    run_subgenome_maf_plots(paths, "all.variant.maf")
    """
}

process variant_subgenome_tajima_plot {
    tag "subgenome Tajima D stacked: ${params.mod}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.{th.tsv,info.tsv}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 4
    memory '8.GB'

    input:
    path tajima_files

    output:
    path "all.variant.tajima_d.subgenome.dist.png", emit: dist
    path "all.variant.tajima_d.subgenome.th.tsv", emit: summary
    path "all.variant.tajima_d.subgenome.info.tsv", emit: info
    path "subgenome_tajima_plot.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    from pathlib import Path
    sys.stdout = open("subgenome_tajima_plot.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.subgenome_maf_tajima import run_subgenome_tajima_plots

    paths = sorted(Path(".").glob("*.Tajima.D"))
    print(f"Stacked Tajima D plots from {len(paths)} tables...")
    run_subgenome_tajima_plots(paths, "all.variant.tajima_d")
    """
}

process variant_circos_ld_map {
    tag "circos rare/common density + LD blocks"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    label 'cpus_16'

    input:
    tuple path(pvars, stageAs: 'density/*'), path(afreqs), path(ld_pvars, stageAs: 'ld_pvar/*'), path(lds)

    output:
    path("all.variant.circos_ld_map.density.info.tsv"), emit: density
    path("all.variant.circos_ld_map.blocks.info.tsv"), emit: blocks
    path("all.variant.circos_ld_map.blocks_1mb.info.tsv"), emit: blocks_1mb
    path("all.variant.circos_ld_map.th.tsv"), emit: th
    path("all.variant.circos_ld_map.png"), emit: plot
    path("circos_ld_map.log"), emit: logs

    script:
    def maf_thr = params.circos_maf_threshold ?: 0.05
    def r2_thr = params.circos_ld_r2_threshold ?: 0.8
    def dprime_q1 = params.circos_ld_dprime_q1 ?: 0.98
    def bin_bp = params.circos_bin_bp ?: 1000000
    def cov_frac = params.circos_ld_block_coverage_frac ?: 0.5
    """
    #!/usr/bin/env python
    import sys
    from pathlib import Path

    sys.stdout = open("circos_ld_map.log", "w", buffering=1)
    sys.stderr = sys.stdout

    from genetics.genomics.variant.circos_ld_map import ana_circos_ld_map

    ids = ["A", "B", "D"]
    pvar_files = [f"density/{sg}_test.plink2.pvar" for sg in ids]
    afreq_files = [f"{sg}.info.afreq" for sg in ids]
    ld_pvar_files = [f"ld_pvar/{sg}_test.plink2.pvar" for sg in ids]
    ld_files = [f"{sg}.info.ld_block.ld" for sg in ids]
    missing = [p for p in pvar_files + afreq_files + ld_pvar_files + ld_files if not Path(p).exists()]
    if missing:
        raise FileNotFoundError(f"Circos LD map missing inputs: {missing}")
    print("Circos rare/common density + LD blocks...", flush=True)
    ana_circos_ld_map(
        pvar_files=pvar_files,
        afreq_files=afreq_files,
        ld_files=ld_files,
        ld_pvar_files=ld_pvar_files,
        output_prefix="all.variant.circos_ld_map",
        maf_threshold=float("${maf_thr}"),
        r2_threshold=float("${r2_thr}"),
        dprime_q1_threshold=float("${dprime_q1}"),
        bin_bp=int("${bin_bp}"),
        coverage_frac=float("${cov_frac}"),
        n_jobs=int("${task.cpus}"),
    )
    """
}

process variant_circos_ld_map_replot {
    tag "circos replot from density + LD-block tables"
    cache false
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple path(density), path(blocks)

    output:
    path("all.variant.circos_ld_map.blocks_1mb.info.tsv"), emit: blocks_1mb
    path("all.variant.circos_ld_map.png"), emit: plot
    path("circos_ld_map_replot.log"), emit: logs

    script:
    def bin_bp = params.circos_bin_bp ?: 1000000
    def cov_frac = params.circos_ld_block_coverage_frac ?: 0.5
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("circos_ld_map_replot.log", "w", buffering=1)
    sys.stderr = sys.stdout

    from genetics.genomics.variant.circos_ld_map import plot_circos_ld_map_from_tables

    print("Circos replot from ${density} + ${blocks} ...", flush=True)
    plot_circos_ld_map_from_tables(
        "${density}",
        "${blocks}",
        "all.variant.circos_ld_map",
        bin_bp=int("${bin_bp}"),
        coverage_frac=float("${cov_frac}"),
    )
    """
}

process variant_vep_atlas {
    tag "VEP atlas composition + window counts"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    label 'cpus_32'
    errorStrategy 'terminate'

    input:
    tuple path(genome_vep, stageAs: 'genome/*'), path(exon_vep, stageAs: 'exon/*'), path(genome_afreq, stageAs: 'genome_afreq/*'), path(genome_acount, stageAs: 'genome_acount/*'), path(genome_pvar, stageAs: 'genome_pvar/*'), path(exon_afreq, stageAs: 'exon_afreq/*'), path(exon_acount, stageAs: 'exon_acount/*'), path(exon_pvar, stageAs: 'exon_pvar/*')

    output:
    path("all.variant.vep_atlas.composition.png"), emit: composition
    path("all.variant.vep_atlas.freq_function.png"), emit: heatmap
    path("all.variant.vep_atlas.circos_windows.png"), emit: circos
    path("all.variant.vep_atlas.mirror.png"), emit: mirror
    path("all.variant.vep_atlas.enrichment.png"), emit: enrichment
    path("all.variant.vep_atlas.windows.violin.png"), emit: violin
    path("all.variant.vep_atlas*.tsv"), emit: tables
    path("vep_atlas.log"), emit: logs

    script:
    def bin_bp = params.vep_atlas_bin_bp ?: 10000000
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("vep_atlas.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.vep_atlas import run_vep_atlas

    ids = ["A", "B", "D"]
    print("VEP atlas: genome common/rare + exon rare predicted LoF...", flush=True)
    run_vep_atlas(
        {sg: f"genome/{sg}.vep.tsv" for sg in ids},
        {sg: f"exon/{sg}.vep.tsv" for sg in ids},
        "all.variant.vep_atlas",
        genome_afreq={sg: f"genome_afreq/{sg}.info.afreq" for sg in ids},
        genome_acount={sg: f"genome_acount/{sg}.assoc.acount" for sg in ids},
        genome_pvar={sg: f"genome_pvar/{sg}_test.plink2.pvar" for sg in ids},
        exon_afreq={sg: f"exon_afreq/{sg}.info.afreq" for sg in ids},
        exon_acount={sg: f"exon_acount/{sg}.assoc.acount" for sg in ids},
        exon_pvar={sg: f"exon_pvar/{sg}_test.plink2.pvar" for sg in ids},
        bin_bp=int("${bin_bp}"),
    )
    """
}

process variant_vep_atlas_circos_replot {
    tag "VEP atlas Circos replot"
    cache false
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    errorStrategy 'terminate'

    input:
    tuple path(windows), path(ld_blocks), path(gene_gff)

    output:
    path("all.variant.vep_atlas.circos_windows.png"), emit: circos
    path("vep_atlas_circos_replot.log"), emit: logs

    script:
    def bin_bp = params.vep_atlas_bin_bp ?: 10000000
    def peri_flank = params.vep_atlas_peri_flank_bp ?: 50000000
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("vep_atlas_circos_replot.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.vep_atlas import plot_vep_atlas_circos

    print("VEP atlas Circos replot from ${windows} + ${ld_blocks} + ${gene_gff} ...", flush=True)
    plot_vep_atlas_circos(
        "${windows}",
        "all.variant.vep_atlas.circos_windows.png",
        bin_bp=int("${bin_bp}"),
        ld_blocks="${ld_blocks}",
        gene_gff="${gene_gff}",
        peri_flank_bp=int("${peri_flank}"),
    )
    """
}

process popdep_density_hist_chr {
    tag "PopDep density hist ${bgz.simpleName}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'test_popdep_chr'
    errorStrategy 'retry'
    maxRetries 1

    input:
    path bgz

    output:
    path("*.rel_mean_sd.counts.npz"), emit: npz

    script:
    def stem = bgz.simpleName.replace('.popdep.txt', '')
    """
    #!/usr/bin/env python
    import sys
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")
    from genetics.genomics.variant.popdep_density_filter import write_chr_histogram
    write_chr_histogram("${bgz}", "${stem}.rel_mean_sd.counts.npz")
    """
}

process popdep_density_merge {
    tag "PopDep density merge"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.{tsv,npz}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    cpus 4
    memory '16.GB'
    errorStrategy 'terminate'

    input:
    path npz_files

    output:
    path("all.variant.popdep_density.counts.npz"), emit: counts
    path("all.variant.popdep_density.masks.npz"), emit: masks
    path("all.variant.popdep_density.trusted_cells.tsv"), emit: cells
    path("all.variant.popdep_density.th.tsv"), emit: thresholds
    path("all.variant.popdep_density.rel_mean_sd.density.png"), emit: plot
    path("popdep_density_merge.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("popdep_density_merge.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")
    from pathlib import Path
    from genetics.genomics.variant.popdep_density_filter import merge_histogram_files, write_density_outputs
    paths = sorted(str(p) for p in Path('.').glob('*.rel_mean_sd.counts.npz'))
    print(f"Merging {len(paths)} chromosome histograms ...", flush=True)
    counts, stats = merge_histogram_files(paths, "all.variant.popdep_density.counts.npz")
    write_density_outputs(counts, stats, "all.variant.popdep_density")
    """
}

process variant_filter_cascade_popdep_density {
    tag "filter cascade popdep density: ${params.mod}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    label 'cpus_16'
    errorStrategy 'terminate'

    input:
    tuple path(masks), path(density_th), path(wap_keep)

    output:
    path("all.variant.popdep_density_cascade.filter_steps_total.tsv"), emit: steps_total
    path("all.variant.popdep_density_cascade.filter_steps_by_subgenome.tsv"), emit: steps_by_sg
    path("all.variant.popdep_density_cascade.filter_steps_genome_estimate.tsv"), emit: estimate
    path("all.variant.popdep_density_cascade.filter.th.tsv"), emit: thresholds
    path("all.variant.popdep_density_cascade.filtered_sites.tsv"), emit: sites, optional: true
    path("all.variant.popdep_density_cascade.maf_rare_common.png"), emit: plot, optional: true
    path("filter_cascade_popdep_density.log"), emit: logs

    script:
    def het_mod = params.filter_cascade_het_mod ?: 'test_common_only'
    def het_base = "${params.output_dir}/${params.job}/process/${het_mod}/sample"
    def alpha = params.filter_cascade_het_alpha ?: '1e-6'
    def thin_rate = params.thin_rate ?: '0.001'
    def process_dir = "${params.output_dir}/${params.job}/process/${params.mod}"
    def popdep_dir = params.popdep_density_dir ?: params.popdep_dir
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("filter_cascade_popdep_density.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from pathlib import Path
    import numpy as np
    from genetics.genomics.variant.popdep import annotate_variants_popdep_from_pvar
    from genetics.genomics.variant.filter_cascade_maf_plot import run_filter_cascade_popdep_density
    from infra.utils.io import load_thresholds

    process_dir = Path("${process_dir}")
    popdep_dir = "${popdep_dir}"
    masks = np.load("${masks}")
    mask_60 = masks["mask_60"]
    mask_80 = masks["mask_80"]
    genome_stats = load_thresholds("${density_th}") or {}
    popdep_tables = {}
    for sg in ("A", "B", "D", "Others"):
        pvar = process_dir / f"{sg}_test.plink2.pvar"
        out = Path(f"{sg}.ntaxa500.popdep.info.tsv")
        print(f"Annotating {pvar} from {popdep_dir} ...", flush=True)
        annotate_variants_popdep_from_pvar(str(pvar), popdep_dir, str(out), use_tabix=True)
        popdep_tables[sg] = out
    het_paths = {sg: f"${het_base}/{sg}.info.het" for sg in ("A", "B", "D", "Others")}
    print("Running popdep-density filter cascade ...", flush=True)
    run_filter_cascade_popdep_density(
        process_dir,
        "all.variant.popdep_density_cascade",
        het_paths=het_paths,
        wap_keep_path="${wap_keep}",
        mask_60=mask_60,
        mask_80=mask_80,
        popdep_tables=popdep_tables,
        genome_density_stats={
            "n_rows": int(float(genome_stats.get("n_rows", 0) or 0)),
            "n_valid": int(float(genome_stats.get("n_valid", 0) or 0)),
            "n_sites_top60": int(float(genome_stats.get("n_sites_top60", 0) or 0)),
            "n_sites_top80": int(float(genome_stats.get("n_sites_top80", 0) or 0)),
        },
        thin_rate=float("${thin_rate}"),
        alpha=float("${alpha}"),
        plink2="plink2",
    )
    """
}

process popdep_cv_reliable_hist_chr {
    tag "PopDep CV hist ${bgz.simpleName}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'test_popdep_chr'
    errorStrategy 'retry'
    maxRetries 1

    input:
    path bgz

    output:
    path("*.rel_mean_cv.counts.npz"), emit: npz

    script:
    def stem = bgz.simpleName.replace('.popdep.txt', '')
    """
    #!/usr/bin/env python
    import sys
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")
    from genetics.genomics.variant.popdep_density_filter import write_chr_mean_cv_histogram
    write_chr_mean_cv_histogram("${bgz}", "${stem}.rel_mean_cv.counts.npz")
    """
}

process popdep_cv_reliable_merge {
    tag "PopDep CV reliable merge"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.{tsv,npz}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    cpus 4
    memory '16.GB'
    errorStrategy 'terminate'

    input:
    path npz_files

    output:
    path("all.site.popdep.cv_reliable.counts.npz"), emit: counts
    path("all.site.popdep.cv_reliable.masks.npz"), emit: masks
    path("all.site.popdep.cv_reliable.trusted_cells.tsv"), emit: cells
    path("all.site.popdep.cv_reliable.site_coverage_curve.tsv"), emit: curve
    path("all.site.popdep.cv_reliable.th.tsv"), emit: thresholds
    path("all.site.popdep.cv_reliable.rel_mean_cv.density.png"), emit: heatmap
    path("all.site.popdep.cv_reliable.site_coverage_curve.png"), emit: curve_plot
    path("popdep_cv_reliable_merge.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("popdep_cv_reliable_merge.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")
    from pathlib import Path
    from genetics.genomics.variant.popdep_density_filter import (
        merge_histogram_files,
        write_cv_reliable_density_outputs,
    )
    paths = sorted(str(p) for p in Path('.').glob('*.rel_mean_cv.counts.npz'))
    print(f"Merging {len(paths)} RelMean-CV histograms ...", flush=True)
    counts, stats = merge_histogram_files(paths, "all.site.popdep.cv_reliable.counts.npz")
    write_cv_reliable_density_outputs(counts, stats, "all.site.popdep.cv_reliable")
    """
}

process popdep_cv_reliable_lib_chr {
    tag "PopDep CV lib ${bgz.simpleName}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'test_popdep_chr'
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/keep60", mode: 'copy', pattern: "*.keep60.bed.gz"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/keep80", mode: 'copy', pattern: "*.keep80.bed.gz"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.lib.info.tsv"
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple path(bgz), path(masks)

    output:
    path("*.keep60.bed.gz"), emit: bed60
    path("*.keep80.bed.gz"), emit: bed80
    path("*.lib.info.tsv"), emit: info

    script:
    def stem = bgz.simpleName.replace('.popdep.txt', '')
    """
    #!/usr/bin/env python
    import sys
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")
    from genetics.genomics.variant.popdep_density_filter import write_chr_cv_reliable_beds
    write_chr_cv_reliable_beds("${bgz}", "${masks}", "${stem}.popdep.cv_reliable")
    """
}

process popdep_cv_reliable_concat {
    tag "PopDep CV concat libraries"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "all.site.popdep.cv_reliable.keep*.bed.gz"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    cpus 2
    memory '8.GB'
    errorStrategy 'terminate'

    input:
    path bed60_files
    path bed80_files
    path info_files

    output:
    path("all.site.popdep.cv_reliable.keep60.bed.gz"), emit: bed60
    path("all.site.popdep.cv_reliable.keep80.bed.gz"), emit: bed80
    path("all.site.popdep.cv_reliable.lib.info.tsv"), emit: info
    path("popdep_cv_reliable_concat.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("popdep_cv_reliable_concat.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")
    from pathlib import Path
    from genetics.genomics.variant.popdep_density_filter import concat_gzip_beds, merge_lib_info_tables
    beds60 = sorted(Path('.').glob('*.keep60.bed.gz'))
    beds80 = sorted(Path('.').glob('*.keep80.bed.gz'))
    infos = sorted(Path('.').glob('*.lib.info.tsv'))
    print(f"Concat keep60={len(beds60)} keep80={len(beds80)} info={len(infos)}", flush=True)
    concat_gzip_beds(beds60, "all.site.popdep.cv_reliable.keep60.bed.gz")
    concat_gzip_beds(beds80, "all.site.popdep.cv_reliable.keep80.bed.gz")
    merge_lib_info_tables(infos, "all.site.popdep.cv_reliable.lib.info.tsv")
    """
}

process variant_filter_cascade_cv_reliable {
    tag "cv reliable cascade: ${params.mod}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    label 'cpus_16'
    errorStrategy 'terminate'

    input:
    val ready

    output:
    path("all.variant.cv_reliable.keep*.filter_steps_total.tsv"), emit: steps
    path("all.variant.cv_reliable.keep*.filter_steps_by_subgenome.tsv"), emit: by_sg
    path("all.variant.cv_reliable.keep*.filter_steps_genome_estimate.tsv"), emit: estimate
    path("all.variant.cv_reliable.keep*.filter.th.tsv"), emit: thresholds
    path("all.variant.cv_reliable.keep*.filtered_sites.tsv"), emit: sites, optional: true
    path("all.variant.cv_reliable.keep*.maf_rare_common.png"), emit: plot, optional: true
    path("filter_cascade_cv_reliable.log"), emit: logs

    script:
    def het_mod = params.filter_cascade_het_mod ?: 'test_common_only'
    def het_base = "${params.output_dir}/${params.job}/process/${het_mod}/sample"
    def wap_mod = params.filter_cascade_wap_mod ?: 'test_wap'
    def keep_name = params.filter_cascade_wap_keep ?: "${wap_mod}.keep"
    def wap_keep = "${params.output_dir}/${params.job}/process/${wap_mod}/${keep_name}"
    def bed_root = params.popdep_cv_reliable_dir
    def alpha = params.filter_cascade_het_alpha ?: '1e-6'
    def thin_rate = params.thin_rate ?: '0.001'
    def process_dir = "${params.output_dir}/${params.job}/process/${params.mod}"
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("filter_cascade_cv_reliable.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from pathlib import Path
    from genetics.genomics.variant.filter_cascade_maf_plot import run_filter_cascade_cv_reliable

    bed_root = Path("${bed_root}")
    het_paths = {sg: f"${het_base}/{sg}.info.het" for sg in ("A", "B", "D", "Others")}
    for pct in (80, 60):
        beds = sorted((bed_root / f"keep{pct}").glob("chr*.bed.gz"))
        print(f"keep{pct}: {len(beds)} chromosome BEDs under {bed_root}", flush=True)
        if len(beds) != 45:
            raise SystemExit(f"expected 45 keep{pct} BEDs, found {len(beds)}")
        run_filter_cascade_cv_reliable(
            "${process_dir}",
            f"all.variant.cv_reliable.keep{pct}",
            bed_paths=beds,
            het_paths=het_paths,
            wap_keep_path="${wap_keep}",
            keep_pct=pct,
            thin_rate=float("${thin_rate}"),
            alpha=float("${alpha}"),
            plink2="plink2",
        )
    """
}

process variant_wap_cv_reliable_filter {
    tag "cv reliable three-step: ${mod}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    maxForks 1
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/info", mode: 'copy', pattern: "all.variant.cv_reliable.keep*.filter_steps_*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/info", mode: 'copy', pattern: "all.variant.cv_reliable.keep*.filtered_sites.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/plots", mode: 'copy', pattern: "all.variant.cv_reliable.keep*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/thresholds", mode: 'copy', pattern: "all.variant.cv_reliable.keep*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${mod}/logs", mode: 'copy', pattern: "*.log"
    label 'cpus_16'
    errorStrategy 'terminate'

    input:
    val mod

    output:
    tuple val(mod), path("all.variant.cv_reliable.keep80.filtered_sites.tsv"), emit: keep80
    tuple val(mod), path("all.variant.cv_reliable.keep60.filtered_sites.tsv"), emit: keep60
    path("all.variant.cv_reliable.keep*.filter_steps_*.tsv"), emit: steps
    path("all.variant.cv_reliable.keep*.png"), emit: plots
    path("all.variant.cv_reliable.keep*.th.tsv"), emit: thresholds
    path("wap_cv_reliable_filter.log"), emit: logs

    script:
    def bed_root = params.popdep_cv_reliable_dir
    def alpha = params.filter_cascade_het_alpha ?: '1e-6'
    def min_mac = params.filter_cascade_min_mac ?: 5
    def process_dir = "${params.output_dir}/${params.job}/process/${mod}"
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("wap_cv_reliable_filter.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from pathlib import Path
    from genetics.genomics.variant.filter_cascade_maf_plot import (
        discover_pfile_subgenomes,
        het_paths_from_pfiles,
        run_cv_reliable_three_step,
    )

    process_dir = Path("${process_dir}")
    bed_root = Path("${bed_root}")
    print(f"three-step filter on {process_dir} (pfiles are not rewritten)", flush=True)
    subgenomes = discover_pfile_subgenomes(process_dir)
    print(f"subgenomes: {', '.join(subgenomes)}", flush=True)
    het_paths = het_paths_from_pfiles(process_dir, subgenomes, Path("het"), plink2="plink2")
    for pct in (80, 60):
        beds = sorted((bed_root / f"keep{pct}").glob("chr*.bed.gz"))
        print(f"keep{pct}: {len(beds)} chromosome BEDs", flush=True)
        if len(beds) != 45:
            raise SystemExit(f"expected 45 keep{pct} BEDs, found {len(beds)}")
        run_cv_reliable_three_step(
            process_dir,
            f"all.variant.cv_reliable.keep{pct}",
            bed_paths=beds,
            keep_pct=pct,
            het_paths=het_paths,
            subgenomes=subgenomes,
            alpha=float("${alpha}"),
            plink2="plink2",
            min_mac=int("${min_mac}"),
        )
    """
}

process variant_vep_atlas_cv_reliable {
    tag "VEP atlas on cv-reliable sites"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "all.variant.vep_atlas.cv_reliable.*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "all.variant.vep_atlas.cv_reliable.*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "all.variant.vep_atlas.cv_reliable.*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    label 'cpus_32'
    errorStrategy 'terminate'

    input:
    tuple path(genome80, stageAs: 'genome.keep80.filtered_sites.tsv'), path(genome60, stageAs: 'genome.keep60.filtered_sites.tsv'), path(exon80, stageAs: 'exon.keep80.filtered_sites.tsv'), path(exon60, stageAs: 'exon.keep60.filtered_sites.tsv'), path(ld_blocks), path(gene_gff)

    output:
    path("all.variant.vep_atlas.cv_reliable.keep*.png"), emit: plots
    path("all.variant.vep_atlas.cv_reliable.keep*.tsv"), emit: tables
    path("all.variant.vep_atlas.cv_reliable.keep*.th.tsv"), emit: thresholds
    path("vep_atlas_cv_reliable.log"), emit: logs

    script:
    def bin_bp = params.vep_atlas_bin_bp ?: 10000000
    def genome_dir = "${params.output_dir}/${params.job}/process/${params.mod}"
    def exon_dir = "${params.output_dir}/${params.job}/process/${params.vep_atlas_exon_mod}"
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("vep_atlas_cv_reliable.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    import pandas as pd
    from genetics.genomics.variant.vep_atlas import run_vep_atlas

    def load_ids(path):
        frame = pd.read_csv(path, sep="\\t", usecols=["ID"], dtype=str)
        return set(frame["ID"].dropna().astype(str))

    ids = ["A", "B", "D"]
    genome_vep = {sg: f"${genome_dir}/variant/{sg}.vep.tsv" for sg in ids}
    exon_vep = {sg: f"${exon_dir}/variant/{sg}.vep.tsv" for sg in ids}
    genome_afreq = {sg: f"${genome_dir}/variant/{sg}.info.afreq" for sg in ids}
    genome_acount = {sg: f"${genome_dir}/variant/{sg}.assoc.acount" for sg in ids}
    genome_pvar = {sg: f"${genome_dir}/{sg}_test.plink2.pvar" for sg in ids}
    exon_afreq = {sg: f"${exon_dir}/variant/{sg}.info.afreq" for sg in ids}
    exon_acount = {sg: f"${exon_dir}/variant/{sg}.assoc.acount" for sg in ids}
    exon_pvar = {sg: f"${exon_dir}/{sg}_test.plink2.pvar" for sg in ids}
    for pct, genome_ids, exon_ids in (
        ("80", "${genome80}", "${exon80}"),
        ("60", "${genome60}", "${exon60}"),
    ):
        prefix = f"all.variant.vep_atlas.cv_reliable.keep{pct}"
        print(f"VEP atlas {prefix}", flush=True)
        run_vep_atlas(
            genome_vep,
            exon_vep,
            prefix,
            genome_afreq=genome_afreq,
            genome_acount=genome_acount,
            genome_pvar=genome_pvar,
            exon_afreq=exon_afreq,
            exon_acount=exon_acount,
            exon_pvar=exon_pvar,
            ld_blocks="${ld_blocks}",
            gene_gff="${gene_gff}",
            bin_bp=int("${bin_bp}"),
            genome_keep_ids=load_ids(genome_ids),
            exon_keep_ids=load_ids(exon_ids),
        )
    """
}


