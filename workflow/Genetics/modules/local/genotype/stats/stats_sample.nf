nextflow.enable.dsl=2

process sample_missing_stats {
    tag "compute missing rate threshold: ${chr}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(smiss)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.smiss_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import ana_sample_missing

    print(f"Processing sample missing rate for ${chr}...")
    ana_sample_missing("${smiss}", "${id}.sample.miss")
    """
}

process sample_coverage_stats {
    tag "compute coverage stats: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path group_file
    tuple val(id), val(chr), path(smiss)
    tuple val(id2), val(chr2), path(tbm)
    tuple val(id3), val(chr3), path(scount)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.coverage_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import (
        coverage_dist,
        reg_missing_coverage,
        reg_coverage_ref_ibs,
        heatmap_ibs_depth_missing,
        sample_coverage_stats_bundle,
    )
    from genetics.germplasm.sample.mosdepth_sg_depth import sg_depth_taxa_bam_map_path

    print(f"Processing coverage for ${chr}...")
    sample_coverage_stats_bundle(
        "${tbm}",
        "${smiss}",
        "${scount}",
        "${id}.sample.cov",
        group_file="${group_file}",
    )

    sg_tbm = sg_depth_taxa_bam_map_path("${params.home_dir}", "${id}")
    print(f"Processing subgenome mosdepth coverage for ${chr} (sg_tbm={sg_tbm})...")
    sample_coverage_stats_bundle(
        str(sg_tbm),
        "${smiss}",
        "${scount}",
        "${id}.sample.sg_cov",
        group_file="${group_file}",
    )
    """
}

process sample_sg_coverage_stats {
    tag "compute sg_cov stats: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path group_file
    tuple val(id), val(chr), path(smiss)
    tuple val(id2), val(chr2), path(scount)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("${chr}.sg_coverage_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import sample_coverage_stats_bundle
    from genetics.germplasm.sample.mosdepth_sg_depth import sg_depth_taxa_bam_map_path

    sg_tbm = sg_depth_taxa_bam_map_path("${params.home_dir}", "${id}")
    print(f"Processing subgenome mosdepth coverage for ${chr} (sg_tbm={sg_tbm})...")
    sample_coverage_stats_bundle(
        str(sg_tbm),
        "${smiss}",
        "${scount}",
        "${id}.sample.sg_cov",
        group_file="${group_file}",
    )
    """
}

process plot_subgenome_gam_ibs_depth_compare {
    tag "A/B/D GAM compare panels"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path ibs_depth_info_files

    output:
    path("ABD.sample.cov.gam_subgenome_summary.info.tsv"), emit: info
    path("ABD.sample.cov.gam_*.panels.png"), emit: plots
    path("ABD.sample.cov.gam_compare.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import glob
    import sys

    sys.stdout = open("ABD.sample.cov.gam_compare.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample.cov import compare_subgenome_ibs_depth_gam

    paths = sorted(glob.glob("*.ibs_depth_miss.info.tsv"))
    print(f"Subgenome GAM compare inputs: {paths}")
    compare_subgenome_ibs_depth_gam(paths, output_prefix="ABD.sample.cov")
    """
}

process plot_subgenome_gam_ibs_depth_compare_sg {
    tag "A/B/D GAM compare panels (sg_cov)"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path ibs_depth_info_files

    output:
    path("ABD.sample.sg_cov.gam_subgenome_summary.info.tsv"), emit: info
    path("ABD.sample.sg_cov.gam_*.panels.png"), emit: plots
    path("ABD.sample.sg_cov.gam_compare.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import glob
    import sys

    sys.stdout = open("ABD.sample.sg_cov.gam_compare.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample.cov import compare_subgenome_ibs_depth_gam

    paths = sorted(glob.glob("*.ibs_depth_miss.info.tsv"))
    print(f"Subgenome GAM compare (sg_cov) inputs: {paths}")
    compare_subgenome_ibs_depth_gam(paths, output_prefix="ABD.sample.sg_cov")
    """
}

process sample_hybrid_thinmiss_commonibs_heatmap {
    tag "hybrid thinmiss+commonibs: ${mod} (${sample_topic})"
    publishDir "${params.output_dir}/${params.job}/stats/${params.thin_stats_mod}/info", mode: 'copy', pattern: "*.hybrid_thinmiss_commonibs.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.thin_stats_mod}/plots", mode: 'copy', pattern: "*.heatmap_thinmiss_commonibs*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.thin_stats_mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(mod), val(sample_topic), val(thin_ibs_info), val(common_ibs_info)

    output:
    path("${mod}.sample.${sample_topic}.hybrid_thinmiss_commonibs.info.tsv"), emit: info
    path("${mod}.sample.${sample_topic}.heatmap_thinmiss_commonibs*.png"), emit: plots
    path("${mod}.sample.${sample_topic}.hybrid_thinmiss_commonibs.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("${mod}.sample.${sample_topic}.hybrid_thinmiss_commonibs.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample.cov import heatmap_thin_miss_common_ibs

    print("Hybrid heatmap: thin=${params.thin_stats_mod}, common=${params.common_stats_mod}, mod=${mod}, topic=${sample_topic}")
    heatmap_thin_miss_common_ibs(
        "${thin_ibs_info}",
        "${common_ibs_info}",
        output_prefix="${mod}.sample.${sample_topic}",
        save_info=True,
    )
    """
}

process sample_gam_residual_outlier_plots {
    tag "GAM residual outliers: ${mod} (${sample_topic})"
    publishDir "${params.output_dir}/${params.job}/stats/${params.stats_mod}/info", mode: 'copy', pattern: "*.gam_residual.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.stats_mod}/plots", mode: 'copy', pattern: "*.gam_residual*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.stats_mod}/logs", mode: 'copy', pattern: "*.gam_residual.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(mod), val(sample_topic), val(stats_mod), val(ibs_depth_miss_info)

    output:
    path("${mod}.sample.${sample_topic}.gam_residual.info.tsv"), emit: info
    path("${mod}.sample.${sample_topic}.gam_residual*.png"), emit: plots
    path("${mod}.sample.${sample_topic}.gam_residual.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("${mod}.sample.${sample_topic}.gam_residual.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample.cov import gam_residual_outlier_diagnostics_from_info

    print("GAM residual outliers: mod=${mod}, topic=${sample_topic}, info=${ibs_depth_miss_info}")
    gam_residual_outlier_diagnostics_from_info(
        "${ibs_depth_miss_info}",
        output_prefix="${mod}.sample.${sample_topic}",
        outlier_frac=${params.gam_residual_outlier_frac},
    )
    """
}

process sample_gam_residual_scatter_extras {
    tag "GAM residual scatter extras: ${mod} (${sample_topic})"
    publishDir "${params.output_dir}/${params.job}/stats/${params.stats_mod}/plots", mode: 'copy', pattern: "*.gam_residual_outliers_vs_*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.stats_mod}/plots", mode: 'copy', pattern: "*.gam_residual.dist.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.stats_mod}/plots", mode: 'copy', pattern: "*.gam_residual.dist.logy.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.stats_mod}/logs", mode: 'copy', pattern: "*.gam_residual_scatter_extras.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(mod), val(sample_topic), val(stats_mod), val(gam_residual_info)

    output:
    path("${mod}.sample.${sample_topic}.gam_residual_outliers_vs_*.png"), emit: plots
    path("${mod}.sample.${sample_topic}.gam_residual.dist.png"), emit: dist
    path("${mod}.sample.${sample_topic}.gam_residual.dist.logy.png"), emit: dist_logy
    path("${mod}.sample.${sample_topic}.gam_residual_scatter_extras.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("${mod}.sample.${sample_topic}.gam_residual_scatter_extras.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample.cov import gam_residual_scatter_extras_from_info

    print("GAM residual scatter extras: mod=${mod}, info=${gam_residual_info}")
    gam_residual_scatter_extras_from_info(
        "${gam_residual_info}",
        output_prefix="${mod}.sample.${sample_topic}",
    )
    """
}

process sample_heterozygosity_stats {
    tag "compute heterozygosity rate: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(scount)
    tuple val(id2), val(chr2), path(smiss)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.het_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import ana_heterozygosity

    print(f"Processing heterozygosity rate for ${chr}...")
    ana_heterozygosity("${scount}", "${smiss}", "${id}.sample.het")
    """
}

process sample_mapping_rate_stats {
    tag "compute mapping rate stats: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path idxstats
    path group_file
    tuple val(id), val(chr), path(smiss)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.mr_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import reg_missing_idx_mapping

    print(f"Processing mapping rate for ${chr}...")
    reg_missing_idx_mapping("${idxstats}", "${smiss}", "${id}.sample.idx.mr", "${group_file}")
    """
}

process sample_ref_ibs_stats {
    tag "compute ref IBS stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path group_file
    path idxstats
    tuple val(id), val(chr), path(scount)
    tuple val(id2), val(chr2), path(smiss)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.ref_ibs_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import safe_ref_ibs_stats

    print(f"Processing reference IBS stats for ${chr}...")
    safe_ref_ibs_stats(
        mapping_file="${idxstats}",
        scount_file="${scount}",
        missing_file="${smiss}",
        group_file="${group_file}",
        output_prefix="${id}.sample.ref_ibs",
    )
    """
}

process sample_king_stats {
    tag "compute king stats: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/thresholds", mode: 'copy', pattern: "*.stats.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(king_file)
    tuple val(id2), val(chr2), path(scount)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.stats.tsv"), emit: stats
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.king_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import ana_king_kinship, ana_derived_het, ana_het_vs_max_kinship

    print(f"Processing KING Kinship for ${chr}...")
    ana_king_kinship("${king_file}", "${id}.sample.king")
    
    print(f"Processing Derived Heterozygosity for ${chr}...")
    ana_derived_het("${king_file}", "${id}.sample.derived_het")
    
    print(f"Processing Het vs Max Kinship for ${chr}...")
    ana_het_vs_max_kinship("${scount}", "${king_file}", "${id}.sample.het_vs_kin")
    """
}

process sample_ibs_stats {
    tag "compute ibs stats: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path group_file
    tuple val(id), val(chr), path(ibs_matrix)
    tuple val(id2), val(chr2), path(ids_file)
    tuple val(id3), val(chr3), path(scount)

    output:
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.png"), emit: plots
    tuple val(id), val(chr), path("*.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.ibs_ana.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample import ana_ibs_matrix, ana_het_vs_max_ibs, ana_ibs_trend

    print(f"Processing IBS Matrix Analysis for ${chr}...")
    # matrix_file, id_file, group_file, output_prefix
    ana_ibs_matrix("${ibs_matrix}", "${ids_file}", "${group_file}", "${id}.sample.ibs")
    
    print(f"Processing Het vs Max IBS for ${chr}...")
    # scount_file, matrix_file, id_file, output_prefix
    ana_het_vs_max_ibs("${scount}", "${ibs_matrix}", "${ids_file}", "${id}.sample.het_vs_ibs")

    print(f"Processing IBS Trend for ${chr}...")
    # matrix_file, id_file, group_file, output_prefix
    ana_ibs_trend("${ibs_matrix}", "${ids_file}", "${group_file}", "${id}.sample.ibs_trend")
    """
}

process sample_germplasm_dedup {
    tag "detect germplasm duplicates: ${id}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), path(smiss)
    path tbm_dir
    path dbone_dir

    output:
    tuple val(id), val(chr), path("${id}.sample.dedup.rm.th.tsv"), emit: th
    tuple val(id), val(chr), path("*.info.tsv"), emit: info
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${chr}.germplasm.duplicates.log", "w")
    sys.stderr = sys.stdout

    from genetics.germplasm.sample import safe_ana_duplication

    print(f"Detecting germplasm duplicates for ${chr}...")
    safe_ana_duplication("${tbm_dir}", "${dbone_dir}", "${id}.sample.dedup")
    """
}

process plink_pca {
    tag "${id}"
    publishDir "${params.output_dir}/${params.job}/stats/pca", mode: 'copy'

    input:
    tuple val(id), path(bed), path(bim), path(fam)

    output:
    tuple val(id), path("${id}.eigenvec"), path("${id}.eigenval"), emit: pca

    script:
    """
    set -euo pipefail
    exec > ${id}.plink_pca.log 2>&1

    plink --bfile ${id} \\
        --pca ${params.pc_num} \\
        --out ${id} \\
        --allow-extra-chr \\
        --chr-set 42 \\
        --threads ${task.cpus}
    """
}

process PLOT_PLINK_PCA {
    tag "${id}"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy'

    input:
    tuple val(id), path(eigenvec), path(eigenval)
    path metadata

    output:
    tuple val(id), path("${id}.pca_plot.pdf"), emit: pca_plot

    script:
    def md_arg = metadata.name != 'NO_METADATA' ? "--metadata ${metadata}" : ""
    """
    set -euo pipefail
    Rscript ${params.src_dir}/r/genetics/plot_pca.r \\
        --eigenvec ${eigenvec} \\
        --eigenval ${eigenval} \\
        ${md_arg} \\
        --output ${id}.pca_plot.pdf
    """
}

process plot_thin_common_sample_missing_compare {
    tag "thin vs common-thin sample missing slopes"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/thin_common_compare/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/thin_common_compare/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/thin_common_compare/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(thin_process_dir), val(compare_process_dir), val(output_prefix), val(left_label), val(right_label)

    output:
    path("${output_prefix}.sample.missing_slope.info.tsv"), emit: info
    path("${output_prefix}.sample.missing_slope.A.png"), emit: missing_slope_a
    path("${output_prefix}.sample.missing_slope.B.png"), emit: missing_slope_b
    path("${output_prefix}.sample.missing_slope.D.png"), emit: missing_slope_d
    path("${output_prefix}.sample.missing_slope.Others.png"), emit: missing_slope_others
    path("${output_prefix}.thin_common_miss_compare.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${output_prefix}.thin_common_miss_compare.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample.miss import plot_thin_common_sample_missing_slope

    group_file = "${params.output_dir}/meta_data/sample_group.txt"
    print("Plotting ${left_label} vs ${right_label} sample missing-rate slopes (A / B / D / Others) ...")
    rows = plot_thin_common_sample_missing_slope(
        "${thin_process_dir}",
        "${compare_process_dir}",
        "${output_prefix}",
        group_file=group_file,
        left_label="${left_label}",
        right_label="${right_label}",
    )
    print(f"Subgenomes plotted: {len(rows)}")
    """
}
