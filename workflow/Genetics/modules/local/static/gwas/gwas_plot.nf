nextflow.enable.dsl=2

process plot_gwas_association {
    tag "plot gwas"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/plots", mode: 'copy', pattern: "*.png"

    input:
    path association
    val source
    val output_prefix

    output:
    path "*.tsv"
    path "*.png"

    script:
    """
    #!/usr/bin/env python
    from genetics.gwas.association_plot import plot_gwas_results
    plot_gwas_results("${association}", "${output_prefix}", source="${source}")
    """
}

process plot_wap_association {
    tag "plot wap assoc: ${output_prefix}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"

    input:
    tuple path(association), val(output_prefix)

    output:
    path "*.tsv", emit: info
    path "*.png", emit: plots

    script:
    """
    #!/usr/bin/env python
    import sys

    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.gwas.association_plot import plot_gwas_results

    plot_gwas_results("${association}", "${output_prefix}", source="plink2")
    """
}

process plot_wap_burden {
    tag "plot wap burden: ${output_prefix}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple path(burden), val(output_prefix)
    path gff

    output:
    path "*.tsv", emit: info
    path "*.png", emit: plots
    path "*.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("${output_prefix}.burden_plot.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.gwas.association_plot import plot_burden_results

    plot_burden_results("${burden}", "${output_prefix}", "${gff}")
    """
}

process plot_gwas {
    tag "${res.baseName}"
    publishDir "${params.outdir}/post", mode: 'copy'


    input:
    path res

    output:
    path "*.manhattan.png", emit: plots
    path "*.qq.png"
    path "*.log"

    script:
    def base = res.baseName
    """
    echo "Plotting Manhattan and QQ for ${res}" > ${base}.log

    Rscript - <<'RSCRIPT' >> ${base}.log 2>&1
    # auto-install qqman if missing (best-effort)
    if (!requireNamespace("qqman", quietly = TRUE)) {
      tryCatch({
        install.packages("qqman", repos = "https://cloud.r-project.org")
      }, error = function(e) cat("qqman install failed:", conditionMessage(e), "\n"))
    }
    suppressPackageStartupMessages(library(qqman))

    infile <- "${res}"
    base <- tools::file_path_sans_ext(basename(infile))

    # Try to detect PLINK *.glm.linear format
    df <- tryCatch({
      read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
    }, error = function(e) NULL)

    ok <- !is.null(df) && all(c("CHR", "BP") %in% names(df)) && any(grepl("^P\$|^P_BOLT|^P\\.|^PVAL\$", names(df)))
    if (isTRUE(ok)) {
      pcol <- intersect(c("P", "P_BOLT", "P.", "PVAL"), names(df))[1]
      # basic cleanup: ensure numeric
      df\$CHR <- as.integer(df\$CHR)
      df\$BP  <- as.integer(df\$BP)
      df\$SNP <- if ("SNP" %in% names(df)) df\$SNP else sprintf("%s:%s", df\$CHR, df\$BP)
      df[[pcol]] <- suppressWarnings(as.numeric(df[[pcol]]))

      png(sprintf("%s.manhattan.png", base), width = 1400, height = 600)
      manhattan(df, chr = "CHR", bp = "BP", p = pcol, snp = "SNP", main = base)
      dev.off()

      png(sprintf("%s.qq.png", base), width = 800, height = 800)
      qq(df[[pcol]], main = paste0("QQ: ", base))
      dev.off()
    } else {
      cat("Unsupported or unreadable result format for:", infile, "\n")
      # create placeholder empty plots to keep pipeline consistent
      png(sprintf("%s.manhattan.png", base)); plot.new(); title(main = paste("No plot for", base)); dev.off()
      png(sprintf("%s.qq.png", base)); plot.new(); title(main = paste("No QQ for", base)); dev.off()
    }
    RSCRIPT
    """
}

