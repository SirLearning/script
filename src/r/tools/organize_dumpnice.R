#!/usr/bin/env Rscript

# Organize and copy R scripts from note/04DumpNice into src/r/dumpnice by category
# Usage: Rscript src/r/tools/organize_dumpnice.R [source_dir] [dest_root]
# Defaults: source_dir = note/04DumpNice, dest_root = src/r/dumpnice

suppressWarnings(suppressMessages({
  # base only
}))

args <- commandArgs(trailingOnly = TRUE)
repo_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "src/r/tools/organize_dumpnice.R"), "..", "..", ".."), mustWork = FALSE)
# Helper for null-coalesce
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a)) a else b

source_dir <- normalizePath(args[1] %||% file.path(repo_root, "note", "04DumpNice"), mustWork = FALSE)
dest_root  <- normalizePath(args[2] %||% file.path(repo_root, "src", "r", "dumpnice"), mustWork = FALSE)

message("Source: ", source_dir)
message("Dest:   ", dest_root)

if (!dir.exists(source_dir)) {
  stop("Source directory does not exist: ", source_dir)
}
if (!dir.exists(dest_root)) {
  dir.create(dest_root, recursive = TRUE, showWarnings = FALSE)
}

# Category regex map (ordered)
category_map <- list(
  vcf        = "(?i)\\bvcf|Vcf\\b",
  haplotype  = "(?i)haplo",
  gwas       = "(?i)gwas|gapit|rmvp|mvp",
  plot       = "(?i)plot|pheatmap|heatmap|venn|density|qmatrix|barcode|map\\b|geneStructure|chrbin",
  fst        = "(?i)\\bfst\\b",
  ibs        = "(?i)\\bibs\\b",
  maf        = "(?i)\\bmaf\\b",
  pi         = "(?i)\\bpi(?![a-z])",
  ld         = "(?i)\\bld|ldjump",
  admixture  = "(?i)admixture|structure|qmatrix",
  enrichment = "(?i)enrich",
  pca        = "(?i)\\bpca\\b",
  rda        = "(?i)\\brda\\b",
  migration  = "(?i)migration",
  permutation= "(?i)permut",
  depth      = "(?i)depth",
  bwa        = "(?i)\\bbwa\\b",
  tree       = "(?i)\\btree\\b|nucdiff",
  statistic  = "(?i)statistic|basic_statistic|seperationsite|tajimasd|nesize|maf-dis",
  qc         = "(?i)qc",
  pipeline   = "(?i)流程|vf_pipe|vmap3|storage|mcmc|bayenv|lfmm|variants_age|method|testgf|39_42-21chr|scripts\\.Rproj"
)

categorize <- function(path) {
  name <- basename(path)
  for (cat in names(category_map)) {
    if (grepl(category_map[[cat]], name, perl = TRUE)) return(cat)
  }
  # also peek into one-level parent dir
  parent <- basename(dirname(path))
  for (cat in names(category_map)) {
    if (grepl(category_map[[cat]], parent, perl = TRUE)) return(cat)
  }
  return("utils")
}

# List all R scripts
all_files <- list.files(source_dir, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
all_files <- all_files[!grepl("(^|/)\\.", all_files)]  # exclude hidden
if (length(all_files) == 0) {
  message("No R files found under ", source_dir)
  quit(status = 0)
}

# Prepare manifest
manifest <- data.frame(
  source = character(0),
  target = character(0),
  category = character(0),
  md5_source = character(0),
  md5_target = character(0),
  stringsAsFactors = FALSE
)

safe_copy <- function(src, dst_dir) {
  if (!dir.exists(dst_dir)) dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
  base <- basename(src)
  target <- file.path(dst_dir, base)
  get_md5 <- function(f) tryCatch(as.character(tools::md5sum(f)), error = function(e) NA_character_)
  src_md5 <- get_md5(src)
  if (file.exists(target)) {
    # If identical, skip
    if (!is.na(src_md5) && identical(src_md5, get_md5(target))) {
      return(list(path = target, md5 = get_md5(target), skipped = TRUE))
    }
    # Else, add suffix
    stem <- sub("\\.[Rr]$", "", base)
    ext <- sub("^.*(\\.[Rr])$", "\\1", base)
    i <- 1
    repeat {
      cand <- file.path(dst_dir, sprintf("%s__%02d%s", stem, i, ext))
      if (!file.exists(cand)) { target <- cand; break }
      i <- i + 1
      if (i > 99) stop("Too many name collisions for ", base)
    }
  }
  ok <- file.copy(src, target, overwrite = FALSE)
  if (!ok) stop("Failed to copy ", src, " -> ", target)
  list(path = target, md5 = get_md5(target), skipped = FALSE)
}

message(sprintf("Found %d R files. Copying...", length(all_files)))

for (src in all_files) {
  rel <- sub(paste0("^", gsub("\\+", "\\\\+", source_dir)), "", src)
  catg <- categorize(src)
  dst_dir <- file.path(dest_root, catg)
  res <- safe_copy(src, dst_dir)
  manifest[nrow(manifest) + 1, ] <- list(src, res$path, catg, tryCatch(as.character(tools::md5sum(src)), error=function(e) NA), res$md5)
}

# Write manifest
man_path <- file.path(dest_root, "_manifest.csv")
utils::write.csv(manifest, man_path, row.names = FALSE)

# Also write a simple summary
summary_txt <- file.path(dest_root, "_summary.txt")
counts <- sort(table(manifest$category), decreasing = TRUE)
cat("DumpNice organization summary\n", file = summary_txt)
cat("Source:", source_dir, "\nDest:", dest_root, "\n\n", file = summary_txt, append = TRUE)
cat("Counts by category:\n", file = summary_txt, append = TRUE)
for (nm in names(counts)) cat(sprintf("- %s: %d\n", nm, counts[[nm]]), file = summary_txt, append = TRUE)

message("Done. Manifest at ", man_path)
