# Load all DumpNice R scripts (categorized) into the current R session
# Usage:
#   source('src/r/load_dumpnice.R')
# Optionally set options(dumpnice.verbose=TRUE) to print sourced files.

load_dumpnice <- function(root = 'src/r/dumpnice') {
  verbose <- isTRUE(getOption('dumpnice.verbose'))
  if (!dir.exists(root)) {
    stop('DumpNice folder not found: ', root)
  }
  files <- list.files(root, pattern = '\\.[Rr]$', recursive = TRUE, full.names = TRUE)
  # Exclude this loader and any files under tools if present
  files <- files[!grepl('load_dumpnice\\.R$|/tools/', files)]
  # Sort by category path then name for determinism
  files <- files[order(dirname(files), basename(files))]
  if (verbose) message('Sourcing ', length(files), ' files from ', root)
  for (f in files) {
    if (verbose) message(' - ', f)
    sys.source(f, envir = .GlobalEnv)
  }
  invisible(files)
}

# Auto-execute if called directly via source()
if (identical(environment(), globalenv())) {
  try(load_dumpnice(), silent = TRUE)
}
