# Plotting Scripts for Genetics Pipeline

This directory contains R scripts for plotting various population genetics statistics.

## Scripts

### `plot_fst_boxplot.R`
Plots boxplots of FST values from multiple comparison files.
**Usage:**
```bash
Rscript plot_fst_boxplot.R -i /path/to/fst_files -o output.pdf
```

### `plot_pi_distribution.R`
Plots distribution (boxplot and density) of Nucleotide Diversity (Pi).
**Usage:**
```bash
Rscript plot_pi_distribution.R -i /path/to/pi_files -o output.pdf
```

### `plot_tajima_d.R`
Plots distribution of Tajima's D values.
**Usage:**
```bash
Rscript plot_tajima_d.R -i /path/to/tajima_files -o output.pdf
```

### `plot_ideogram_density.R`
Plots feature density on chromosomes using `RIdeogram`.
**Usage:**
```bash
Rscript plot_ideogram_density.R -k karyotype.txt -d density.txt -o output.svg
```
