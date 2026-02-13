# Genetics Analysis Pipeline & Script Library

This project contains a complete genetics data analysis pipeline (based on Nextflow) and accompanying Python/R/Java script libraries. The project integrates whole genome sequencing data processing, variant calling, population genetics analysis, and GWAS features.

## 📂 Project Structure

```
.
├── environment_*.yml       # Conda environment configuration files
├── setup.py               # Python package installation script
├── README.md              # Project documentation
├── src/                   # Source code directory
│   ├── python/            # Python core library (python_script)
│   │   ├── genetics/      # Genetics analysis modules (genomics, gwas, phenotype, etc.)
│   │   ├── infra/         # Infrastructure and tools (server, stats, utils)
│   │   └── WeaTE/         # Transposon analysis module
│   ├── r/                 # R statistics and plotting scripts
│   └── java/              # Java tools (e.g., BamHeader, TaxaBamMap)
├── workflow/              # Nextflow workflows
│   ├── Genetics/          # Core genetics analysis pipelines
│   │   ├── genotype/      # Genotype processing (alignment, calling, statistics)
│   │   ├── static/        # Static analysis (GWAS, phenotype)
│   │   ├── dynamic/       # Dynamic/Population analysis (Kinship, XP-CLR)
│   │   └── main.nf        # Main entry script
│   └── ALiYun/            # Alibaba Cloud WDL workflow backup
└── note/                  # Analysis notes (Jupyter Notebooks)
```

## 🛠️ Environment Setup

This project relies on Conda for environment management and provides multiple dedicated environments for different tasks.

### 1. Create Conda Environments

Create the appropriate environment based on your task requirements:

*   **`run`**: Workflow execution environment (Nextflow, Screen)
    ```bash
    conda env create -f environment_run.yml
    ```
*   **`stats`**: Main statistical analysis environment (Python 3.12, Hail, Plink, Samtools, Bcftools)
    *   *Used for script development and most Python analysis tasks*
    ```bash
    conda env create -f environment_stats.yml
    ```
*   **`stats_r`**: R language statistics and plotting environment (R 4.3.1, Tidyverse, BioConductor)
    ```bash
    conda env create -f environment_stats_r.yml
    ```
*   **`tiger`**: Compute-intensive task environment (BWA-MEM2, Samtools)
    *   *Used for alignment and variant calling*
    ```bash
    conda env create -f environment_tiger.yml
    ```
*   **`dbone`**: Database and basic operations environment
    ```bash
    conda env create -f environment_dbone.yml
    ```

### 2. Install Python Package

To enable calls to code under `src/python` within the environment, you need to install this project in editable mode in the **stats** environment (or other environments requiring Python scripts):

```bash
conda activate stats
pip install -e .
```
This will install the `python_script` package and its dependencies (pandas, numpy, scipy, seaborn, hail, etc.).

## 🚀 Workflows

The main workflow scripts are located in the `workflow/Genetics` directory.

### Genotype
Path: `workflow/Genetics/genotype/`

*   **`align.nf`**: Sequence alignment workflow. Includes FASTQ QC, BWA-MEM2 alignment, Samtools sort/dedup/index, and Mosdepth depth calculation.
    *   *New Feature*: Includes smart USB file distribution (`cp_based_on_usb_size`), which can copy files in parallel using multiple threads based on available disk space.
*   **`stats.nf`**: Genotype statistics workflow. Integrates Python scripts to compute and plot metrics such as missing rates, heterozygosity, IBS, King kinship, etc.
*   **`caller.nf`**: Variant calling workflow.
*   **`hail.nf`**: Distributed data processing related to Hail.

### Static (Association Analysis)
Path: `workflow/Genetics/static/`

*   **`gwas/`**: Genome-wide association study workflow.
*   **`phenotype/`**: Phenotype data processing and statistics.

### Dynamic (Population/Evolution Analysis)
Path: `workflow/Genetics/dynamic/`

*   **`kinship.nf`**: Kinship analysis.
*   **`xp_clr.nf`**: Selective sweep analysis.

## 📦 Source Code Library (Src)

### Python (`src/python`)
Core logic is encapsulated in the `python_script` package:
*   **`genetics.genomics`**: Genomics analysis (Sample QC, Variant QC, IBS, Kinship, etc.).
*   **`infra.server.cp`**: Contains smart file copying logic (`run_copy_process`) for use by Nextflow pipelines.
*   **`infra.utils.graph`**: General plotting utility library.
*   **`WeaTE`**: Wheat transposon analysis specialized module.

### R (`src/r`)
Contains R scripts for GWAS result visualization (Manhattan/QQ plots), PCA visualization, etc.

### Java (`src/java`)
Contains utility classes for handling BAM headers (`BamHeader.java`) and Taxa mapping (`TaxaBamMap.java`).

## 📝 Usage Examples

### Run Alignment Workflow
```bash
conda activate run
nextflow run workflow/Genetics/genotype/align.nf \
    --fastq_dir ./raw_data \
    --reference ref.fa \
    -profile tiger
```

### Run Statistics Workflow
```bash
conda activate run
nextflow run workflow/Genetics/genotype/stats.nf \
    --output_dir ./results \
    -profile stats
```
*(Note: Please adjust parameters according to the specific `nextflow.config` file)*
