# FastCall2 High-throughput Nextflow Pipeline

This pipeline implements a high-throughput version of the FastCall2 variant calling workflow using Nextflow DSL2.

## Overview

FastCall2 is a variant calling pipeline that consists of three main steps:
1. **disc**: Discovery of potential variant sites
2. **blib**: Building of variant libraries
3. **scan**: Scanning and genotyping variants

This Nextflow implementation parallelizes the process by chromosome, allowing for efficient processing of large datasets.

## Features

- **Parallel processing**: Each chromosome is processed independently
- **Scalable**: Can run on local machines, clusters, or cloud environments
- **Configurable**: All FastCall2 parameters can be customized
- **Resume capability**: Can resume interrupted runs
- **Comprehensive reporting**: Generates execution reports and summaries

## Requirements

### Software
- Nextflow (>= 20.04.0)
- Java (for TIGER/FastCall2)
- samtools (version 1.18 recommended)
- bcftools (optional, for VCF merging)

### Input Files
- Reference genome (FASTA format)
- Taxa-BAM mapping file (created using Workshop.jar)
- TIGER jar file containing FastCall2
- BAM files with proper indexing

## Quick Start

1. **Prepare your input files**:
   ```bash
   # Create taxaBamMap.txt using Workshop.jar
   java -Xmx100g -jar Workshop.jar -d /path/to/depth/files -b /path/to/bam/files -o taxaBamMap.txt
   ```

2. **Run the pipeline**:
   ```bash
   nextflow run runFastCall2.nf \
     --reference genome.fa \
     --taxaBamMap taxaBamMap.txt \
     --tiger_jar TIGER.jar \
     --samtools_path /path/to/samtools
   ```

3. **View results**:
   Results will be saved in the `fastcall2_output` directory with the following structure:
   ```
   fastcall2_output/
   ├── disc/          # Discovery files (.ing)
   ├── blib/          # Library files (.lib.gz)
   ├── scan/          # VCF files
   ├── final/         # Merged results and summary
   └── pipeline_info/ # Execution reports
   ```

## Parameters

### Required Parameters
- `--reference`: Reference genome file (FASTA)
- `--taxaBamMap`: Taxa-BAM mapping file
- `--tiger_jar`: Path to TIGER jar file
- `--samtools_path`: Path to samtools executable

### Optional Parameters
- `--output_dir`: Output directory (default: fastcall2_output)
- `--threads`: Number of threads (default: 32)
- `--memory`: Memory allocation (default: 100g)
- `--chromosomes`: List of chromosomes to process

### FastCall2 Disc Parameters
- `--disc_min_depth`: Minimum depth (default: 30)
- `--disc_min_qual`: Minimum quality (default: 20)
- `--disc_min_snp_count`: Minimum SNP count (default: 2)
- `--disc_min_allele_freq`: Minimum allele frequency (default: 0.2)
- `--disc_min_coverage`: Minimum coverage (default: 3)
- `--disc_max_missing`: Maximum missing rate (default: 0.8)
- `--disc_min_het_freq`: Minimum heterozygous frequency (default: 0.35)
- `--disc_max_het_freq`: Maximum heterozygous frequency (default: 0.2)

### FastCall2 Scan Parameters
- `--scan_min_depth`: Minimum depth (default: 30)
- `--scan_min_qual`: Minimum quality (default: 20)
- `--scan_p_value`: P-value threshold (default: 0.05)

## Execution Profiles

The pipeline includes several execution profiles:

### Standard (Local)
```bash
nextflow run runFastCall2.nf --reference genome.fa ... -profile standard
```

### High Performance Computing (SLURM)
```bash
nextflow run runFastCall2.nf --reference genome.fa ... -profile hpc
```

### Test (Reduced resources)
```bash
nextflow run runFastCall2.nf --reference genome.fa ... -profile test
```

## Advanced Usage

### Custom Chromosome List
```bash
nextflow run runFastCall2.nf \
  --chromosomes "1A,1B,1D,2A,2B,2D" \
  --reference genome.fa \
  --taxaBamMap taxaBamMap.txt \
  --tiger_jar TIGER.jar \
  --samtools_path samtools
```

### Resume Failed Run
```bash
nextflow run runFastCall2.nf \
  --reference genome.fa \
  --taxaBamMap taxaBamMap.txt \
  --tiger_jar TIGER.jar \
  --samtools_path samtools \
  -resume
```

### Custom Parameters
```bash
nextflow run runFastCall2.nf \
  --reference genome.fa \
  --taxaBamMap taxaBamMap.txt \
  --tiger_jar TIGER.jar \
  --samtools_path samtools \
  --disc_min_depth 20 \
  --scan_p_value 0.01 \
  --threads 48 \
  --memory "200g"
```

## Output Structure

```
fastcall2_output/
├── disc/
│   ├── 1A_1_***.ing
│   ├── 1B_1_***.ing
│   └── ...
├── blib/
│   ├── 1A_1_***.lib.gz
│   ├── 1B_1_***.lib.gz
│   └── ...
├── scan/
│   ├── *.vcf or *.vcf.gz files
│   └── scan_*.log
├── final/
│   ├── merged_variants.vcf.gz
│   └── summary_stats.txt
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    ├── trace.txt
    └── dag.dot
```

## Troubleshooting

### Common Issues

1. **Java memory errors**: Increase `--memory` parameter
2. **Missing BAM index files**: Ensure all BAM files are indexed
3. **Chromosome naming**: Verify chromosome names match between reference and BAM files
4. **Permission errors**: Check file permissions and paths

### Performance Optimization

1. **For large datasets**: Use the HPC profile
2. **For testing**: Use the test profile with reduced chromosomes
3. **Memory allocation**: Adjust based on your system capabilities
4. **Thread count**: Set to match your CPU cores

## Citation

If you use this pipeline, please cite:
- FastCall2/TIGER: [Original FastCall2 publication]
- Nextflow: Di Tommaso, P. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316–319 (2017).

## Support

For issues related to:
- FastCall2 parameters: Check the TIGER/FastCall2 documentation
- Nextflow execution: Check Nextflow documentation
- Pipeline bugs: Submit an issue to this repository
