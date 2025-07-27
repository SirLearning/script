#!/bin/bash

# FastCall2 Pipeline Execution Script
# This script demonstrates how to run the FastCall2 Nextflow pipeline

# Set your paths here
REFERENCE_GENOME="/path/to/your/reference/genome.fa"
TAXA_BAM_MAP="/path/to/your/taxaBamMap.txt"
TIGER_JAR="/path/to/TIGER.jar"
SAMTOOLS_PATH="/path/to/samtools"
OUTPUT_DIR="fastcall2_results"

# Basic run command
echo "Running FastCall2 pipeline..."
nextflow run runFastCall2.nf \
    --reference $REFERENCE_GENOME \
    --taxaBamMap $TAXA_BAM_MAP \
    --tiger_jar $TIGER_JAR \
    --samtools_path $SAMTOOLS_PATH \
    --output_dir $OUTPUT_DIR \
    --threads 32 \
    --memory "100g"

# Alternative: Run with custom parameters
echo "Running with custom parameters..."
nextflow run runFastCall2.nf \
    --reference $REFERENCE_GENOME \
    --taxaBamMap $TAXA_BAM_MAP \
    --tiger_jar $TIGER_JAR \
    --samtools_path $SAMTOOLS_PATH \
    --output_dir $OUTPUT_DIR \
    --disc_min_depth 20 \
    --disc_min_qual 15 \
    --scan_p_value 0.01 \
    --chromosomes "1A,1B,1D,2A,2B,2D" \
    -profile hpc

# Test run with reduced dataset
echo "Running test with reduced dataset..."
nextflow run runFastCall2.nf \
    --reference $REFERENCE_GENOME \
    --taxaBamMap $TAXA_BAM_MAP \
    --tiger_jar $TIGER_JAR \
    --samtools_path $SAMTOOLS_PATH \
    --output_dir test_results \
    -profile test

# Resume a failed run
echo "Resuming a previous run..."
nextflow run runFastCall2.nf \
    --reference $REFERENCE_GENOME \
    --taxaBamMap $TAXA_BAM_MAP \
    --tiger_jar $TIGER_JAR \
    --samtools_path $SAMTOOLS_PATH \
    --output_dir $OUTPUT_DIR \
    -resume

# Show help
echo "Showing help message..."
nextflow run runFastCall2.nf --help
