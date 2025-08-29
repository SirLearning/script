#!/bin/bash

# Vmap4 Data Processing Pipeline Runner
# Based on methods from 00seq_data_process.ipynb
# This script demonstrates how to run the complete workflow

set -e  # Exit on any error

echo "======================================"
echo "Vmap4 Data Processing Pipeline"
echo "======================================"
echo "Date: $(date)"
echo ""

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SCRIPT_DIR}/.."
NEXTFLOW_CONFIG="${SCRIPT_DIR}/nextflow.config"

# Default parameters (can be overridden by command line arguments)
WORKFLOW_TYPE="alignment"  # Options: data_processing, alignment, depth_qc, fastcall2, performance
OUTPUT_DIR="./vmap4_output"
THREADS=20
MEMORY="8G"
PROFILE="standard"

# Input files (modify these paths according to your setup)
SRA_LIST=""
FASTQ_DIR=""
REFERENCE=""
BAM_DIR=""
BAM_LIST=""
TAXA_BAM_MAP=""
TIGER_JAR=""
SAMTOOLS_PATH=""

# Function to display help
show_help() {
    cat << EOF
Vmap4 Data Processing Pipeline Runner

Usage: $0 [OPTIONS]

WORKFLOW OPTIONS:
    -w, --workflow WORKFLOW     Workflow type (data_processing, alignment, depth_qc, fastcall2, performance)
    -o, --output OUTPUT_DIR     Output directory (default: ./vmap4_output)
    -t, --threads THREADS       Number of threads (default: 20)
    -m, --memory MEMORY         Memory allocation (default: 8G)
    -p, --profile PROFILE       Nextflow profile (default: standard)

INPUT FILES:
    --sra-list FILE            SRA accession list (for data_processing)
    --fastq-dir DIR            FASTQ directory (for alignment)
    --reference FILE           Reference genome file
    --bam-dir DIR              BAM directory (for depth_qc)
    --bam-list FILE            BAM list file (for depth_qc)
    --taxa-bam-map FILE        Taxa-BAM mapping file (for fastcall2)
    --tiger-jar FILE           TIGER jar file (for fastcall2)
    --samtools-path PATH       Samtools executable path (for fastcall2)

EXAMPLES:
    # Data processing workflow
    $0 -w data_processing --sra-list sra_accessions.txt

    # Alignment workflow
    $0 -w alignment --fastq-dir ./fastq --reference genome.fa --bam-list samples.txt

    # Depth and QC workflow
    $0 -w depth_qc --bam-dir ./bam --bam-list bam_files.txt

    # FastCall2 variant calling
    $0 -w fastcall2 --reference genome.fa --taxa-bam-map taxa.txt --tiger-jar TIGER.jar --samtools-path samtools

    # Performance analysis
    $0 -w performance --reference genome.fa --taxa-bam-map taxa.txt --tiger-jar TIGER.jar --samtools-path samtools

EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -w|--workflow)
            WORKFLOW_TYPE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -m|--memory)
            MEMORY="$2"
            shift 2
            ;;
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        --sra-list)
            SRA_LIST="$2"
            shift 2
            ;;
        --fastq-dir)
            FASTQ_DIR="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --bam-dir)
            BAM_DIR="$2"
            shift 2
            ;;
        --bam-list)
            BAM_LIST="$2"
            shift 2
            ;;
        --taxa-bam-map)
            TAXA_BAM_MAP="$2"
            shift 2
            ;;
        --tiger-jar)
            TIGER_JAR="$2"
            shift 2
            ;;
        --samtools-path)
            SAMTOOLS_PATH="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Validate workflow type
case $WORKFLOW_TYPE in
    data_processing|alignment|depth_qc|fastcall2|performance)
        ;;
    *)
        echo "Error: Invalid workflow type '$WORKFLOW_TYPE'"
        echo "Valid options: data_processing, alignment, depth_qc, fastcall2, performance"
        exit 1
        ;;
esac

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Configuration:"
echo "  Workflow type: $WORKFLOW_TYPE"
echo "  Output directory: $OUTPUT_DIR"
echo "  Threads: $THREADS"
echo "  Memory: $MEMORY"
echo "  Profile: $PROFILE"
echo ""

# Set Nextflow command based on workflow type
case $WORKFLOW_TYPE in
    data_processing)
        WORKFLOW_FILE="${BASE_DIR}/DataProcess/overall/data_processing.nf"
        NEXTFLOW_PARAMS="--output_dir $OUTPUT_DIR --threads $THREADS --memory $MEMORY"
        
        if [[ -n "$SRA_LIST" ]]; then
            NEXTFLOW_PARAMS="$NEXTFLOW_PARAMS --sra_list $SRA_LIST"
        else
            echo "Error: SRA list is required for data processing workflow"
            exit 1
        fi
        ;;
        
    alignment)
        WORKFLOW_FILE="${BASE_DIR}/DataProcess/alignment/sequence_alignment.nf"
        NEXTFLOW_PARAMS="--output_dir $OUTPUT_DIR --threads $THREADS --memory $MEMORY"
        
        if [[ -n "$FASTQ_DIR" && -n "$REFERENCE" && -n "$BAM_LIST" ]]; then
            NEXTFLOW_PARAMS="$NEXTFLOW_PARAMS --fastq_dir $FASTQ_DIR --reference $REFERENCE --fq_list $BAM_LIST --sample_list $BAM_LIST"
        else
            echo "Error: FASTQ directory, reference genome, and sample list are required for alignment workflow"
            exit 1
        fi
        ;;
        
    depth_qc)
        WORKFLOW_FILE="${BASE_DIR}/DataProcess/overall/depth_qc.nf"
        NEXTFLOW_PARAMS="--output_dir $OUTPUT_DIR --threads $THREADS --memory $MEMORY"
        
        if [[ -n "$BAM_DIR" && -n "$BAM_LIST" ]]; then
            NEXTFLOW_PARAMS="$NEXTFLOW_PARAMS --bam_dir $BAM_DIR --bam_list $BAM_LIST"
        else
            echo "Error: BAM directory and BAM list are required for depth QC workflow"
            exit 1
        fi
        ;;
        
    fastcall2)
        WORKFLOW_FILE="${BASE_DIR}/DataProcess/calling/run/runFastCall2.nf"
        NEXTFLOW_PARAMS="--output_dir $OUTPUT_DIR --threads $THREADS --memory $MEMORY"
        
        if [[ -n "$REFERENCE" && -n "$TAXA_BAM_MAP" && -n "$TIGER_JAR" && -n "$SAMTOOLS_PATH" ]]; then
            NEXTFLOW_PARAMS="$NEXTFLOW_PARAMS --reference $REFERENCE --taxaBamMap $TAXA_BAM_MAP --tiger_jar $TIGER_JAR --samtools_path $SAMTOOLS_PATH"
        else
            echo "Error: Reference, taxa-BAM map, TIGER jar, and samtools path are required for FastCall2 workflow"
            exit 1
        fi
        ;;
        
    performance)
        WORKFLOW_FILE="${BASE_DIR}/DataProcess/calling/perf/fastcall.nf"
        NEXTFLOW_PARAMS="--output_dir $OUTPUT_DIR --threads $THREADS --memory $MEMORY --benchmark_mode true"
        
        if [[ -n "$REFERENCE" && -n "$TAXA_BAM_MAP" && -n "$TIGER_JAR" && -n "$SAMTOOLS_PATH" ]]; then
            NEXTFLOW_PARAMS="$NEXTFLOW_PARAMS --reference $REFERENCE --taxaBamMap $TAXA_BAM_MAP --tiger_jar $TIGER_JAR --samtools_path $SAMTOOLS_PATH"
        else
            echo "Error: Reference, taxa-BAM map, TIGER jar, and samtools path are required for performance workflow"
            exit 1
        fi
        ;;
esac

# Check if workflow file exists
if [[ ! -f "$WORKFLOW_FILE" ]]; then
    echo "Error: Workflow file not found: $WORKFLOW_FILE"
    exit 1
fi

# Check if Nextflow is available
if ! command -v nextflow &> /dev/null; then
    echo "Error: Nextflow is not installed or not in PATH"
    echo "Please install Nextflow: https://www.nextflow.io/docs/latest/getstarted.html"
    exit 1
fi

echo "Starting workflow..."
echo "Workflow file: $WORKFLOW_FILE"
echo "Parameters: $NEXTFLOW_PARAMS"
echo ""

# Construct and run Nextflow command
NEXTFLOW_CMD="nextflow run $WORKFLOW_FILE"

if [[ -f "$NEXTFLOW_CONFIG" ]]; then
    NEXTFLOW_CMD="$NEXTFLOW_CMD -c $NEXTFLOW_CONFIG"
fi

NEXTFLOW_CMD="$NEXTFLOW_CMD -profile $PROFILE $NEXTFLOW_PARAMS"

echo "Executing: $NEXTFLOW_CMD"
echo ""

# Run the workflow
eval $NEXTFLOW_CMD

# Check exit status
if [[ $? -eq 0 ]]; then
    echo ""
    echo "======================================"
    echo "Workflow completed successfully!"
    echo "======================================"
    echo "Output directory: $OUTPUT_DIR"
    echo "Check the results and reports in the output directory."
    
    # Show output summary
    if [[ -d "$OUTPUT_DIR" ]]; then
        echo ""
        echo "Output structure:"
        find "$OUTPUT_DIR" -type d -maxdepth 2 | sort
    fi
else
    echo ""
    echo "======================================"
    echo "Workflow failed!"
    echo "======================================"
    echo "Check the error messages above and the Nextflow log files."
    exit 1
fi
