@echo off
REM FastCall2 Pipeline Execution Script for Windows
REM This script demonstrates how to run the FastCall2 Nextflow pipeline

REM Set your paths here (modify these paths according to your system)
set REFERENCE_GENOME=D:\path\to\your\reference\genome.fa
set TAXA_BAM_MAP=D:\path\to\your\taxaBamMap.txt
set TIGER_JAR=D:\path\to\TIGER.jar
set SAMTOOLS_PATH=D:\path\to\samtools.exe
set OUTPUT_DIR=fastcall2_results

REM Check if Nextflow is installed
where nextflow >nul 2>nul
if %errorlevel% neq 0 (
    echo Error: Nextflow is not installed or not in PATH
    echo Please install Nextflow from https://www.nextflow.io/
    pause
    exit /b 1
)

echo ========================================
echo FastCall2 Pipeline - Windows Runner
echo ========================================

REM Basic run command
echo Running FastCall2 pipeline...
nextflow run runFastCall2.nf ^
    --reference %REFERENCE_GENOME% ^
    --taxaBamMap %TAXA_BAM_MAP% ^
    --tiger_jar %TIGER_JAR% ^
    --samtools_path %SAMTOOLS_PATH% ^
    --output_dir %OUTPUT_DIR% ^
    --threads 32 ^
    --memory "100g"

if %errorlevel% neq 0 (
    echo Pipeline execution failed!
    pause
    exit /b 1
)

echo Pipeline completed successfully!
echo Results can be found in: %OUTPUT_DIR%
pause
