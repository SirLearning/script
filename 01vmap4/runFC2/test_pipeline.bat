@echo off
REM FastCall2 Pipeline Test Script for Windows
REM This script tests the Nextflow pipeline syntax and basic functionality

echo =========================================
echo FastCall2 Pipeline Test Suite - Windows
echo =========================================

REM Test 1: Check if Nextflow is available
echo Test 1: Checking Nextflow installation...
where nextflow >nul 2>nul
if %errorlevel% equ 0 (
    echo [OK] Nextflow found
    REM Test syntax
    nextflow run runFastCall2.nf --help >nul 2>nul
    if %errorlevel% equ 0 (
        echo [OK] Nextflow syntax check passed
    ) else (
        echo [ERROR] Nextflow syntax check failed
        goto :error
    )
) else (
    echo [WARNING] Nextflow not found, please install from https://www.nextflow.io/
)

REM Test 2: Check required files
echo.
echo Test 2: Checking required files...
set "missing_files="
if not exist "runFastCall2.nf" (
    echo [ERROR] runFastCall2.nf missing
    set "missing_files=1"
) else (
    echo [OK] runFastCall2.nf exists
)

if not exist "nextflow.config" (
    echo [ERROR] nextflow.config missing
    set "missing_files=1"
) else (
    echo [OK] nextflow.config exists
)

if not exist "README.md" (
    echo [ERROR] README.md missing
    set "missing_files=1"
) else (
    echo [OK] README.md exists
)

if defined missing_files (
    echo [ERROR] Some required files are missing
    goto :error
) else (
    echo [OK] All required files present
)

REM Test 3: Check help message
echo.
echo Test 3: Checking help message...
nextflow run runFastCall2.nf --help > help_test.tmp 2>&1
findstr "FastCall2 High-throughput Pipeline" help_test.tmp >nul
if %errorlevel% equ 0 (
    echo [OK] Help message working correctly
) else (
    echo [ERROR] Help message not working
)
del help_test.tmp >nul 2>nul

echo.
echo =========================================
echo Test Summary Complete
echo =========================================
echo.
echo To run the actual pipeline, modify the paths in run_example.bat
echo and execute it, or use:
echo.
echo nextflow run runFastCall2.nf ^
echo   --reference C:\path\to\genome.fa ^
echo   --taxaBamMap C:\path\to\taxaBamMap.txt ^
echo   --tiger_jar C:\path\to\TIGER.jar ^
echo   --samtools_path C:\path\to\samtools.exe
echo.
pause
exit /b 0

:error
echo.
echo [ERROR] Tests failed. Please check the pipeline setup.
pause
exit /b 1
