#!/bin/bash

# FastCall2 Pipeline Test Script
# This script tests the Nextflow pipeline syntax and basic functionality

echo "========================================="
echo "FastCall2 Pipeline Test Suite"
echo "========================================="

# Test 1: Check Nextflow syntax
echo "Test 1: Checking Nextflow syntax..."
if command -v nextflow &> /dev/null; then
    nextflow run runFastCall2.nf --help > /dev/null 2>&1
    if [ $? -eq 0 ]; then
        echo "✓ Nextflow syntax check passed"
    else
        echo "✗ Nextflow syntax check failed"
        exit 1
    fi
else
    echo "⚠ Nextflow not found, skipping syntax check"
fi

# Test 2: Check required parameters validation
echo ""
echo "Test 2: Checking parameter validation..."
nextflow run runFastCall2.nf > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "✓ Parameter validation working correctly"
else
    echo "✗ Parameter validation not working"
fi

# Test 3: Check help message
echo ""
echo "Test 3: Checking help message..."
nextflow run runFastCall2.nf --help > help_output.txt 2>&1
if grep -q "FastCall2 High-throughput Pipeline" help_output.txt; then
    echo "✓ Help message displayed correctly"
    rm help_output.txt
else
    echo "✗ Help message not working"
    rm help_output.txt
fi

# Test 4: Check if all required files exist
echo ""
echo "Test 4: Checking required files..."
required_files=("runFastCall2.nf" "nextflow.config" "README.md")
all_exist=true

for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✓ $file exists"
    else
        echo "✗ $file missing"
        all_exist=false
    fi
done

if [ "$all_exist" = true ]; then
    echo "✓ All required files present"
else
    echo "✗ Some required files are missing"
fi

# Test 5: Check configuration file syntax
echo ""
echo "Test 5: Checking configuration file..."
if [ -f "nextflow.config" ]; then
    # Basic syntax check for the config file
    if grep -q "nextflow.enable.dsl" nextflow.config; then
        echo "✓ Configuration file syntax appears correct"
    else
        echo "⚠ Configuration file may have issues"
    fi
else
    echo "✗ Configuration file missing"
fi

echo ""
echo "========================================="
echo "Test Summary Complete"
echo "========================================="
echo ""
echo "To run the actual pipeline, use:"
echo "nextflow run runFastCall2.nf \\"
echo "  --reference /path/to/genome.fa \\"
echo "  --taxaBamMap /path/to/taxaBamMap.txt \\"
echo "  --tiger_jar /path/to/TIGER.jar \\"
echo "  --samtools_path /path/to/samtools"
