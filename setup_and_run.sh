#!/bin/bash

echo "=========================================="
echo "APOE Analysis - Complete Setup"
echo "=========================================="
echo ""

# Check if tools are installed
echo "Checking for required tools..."
echo ""

TOOLS_MISSING=false

if ! command -v bwa &> /dev/null; then
    echo "✗ BWA not installed"
    TOOLS_MISSING=true
else
    echo "✓ BWA installed"
fi

if ! command -v samtools &> /dev/null; then
    echo "✗ SAMtools not installed"
    TOOLS_MISSING=true
else
    echo "✓ SAMtools installed"
fi

if ! command -v bcftools &> /dev/null; then
    echo "✗ BCFtools not installed"
    TOOLS_MISSING=true
else
    echo "✓ BCFtools installed"
fi

if ! command -v python3 &> /dev/null; then
    echo "✗ Python3 not installed"
    TOOLS_MISSING=true
else
    echo "✓ Python3 installed"
fi

if [ "$TOOLS_MISSING" = true ]; then
    echo ""
    echo "=========================================="
    echo "Installing Missing Tools"
    echo "=========================================="
    echo ""
    echo "This will install: BWA, SAMtools, BCFtools, wget"
    echo "You may be prompted for your password..."
    echo ""
    
    # Update package lists
    sudo apt update
    
    # Install tools
    sudo apt install -y bwa samtools bcftools wget
    
    echo ""
    echo "✓ Installation complete!"
fi

echo ""
echo "=========================================="
echo "All Required Tools are Installed!"
echo "=========================================="
echo ""

# Navigate to APOE directory
cd /mnt/d/APOE

# Make pipeline executable
chmod +x apoe_pipeline.sh

# Check if FASTQ files exist
if [ ! -f "HRR024685_f1.fq.gz" ]; then
    echo "ERROR: FASTQ file HRR024685_f1.fq.gz not found!"
    exit 1
fi

if [ ! -f "HRR024685_r2.fq.gz" ]; then
    echo "ERROR: FASTQ file HRR024685_r2.fq.gz not found!"
    exit 1
fi

echo "✓ FASTQ files found"
echo ""

echo "=========================================="
echo "Starting APOE Genotyping Pipeline"
echo "=========================================="
echo ""
echo "Sample: HRR024685"
echo "Working directory: /mnt/d/APOE"
echo "Expected runtime: 30-120 minutes"
echo ""
echo "The pipeline will:"
echo "  1. Download reference genome (~3 GB)"
echo "  2. Index reference genome"
echo "  3. Align FASTQ reads to reference"
echo "  4. Extract APOE gene region"
echo "  5. Call variants"
echo "  6. Determine APOE genotype"
echo "  7. Generate clinical report"
echo ""
echo "Starting in 3 seconds..."
echo "Press Ctrl+C to cancel..."
sleep 3

echo ""
echo "=========================================="
echo "Pipeline Starting..."
echo "=========================================="
echo ""

# Run the pipeline
./apoe_pipeline.sh

EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo "=========================================="
    echo "✓ ANALYSIS COMPLETE!"
    echo "=========================================="
    echo ""
    echo "Your APOE genotype report is ready:"
    echo "  /mnt/d/APOE/apoe_analysis/results/HRR024685_apoe_report.txt"
    echo ""
    echo "View it with:"
    echo "  cat apoe_analysis/results/HRR024685_apoe_report.txt"
    echo ""
    echo "Or in Windows:"
    echo "  D:\\APOE\\apoe_analysis\\results\\HRR024685_apoe_report.txt"
    echo ""
else
    echo "=========================================="
    echo "✗ Pipeline failed with exit code: $EXIT_CODE"
    echo "=========================================="
    echo ""
    echo "Check the error messages above for details."
fi

