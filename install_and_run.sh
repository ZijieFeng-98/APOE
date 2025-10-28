#!/bin/bash

echo "=========================================="
echo "APOE Analysis - Quick Install & Run"
echo "=========================================="
echo ""

# Install tools
echo "Installing tools (you'll need to enter your password)..."
sudo apt update
sudo apt install -y bwa samtools bcftools wget

echo ""
echo "✓ Tools installed!"
echo ""

# Navigate to directory
cd /mnt/d/APOE

# Run pipeline
echo "Starting APOE pipeline..."
bash apoe_pipeline.sh

