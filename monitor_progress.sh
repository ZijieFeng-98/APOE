#!/bin/bash

# Progress monitoring script

LOG_FILE="/mnt/d/APOE/pipeline_output.log"
PROGRESS_FILE="/mnt/d/APOE/CURRENT_STATUS.txt"

while true; do
    clear
    echo "=========================================="
    echo "APOE Pipeline - Live Progress Monitor"
    echo "Time: $(date)"
    echo "=========================================="
    echo ""
    
    # Check if pipeline is running
    if ps aux | grep -q "[a]poe_pipeline.sh"; then
        echo "✓ Pipeline Status: RUNNING"
    else
        echo "⚠ Pipeline Status: NOT RUNNING"
    fi
    echo ""
    
    # Check completed steps
    echo "Progress Checklist:"
    echo ""
    
    if grep -q "Reference genome already exists\|Reference genome downloaded" "$LOG_FILE" 2>/dev/null; then
        echo "✓ Reference genome ready"
    else
        echo "⏳ Downloading reference genome..."
    fi
    
    if grep -q "BWA index created\|BWA index already exists" "$LOG_FILE" 2>/dev/null; then
        echo "✓ BWA indexing complete"
    elif grep -q "Indexing reference genome for BWA" "$LOG_FILE" 2>/dev/null; then
        echo "⏳ BWA indexing in progress..."
        tail -3 "$LOG_FILE" 2>/dev/null | grep -E "BWTInc|bwa_index"
    fi
    
    if grep -q "SAMtools index already exists\|SAMtools index created" "$LOG_FILE" 2>/dev/null; then
        echo "✓ SAMtools indexing complete"
    fi
    
    if grep -q "Alignment complete" "$LOG_FILE" 2>/dev/null; then
        echo "✓ Read alignment complete"
    elif grep -q "Aligning reads" "$LOG_FILE" 2>/dev/null; then
        echo "⏳ Aligning reads (this takes longest)..."
    fi
    
    if grep -q "BAM file sorted and indexed" "$LOG_FILE" 2>/dev/null; then
        echo "✓ BAM sorting complete"
    fi
    
    if grep -q "APOE region extracted" "$LOG_FILE" 2>/dev/null; then
        echo "✓ APOE region extracted"
    fi
    
    if grep -q "Variant calling complete" "$LOG_FILE" 2>/dev/null; then
        echo "✓ Variant calling complete"
    fi
    
    if [ -f "/mnt/d/APOE/apoe_analysis/results/HRR024685_apoe_report.txt" ]; then
        echo "✓ APOE genotype report generated"
        echo ""
        echo "=========================================="
        echo "✓✓✓ PIPELINE COMPLETED SUCCESSFULLY! ✓✓✓"
        echo "=========================================="
        echo ""
        echo "Your results: /mnt/d/APOE/apoe_analysis/results/HRR024685_apoe_report.txt"
        break
    fi
    
    echo ""
    echo "Last 5 log lines:"
    echo "------------------------------------------"
    tail -5 "$LOG_FILE" 2>/dev/null
    echo "------------------------------------------"
    echo ""
    echo "Refreshing in 30 seconds... (Ctrl+C to exit)"
    
    sleep 30
done

