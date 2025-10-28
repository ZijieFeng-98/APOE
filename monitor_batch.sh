#!/bin/bash

# Real-time batch analysis monitor

OUTPUT_BASE="/mnt/d/APOE/batch_analysis"
LOG_FILE="/mnt/d/APOE/batch_analysis.log"

while true; do
    clear
    echo "=========================================="
    echo "BATCH APOE ANALYSIS - LIVE MONITOR"
    echo "Time: $(date)"
    echo "=========================================="
    echo ""
    
    # Check if pipeline is running
    if ps aux | grep -q "[b]atch_apoe_pipeline"; then
        echo "✓ Pipeline Status: RUNNING"
    else
        echo "⚠ Pipeline Status: NOT RUNNING"
    fi
    echo ""
    
    # Current sample
    CURRENT=$(grep "Processing:" "$LOG_FILE" 2>/dev/null | tail -1 | awk '{print $3}')
    if [ ! -z "$CURRENT" ]; then
        echo "Current Sample: $CURRENT"
        echo ""
    fi
    
    # Progress
    COMPLETED=$(ls "$OUTPUT_BASE"/results/*_summary.txt 2>/dev/null | wc -l)
    echo "Progress: $COMPLETED / 26 samples completed"
    echo ""
    
    # Success/Fail counts
    if [ -f "$LOG_FILE" ]; then
        SUCCESS=$(grep -c "SUCCESS" "$LOG_FILE" 2>/dev/null || echo "0")
        FAILED=$(grep -c "ERROR" "$LOG_FILE" 2>/dev/null || echo "0")
        echo "Success: $SUCCESS | Failed: $FAILED"
        echo ""
    fi
    
    # Latest completions
    echo "Recently Completed Samples:"
    echo "------------------------------------------"
    if [ -d "$OUTPUT_BASE/results" ]; then
        ls -t "$OUTPUT_BASE"/results/*_summary.txt 2>/dev/null | head -5 | while read file; do
            SAMPLE=$(basename "$file" | sed 's/_summary.txt//')
            GENOTYPE=$(grep "APOE Genotype:" "$file" | awk '{print $3}')
            echo "  $SAMPLE: $GENOTYPE"
        done
    fi
    echo ""
    
    # Latest log entries
    echo "Latest Activity:"
    echo "------------------------------------------"
    tail -10 "$LOG_FILE" 2>/dev/null | grep -E "\[.*\]" | tail -5
    echo ""
    
    echo "=========================================="
    echo "Press Ctrl+C to exit monitor"
    echo "=========================================="
    
    sleep 10
done

