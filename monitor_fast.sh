#!/bin/bash

# Monitor the fast pipeline progress

while true; do
    clear
    echo "=============================================="
    echo "  APOE HIGH-PERFORMANCE PIPELINE MONITOR"
    echo "=============================================="
    echo "Time: $(date)"
    echo ""
    
    # Check if pipeline is running
    if pgrep -f "fast_pipeline.sh" > /dev/null; then
        echo "Status: 🟢 RUNNING"
    else
        echo "Status: 🔴 STOPPED or COMPLETED"
    fi
    echo ""
    
    # Count completed samples
    COMPLETED=$(ls /mnt/d/APOE/fast_analysis/results/*.txt 2>/dev/null | wc -l)
    echo "Completed: $COMPLETED / 24"
    echo ""
    
    # Show current process
    echo "Current processes:"
    ps aux | grep -E 'bwa|samtools|bcftools' | grep -v grep | awk '{printf "  %-10s %5s%% CPU %5s%% MEM %s\n", $11, $3, $4, substr($0, index($0,$12))}' | head -10
    echo ""
    
    # Show recent completions
    if [ $COMPLETED -gt 0 ]; then
        echo "Recent completions:"
        ls -t /mnt/d/APOE/fast_analysis/results/*.txt 2>/dev/null | head -5 | while read file; do
            SAMPLE=$(basename "$file" _result.txt)
            GT=$(grep "APOE GENOTYPE:" "$file" | cut -d: -f2 | xargs 2>/dev/null)
            echo "  ✓ $SAMPLE: $GT"
        done
        echo ""
    fi
    
    # Show current sample being processed
    CURRENT_LOG=$(ls -t /mnt/d/APOE/fast_analysis/logs/*.log 2>/dev/null | head -1)
    if [ -n "$CURRENT_LOG" ]; then
        CURRENT_SAMPLE=$(basename "$CURRENT_LOG" .log)
        echo "Currently processing: $CURRENT_SAMPLE"
        echo "Last log entries:"
        tail -3 "$CURRENT_LOG" 2>/dev/null | sed 's/^/  /'
    fi
    
    echo ""
    echo "=============================================="
    echo "Press Ctrl+C to exit monitor (pipeline continues)"
    echo "=============================================="
    
    sleep 30
done


