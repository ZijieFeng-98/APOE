#!/bin/bash

# Live monitoring with automatic updates

OUTPUT_DIR="/mnt/d/APOE/final_analysis"
TOTAL_SAMPLES=24
START_TIME=$(date +%s)

while true; do
    clear
    echo "================================================================================"
    echo "                    APOE GENOTYPING - LIVE MONITOR"
    echo "================================================================================"
    echo ""
    echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Uptime: $((( $(date +%s) - START_TIME ) / 60)) minutes"
    echo ""
    
    # Check if pipeline is running
    if pgrep -f "bulletproof_pipeline" > /dev/null; then
        echo "Pipeline Status: 🟢 RUNNING"
    else
        echo "Pipeline Status: 🟡 COMPLETED OR STOPPED"
    fi
    echo ""
    
    # Count completed samples
    COMPLETED=$(ls "$OUTPUT_DIR/results/"*.txt 2>/dev/null | wc -l)
    PERCENT=$((COMPLETED * 100 / TOTAL_SAMPLES))
    
    echo "Progress: $COMPLETED / $TOTAL_SAMPLES ($PERCENT%)"
    
    # Progress bar
    BARS=$((COMPLETED * 50 / TOTAL_SAMPLES))
    printf "["
    for ((i=0; i<BARS; i++)); do printf "█"; done
    for ((i=BARS; i<50; i++)); do printf "░"; done
    printf "]\n"
    echo ""
    
    # Show currently processing
    CURRENT_PROCESS=$(ps aux | grep 'bwa.*HRR' | grep -v grep | awk '{print $NF}' | grep -o 'HRR[0-9]*' | head -1)
    if [ -n "$CURRENT_PROCESS" ]; then
        echo "Currently Processing: $CURRENT_PROCESS"
        CURRENT_LOG="$OUTPUT_DIR/logs/${CURRENT_PROCESS}.log"
        if [ -f "$CURRENT_LOG" ]; then
            echo "Latest: $(tail -1 "$CURRENT_LOG" 2>/dev/null)"
        fi
    else
        echo "Currently Processing: [Checking...]"
    fi
    echo ""
    
    # Show genotype distribution
    if [ $COMPLETED -gt 0 ]; then
        echo "Genotype Distribution:"
        echo "────────────────────────────────────────────────────────────────────────────────"
        for gt in "ε2/ε2" "ε2/ε3" "ε2/ε4" "ε3/ε3" "ε3/ε4" "ε4/ε4"; do
            COUNT=$(grep -l "APOE GENOTYPE: $gt" "$OUTPUT_DIR/results/"*.txt 2>/dev/null | wc -l)
            if [ $COUNT -gt 0 ]; then
                printf "  %-8s: %2d patients\n" "$gt" "$COUNT"
            fi
        done
        echo ""
        
        # Show last 5 completed
        echo "Recently Completed:"
        echo "────────────────────────────────────────────────────────────────────────────────"
        ls -t "$OUTPUT_DIR/results/"*.txt 2>/dev/null | head -5 | while read f; do
            SAMPLE=$(basename "$f" _result.txt)
            GT=$(grep "APOE GENOTYPE:" "$f" | cut -d: -f2 | xargs)
            RISK=$(grep "Alzheimer's Disease Risk:" "$f" | cut -d: -f2 | xargs | cut -d'(' -f1)
            printf "  ✓ %-12s %8s  %s\n" "$SAMPLE" "$GT" "$RISK"
        done
    fi
    
    echo ""
    echo "================================================================================"
    echo "Press Ctrl+C to exit monitor (pipeline will continue running)"
    echo "Results: $OUTPUT_DIR/results/ | Logs: $OUTPUT_DIR/logs/"
    echo "================================================================================"
    
    # Exit if all complete
    if [ $COMPLETED -ge $TOTAL_SAMPLES ]; then
        echo ""
        echo "🎉 ALL SAMPLES COMPLETED! 🎉"
        echo ""
        echo "View cohort summary: cat $OUTPUT_DIR/COHORT_SUMMARY_REPORT.txt"
        break
    fi
    
    sleep 30
done


