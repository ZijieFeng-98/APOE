#!/bin/bash

# Automatic monitoring and reporting system

OUTPUT_DIR="/mnt/d/APOE/final_analysis"
REPORT_FILE="/mnt/d/APOE/PROGRESS_UPDATES.txt"
TOTAL=24
LAST_COUNT=0

echo "=== APOE PIPELINE AUTO-MONITOR ===" > "$REPORT_FILE"
echo "Started: $(date)" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"

while true; do
    COMPLETED=$(ls "$OUTPUT_DIR/results/"*.txt 2>/dev/null | wc -l)
    
    # Report new completions
    if [ $COMPLETED -gt $LAST_COUNT ]; then
        echo "[$(date '+%H:%M:%S')] NEW COMPLETION! Total: $COMPLETED/$TOTAL" | tee -a "$REPORT_FILE"
        
        # Get the newest result
        NEWEST=$(ls -t "$OUTPUT_DIR/results/"*.txt 2>/dev/null | head -1)
        if [ -f "$NEWEST" ]; then
            SAMPLE=$(basename "$NEWEST" _result.txt)
            GT=$(grep "APOE GENOTYPE:" "$NEWEST" | cut -d: -f2 | xargs)
            RISK=$(grep "Alzheimer's Disease Risk:" "$NEWEST" | cut -d: -f2 | xargs)
            COV=$(grep "APOE region" "$NEWEST" | awk -F': ' '{print $NF}')
            
            echo "  Sample: $SAMPLE" | tee -a "$REPORT_FILE"
            echo "  Genotype: $GT" | tee -a "$REPORT_FILE"
            echo "  Risk: $RISK" | tee -a "$REPORT_FILE"
            echo "  Coverage: $COV" | tee -a "$REPORT_FILE"
            echo "" | tee -a "$REPORT_FILE"
        fi
        
        LAST_COUNT=$COMPLETED
    fi
    
    # Check if all complete
    if [ $COMPLETED -ge $TOTAL ]; then
        echo "" | tee -a "$REPORT_FILE"
        echo "========================================" | tee -a "$REPORT_FILE"
        echo "🎉 ALL $TOTAL SAMPLES COMPLETED! 🎉" | tee -a "$REPORT_FILE"
        echo "========================================" | tee -a "$REPORT_FILE"
        echo "Completed: $(date)" | tee -a "$REPORT_FILE"
        echo "" | tee -a "$REPORT_FILE"
        
        # Generate final summary
        echo "FINAL GENOTYPE DISTRIBUTION:" | tee -a "$REPORT_FILE"
        echo "----------------------------" | tee -a "$REPORT_FILE"
        for gt in "ε2/ε2" "ε2/ε3" "ε2/ε4" "ε3/ε3" "ε3/ε4" "ε4/ε4"; do
            COUNT=$(grep -l "APOE GENOTYPE: $gt" "$OUTPUT_DIR/results/"*.txt 2>/dev/null | wc -l)
            if [ $COUNT -gt 0 ]; then
                echo "  $gt: $COUNT patients" | tee -a "$REPORT_FILE"
            fi
        done
        
        echo "" | tee -a "$REPORT_FILE"
        echo "Full report: $OUTPUT_DIR/COHORT_SUMMARY_REPORT.txt" | tee -a "$REPORT_FILE"
        
        break
    fi
    
    # Update every 2 minutes
    sleep 120
done


