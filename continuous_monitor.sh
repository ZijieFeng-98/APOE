#!/bin/bash

OUTPUT="/mnt/d/APOE/working_analysis"
LOG_FILE="/mnt/d/APOE/COMPLETION_LOG.txt"

echo "=== CONTINUOUS MONITORING STARTED ===" | tee "$LOG_FILE"
echo "Time: $(date)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

LAST_COMPLETED=0

while true; do
    COMPLETED=$(ls "$OUTPUT/results/"*.txt 2>/dev/null | wc -l)
    
    if [ $COMPLETED -gt $LAST_COMPLETED ]; then
        echo "======================================" | tee -a "$LOG_FILE"
        echo "🎉 NEW SAMPLE COMPLETED!" | tee -a "$LOG_FILE"
        echo "Time: $(date)" | tee -a "$LOG_FILE"
        echo "Total completed: $COMPLETED/24" | tee -a "$LOG_FILE"
        echo "======================================" | tee -a "$LOG_FILE"
        
        # Show latest result
        LATEST=$(ls -t "$OUTPUT/results/"*.txt | head -1)
        cat "$LATEST" | tee -a "$LOG_FILE"
        echo "" | tee -a "$LOG_FILE"
        
        LAST_COMPLETED=$COMPLETED
    fi
    
    # Show current sample being processed every 5 minutes
    CURRENT=$(ps aux | grep 'bwa.*HRR' | grep -v grep | awk '{print $NF}' | grep -o 'HRR[0-9]*' | head -1)
    if [ -n "$CURRENT" ] && [ $(($(date +%s) % 300)) -lt 30 ]; then
        echo "[$(date '+%H:%M:%S')] Processing: $CURRENT | Completed: $COMPLETED/24"
    fi
    
    # Exit when all complete
    if [ $COMPLETED -ge 24 ]; then
        echo "" | tee -a "$LOG_FILE"
        echo "========================================" | tee -a "$LOG_FILE"
        echo "🎉 ALL 24 SAMPLES COMPLETED!" | tee -a "$LOG_FILE"
        echo "Finished: $(date)" | tee -a "$LOG_FILE"
        echo "========================================" | tee -a "$LOG_FILE"
        break
    fi
    
    sleep 30
done


