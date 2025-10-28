#!/bin/bash

################################################################################
# Hourly Monitoring & Auto-Recovery Script
# Checks batch analysis progress and fixes issues automatically
################################################################################

LOG_FILE="/mnt/d/APOE/batch_analysis.log"
MONITOR_LOG="/mnt/d/APOE/hourly_monitor.log"
STATUS_FILE="/mnt/d/APOE/HOURLY_STATUS.txt"
RESULTS_DIR="/mnt/d/APOE/batch_analysis/results"
PIPELINE_SCRIPT="/mnt/d/APOE/batch_apoe_pipeline.sh"

echo "========================================" | tee -a $MONITOR_LOG
echo "Hourly Monitor Started: $(date)" | tee -a $MONITOR_LOG
echo "========================================" | tee -a $MONITOR_LOG
echo ""

while true; do
    echo "[$(date)] Running hourly check..." | tee -a $MONITOR_LOG
    
    # Check if pipeline is running
    if ps aux | grep -q "[b]atch_apoe_pipeline"; then
        PIPELINE_STATUS="RUNNING"
        echo "[$(date)] ✓ Pipeline is running" | tee -a $MONITOR_LOG
    else
        PIPELINE_STATUS="STOPPED"
        echo "[$(date)] ✗ Pipeline has stopped!" | tee -a $MONITOR_LOG
    fi
    
    # Count completed samples
    COMPLETED=$(ls "$RESULTS_DIR"/*.txt 2>/dev/null | wc -l)
    echo "[$(date)] Completed: $COMPLETED/26 samples" | tee -a $MONITOR_LOG
    
    # Check for errors
    ERRORS=$(grep -c "ERROR\|truncated\|Aborting" "$LOG_FILE" 2>/dev/null || echo "0")
    echo "[$(date)] Errors detected: $ERRORS" | tee -a $MONITOR_LOG
    
    # Get current sample
    CURRENT=$(grep "Processing:" "$LOG_FILE" 2>/dev/null | tail -1 | awk '{print $3}')
    if [ ! -z "$CURRENT" ]; then
        echo "[$(date)] Current sample: $CURRENT" | tee -a $MONITOR_LOG
    fi
    
    # Check if pipeline stopped but not finished
    if [ "$PIPELINE_STATUS" = "STOPPED" ] && [ "$COMPLETED" -lt 26 ]; then
        echo "[$(date)] ⚠️ Pipeline stopped before completion!" | tee -a $MONITOR_LOG
        echo "[$(date)] Attempting to restart..." | tee -a $MONITOR_LOG
        
        cd /mnt/d/APOE
        nohup bash batch_apoe_pipeline.sh >> batch_analysis.log 2>&1 &
        
        echo "[$(date)] ✓ Pipeline restarted" | tee -a $MONITOR_LOG
    fi
    
    # Check for stuck processes (no progress in 2 hours)
    if [ -f "$LOG_FILE" ]; then
        LAST_ACTIVITY=$(stat -c %Y "$LOG_FILE")
        CURRENT_TIME=$(date +%s)
        TIME_DIFF=$((CURRENT_TIME - LAST_ACTIVITY))
        
        if [ $TIME_DIFF -gt 7200 ]; then  # 2 hours
            echo "[$(date)] ⚠️ No activity for 2 hours - pipeline may be stuck" | tee -a $MONITOR_LOG
            echo "[$(date)] Killing stuck processes..." | tee -a $MONITOR_LOG
            pkill -f batch_apoe_pipeline
            pkill -f "bwa mem"
            sleep 5
            
            echo "[$(date)] Restarting pipeline..." | tee -a $MONITOR_LOG
            cd /mnt/d/APOE
            nohup bash batch_apoe_pipeline.sh >> batch_analysis.log 2>&1 &
            echo "[$(date)] ✓ Pipeline restarted after stuck detection" | tee -a $MONITOR_LOG
        fi
    fi
    
    # Generate status report
    cat > "$STATUS_FILE" << EOF
========================================
BATCH ANALYSIS - HOURLY STATUS REPORT
========================================
Last Check: $(date)

Pipeline Status: $PIPELINE_STATUS
Completed: $COMPLETED / 26 samples
Errors: $ERRORS
Current Sample: ${CURRENT:-"Unknown"}

Progress: $(awk "BEGIN {printf \"%.1f\", ($COMPLETED/26)*100}")%

Recent Activity:
$(tail -5 "$LOG_FILE" 2>/dev/null | grep -E "\[.*\]" | tail -3)

Next Check: $(date -d '+1 hour' 2>/dev/null || date -v +1H 2>/dev/null || echo "In 1 hour")
========================================
EOF
    
    cat "$STATUS_FILE" | tee -a $MONITOR_LOG
    
    # Check if complete
    if [ "$COMPLETED" -eq 26 ]; then
        echo "" | tee -a $MONITOR_LOG
        echo "========================================" | tee -a $MONITOR_LOG
        echo "✓✓✓ ALL 26 SAMPLES COMPLETED! ✓✓✓" | tee -a $MONITOR_LOG
        echo "========================================" | tee -a $MONITOR_LOG
        break
    fi
    
    echo "[$(date)] Next check in 1 hour..." | tee -a $MONITOR_LOG
    echo "" | tee -a $MONITOR_LOG
    
    # Wait 1 hour
    sleep 3600
done

echo ""
echo "Monitoring complete. All samples processed."


