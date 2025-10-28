#!/bin/bash

# Automatic monitoring and recovery script

RESULT="/mnt/d/APOE/apoe_analysis/results/HRR024685_apoe_report.txt"
COMPLETION_SCRIPT="/mnt/d/APOE/complete_pipeline.sh"
MONITOR_LOG="/mnt/d/APOE/monitor.log"

echo "Auto-monitor started at $(date)" | tee $MONITOR_LOG

while true; do
    # Check if result exists
    if [ -f "$RESULT" ]; then
        echo "[$(date)] SUCCESS! Pipeline complete!" | tee -a $MONITOR_LOG
        echo "Your APOE report is ready at: $RESULT" | tee -a $MONITOR_LOG
        exit 0
    fi
    
    # Check if completion script is running
    if ! ps aux | grep -q "[c]omplete_pipeline.sh"; then
        echo "[$(date)] Starting completion script..." | tee -a $MONITOR_LOG
        cd /mnt/d/APOE
        bash $COMPLETION_SCRIPT &
        sleep 30
    else
        echo "[$(date)] Completion script is running..." | tee -a $MONITOR_LOG
    fi
    
    # Wait before next check
    sleep 20
done


