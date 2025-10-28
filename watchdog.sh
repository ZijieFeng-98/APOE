#!/bin/bash

# Watchdog script to monitor and restart the pipeline if needed

LOG_FILE="/mnt/d/APOE/pipeline_output.log"
RESULT_FILE="/mnt/d/APOE/apoe_analysis/results/HRR024685_apoe_report.txt"
WATCHDOG_LOG="/mnt/d/APOE/watchdog.log"

echo "========================================" | tee -a $WATCHDOG_LOG
echo "Watchdog started at $(date)" | tee -a $WATCHDOG_LOG
echo "========================================" | tee -a $WATCHDOG_LOG

while true; do
    # Check if result file exists (pipeline completed successfully)
    if [ -f "$RESULT_FILE" ]; then
        echo "[$(date)] SUCCESS! Pipeline completed. Result file found." | tee -a $WATCHDOG_LOG
        echo "Your APOE genotype report is ready at: $RESULT_FILE" | tee -a $WATCHDOG_LOG
        exit 0
    fi
    
    # Check if pipeline is running
    if ps aux | grep -q "[a]poe_pipeline.sh"; then
        echo "[$(date)] Pipeline is running..." | tee -a $WATCHDOG_LOG
    else
        # Pipeline not running, check if it completed or failed
        if [ -f "$RESULT_FILE" ]; then
            echo "[$(date)] Pipeline completed successfully!" | tee -a $WATCHDOG_LOG
            exit 0
        else
            echo "[$(date)] Pipeline stopped! Restarting..." | tee -a $WATCHDOG_LOG
            cd /mnt/d/APOE
            nohup bash apoe_pipeline.sh >> $LOG_FILE 2>&1 &
            sleep 10
        fi
    fi
    
    # Wait before next check
    sleep 60
done

