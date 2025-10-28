#!/bin/bash

echo "=== APOE PIPELINE STATUS ==="
echo "Time: $(date)"
echo ""

COMPLETED=$(ls /mnt/d/APOE/working_analysis/results/*.txt 2>/dev/null | wc -l)
echo "Completed: $COMPLETED / 24"

CURRENT=$(ps aux | grep 'bwa.*HRR' | grep -v grep | grep -o 'HRR[0-9]*' | head -1)
if [ -n "$CURRENT" ]; then
    echo "Currently processing: $CURRENT"
else
    echo "No active BWA process"
fi

echo ""

if [ $COMPLETED -gt 0 ]; then
    echo "Latest completions:"
    ls -t /mnt/d/APOE/working_analysis/results/*.txt 2>/dev/null | head -3 | while read f; do
        SAMPLE=$(basename "$f" _result.txt)
        GT=$(grep "APOE Genotype:" "$f" | cut -d: -f2 | xargs)
        echo "  - $SAMPLE: $GT"
    done
fi


