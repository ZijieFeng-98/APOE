#!/bin/bash

echo "Starting monitor... checking every 3 minutes"
echo "Pipeline started at 21:27:53"
echo ""

for i in {1..20}; do
    sleep 180
    
    COMPLETED=$(ls /mnt/d/APOE/ultra_fast/results/*.txt 2>/dev/null | wc -l)
    TIMESTAMP=$(date '+%H:%M:%S')
    
    echo "[$TIMESTAMP] Completed: $COMPLETED / 24"
    
    if [ $COMPLETED -gt 0 ]; then
        echo ""
        echo "✓ Samples completed!"
        for f in /mnt/d/APOE/ultra_fast/results/*.txt; do
            SAMPLE=$(basename "$f" _result.txt)
            GENOTYPE=$(grep "APOE Genotype:" "$f" | cut -d: -f2 | xargs)
            echo "  $SAMPLE: $GENOTYPE"
        done
    fi
    
    if [ $COMPLETED -ge 24 ]; then
        echo ""
        echo "🎉 ALL SAMPLES COMPLETE!"
        break
    fi
    
    echo ""
done


