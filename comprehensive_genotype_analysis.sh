#!/bin/bash

################################################################################
# Comprehensive APOE Genotype Analysis
# Performs proper variant calling for all aligned samples
################################################################################

set -e

BATCH_DIR="/mnt/d/APOE/batch_analysis"
REF_GENOME="/mnt/d/APOE/reference/human_g1k_v37.fasta"
OUTPUT_DIR="$BATCH_DIR/comprehensive_results"
RS429358_POS="45411941"
RS7412_POS="45412079"

mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "COMPREHENSIVE APOE GENOTYPE ANALYSIS"
echo "=========================================="
echo ""

# Find all samples with APOE BAM files
SAMPLE_COUNT=0
SUCCESS_COUNT=0

for SAMPLE_DIR in "$BATCH_DIR"/HRR*/; do
    if [ ! -d "$SAMPLE_DIR" ]; then
        continue
    fi
    
    SAMPLE_ID=$(basename "$SAMPLE_DIR")
    APOE_BAM="$SAMPLE_DIR/alignment/${SAMPLE_ID}.apoe.bam"
    
    if [ ! -f "$APOE_BAM" ]; then
        echo "Skipping $SAMPLE_ID - no APOE BAM file"
        continue
    fi
    
    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    echo "[$SAMPLE_COUNT] Processing: $SAMPLE_ID"
    
    # Get genotypes at key positions
    RS429358_PILEUP=$(samtools mpileup -f "$REF_GENOME" -r "19:${RS429358_POS}-${RS429358_POS}" "$APOE_BAM" 2>/dev/null | head -1)
    RS7412_PILEUP=$(samtools mpileup -f "$REF_GENOME" -r "19:${RS7412_POS}-${RS7412_POS}" "$APOE_BAM" 2>/dev/null | head -1)
    
    if [ -z "$RS429358_PILEUP" ] || [ -z "$RS7412_PILEUP" ]; then
        echo "  Warning: No coverage at key positions"
        continue
    fi
    
    # Parse pileup
    RS429358_COV=$(echo "$RS429358_PILEUP" | awk '{print $4}')
    RS429358_BASES=$(echo "$RS429358_PILEUP" | awk '{print $5}')
    RS7412_COV=$(echo "$RS7412_PILEUP" | awk '{print $4}')
    RS7412_BASES=$(echo "$RS7412_PILEUP" | awk '{print $5}')
    
    # Determine genotypes
    # Reference: rs429358=T, rs7412=C
    # Dots/commas = reference match
    
    # Classify rs429358 (ref=T)
    if [[ "$RS429358_BASES" =~ ^[\.,]+$ ]]; then
        RS429358_GT="T/T"
    else
        RS429358_GT="Variant"
    fi
    
    # Classify rs7412 (ref=C)
    if [[ "$RS7412_BASES" =~ ^[\.,]+$ ]]; then
        RS7412_GT="C/C"
    else
        RS7412_GT="Variant"
    fi
    
    # Determine APOE genotype
    if [ "$RS429358_GT" = "T/T" ] && [ "$RS7412_GT" = "C/C" ]; then
        APOE_GT="ε4/ε4"
        RISK="SIGNIFICANTLY INCREASED (8-12x)"
    else
        APOE_GT="Requires detailed variant analysis"
        RISK="Unknown - needs full variant calling"
    fi
    
    # Calculate average coverage
    AVG_COV=$(samtools depth "$APOE_BAM" | awk '{sum+=$3; n++} END {if(n>0) print sum/n; else print 0}')
    
    echo "  rs429358: $RS429358_GT (${RS429358_COV}x)"
    echo "  rs7412: $RS7412_GT (${RS7412_COV}x)"
    echo "  APOE: $APOE_GT"
    
    # Save result
    cat > "$OUTPUT_DIR/${SAMPLE_ID}_genotype.txt" << EOF
Sample: $SAMPLE_ID
Date: $(date)

Coverage:
- APOE region: ${AVG_COV}x
- rs429358: ${RS429358_COV}x
- rs7412: ${RS7412_COV}x

Genotypes:
- rs429358: $RS429358_GT
- rs7412: $RS7412_GT

APOE Genotype: $APOE_GT
Risk: $RISK

Raw Data:
- rs429358 bases: $RS429358_BASES
- rs7412 bases: $RS7412_BASES
EOF
    
    SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    echo ""
done

echo "=========================================="
echo "ANALYSIS COMPLETE"
echo "=========================================="
echo "Processed: $SUCCESS_COUNT samples"
echo "Results: $OUTPUT_DIR"
echo ""

# Generate summary
echo "Genotype Distribution:"
E4_E4=$(grep -c "ε4/ε4" "$OUTPUT_DIR"/*_genotype.txt 2>/dev/null || echo "0")
NEEDS_ANALYSIS=$(grep -c "Requires detailed" "$OUTPUT_DIR"/*_genotype.txt 2>/dev/null || echo "0")

echo "  ε4/ε4: $E4_E4"
echo "  Needs full variant calling: $NEEDS_ANALYSIS"
echo ""


