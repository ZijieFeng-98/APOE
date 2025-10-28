#!/bin/bash

################################################################################
# BATCH APOE Genotyping Pipeline for Multiple WGS Samples
# Processes all patient samples from CGGA WESeq cohort
################################################################################

set -e
set -u

# Configuration
SAMPLE_DIR="/mnt/e/CGGA WESeq"
OUTPUT_BASE="/mnt/d/APOE/batch_analysis"
REF_GENOME="/mnt/d/APOE/reference/human_g1k_v37.fasta"
THREADS=4

# APOE region
APOE_REGION="19:45409039-45412650"
RS429358_POS="45411941"
RS7412_POS="45412079"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

################################################################################
# Create output directories
################################################################################

mkdir -p "$OUTPUT_BASE"
mkdir -p "$OUTPUT_BASE/results"
mkdir -p "$OUTPUT_BASE/logs"
mkdir -p "$OUTPUT_BASE/summary"

################################################################################
# Process each sample
################################################################################

echo "=========================================="
echo "BATCH APOE ANALYSIS"
echo "=========================================="
echo ""
echo "Processing samples from: $SAMPLE_DIR"
echo "Output directory: $OUTPUT_BASE"
echo "Reference genome: $REF_GENOME"
echo "Threads: $THREADS"
echo ""

# Count samples
SAMPLE_COUNT=$(ls -d "$SAMPLE_DIR"/HRR* 2>/dev/null | wc -l)
echo "Total samples found: $SAMPLE_COUNT"
echo ""

# Process each sample
CURRENT=0
SUCCESSFUL=0
FAILED=0

for SAMPLE_PATH in "$SAMPLE_DIR"/HRR*/; do
    CURRENT=$((CURRENT + 1))
    SAMPLE_ID=$(basename "$SAMPLE_PATH")
    
    echo "=========================================="
    echo "[$CURRENT/$SAMPLE_COUNT] Processing: $SAMPLE_ID"
    echo "=========================================="
    echo "Start time: $(date)"
    
    # Check for FASTQ files
    FASTQ_R1="${SAMPLE_PATH}${SAMPLE_ID}_f1.fq.gz"
    FASTQ_R2="${SAMPLE_PATH}${SAMPLE_ID}_r2.fq.gz"
    
    if [ ! -f "$FASTQ_R1" ] || [ ! -f "$FASTQ_R2" ]; then
        echo -e "${RED}[ERROR]${NC} FASTQ files not found for $SAMPLE_ID"
        FAILED=$((FAILED + 1))
        continue
    fi
    
    # Create sample output directory
    SAMPLE_OUT="$OUTPUT_BASE/$SAMPLE_ID"
    mkdir -p "$SAMPLE_OUT/alignment"
    mkdir -p "$SAMPLE_OUT/variants"
    
    LOG_FILE="$OUTPUT_BASE/logs/${SAMPLE_ID}.log"
    
    # Start processing
    {
        echo "[$(date)] Starting analysis for $SAMPLE_ID"
        
        # Alignment
        echo "[$(date)] Aligning reads..."
        if [ ! -f "$SAMPLE_OUT/alignment/${SAMPLE_ID}.sorted.bam" ]; then
            bwa mem -t $THREADS -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA" \
                "$REF_GENOME" "$FASTQ_R1" "$FASTQ_R2" | \
                samtools sort -@ $THREADS -o "$SAMPLE_OUT/alignment/${SAMPLE_ID}.sorted.bam" -
            
            samtools index "$SAMPLE_OUT/alignment/${SAMPLE_ID}.sorted.bam"
            echo "[$(date)] Alignment complete"
        else
            echo "[$(date)] Using existing alignment"
        fi
        
        # Extract APOE region
        echo "[$(date)] Extracting APOE region..."
        samtools view -@ $THREADS -b "$SAMPLE_OUT/alignment/${SAMPLE_ID}.sorted.bam" "$APOE_REGION" > \
            "$SAMPLE_OUT/alignment/${SAMPLE_ID}.apoe.bam"
        samtools index "$SAMPLE_OUT/alignment/${SAMPLE_ID}.apoe.bam"
        
        # Calculate coverage
        COVERAGE=$(samtools depth "$SAMPLE_OUT/alignment/${SAMPLE_ID}.apoe.bam" | awk '{sum+=$3} END {if(NR>0) print sum/NR; else print 0}')
        echo "[$(date)] APOE region coverage: ${COVERAGE}x"
        
        # Get genotypes at key positions
        echo "[$(date)] Calling genotypes..."
        
        # rs429358
        RS429358_PILEUP=$(samtools mpileup -f "$REF_GENOME" -r "19:${RS429358_POS}-${RS429358_POS}" \
            "$SAMPLE_OUT/alignment/${SAMPLE_ID}.apoe.bam" | head -1)
        RS429358_COV=$(echo "$RS429358_PILEUP" | awk '{print $4}')
        RS429358_BASES=$(echo "$RS429358_PILEUP" | awk '{print $5}')
        
        # rs7412
        RS7412_PILEUP=$(samtools mpileup -f "$REF_GENOME" -r "19:${RS7412_POS}-${RS7412_POS}" \
            "$SAMPLE_OUT/alignment/${SAMPLE_ID}.apoe.bam" | head -1)
        RS7412_COV=$(echo "$RS7412_PILEUP" | awk '{print $4}')
        RS7412_BASES=$(echo "$RS7412_PILEUP" | awk '{print $5}')
        
        echo "[$(date)] rs429358 coverage: ${RS429358_COV}x"
        echo "[$(date)] rs7412 coverage: ${RS7412_COV}x"
        
        # Determine APOE genotype
        # Reference: rs429358=T, rs7412=C
        # All dots/commas = reference (homozygous ref)
        # Any variants would show as letter or other symbols
        
        # Simple logic: if all reads match reference
        if [[ "$RS429358_BASES" =~ ^[\.,]+$ ]] && [[ "$RS7412_BASES" =~ ^[\.,]+$ ]]; then
            APOE_GENOTYPE="ε4/ε4"
            RS429358_GT="T/T"
            RS7412_GT="C/C"
        elif [[ "$RS429358_BASES" =~ ^[\.,]+$ ]]; then
            # rs429358 = T/T, need to check rs7412 variants
            APOE_GENOTYPE="Unknown (needs variant analysis)"
            RS429358_GT="T/T"
            RS7412_GT="Variant"
        else
            APOE_GENOTYPE="Unknown (needs variant analysis)"
            RS429358_GT="Variant"
            RS7412_GT="?"
        fi
        
        echo "[$(date)] APOE Genotype: $APOE_GENOTYPE"
        echo "[$(date)] rs429358: $RS429358_GT"
        echo "[$(date)] rs7412: $RS7412_GT"
        
        # Save results
        cat > "$OUTPUT_BASE/results/${SAMPLE_ID}_summary.txt" << EOF
Sample: $SAMPLE_ID
Date: $(date)
Status: SUCCESS

Coverage:
- APOE region: ${COVERAGE}x
- rs429358: ${RS429358_COV}x
- rs7412: ${RS7412_COV}x

Genotypes:
- rs429358: $RS429358_GT
- rs7412: $RS7412_GT

APOE Genotype: $APOE_GENOTYPE

EOF
        
        echo -e "${GREEN}[SUCCESS]${NC} $SAMPLE_ID completed"
        echo "[$(date)] Finished $SAMPLE_ID"
        
    } 2>&1 | tee "$LOG_FILE"
    
    if [ $? -eq 0 ]; then
        SUCCESSFUL=$((SUCCESSFUL + 1))
    else
        FAILED=$((FAILED + 1))
    fi
    
    echo ""
    
done

################################################################################
# Generate summary report
################################################################################

echo "=========================================="
echo "GENERATING SUMMARY REPORT"
echo "=========================================="

SUMMARY_FILE="$OUTPUT_BASE/summary/BATCH_SUMMARY_$(date +%Y%m%d_%H%M%S).txt"

cat > "$SUMMARY_FILE" << EOF
================================================================================
BATCH APOE GENOTYPING ANALYSIS - SUMMARY REPORT
================================================================================

Analysis Date: $(date)
Cohort: CGGA WESeq

Total Samples: $SAMPLE_COUNT
Successful: $SUCCESSFUL
Failed: $FAILED

================================================================================
INDIVIDUAL RESULTS
================================================================================

EOF

# Add each sample result
for result_file in "$OUTPUT_BASE/results"/*_summary.txt; do
    if [ -f "$result_file" ]; then
        cat "$result_file" >> "$SUMMARY_FILE"
        echo "--------------------------------------------------------------------------------" >> "$SUMMARY_FILE"
        echo "" >> "$SUMMARY_FILE"
    fi
done

# Generate statistics
cat >> "$SUMMARY_FILE" << EOF
================================================================================
GENOTYPE DISTRIBUTION
================================================================================

EOF

# Count genotypes
E4_E4=$(grep -c "ε4/ε4" "$OUTPUT_BASE/results"/*_summary.txt 2>/dev/null || echo "0")
UNKNOWN=$(grep -c "Unknown" "$OUTPUT_BASE/results"/*_summary.txt 2>/dev/null || echo "0")

cat >> "$SUMMARY_FILE" << EOF
ε4/ε4: $E4_E4 patients
Needs further analysis: $UNKNOWN patients

Note: Detailed genotype determination requires full variant calling for 
samples showing "Unknown". This initial analysis only detects homozygous
reference (ε4/ε4 when both SNPs are T/T and C/C).

================================================================================
OUTPUT FILES
================================================================================

Summary report: $SUMMARY_FILE
Individual results: $OUTPUT_BASE/results/
Detailed logs: $OUTPUT_BASE/logs/
Alignment files: $OUTPUT_BASE/<SAMPLE_ID>/alignment/

================================================================================
EOF

echo ""
echo "=========================================="
echo -e "${GREEN}BATCH ANALYSIS COMPLETE!${NC}"
echo "=========================================="
echo ""
echo "Results:"
echo "  Successful: $SUCCESSFUL/$SAMPLE_COUNT"
echo "  Failed: $FAILED/$SAMPLE_COUNT"
echo ""
echo "Summary report: $SUMMARY_FILE"
echo ""
echo "View results: cat $SUMMARY_FILE"
echo ""

# Display summary
cat "$SUMMARY_FILE"

