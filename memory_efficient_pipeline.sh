#!/bin/bash

################################################################################
# Memory-Efficient APOE Genotyping Pipeline
# Optimized for low-memory systems (uses 2 threads, streaming sorting)
################################################################################

set -e
set -u

SAMPLE_DIR="/mnt/e/CGGA WESeq"
OUTPUT_BASE="/mnt/d/APOE/reanalysis"
REF_GENOME="/mnt/d/APOE/reference/human_g1k_v37.fasta"
THREADS=2  # Reduced for memory efficiency
APOE_REGION="19:45409039-45412650"
RS429358_POS="45411941"
RS7412_POS="45412079"

mkdir -p "$OUTPUT_BASE"
mkdir -p "$OUTPUT_BASE/results"
mkdir -p "$OUTPUT_BASE/logs"

# Color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

# Failed samples from previous run (skip HRR024685 and HRR024697 - already done)
FAILED_SAMPLES=(
    "HRR024686" "HRR024687" "HRR024688" "HRR024689" "HRR024690"
    "HRR024691" "HRR024692" "HRR024693" "HRR024694" "HRR024695"
    "HRR024696" "HRR024698" "HRR024699" "HRR024700" "HRR024701"
    "HRR024702" "HRR024703" "HRR024704" "HRR024705" "HRR024706"
    "HRR024707" "HRR024708" "HRR024709" "HRR024710"
)

TOTAL=${#FAILED_SAMPLES[@]}
CURRENT=0
SUCCESS=0
FAILED=0

echo "=========================================="
echo "MEMORY-EFFICIENT APOE RE-ANALYSIS"
echo "=========================================="
echo ""
echo "Samples to process: $TOTAL"
echo "Threads: $THREADS (memory-efficient)"
echo ""

for SAMPLE_ID in "${FAILED_SAMPLES[@]}"; do
    CURRENT=$((CURRENT + 1))
    
    echo "=========================================="
    echo "[$CURRENT/$TOTAL] Processing: $SAMPLE_ID"
    echo "=========================================="
    echo "Start: $(date)"
    
    # Input files
    FASTQ_R1="${SAMPLE_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_f1.fq.gz"
    FASTQ_R2="${SAMPLE_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_r2.fq.gz"
    
    if [ ! -f "$FASTQ_R1" ] || [ ! -f "$FASTQ_R2" ]; then
        echo -e "${RED}[ERROR]${NC} FASTQ files not found"
        FAILED=$((FAILED + 1))
        continue
    fi
    
    # Output directory
    SAMPLE_OUT="$OUTPUT_BASE/$SAMPLE_ID"
    mkdir -p "$SAMPLE_OUT"
    
    LOG_FILE="$OUTPUT_BASE/logs/${SAMPLE_ID}.log"
    
    {
        echo "[$(date)] Starting $SAMPLE_ID"
        
        # Align and sort in one go (streaming, memory efficient)
        echo "[$(date)] Aligning and sorting..."
        bwa mem -t $THREADS -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA" \
            "$REF_GENOME" "$FASTQ_R1" "$FASTQ_R2" 2>&1 | \
            samtools view -@ $THREADS -b - | \
            samtools sort -@ $THREADS -m 500M -o "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" - 2>&1
        
        if [ ! -f "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" ]; then
            echo "[$(date)] ERROR: Sorted BAM not created"
            exit 1
        fi
        
        echo "[$(date)] Indexing BAM..."
        samtools index "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" 2>&1
        
        # Extract APOE region
        echo "[$(date)] Extracting APOE region..."
        samtools view -@ $THREADS -b "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" "$APOE_REGION" > \
            "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>&1
        samtools index "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>&1
        
        # Calculate coverage
        AVG_COV=$(samtools depth "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" | \
            awk '{sum+=$3; n++} END {if(n>0) print sum/n; else print 0}')
        echo "[$(date)] APOE coverage: ${AVG_COV}x"
        
        # Get genotypes at key SNPs
        echo "[$(date)] Calling genotypes..."
        RS429358_PILEUP=$(samtools mpileup -f "$REF_GENOME" \
            -r "19:${RS429358_POS}-${RS429358_POS}" \
            "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>/dev/null | head -1)
        
        RS7412_PILEUP=$(samtools mpileup -f "$REF_GENOME" \
            -r "19:${RS7412_POS}-${RS7412_POS}" \
            "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>/dev/null | head -1)
        
        if [ -z "$RS429358_PILEUP" ] || [ -z "$RS7412_PILEUP" ]; then
            echo "[$(date)] ERROR: No coverage at key positions"
            exit 1
        fi
        
        # Parse pileup results
        RS429358_COV=$(echo "$RS429358_PILEUP" | awk '{print $4}')
        RS429358_BASES=$(echo "$RS429358_PILEUP" | awk '{print $5}')
        RS7412_COV=$(echo "$RS7412_PILEUP" | awk '{print $4}')
        RS7412_BASES=$(echo "$RS7412_PILEUP" | awk '{print $5}')
        
        echo "[$(date)] rs429358: ${RS429358_COV}x coverage"
        echo "[$(date)] rs7412: ${RS7412_COV}x coverage"
        
        # Determine genotypes (Reference: rs429358=T, rs7412=C)
        # Dots and commas = match reference
        if [[ "$RS429358_BASES" =~ ^[\.,]+$ ]]; then
            RS429358_GT="T/T"
        else
            RS429358_GT="Has variants"
        fi
        
        if [[ "$RS7412_BASES" =~ ^[\.,]+$ ]]; then
            RS7412_GT="C/C"
        else
            RS7412_GT="Has variants"
        fi
        
        # Determine APOE genotype
        if [ "$RS429358_GT" = "T/T" ] && [ "$RS7412_GT" = "C/C" ]; then
            APOE_GT="ε4/ε4"
            RISK="SIGNIFICANTLY INCREASED RISK (8-12x)"
        else
            APOE_GT="Requires full variant analysis"
            RISK="Unknown - needs detailed variant calling"
        fi
        
        echo "[$(date)] APOE Genotype: $APOE_GT"
        
        # Save result
        cat > "$OUTPUT_BASE/results/${SAMPLE_ID}_result.txt" << EOF
Sample: $SAMPLE_ID
Date: $(date)
Status: SUCCESS

Coverage:
- APOE region: ${AVG_COV}x
- rs429358: ${RS429358_COV}x
- rs7412: ${RS7412_COV}x

Raw Genotype Data:
- rs429358 bases: $RS429358_BASES
- rs7412 bases: $RS7412_BASES

Genotypes:
- rs429358: $RS429358_GT
- rs7412: $RS7412_GT

APOE Genotype: $APOE_GT
Risk Assessment: $RISK

Files:
- Full alignment: ${SAMPLE_OUT}/${SAMPLE_ID}.sorted.bam
- APOE region: ${SAMPLE_OUT}/${SAMPLE_ID}.apoe.bam
EOF
        
        echo -e "${GREEN}[SUCCESS]${NC} $SAMPLE_ID completed"
        echo "[$(date)] Finished $SAMPLE_ID"
        
    } 2>&1 | tee "$LOG_FILE"
    
    if [ $? -eq 0 ]; then
        SUCCESS=$((SUCCESS + 1))
        echo -e "${GREEN}✓ $SAMPLE_ID: SUCCESS${NC}"
    else
        FAILED=$((FAILED + 1))
        echo -e "${RED}✗ $SAMPLE_ID: FAILED${NC}"
    fi
    
    echo ""
    
    # Clean up large sorted BAM to save space (keep APOE region only)
    if [ -f "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" ]; then
        rm -f "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam"
        rm -f "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam.bai"
        echo "Cleaned up full BAM (kept APOE region only)"
    fi
    
    echo ""
done

echo "=========================================="
echo "RE-ANALYSIS COMPLETE"
echo "=========================================="
echo "Total: $TOTAL"
echo -e "${GREEN}Success: $SUCCESS${NC}"
echo -e "${RED}Failed: $FAILED${NC}"
echo ""
echo "Results: $OUTPUT_BASE/results/"
echo ""


