#!/bin/bash

################################################################################
# ULTRA-FAST APOE PIPELINE - Step-by-step processing for reliability
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit if any command in pipeline fails

SAMPLE_DIR="/mnt/e/CGGA WESeq"
OUTPUT_BASE="/mnt/d/APOE/ultra_fast"
REF_GENOME="/mnt/d/APOE/reference/human_g1k_v37.fasta"
THREADS=12
APOE_REGION="19:45409039-45412650"
RS429358_POS="45411941"
RS7412_POS="45412079"

mkdir -p "$OUTPUT_BASE"/{results,logs}

GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

SAMPLES=(
    "HRR024686" "HRR024687" "HRR024688" "HRR024689" "HRR024690"
    "HRR024691" "HRR024692" "HRR024693" "HRR024694" "HRR024695"
    "HRR024696" "HRR024698" "HRR024699" "HRR024700" "HRR024701"
    "HRR024702" "HRR024703" "HRR024704" "HRR024705" "HRR024706"
    "HRR024707" "HRR024708" "HRR024709" "HRR024710"
)

TOTAL=${#SAMPLES[@]}
SUCCESS=0
FAILED=0

echo "=========================================="
echo -e "${YELLOW}ULTRA-FAST APOE PIPELINE${NC}"
echo "=========================================="
echo "Samples: $TOTAL"
echo "Threads: $THREADS"
echo "Start: $(date)"
echo "=========================================="

for i in "${!SAMPLES[@]}"; do
    SAMPLE_ID="${SAMPLES[$i]}"
    NUM=$((i + 1))
    
    echo ""
    echo -e "${BLUE}==========================================${NC}"
    echo -e "${BLUE}[$NUM/$TOTAL] $SAMPLE_ID${NC}"
    echo -e "${BLUE}==========================================${NC}"
    
    FASTQ_R1="${SAMPLE_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_f1.fq.gz"
    FASTQ_R2="${SAMPLE_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_r2.fq.gz"
    
    if [[ ! -f "$FASTQ_R1" ]] || [[ ! -f "$FASTQ_R2" ]]; then
        echo -e "${RED}✗ FASTQ files missing${NC}"
        FAILED=$((FAILED + 1))
        continue
    fi
    
    SAMPLE_OUT="$OUTPUT_BASE/$SAMPLE_ID"
    mkdir -p "$SAMPLE_OUT"
    LOG="$OUTPUT_BASE/logs/${SAMPLE_ID}.log"
    
    {
        set +e  # Don't exit on error in subshell
        
        echo "[$(date '+%H:%M:%S')] Processing $SAMPLE_ID"
        
        # Step 1: Align with BWA
        echo "[$(date '+%H:%M:%S')] Step 1/6: BWA alignment..."
        bwa mem -t $THREADS \
            -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA" \
            "$REF_GENOME" "$FASTQ_R1" "$FASTQ_R2" \
            > "$SAMPLE_OUT/${SAMPLE_ID}.sam" 2>&1
        
        if [[ ! -s "$SAMPLE_OUT/${SAMPLE_ID}.sam" ]]; then
            echo "ERROR: SAM file not created or empty"
            exit 1
        fi
        
        #Step 2: Convert to BAM
        echo "[$(date '+%H:%M:%S')] Step 2/6: Converting to BAM..."
        samtools view -@ $THREADS -b "$SAMPLE_OUT/${SAMPLE_ID}.sam" \
            > "$SAMPLE_OUT/${SAMPLE_ID}.bam" 2>&1
        rm -f "$SAMPLE_OUT/${SAMPLE_ID}.sam"  # Save space
        
        # Step 3: Sort BAM
        echo "[$(date '+%H:%M:%S')] Step 3/6: Sorting BAM..."
        samtools sort -@ $THREADS -m 2G \
            -o "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" \
            "$SAMPLE_OUT/${SAMPLE_ID}.bam" 2>&1
        rm -f "$SAMPLE_OUT/${SAMPLE_ID}.bam"  # Save space
        
        if [[ ! -s "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" ]]; then
            echo "ERROR: Sorted BAM not created"
            exit 1
        fi
        
        # Step 4: Index BAM
        echo "[$(date '+%H:%M:%S')] Step 4/6: Indexing BAM..."
        samtools index -@ $THREADS "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" 2>&1
        
        # Step 5: Extract APOE region
        echo "[$(date '+%H:%M:%S')] Step 5/6: Extracting APOE..."
        samtools view -@ $THREADS -b \
            "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" "$APOE_REGION" \
            > "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>&1
        samtools index "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>&1
        
        # Step 6: Genotype
        echo "[$(date '+%H:%M:%S')] Step 6/6: Genotyping..."
        
        # Get coverage
        AVG_COV=$(samtools depth "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" | \
            awk '{sum+=$3; n++} END {if(n>0) printf "%.2f", sum/n; else print "0"}')
        
        # Get genotypes
        RS429358_DATA=$(samtools mpileup -f "$REF_GENOME" \
            -r "19:${RS429358_POS}-${RS429358_POS}" \
            "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>/dev/null | head -1)
        
        RS7412_DATA=$(samtools mpileup -f "$REF_GENOME" \
            -r "19:${RS7412_POS}-${RS7412_POS}" \
            "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>/dev/null | head -1)
        
        RS429358_COV=$(echo "$RS429358_DATA" | awk '{print $4}')
        RS429358_BASES=$(echo "$RS429358_DATA" | awk '{print $5}')
        RS7412_COV=$(echo "$RS7412_DATA" | awk '{print $4}')
        RS7412_BASES=$(echo "$RS7412_DATA" | awk '{print $5}')
        
        # Determine genotypes (Ref: rs429358=T, rs7412=C)
        # Dots/commas = reference
        # Letters = variants
        
        # rs429358 (T->C mutation = ε4)
        if [[ "$RS429358_BASES" =~ ^[\.,]+$ ]]; then
            RS429358_GT="T/T"
        elif echo "$RS429358_BASES" | grep -qi "c"; then
            C_COUNT=$(echo "$RS429358_BASES" | grep -o "[Cc]" | wc -l)
            REF_COUNT=$(echo "$RS429358_BASES" | grep -o "[\.,]" | wc -l)
            if (( C_COUNT > 0 && REF_COUNT == 0 )); then
                RS429358_GT="C/C"
            elif (( C_COUNT > 0 )); then
                RS429358_GT="T/C"
            else
                RS429358_GT="T/T"
            fi
        else
            RS429358_GT="T/T"
        fi
        
        # rs7412 (C->T mutation = ε2)
        if [[ "$RS7412_BASES" =~ ^[\.,]+$ ]]; then
            RS7412_GT="C/C"
        elif echo "$RS7412_BASES" | grep -qi "t"; then
            T_COUNT=$(echo "$RS7412_BASES" | grep -o "[Tt]" | wc -l)
            REF_COUNT=$(echo "$RS7412_BASES" | grep -o "[\.,]" | wc -l)
            if (( T_COUNT > 0 && REF_COUNT == 0 )); then
                RS7412_GT="T/T"
            elif (( T_COUNT > 0 )); then
                RS7412_GT="C/T"
            else
                RS7412_GT="C/C"
            fi
        else
            RS7412_GT="C/C"
        fi
        
        # Determine APOE genotype
        if [[ "$RS429358_GT" == "T/T" ]] && [[ "$RS7412_GT" == "C/C" ]]; then
            APOE_GT="ε3/ε3"
            RISK="Normal (1x)"
        elif [[ "$RS429358_GT" == "C/C" ]] && [[ "$RS7412_GT" == "C/C" ]]; then
            APOE_GT="ε4/ε4"
            RISK="VERY HIGH (12-15x)"
        elif [[ "$RS429358_GT" == "T/T" ]] && [[ "$RS7412_GT" == "T/T" ]]; then
            APOE_GT="ε2/ε2"
            RISK="Protective (0.5x)"
        elif [[ "$RS429358_GT" == "T/C" ]] && [[ "$RS7412_GT" == "C/C" ]]; then
            APOE_GT="ε3/ε4"
            RISK="Increased (3-4x)"
        elif [[ "$RS429358_GT" == "T/T" ]] && [[ "$RS7412_GT" == "C/T" ]]; then
            APOE_GT="ε2/ε3"
            RISK="Protective (0.7x)"
        elif [[ "$RS429358_GT" == "T/C" ]] && [[ "$RS7412_GT" == "C/T" ]]; then
            APOE_GT="ε2/ε4"
            RISK="Normal (1x)"
        else
            APOE_GT="Undetermined"
            RISK="Unknown"
        fi
        
        # Save result
        cat > "$OUTPUT_BASE/results/${SAMPLE_ID}_result.txt" << EOF
Sample: $SAMPLE_ID
Date: $(date)
Status: SUCCESS

Coverage:
- APOE region: ${AVG_COV}x
- rs429358: ${RS429358_COV}x  
- rs7412: ${RS7412_COV}x

Genotypes:
- rs429358: $RS429358_GT
- rs7412: $RS7412_GT

APOE Genotype: $APOE_GT
Alzheimer's Risk: $RISK
EOF
        
        echo "[$(date '+%H:%M:%S')] ✓ COMPLETE: $APOE_GT"
        
        # Cleanup to save space
        rm -f "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam"
        rm -f "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam.bai"
        
    } > "$LOG" 2>&1
    
    if [[ $? -eq 0 ]]; then
        SUCCESS=$((SUCCESS + 1))
        RESULT=$(grep "APOE Genotype:" "$OUTPUT_BASE/results/${SAMPLE_ID}_result.txt" | cut -d: -f2 | xargs)
        echo -e "${GREEN}✓ $RESULT${NC}"
    else
        FAILED=$((FAILED + 1))
        echo -e "${RED}✗ FAILED - see $LOG${NC}"
    fi
done

echo ""
echo -e "${BLUE}=========================================="
echo "COMPLETE!"
echo "==========================================${NC}"
echo -e "${GREEN}Success: $SUCCESS${NC} / ${RED}Failed: $FAILED${NC} / Total: $TOTAL"
echo ""
echo "Results: $OUTPUT_BASE/results/"


