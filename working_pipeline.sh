#!/bin/bash

################################################################################
# WORKING APOE PIPELINE - Write to disk, no fragile pipes
################################################################################

set -u

SAMPLE_DIR="/mnt/e/CGGA WESeq"
OUTPUT_BASE="/mnt/d/APOE/working_analysis"
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
echo -e "${YELLOW}WORKING APOE PIPELINE${NC}"
echo "=========================================="
echo "Samples: $TOTAL | Threads: $THREADS"
echo "Strategy: Write to disk (no pipes)"
echo "Start: $(date)"
echo "=========================================="

for i in "${!SAMPLES[@]}"; do
    SAMPLE_ID="${SAMPLES[$i]}"
    NUM=$((i + 1))
    
    echo ""
    echo -e "${BLUE}[$NUM/$TOTAL] $SAMPLE_ID - $(date '+%H:%M:%S')${NC}"
    
    FASTQ_R1="${SAMPLE_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_f1.fq.gz"
    FASTQ_R2="${SAMPLE_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_r2.fq.gz"
    
    if [[ ! -f "$FASTQ_R1" ]] || [[ ! -f "$FASTQ_R2" ]]; then
        echo -e "${RED}✗ FASTQ not found${NC}"
        ((FAILED++))
        continue
    fi
    
    SAMPLE_OUT="$OUTPUT_BASE/$SAMPLE_ID"
    mkdir -p "$SAMPLE_OUT"
    LOG="$OUTPUT_BASE/logs/${SAMPLE_ID}.log"
    
    (
        set -e  # Exit on any error
        
        echo "[$(date '+%H:%M:%S')] START: $SAMPLE_ID"
        
        # Step 1: BWA alignment to SAM file
        echo "[$(date '+%H:%M:%S')] Step 1: BWA alignment to SAM..."
        bwa mem -t $THREADS \
            -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA" \
            "$REF_GENOME" "$FASTQ_R1" "$FASTQ_R2" \
            -o "$SAMPLE_OUT/${SAMPLE_ID}.sam" 2>&1
        
        if [[ ! -s "$SAMPLE_OUT/${SAMPLE_ID}.sam" ]]; then
            echo "ERROR: SAM file not created"
            exit 1
        fi
        
        echo "[$(date '+%H:%M:%S')] Step 2: SAM to BAM conversion..."
        samtools view -@ $THREADS -b "$SAMPLE_OUT/${SAMPLE_ID}.sam" \
            -o "$SAMPLE_OUT/${SAMPLE_ID}.bam" 2>&1
        rm -f "$SAMPLE_OUT/${SAMPLE_ID}.sam"  # Free space
        
        echo "[$(date '+%H:%M:%S')] Step 3: Sorting BAM..."
        samtools sort -@ 4 -m 1G \
            "$SAMPLE_OUT/${SAMPLE_ID}.bam" \
            -o "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" 2>&1
        rm -f "$SAMPLE_OUT/${SAMPLE_ID}.bam"  # Free space
        
        echo "[$(date '+%H:%M:%S')] Step 4: Indexing..."
        samtools index -@ $THREADS "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" 2>&1
        
        echo "[$(date '+%H:%M:%S')] Step 5: Extracting APOE region..."
        samtools view -@ $THREADS -b \
            "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" "$APOE_REGION" \
            -o "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>&1
        samtools index "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>&1
        
        echo "[$(date '+%H:%M:%S')] Step 6: Genotyping..."
        
        # Coverage
        AVG_COV=$(samtools depth "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>/dev/null | \
            awk '{sum+=$3; n++} END {if(n>0) printf "%.2f", sum/n; else print "0"}')
        
        # Genotype SNPs
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
        
        # Determine genotypes
        # rs429358: Ref=T, Alt=C (ε4)
        # rs7412: Ref=C, Alt=T (ε2)
        
        C_429=$(echo "$RS429358_BASES" | grep -o "[Cc]" | wc -l)
        REF_429=$(echo "$RS429358_BASES" | grep -o "[\.,]" | wc -l)
        TOT_429=$((C_429 + REF_429))
        
        if [[ $TOT_429 -eq 0 ]]; then
            RS429358_GT="No coverage"
        elif [[ $C_429 -eq $TOT_429 ]]; then
            RS429358_GT="C/C"
        elif [[ $C_429 -gt 0 ]]; then
            RS429358_GT="T/C"
        else
            RS429358_GT="T/T"
        fi
        
        T_7412=$(echo "$RS7412_BASES" | grep -o "[Tt]" | wc -l)
        REF_7412=$(echo "$RS7412_BASES" | grep -o "[\.,]" | wc -l)
        TOT_7412=$((T_7412 + REF_7412))
        
        if [[ $TOT_7412 -eq 0 ]]; then
            RS7412_GT="No coverage"
        elif [[ $T_7412 -eq $TOT_7412 ]]; then
            RS7412_GT="T/T"
        elif [[ $T_7412 -gt 0 ]]; then
            RS7412_GT="C/T"
        else
            RS7412_GT="C/C"
        fi
        
        # APOE genotype
        if [[ "$RS429358_GT" == "T/T" ]] && [[ "$RS7412_GT" == "C/C" ]]; then
            APOE_GT="ε3/ε3"; RISK="Baseline (1.0x)"
        elif [[ "$RS429358_GT" == "C/C" ]] && [[ "$RS7412_GT" == "C/C" ]]; then
            APOE_GT="ε4/ε4"; RISK="VERY HIGH (12-15x)"
        elif [[ "$RS429358_GT" == "T/T" ]] && [[ "$RS7412_GT" == "T/T" ]]; then
            APOE_GT="ε2/ε2"; RISK="Protective (0.4x)"
        elif [[ "$RS429358_GT" == "T/C" ]] && [[ "$RS7412_GT" == "C/C" ]]; then
            APOE_GT="ε3/ε4"; RISK="Increased (3-4x)"
        elif [[ "$RS429358_GT" == "T/T" ]] && [[ "$RS7412_GT" == "C/T" ]]; then
            APOE_GT="ε2/ε3"; RISK="Protective (0.7x)"
        elif [[ "$RS429358_GT" == "T/C" ]] && [[ "$RS7412_GT" == "C/T" ]]; then
            APOE_GT="ε2/ε4"; RISK="Baseline (1.0x)"
        else
            APOE_GT="Undetermined"; RISK="Cannot determine"
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
- rs429358 (Arg112Cys): $RS429358_GT
- rs7412 (Cys158Arg): $RS7412_GT

APOE Genotype: $APOE_GT
Alzheimer's Risk: $RISK

Raw Data:
- rs429358 bases: $RS429358_BASES
- rs7412 bases: $RS7412_BASES
EOF
        
        echo "[$(date '+%H:%M:%S')] ✓ SUCCESS: $APOE_GT ($RISK)"
        
        # Cleanup
        rm -f "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam"
        rm -f "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam.bai"
        
    ) > "$LOG" 2>&1
    
    if [[ $? -eq 0 ]] && [[ -f "$OUTPUT_BASE/results/${SAMPLE_ID}_result.txt" ]]; then
        ((SUCCESS++))
        RESULT=$(grep "APOE Genotype:" "$OUTPUT_BASE/results/${SAMPLE_ID}_result.txt" | cut -d: -f2 | xargs)
        echo -e "${GREEN}✓ $RESULT${NC}"
    else
        ((FAILED++))
        echo -e "${RED}✗ FAILED${NC}"
    fi
    
    echo -e "${BLUE}Progress: $SUCCESS/$TOTAL success, $FAILED failed${NC}"
done

echo ""
echo "=========================================="
echo "COMPLETE!"
echo "=========================================="
echo "Success: $SUCCESS | Failed: $FAILED"

