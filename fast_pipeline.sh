#!/bin/bash

################################################################################
# HIGH-PERFORMANCE APOE GENOTYPING PIPELINE
# Uses ALL available CPU cores for maximum speed
################################################################################

set -e
set -u

SAMPLE_DIR="/mnt/e/CGGA WESeq"
OUTPUT_BASE="/mnt/d/APOE/fast_analysis"
REF_GENOME="/mnt/d/APOE/reference/human_g1k_v37.fasta"
THREADS=12  # MAXIMUM THREADS FOR SPEED
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
YELLOW='\033[1;33m'
NC='\033[0m'

# All samples except HRR024685 and HRR024697 (already done)
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
echo -e "${YELLOW}HIGH-PERFORMANCE APOE ANALYSIS${NC}"
echo "=========================================="
echo "Samples: $TOTAL"
echo "Threads per sample: $THREADS"
echo "Mode: MAXIMUM SPEED"
echo "Start: $(date)"
echo "=========================================="
echo ""

for i in "${!SAMPLES[@]}"; do
    SAMPLE_ID="${SAMPLES[$i]}"
    NUM=$((i + 1))
    
    echo -e "${BLUE}===========================================${NC}"
    echo -e "${BLUE}[$NUM/$TOTAL] $SAMPLE_ID - $(date)${NC}"
    echo -e "${BLUE}===========================================${NC}"
    
    FASTQ_R1="${SAMPLE_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_f1.fq.gz"
    FASTQ_R2="${SAMPLE_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_r2.fq.gz"
    
    if [ ! -f "$FASTQ_R1" ] || [ ! -f "$FASTQ_R2" ]; then
        echo -e "${RED}✗ FASTQ files not found${NC}"
        FAILED=$((FAILED + 1))
        continue
    fi
    
    SAMPLE_OUT="$OUTPUT_BASE/$SAMPLE_ID"
    mkdir -p "$SAMPLE_OUT"
    LOG_FILE="$OUTPUT_BASE/logs/${SAMPLE_ID}.log"
    
    (
        echo "▶ Aligning with BWA (${THREADS} threads)..."
        bwa mem -t $THREADS -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA" \
            "$REF_GENOME" "$FASTQ_R1" "$FASTQ_R2" 2>&1 | \
            samtools view -@ $THREADS -b - | \
            samtools sort -@ $THREADS -m 2G -o "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" - 2>&1
        
        if [ ! -f "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" ] || [ ! -s "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" ]; then
            echo "ERROR: BAM creation failed"
            exit 1
        fi
        
        echo "▶ Indexing..."
        samtools index -@ $THREADS "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" 2>&1
        
        echo "▶ Extracting APOE region..."
        samtools view -@ $THREADS -b "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" "$APOE_REGION" > \
            "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>&1
        samtools index "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>&1
        
        echo "▶ Calling variants with BCFtools..."
        bcftools mpileup -f "$REF_GENOME" -r "$APOE_REGION" \
            -Ou "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>/dev/null | \
            bcftools call -mv -Ov -o "$SAMPLE_OUT/${SAMPLE_ID}.apoe.vcf" 2>&1
        
        # Get coverage
        AVG_COV=$(samtools depth "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" | \
            awk '{sum+=$3; n++} END {if(n>0) printf "%.2f", sum/n; else print "0"}')
        
        # Get genotypes at key SNPs
        RS429358_PILEUP=$(samtools mpileup -f "$REF_GENOME" \
            -r "19:${RS429358_POS}-${RS429358_POS}" \
            "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>/dev/null | head -1)
        
        RS7412_PILEUP=$(samtools mpileup -f "$REF_GENOME" \
            -r "19:${RS7412_POS}-${RS7412_POS}" \
            "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>/dev/null | head -1)
        
        RS429358_COV=$(echo "$RS429358_PILEUP" | awk '{print $4}')
        RS429358_BASES=$(echo "$RS429358_PILEUP" | awk '{print $5}')
        RS7412_COV=$(echo "$RS7412_PILEUP" | awk '{print $4}')
        RS7412_BASES=$(echo "$RS7412_PILEUP" | awk '{print $5}')
        
        # Try to extract from VCF
        RS429358_VCF=$(grep "^19\s*${RS429358_POS}" "$SAMPLE_OUT/${SAMPLE_ID}.apoe.vcf" 2>/dev/null || echo "")
        RS7412_VCF=$(grep "^19\s*${RS7412_POS}" "$SAMPLE_OUT/${SAMPLE_ID}.apoe.vcf" 2>/dev/null || echo "")
        
        # Determine genotypes
        # Reference: rs429358=T, rs7412=C
        # ε2: rs429358=T/T, rs7412=T/T (C->T mutation)
        # ε3: rs429358=T/T, rs7412=C/C (reference)
        # ε4: rs429358=C/C, rs7412=C/C (T->C mutation)
        
        if [ -n "$RS429358_VCF" ]; then
            RS429358_ALT=$(echo "$RS429358_VCF" | awk '{print $5}')
            RS429358_GT_INFO=$(echo "$RS429358_VCF" | awk '{print $NF}')
            if [[ "$RS429358_GT_INFO" =~ ^1/1 ]]; then
                RS429358_GT="C/C"  # Homozygous alt (ε4)
            elif [[ "$RS429358_GT_INFO" =~ ^0/1 ]] || [[ "$RS429358_GT_INFO" =~ ^1/0 ]]; then
                RS429358_GT="T/C"  # Het
            else
                RS429358_GT="T/T"  # Ref (ε3)
            fi
        elif [[ -n "$RS429358_BASES" ]] && [[ "$RS429358_BASES" =~ ^[\.,]+$ ]]; then
            RS429358_GT="T/T"  # Reference
        elif [[ -n "$RS429358_BASES" ]]; then
            # Count C's in bases
            C_COUNT=$(echo "$RS429358_BASES" | grep -o "[Cc]" | wc -l)
            REF_COUNT=$(echo "$RS429358_BASES" | grep -o "[\.,]" | wc -l)
            if [ $C_COUNT -gt 0 ] && [ $REF_COUNT -eq 0 ]; then
                RS429358_GT="C/C"
            elif [ $C_COUNT -gt 0 ]; then
                RS429358_GT="T/C"
            else
                RS429358_GT="T/T"
            fi
        else
            RS429358_GT="No coverage"
        fi
        
        if [ -n "$RS7412_VCF" ]; then
            RS7412_ALT=$(echo "$RS7412_VCF" | awk '{print $5}')
            RS7412_GT_INFO=$(echo "$RS7412_VCF" | awk '{print $NF}')
            if [[ "$RS7412_GT_INFO" =~ ^1/1 ]]; then
                RS7412_GT="T/T"  # Homozygous alt (ε2)
            elif [[ "$RS7412_GT_INFO" =~ ^0/1 ]] || [[ "$RS7412_GT_INFO" =~ ^1/0 ]]; then
                RS7412_GT="C/T"  # Het
            else
                RS7412_GT="C/C"  # Ref
            fi
        elif [[ -n "$RS7412_BASES" ]] && [[ "$RS7412_BASES" =~ ^[\.,]+$ ]]; then
            RS7412_GT="C/C"  # Reference
        elif [[ -n "$RS7412_BASES" ]]; then
            # Count T's in bases
            T_COUNT=$(echo "$RS7412_BASES" | grep -o "[Tt]" | wc -l)
            REF_COUNT=$(echo "$RS7412_BASES" | grep -o "[\.,]" | wc -l)
            if [ $T_COUNT -gt 0 ] && [ $REF_COUNT -eq 0 ]; then
                RS7412_GT="T/T"
            elif [ $T_COUNT -gt 0 ]; then
                RS7412_GT="C/T"
            else
                RS7412_GT="C/C"
            fi
        else
            RS7412_GT="No coverage"
        fi
        
        # Determine APOE genotype
        if [ "$RS429358_GT" = "T/T" ] && [ "$RS7412_GT" = "C/C" ]; then
            APOE_GT="ε3/ε3"
            RISK="Normal/Baseline (1x)"
        elif [ "$RS429358_GT" = "C/C" ] && [ "$RS7412_GT" = "C/C" ]; then
            APOE_GT="ε4/ε4"
            RISK="VERY HIGH RISK (12-15x)"
        elif [ "$RS429358_GT" = "T/T" ] && [ "$RS7412_GT" = "T/T" ]; then
            APOE_GT="ε2/ε2"
            RISK="Protective (0.5x)"
        elif [ "$RS429358_GT" = "T/C" ] && [ "$RS7412_GT" = "C/C" ]; then
            APOE_GT="ε3/ε4"
            RISK="Increased risk (3-4x)"
        elif [ "$RS429358_GT" = "T/T" ] && [ "$RS7412_GT" = "C/T" ]; then
            APOE_GT="ε2/ε3"
            RISK="Slightly protective (0.7x)"
        elif [ "$RS429358_GT" = "T/C" ] && [ "$RS7412_GT" = "C/T" ]; then
            APOE_GT="ε2/ε4"
            RISK="Normal (1x)"
        else
            APOE_GT="Unknown (${RS429358_GT} / ${RS7412_GT})"
            RISK="Cannot determine"
        fi
        
        # Save result
        cat > "$OUTPUT_BASE/results/${SAMPLE_ID}_result.txt" << EOF
================================================================================
APOE GENOTYPING REPORT
================================================================================
Sample ID: $SAMPLE_ID
Analysis Date: $(date)
Status: SUCCESS

COVERAGE STATISTICS
--------------------------------------------------------------------------------
Average APOE region coverage: ${AVG_COV}x
rs429358 coverage: ${RS429358_COV}x
rs7412 coverage: ${RS7412_COV}x

GENOTYPE RESULTS
--------------------------------------------------------------------------------
rs429358 (Arg112Cys): $RS429358_GT
rs7412 (Cys158Arg): $RS7412_GT

APOE GENOTYPE: $APOE_GT
ALZHEIMER'S RISK: $RISK

RAW DATA
--------------------------------------------------------------------------------
rs429358 pileup: $RS429358_BASES
rs7412 pileup: $RS7412_BASES

FILES GENERATED
--------------------------------------------------------------------------------
Full alignment: ${SAMPLE_OUT}/${SAMPLE_ID}.sorted.bam
APOE region: ${SAMPLE_OUT}/${SAMPLE_ID}.apoe.bam
Variants: ${SAMPLE_OUT}/${SAMPLE_ID}.apoe.vcf
================================================================================
EOF
        
        echo "✓ SUCCESS: $APOE_GT"
        
    ) > "$LOG_FILE" 2>&1
    
    if [ $? -eq 0 ]; then
        SUCCESS=$((SUCCESS + 1))
        RESULT=$(grep "APOE GENOTYPE:" "$OUTPUT_BASE/results/${SAMPLE_ID}_result.txt" | cut -d: -f2 | xargs)
        echo -e "${GREEN}✓ SUCCESS: $RESULT${NC}"
    else
        FAILED=$((FAILED + 1))
        echo -e "${RED}✗ FAILED - check log: $LOG_FILE${NC}"
    fi
    
    # Clean up full BAM to save space
    if [ -f "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" ]; then
        rm -f "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam"
        rm -f "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam.bai"
    fi
    
    echo ""
done

echo -e "${BLUE}=========================================="
echo "ANALYSIS COMPLETE!"
echo "==========================================${NC}"
echo "Total: $TOTAL"
echo -e "${GREEN}Success: $SUCCESS${NC}"
echo -e "${RED}Failed: $FAILED${NC}"
echo ""
echo "Results directory: $OUTPUT_BASE/results/"
echo ""

# Generate summary
cat > "$OUTPUT_BASE/COHORT_SUMMARY.txt" << EOF
================================================================================
CGGA WESeq COHORT - APOE GENOTYPING SUMMARY
================================================================================
Analysis Date: $(date)
Total Samples: $TOTAL
Successful: $SUCCESS
Failed: $FAILED

GENOTYPE DISTRIBUTION
================================================================================
EOF

for gt in "ε2/ε2" "ε2/ε3" "ε2/ε4" "ε3/ε3" "ε3/ε4" "ε4/ε4"; do
    COUNT=$(grep -l "APOE GENOTYPE: $gt" "$OUTPUT_BASE/results/"*.txt 2>/dev/null | wc -l)
    if [ $COUNT -gt 0 ]; then
        echo "$gt: $COUNT patients" >> "$OUTPUT_BASE/COHORT_SUMMARY.txt"
    fi
done

echo "" >> "$OUTPUT_BASE/COHORT_SUMMARY.txt"
echo "INDIVIDUAL RESULTS" >> "$OUTPUT_BASE/COHORT_SUMMARY.txt"
echo "================================================================================" >> "$OUTPUT_BASE/COHORT_SUMMARY.txt"

for result_file in "$OUTPUT_BASE/results/"*.txt; do
    if [ -f "$result_file" ]; then
        SAMPLE=$(basename "$result_file" _result.txt)
        GT=$(grep "APOE GENOTYPE:" "$result_file" | cut -d: -f2 | xargs)
        RISK=$(grep "ALZHEIMER'S RISK:" "$result_file" | cut -d: -f2 | xargs)
        echo "$SAMPLE: $GT - $RISK" >> "$OUTPUT_BASE/COHORT_SUMMARY.txt"
    fi
done

echo ""
echo "Summary report: $OUTPUT_BASE/COHORT_SUMMARY.txt"
cat "$OUTPUT_BASE/COHORT_SUMMARY.txt"


