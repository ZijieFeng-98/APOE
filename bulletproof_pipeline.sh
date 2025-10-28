#!/bin/bash

################################################################################
# BULLETPROOF APOE PIPELINE - Direct to BAM, robust error handling
################################################################################

set -u

SAMPLE_DIR="/mnt/e/CGGA WESeq"
OUTPUT_BASE="/mnt/d/APOE/final_analysis"
REF_GENOME="/mnt/d/APOE/reference/human_g1k_v37.fasta"
THREADS=12
APOE_REGION="19:45409039-45412650"
RS429358_POS="45411941"
RS7412_POS="45412079"

mkdir -p "$OUTPUT_BASE"/{results,logs,temp}

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
echo -e "${YELLOW}BULLETPROOF APOE PIPELINE${NC}"
echo "=========================================="
echo "Samples: $TOTAL | Threads: $THREADS"
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
        FAILED=$((FAILED + 1))
        continue
    fi
    
    SAMPLE_OUT="$OUTPUT_BASE/$SAMPLE_ID"
    mkdir -p "$SAMPLE_OUT"
    LOG="$OUTPUT_BASE/logs/${SAMPLE_ID}.log"
    
    (
        echo "[$(date '+%H:%M:%S')] START: $SAMPLE_ID"
        
        # Direct BWA to sorted BAM (streaming pipeline)
        echo "[$(date '+%H:%M:%S')] Aligning and sorting..."
        
        bwa mem -t $THREADS \
            -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA" \
            "$REF_GENOME" "$FASTQ_R1" "$FASTQ_R2" 2>&1 | \
        samtools view -@ 4 -b -u - 2>&1 | \
        samtools sort -@ 4 -m 4G -T "$OUTPUT_BASE/temp/${SAMPLE_ID}" \
            -o "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" - 2>&1
        
        if [[ $? -ne 0 ]] || [[ ! -s "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" ]]; then
            echo "ERROR: Pipeline failed"
            exit 1
        fi
        
        echo "[$(date '+%H:%M:%S')] Indexing..."
        samtools index -@ $THREADS "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" 2>&1
        
        echo "[$(date '+%H:%M:%S')] Extracting APOE..."
        samtools view -@ $THREADS -b \
            "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam" "$APOE_REGION" \
            > "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>&1
        
        if [[ ! -s "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" ]]; then
            echo "ERROR: APOE extraction failed"
            exit 1
        fi
        
        samtools index "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>&1
        
        echo "[$(date '+%H:%M:%S')] Genotyping..."
        
        # Coverage
        AVG_COV=$(samtools depth "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>/dev/null | \
            awk '{sum+=$3; n++} END {if(n>0) printf "%.2f", sum/n; else print "0"}')
        
        # Get genotypes at SNP positions
        RS429358_DATA=$(samtools mpileup -f "$REF_GENOME" \
            -r "19:${RS429358_POS}-${RS429358_POS}" \
            "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>/dev/null | head -1)
        
        RS7412_DATA=$(samtools mpileup -f "$REF_GENOME" \
            -r "19:${RS7412_POS}-${RS7412_POS}" \
            "$SAMPLE_OUT/${SAMPLE_ID}.apoe.bam" 2>/dev/null | head -1)
        
        if [[ -z "$RS429358_DATA" ]] || [[ -z "$RS7412_DATA" ]]; then
            echo "ERROR: No coverage at SNP positions"
            exit 1
        fi
        
        RS429358_COV=$(echo "$RS429358_DATA" | awk '{print $4}')
        RS429358_BASES=$(echo "$RS429358_DATA" | awk '{print $5}')
        RS7412_COV=$(echo "$RS7412_DATA" | awk '{print $4}')
        RS7412_BASES=$(echo "$RS7412_DATA" | awk '{print $5}')
        
        # Genotype determination
        # rs429358: Reference=T, Alt=C (ε4 allele)
        # rs7412: Reference=C, Alt=T (ε2 allele)
        
        # Count alleles in rs429358
        C_429358=$(echo "$RS429358_BASES" | grep -o "[Cc]" | wc -l)
        REF_429358=$(echo "$RS429358_BASES" | grep -o "[\.,]" | wc -l)
        TOTAL_429358=$((C_429358 + REF_429358))
        
        if (( TOTAL_429358 == 0 )); then
            RS429358_GT="No coverage"
        elif (( C_429358 == TOTAL_429358 )); then
            RS429358_GT="C/C"
        elif (( C_429358 > 0 )); then
            RS429358_GT="T/C"
        else
            RS429358_GT="T/T"
        fi
        
        # Count alleles in rs7412
        T_7412=$(echo "$RS7412_BASES" | grep -o "[Tt]" | wc -l)
        REF_7412=$(echo "$RS7412_BASES" | grep -o "[\.,]" | wc -l)
        TOTAL_7412=$((T_7412 + REF_7412))
        
        if (( TOTAL_7412 == 0 )); then
            RS7412_GT="No coverage"
        elif (( T_7412 == TOTAL_7412 )); then
            RS7412_GT="T/T"
        elif (( T_7412 > 0 )); then
            RS7412_GT="C/T"
        else
            RS7412_GT="C/C"
        fi
        
        # Determine APOE allele based on haplotype
        # ε2 = T at rs429358, T at rs7412 (112C, 158C in protein)
        # ε3 = T at rs429358, C at rs7412 (112C, 158R in protein) - REFERENCE
        # ε4 = C at rs429358, C at rs7412 (112R, 158R in protein)
        
        if [[ "$RS429358_GT" == "T/T" ]] && [[ "$RS7412_GT" == "C/C" ]]; then
            APOE_GT="ε3/ε3"
            RISK="Baseline (1.0x)"
        elif [[ "$RS429358_GT" == "C/C" ]] && [[ "$RS7412_GT" == "C/C" ]]; then
            APOE_GT="ε4/ε4"
            RISK="VERY HIGH (12-15x)"
        elif [[ "$RS429358_GT" == "T/T" ]] && [[ "$RS7412_GT" == "T/T" ]]; then
            APOE_GT="ε2/ε2"
            RISK="Protective (0.4x)"
        elif [[ "$RS429358_GT" == "T/C" ]] && [[ "$RS7412_GT" == "C/C" ]]; then
            APOE_GT="ε3/ε4"
            RISK="Increased (3-4x)"
        elif [[ "$RS429358_GT" == "T/T" ]] && [[ "$RS7412_GT" == "C/T" ]]; then
            APOE_GT="ε2/ε3"
            RISK="Slightly protective (0.7x)"
        elif [[ "$RS429358_GT" == "T/C" ]] && [[ "$RS7412_GT" == "C/T" ]]; then
            APOE_GT="ε2/ε4"
            RISK="Baseline (1.0x)"
        else
            APOE_GT="Undetermined (${RS429358_GT}, ${RS7412_GT})"
            RISK="Cannot determine"
        fi
        
        # Save result
        cat > "$OUTPUT_BASE/results/${SAMPLE_ID}_result.txt" << EOF
================================================================================
APOE GENOTYPING REPORT - CGGA WESeq Cohort
================================================================================
Sample ID: $SAMPLE_ID
Analysis Date: $(date)
Pipeline: Bulletproof v1.0

COVERAGE
--------------------------------------------------------------------------------
APOE region (chr19:45409039-45412650): ${AVG_COV}x
rs429358 (chr19:45411941): ${RS429358_COV}x
rs7412 (chr19:45412079): ${RS7412_COV}x

GENOTYPES
--------------------------------------------------------------------------------
rs429358 (Arg112Cys): $RS429358_GT
rs7412 (Cys158Arg): $RS7412_GT

APOE GENOTYPE: $APOE_GT
Alzheimer's Disease Risk: $RISK

RAW DATA
--------------------------------------------------------------------------------
rs429358 bases: $RS429358_BASES
rs7412 bases: $RS7412_BASES

FILES
--------------------------------------------------------------------------------
Full BAM: ${SAMPLE_OUT}/${SAMPLE_ID}.sorted.bam
APOE BAM: ${SAMPLE_OUT}/${SAMPLE_ID}.apoe.bam
================================================================================
EOF
        
        echo "[$(date '+%H:%M:%S')] ✓ COMPLETE: $APOE_GT ($RISK)"
        
        # Cleanup full genome BAM to save space
        rm -f "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam"
        rm -f "$SAMPLE_OUT/${SAMPLE_ID}.sorted.bam.bai"
        
        exit 0
        
    ) > "$LOG" 2>&1
    
    if [[ $? -eq 0 ]] && [[ -f "$OUTPUT_BASE/results/${SAMPLE_ID}_result.txt" ]]; then
        SUCCESS=$((SUCCESS + 1))
        RESULT=$(grep "APOE GENOTYPE:" "$OUTPUT_BASE/results/${SAMPLE_ID}_result.txt" | cut -d: -f2 | xargs)
        echo -e "${GREEN}✓ $RESULT${NC}"
    else
        FAILED=$((FAILED + 1))
        echo -e "${RED}✗ FAILED (see $LOG)${NC}"
    fi
    
    # Progress update
    echo -e "${BLUE}Progress: $SUCCESS success, $FAILED failed, $((TOTAL - NUM)) remaining${NC}"
done

echo ""
echo -e "${BLUE}=========================================="
echo "ANALYSIS COMPLETE!"
echo "==========================================${NC}"
echo -e "Success: ${GREEN}$SUCCESS${NC} | Failed: ${RED}$FAILED${NC} | Total: $TOTAL"
echo ""

# Generate cohort summary
echo "Generating cohort summary..."

cat > "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt" << EOF
================================================================================
CGGA WESeq COHORT - APOE GENOTYPING ANALYSIS
================================================================================
Analysis Date: $(date)
Total Samples Analyzed: $TOTAL
Successful: $SUCCESS
Failed: $FAILED
Pipeline: Bulletproof v1.0 (12 threads, streaming BWA→BAM)

================================================================================
GENOTYPE DISTRIBUTION
================================================================================
EOF

for genotype in "ε2/ε2" "ε2/ε3" "ε2/ε4" "ε3/ε3" "ε3/ε4" "ε4/ε4"; do
    COUNT=$(grep -l "APOE GENOTYPE: $genotype" "$OUTPUT_BASE/results/"*.txt 2>/dev/null | wc -l)
    if (( COUNT > 0 )); then
        PERCENT=$(awk "BEGIN {printf \"%.1f\", ($COUNT/$SUCCESS)*100}")
        echo "$genotype: $COUNT patients (${PERCENT}%)" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"
    fi
done

echo "" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"
echo "RISK STRATIFICATION" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"
echo "================================================================================" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"

HIGH_RISK=$(grep -l "ε4/ε4" "$OUTPUT_BASE/results/"*.txt 2>/dev/null | wc -l)
INCREASED_RISK=$(grep -l "ε3/ε4" "$OUTPUT_BASE/results/"*.txt 2>/dev/null | wc -l)
BASELINE=$(grep -E "ε3/ε3|ε2/ε4" "$OUTPUT_BASE/results/"*.txt 2>/dev/null | grep -c "APOE GENOTYPE:")
PROTECTIVE=$(grep -E "ε2/ε2|ε2/ε3" "$OUTPUT_BASE/results/"*.txt 2>/dev/null | grep -c "APOE GENOTYPE:")

echo "Very High Risk (ε4/ε4): $HIGH_RISK" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"
echo "Increased Risk (ε3/ε4): $INCREASED_RISK" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"
echo "Baseline Risk (ε3/ε3, ε2/ε4): $BASELINE" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"
echo "Protective (ε2/ε2, ε2/ε3): $PROTECTIVE" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"

echo "" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"
echo "INDIVIDUAL RESULTS" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"
echo "================================================================================" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"

for result_file in "$OUTPUT_BASE/results/"*_result.txt; do
    if [[ -f "$result_file" ]]; then
        SAMPLE=$(basename "$result_file" _result.txt)
        GENOTYPE=$(grep "APOE GENOTYPE:" "$result_file" | cut -d: -f2 | xargs)
        RISK=$(grep "Alzheimer's Disease Risk:" "$result_file" | cut -d: -f2 | xargs)
        COV=$(grep "APOE region" "$result_file" | awk '{print $NF}')
        echo "$SAMPLE: $GENOTYPE - $RISK (Coverage: $COV)" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"
    fi
done

echo "" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"
echo "================================================================================" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"
echo "Results directory: $OUTPUT_BASE/results/" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"
echo "Logs directory: $OUTPUT_BASE/logs/" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"
echo "================================================================================" >> "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"

echo ""
echo "Cohort summary: $OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"
cat "$OUTPUT_BASE/COHORT_SUMMARY_REPORT.txt"


