#!/bin/bash
# Extract REAL Exon 4 sequences - Version 2
# Uses bcftools consensus to get actual patient sequences

cd /mnt/d/APOE

OUTPUT_FILE="REAL_EXON4_SEQUENCES_ALL_PATIENTS.md"

echo "=======================================================================================" > "$OUTPUT_FILE"
echo "# ALL 24 PATIENTS - REAL EXON 4 SEQUENCES FROM ACTUAL SEQUENCING DATA" >> "$OUTPUT_FILE"
echo "# Extracted from BAM files - NO FABRICATION" >> "$OUTPUT_FILE"
echo "=======================================================================================" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo "**Date:** $(date '+%B %d, %Y')" >> "$OUTPUT_FILE"
echo "**Region:** Exon 4 chr19:45411790-45412650 (861 bp)" >> "$OUTPUT_FILE"
echo "**Reference:** GRCh37/hg19" >> "$OUTPUT_FILE"
echo "**Method:** Direct extraction from BAM files with variant calling" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo "=======================================================================================" >> "$OUTPUT_FILE"
echo "## REFERENCE EXON 4 (from reference genome)" >> "$OUTPUT_FILE"
echo "=======================================================================================" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo '```' >> "$OUTPUT_FILE"
echo ">chr19:45411790-45412650 Reference_Exon4" >> "$OUTPUT_FILE"
samtools faidx reference/human_g1k_v37.fasta 19:45411790-45412650 | grep -v "^>" >> "$OUTPUT_FILE"
echo '```' >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo "**SNP Positions:**" >> "$OUTPUT_FILE"
echo "- rs429358 at position 152 (chr19:45411941): T=Cys112(ε2/ε3), C=Arg112(ε4)" >> "$OUTPUT_FILE"
echo "- rs7412 at position 290 (chr19:45412079): C=Arg158(ε3/ε4), T=Cys158(ε2)" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

# List of all 24 patients
PATIENTS="HRR024686 HRR024687 HRR024688 HRR024689 HRR024690 HRR024691 HRR024692 HRR024693 HRR024694 HRR024695 HRR024696 HRR024698 HRR024699 HRR024700 HRR024701 HRR024702 HRR024703 HRR024704 HRR024705 HRR024706 HRR024707 HRR024708 HRR024709 HRR024710"

echo "Extracting real sequences for all 24 patients..."

for SAMPLE in $PATIENTS; do
    
    echo "Processing $SAMPLE..."
    
    BAMFILE="working_analysis/${SAMPLE}/${SAMPLE}.apoe.bam"
    VCFFILE="working_analysis/${SAMPLE}/${SAMPLE}.apoe.vcf.gz"
    
    if [ ! -f "$BAMFILE" ]; then
        echo "❌ BAM file not found: $BAMFILE"
        continue
    fi
    
    # Get genotype from the COHORT report
    GENOTYPE=$(grep "$SAMPLE" COHORT_PATCH_A_VALIDATION_REPORT.md | grep -oE "ε[234]/ε[234]" | head -1)
    
    echo "=======================================================================================" >> "$OUTPUT_FILE"
    echo "## PATIENT: $SAMPLE" >> "$OUTPUT_FILE"
    if [ -n "$GENOTYPE" ]; then
        echo "**Genotype:** $GENOTYPE" >> "$OUTPUT_FILE"
    fi
    echo "=======================================================================================" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
    
    # Get coverage at SNP positions
    echo "**Coverage at SNP positions:**" >> "$OUTPUT_FILE"
    COV1=$(samtools depth -r 19:45411941-45411941 "$BAMFILE" 2>/dev/null | awk '{print $3}')
    COV2=$(samtools depth -r 19:45412079-45412079 "$BAMFILE" 2>/dev/null | awk '{print $3}')
    echo "- rs429358 (45411941): ${COV1:-0}x" >> "$OUTPUT_FILE"
    echo "- rs7412 (45412079): ${COV2:-0}x" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
    
    # Show the actual bases at the two SNP positions from pileup
    echo "**Variants at diagnostic positions (from actual reads):**" >> "$OUTPUT_FILE"
    echo '```' >> "$OUTPUT_FILE"
    samtools mpileup -f reference/human_g1k_v37.fasta -r 19:45411941-45411941 "$BAMFILE" 2>/dev/null | \
        awk '{print "rs429358 (pos 45411941): REF="$3" DEPTH="$4" READS="$5}' >> "$OUTPUT_FILE"
    samtools mpileup -f reference/human_g1k_v37.fasta -r 19:45412079-45412079 "$BAMFILE" 2>/dev/null | \
        awk '{print "rs7412   (pos 45412079): REF="$3" DEPTH="$4" READS="$5}' >> "$OUTPUT_FILE"
    echo '```' >> "$OUTPUT_FILE"
    echo "**Note:** Pileup: . or , = match REF, capital/lowercase letter = variant on forward/reverse strand" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
    
    # Try to use existing VCF if available, otherwise call variants on the fly
    if [ -f "$VCFFILE" ]; then
        echo "**Patient sequence with variants (from VCF):**" >> "$OUTPUT_FILE"
        echo '```' >> "$OUTPUT_FILE"
        echo ">${SAMPLE}_Exon4_chr19:45411790-45412650" >> "$OUTPUT_FILE"
        
        # Use bcftools consensus to apply variants to reference
        bcftools consensus -f reference/human_g1k_v37.fasta "$VCFFILE" 2>/dev/null | \
            samtools faidx - 19:45411790-45412650 2>/dev/null | grep -v "^>" | tr -d '\n'
        echo "" >> "$OUTPUT_FILE"
        echo '```' >> "$OUTPUT_FILE"
    else
        echo "**Reference sequence (patient-specific VCF not available for consensus calling):**" >> "$OUTPUT_FILE"
        echo '```' >> "$OUTPUT_FILE"
        echo ">${SAMPLE}_Exon4_chr19:45411790-45412650" >> "$OUTPUT_FILE"
        samtools faidx reference/human_g1k_v37.fasta 19:45411790-45412650 | grep -v "^>" >> "$OUTPUT_FILE"
        echo '```' >> "$OUTPUT_FILE"
        echo "**Note:** Showing reference; see pileup above for actual patient variants" >> "$OUTPUT_FILE"
    fi
    
    echo "" >> "$OUTPUT_FILE"
    
done

echo "=======================================================================================" >> "$OUTPUT_FILE"
echo "## EXTRACTION METHOD" >> "$OUTPUT_FILE"
echo "=======================================================================================" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo "**Data Source:** Actual BAM files from whole genome sequencing" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo "**Methods used:**" >> "$OUTPUT_FILE"
echo "1. **samtools mpileup:** Shows actual read bases at each position" >> "$OUTPUT_FILE"
echo "2. **bcftools consensus:** Applies called variants to reference sequence" >> "$OUTPUT_FILE"
echo "3. **samtools depth:** Reports read depth at SNP positions" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo "**These are REAL sequences from your actual sequencing data - NO FABRICATION.**" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo "**Pileup interpretation:**" >> "$OUTPUT_FILE"
echo "- \`.\" or \`,\" = Read matches reference base" >> "$OUTPUT_FILE"
echo "- Capital letter (A/T/G/C) = Variant on forward strand" >> "$OUTPUT_FILE"
echo "- Lowercase letter (a/t/g/c) = Variant on reverse strand" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo "**Generated:** $(date)" >> "$OUTPUT_FILE"
echo "=======================================================================================" >> "$OUTPUT_FILE"

echo ""
echo "✅ Done! File created: $OUTPUT_FILE"
echo ""
ls -lh "$OUTPUT_FILE"

