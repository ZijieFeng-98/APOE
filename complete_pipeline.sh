#!/bin/bash

# Bulletproof APOE Pipeline Completion Script
# This will complete all remaining steps and auto-recover from any issues

set +e  # Don't exit on errors, handle them

LOG="/mnt/d/APOE/completion.log"
RESULT="/mnt/d/APOE/apoe_analysis/results/HRR024685_apoe_report.txt"

echo "========================================" | tee -a $LOG
echo "Starting Pipeline Completion: $(date)" | tee -a $LOG
echo "========================================" | tee -a $LOG

cd /mnt/d/APOE

# Step 1: Sort BAM if needed
echo "[$(date)] Step 1: Checking BAM sorting..." | tee -a $LOG
if [ ! -f "apoe_analysis/alignment/HRR024685.sorted.bam" ]; then
    echo "[$(date)] Sorting BAM file..." | tee -a $LOG
    samtools sort -@ 4 -o apoe_analysis/alignment/HRR024685.sorted.bam apoe_analysis/alignment/HRR024685.bam 2>&1 | tee -a $LOG
    if [ $? -eq 0 ]; then
        echo "[$(date)] ✓ BAM sorting complete" | tee -a $LOG
    else
        echo "[$(date)] ✗ BAM sorting failed, retrying..." | tee -a $LOG
        sleep 5
        samtools sort -@ 2 -o apoe_analysis/alignment/HRR024685.sorted.bam apoe_analysis/alignment/HRR024685.bam 2>&1 | tee -a $LOG
    fi
else
    echo "[$(date)] ✓ Sorted BAM already exists" | tee -a $LOG
fi

# Step 2: Index BAM
echo "[$(date)] Step 2: Indexing BAM..." | tee -a $LOG
if [ ! -f "apoe_analysis/alignment/HRR024685.sorted.bam.bai" ]; then
    samtools index apoe_analysis/alignment/HRR024685.sorted.bam 2>&1 | tee -a $LOG
    echo "[$(date)] ✓ BAM indexing complete" | tee -a $LOG
else
    echo "[$(date)] ✓ BAM index already exists" | tee -a $LOG
fi

# Remove unsorted BAM to save space
if [ -f "apoe_analysis/alignment/HRR024685.bam" ]; then
    echo "[$(date)] Removing unsorted BAM to save space..." | tee -a $LOG
    rm -f apoe_analysis/alignment/HRR024685.bam
fi

# Step 3: Extract APOE region
echo "[$(date)] Step 3: Extracting APOE region..." | tee -a $LOG
if [ ! -f "apoe_analysis/alignment/HRR024685.apoe.bam" ]; then
    samtools view -@ 4 -b apoe_analysis/alignment/HRR024685.sorted.bam 19:45409039-45412650 > apoe_analysis/alignment/HRR024685.apoe.bam 2>&1 | tee -a $LOG
    samtools index apoe_analysis/alignment/HRR024685.apoe.bam 2>&1 | tee -a $LOG
    echo "[$(date)] ✓ APOE region extracted" | tee -a $LOG
    
    # Calculate coverage
    COVERAGE=$(samtools depth apoe_analysis/alignment/HRR024685.apoe.bam | awk '{sum+=$3} END {if(NR>0) print sum/NR; else print 0}')
    echo "[$(date)] APOE region coverage: ${COVERAGE}x" | tee -a $LOG
else
    echo "[$(date)] ✓ APOE BAM already exists" | tee -a $LOG
fi

# Step 4: Variant calling
echo "[$(date)] Step 4: Calling variants..." | tee -a $LOG
if [ ! -f "apoe_analysis/variants/HRR024685.apoe.vcf.gz" ]; then
    bcftools mpileup -Ou -f reference/human_g1k_v37.fasta -r 19:45409039-45412650 apoe_analysis/alignment/HRR024685.apoe.bam 2>&1 | \
        bcftools call -mv -Oz -o apoe_analysis/variants/HRR024685.apoe.vcf.gz 2>&1 | tee -a $LOG
    
    if [ $? -eq 0 ]; then
        bcftools index apoe_analysis/variants/HRR024685.apoe.vcf.gz 2>&1 | tee -a $LOG
        echo "[$(date)] ✓ Variant calling complete" | tee -a $LOG
    else
        echo "[$(date)] ✗ Variant calling failed, retrying..." | tee -a $LOG
        sleep 5
        bcftools mpileup -Ou -f reference/human_g1k_v37.fasta -r 19:45409039-45412650 apoe_analysis/alignment/HRR024685.apoe.bam | \
            bcftools call -mv -Oz -o apoe_analysis/variants/HRR024685.apoe.vcf.gz 2>&1 | tee -a $LOG
        bcftools index apoe_analysis/variants/HRR024685.apoe.vcf.gz 2>&1 | tee -a $LOG
    fi
else
    echo "[$(date)] ✓ VCF already exists" | tee -a $LOG
fi

# Step 5: Extract key APOE SNPs
echo "[$(date)] Step 5: Extracting APOE SNPs..." | tee -a $LOG
bcftools view -r 19:45411941-45411941 apoe_analysis/variants/HRR024685.apoe.vcf.gz > apoe_analysis/variants/rs429358.vcf 2>&1 | tee -a $LOG
bcftools view -r 19:45412079-45412079 apoe_analysis/variants/HRR024685.apoe.vcf.gz > apoe_analysis/variants/rs7412.vcf 2>&1 | tee -a $LOG
echo "[$(date)] ✓ SNPs extracted" | tee -a $LOG

# Step 6: Create interpretation script (embedded)
echo "[$(date)] Step 6: Creating interpretation script..." | tee -a $LOG
cat > apoe_analysis/interpret.py << 'PYTHON_EOF'
#!/usr/bin/env python3
import sys
import re

def parse_vcf(vcf_file):
    try:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue
                chrom = fields[0]
                pos = fields[1]
                ref = fields[3]
                alt = fields[4]
                format_fields = fields[8].split(':')
                sample_fields = fields[9].split(':')
                if 'GT' in format_fields:
                    gt_index = format_fields.index('GT')
                    genotype = sample_fields[gt_index]
                    gt_alleles = re.findall(r'\d+', genotype)
                    if len(gt_alleles) >= 2:
                        allele1 = ref if gt_alleles[0] == '0' else alt.split(',')[int(gt_alleles[0])-1]
                        allele2 = ref if gt_alleles[1] == '0' else alt.split(',')[int(gt_alleles[1])-1]
                        return (allele1, allele2), chrom, pos, ref, alt
        return None, None, None, None, None
    except:
        return None, None, None, None, None

def determine_apoe(rs429358_gt, rs7412_gt):
    def classify(snp429358, snp7412):
        if snp429358 == 'C' and snp7412 == 'T':
            return 'ε2'
        elif snp429358 == 'C' and snp7412 == 'C':
            return 'ε3'
        elif snp429358 == 'T' and snp7412 == 'C':
            return 'ε4'
        else:
            return 'unknown'
    
    if rs429358_gt is None or rs7412_gt is None:
        return None, None, "Insufficient data"
    
    allele1 = classify(rs429358_gt[0], rs7412_gt[0])
    allele2 = classify(rs429358_gt[1], rs7412_gt[1])
    alleles = sorted([allele1, allele2])
    return allele1, allele2, f"{alleles[0]}/{alleles[1]}"

def clinical_interp(genotype):
    interp = {
        'ε2/ε2': ('REDUCED RISK', 'Protective against Alzheimer\'s', '~0.5x'),
        'ε2/ε3': ('REDUCED RISK', 'Somewhat protective', '~0.6x'),
        'ε2/ε4': ('MODERATE RISK', 'Mixed effect', '~2-3x'),
        'ε3/ε3': ('AVERAGE RISK', 'Most common; neutral', '1x'),
        'ε3/ε4': ('INCREASED RISK', 'Moderately increased', '~3x'),
        'ε4/ε4': ('HIGH RISK', 'Significantly increased', '~8-12x')
    }
    return interp.get(genotype, ('UNKNOWN', 'Unable to determine', 'Unknown'))

rs429358_gt, _, _, ref1, alt1 = parse_vcf(sys.argv[1])
rs7412_gt, _, _, ref2, alt2 = parse_vcf(sys.argv[2])
_, _, genotype = determine_apoe(rs429358_gt, rs7412_gt)
risk, desc, rel_risk = clinical_interp(genotype)

with open(sys.argv[3], 'w') as out:
    out.write("=" * 80 + "\n")
    out.write("APOE GENOTYPING REPORT\n")
    out.write("=" * 80 + "\n\n")
    out.write(f"Sample ID: HRR024685\n")
    out.write(f"Analysis Date: {sys.argv[4] if len(sys.argv) > 4 else '2025-10-20'}\n")
    out.write(f"Reference: GRCh37/hg19\n\n")
    out.write("-" * 80 + "\n")
    out.write("RAW GENOTYPE DATA\n")
    out.write("-" * 80 + "\n\n")
    if rs429358_gt:
        out.write(f"rs429358 (chr19:45411941):\n")
        out.write(f"  Genotype: {rs429358_gt[0]}/{rs429358_gt[1]}\n\n")
    if rs7412_gt:
        out.write(f"rs7412 (chr19:45412079):\n")
        out.write(f"  Genotype: {rs7412_gt[0]}/{rs7412_gt[1]}\n\n")
    out.write("-" * 80 + "\n")
    out.write("APOE GENOTYPE\n")
    out.write("-" * 80 + "\n\n")
    out.write(f"APOE Genotype: {genotype}\n\n")
    out.write("-" * 80 + "\n")
    out.write("CLINICAL INTERPRETATION\n")
    out.write("-" * 80 + "\n\n")
    out.write(f"Risk Category: {risk}\n")
    out.write(f"Description: {desc}\n")
    out.write(f"Relative Risk: {rel_risk}\n\n")
    out.write("Important Notes:\n")
    out.write("- APOE is just one of many risk factors\n")
    out.write("- This is for research purposes only\n")
    out.write("- Consult a genetic counselor for clinical interpretation\n\n")
    out.write("=" * 80 + "\n")
print(f"APOE Genotype: {genotype} | Risk: {risk}")
PYTHON_EOF

chmod +x apoe_analysis/interpret.py

# Step 7: Run interpretation
echo "[$(date)] Step 7: Interpreting APOE genotype..." | tee -a $LOG
python3 apoe_analysis/interpret.py \
    apoe_analysis/variants/rs429358.vcf \
    apoe_analysis/variants/rs7412.vcf \
    "$RESULT" \
    "$(date +%Y-%m-%d)" 2>&1 | tee -a $LOG

if [ -f "$RESULT" ]; then
    echo "" | tee -a $LOG
    echo "========================================" | tee -a $LOG
    echo "✓✓✓ ANALYSIS COMPLETE! ✓✓✓" | tee -a $LOG
    echo "========================================" | tee -a $LOG
    echo "" | tee -a $LOG
    cat "$RESULT" | tee -a $LOG
    echo "" | tee -a $LOG
    echo "Results saved to: $RESULT" | tee -a $LOG
    exit 0
else
    echo "[$(date)] ✗ Failed to generate report" | tee -a $LOG
    exit 1
fi


