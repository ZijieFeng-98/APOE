#!/bin/bash

################################################################################
# APOE Genotyping Pipeline for WGS Data
# 
# This script processes FASTQ files from whole genome sequencing to determine
# APOE genotype (ε2, ε3, or ε4) based on two key SNPs:
# - rs429358 (chr19:45411941, C>T)
# - rs7412 (chr19:45412079, C>T)
#
# Author: APOE Analysis Pipeline
# Date: 2025-10-20
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

################################################################################
# Configuration
################################################################################

# Input files (update these paths if needed)
FASTQ_R1="/mnt/d/APOE/HRR024685_f1.fq.gz"
FASTQ_R2="/mnt/d/APOE/HRR024685_r2.fq.gz"
SAMPLE_ID="HRR024685"

# Output directory
OUTPUT_DIR="/mnt/d/APOE/apoe_analysis"

# Reference genome
REF_DIR="/mnt/d/APOE/reference"
REF_GENOME="${REF_DIR}/human_g1k_v37.fasta"
REF_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"

# APOE gene region (GRCh37/hg19)
APOE_CHR="19"
APOE_START="45409039"
APOE_END="45412650"
APOE_REGION="${APOE_CHR}:${APOE_START}-${APOE_END}"

# Key APOE SNPs (GRCh37 coordinates)
RS429358_POS="45411941"  # C>T
RS7412_POS="45412079"    # C>T

# System settings
THREADS=4  # Adjust based on your CPU cores

################################################################################
# Helper Functions
################################################################################

log_info() {
    echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
    exit 1
}

check_tool() {
    if ! command -v $1 &> /dev/null; then
        log_error "$1 is not installed. Please install it first (see INSTALLATION_GUIDE.md)"
    fi
    log_success "$1 is installed: $(which $1)"
}

################################################################################
# Pre-flight Checks
################################################################################

log_info "Starting APOE Genotyping Pipeline"
log_info "Sample ID: ${SAMPLE_ID}"
echo ""

# Check required tools
log_info "Checking required tools..."
check_tool bwa
check_tool samtools
check_tool bcftools
check_tool python3

if command -v fastqc &> /dev/null; then
    log_success "fastqc is installed (optional)"
    RUN_FASTQC=true
else
    log_warning "fastqc not found (optional, skipping QC)"
    RUN_FASTQC=false
fi

echo ""

# Check input files
log_info "Checking input files..."
if [ ! -f "$FASTQ_R1" ]; then
    log_error "FASTQ R1 file not found: $FASTQ_R1"
fi
log_success "Found R1: $FASTQ_R1"

if [ ! -f "$FASTQ_R2" ]; then
    log_error "FASTQ R2 file not found: $FASTQ_R2"
fi
log_success "Found R2: $FASTQ_R2"

echo ""

################################################################################
# Create Output Directories
################################################################################

log_info "Creating output directories..."
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/qc"
mkdir -p "${OUTPUT_DIR}/alignment"
mkdir -p "${OUTPUT_DIR}/variants"
mkdir -p "${OUTPUT_DIR}/results"
mkdir -p "${REF_DIR}"
log_success "Directories created"

echo ""

################################################################################
# Download and Index Reference Genome
################################################################################

if [ ! -f "$REF_GENOME" ]; then
    log_info "Downloading reference genome (this may take a while)..."
    log_info "Source: ${REF_URL}"
    
    wget -O "${REF_GENOME}.gz" "${REF_URL}" || \
        log_error "Failed to download reference genome"
    
    log_info "Decompressing reference genome..."
    gunzip "${REF_GENOME}.gz"
    
    log_success "Reference genome downloaded"
else
    log_success "Reference genome already exists"
fi

# Index reference genome for BWA
if [ ! -f "${REF_GENOME}.bwt" ]; then
    log_info "Indexing reference genome for BWA (this will take time)..."
    bwa index "$REF_GENOME" || log_error "BWA indexing failed"
    log_success "BWA index created"
else
    log_success "BWA index already exists"
fi

# Index reference genome for SAMtools
if [ ! -f "${REF_GENOME}.fai" ]; then
    log_info "Indexing reference genome for SAMtools..."
    samtools faidx "$REF_GENOME" || log_error "SAMtools indexing failed"
    log_success "SAMtools index created"
else
    log_success "SAMtools index already exists"
fi

echo ""

################################################################################
# Quality Control (Optional)
################################################################################

if [ "$RUN_FASTQC" = true ]; then
    log_info "Running FastQC quality control..."
    fastqc -o "${OUTPUT_DIR}/qc" -t $THREADS "$FASTQ_R1" "$FASTQ_R2" || \
        log_warning "FastQC completed with warnings"
    log_success "Quality control complete. Check ${OUTPUT_DIR}/qc/ for reports"
    echo ""
fi

################################################################################
# Read Alignment
################################################################################

BAM_FILE="${OUTPUT_DIR}/alignment/${SAMPLE_ID}.bam"
SORTED_BAM="${OUTPUT_DIR}/alignment/${SAMPLE_ID}.sorted.bam"

if [ ! -f "$SORTED_BAM" ]; then
    log_info "Aligning reads to reference genome (this will take time)..."
    log_info "Using ${THREADS} threads"
    
    bwa mem -t $THREADS -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA" \
        "$REF_GENOME" "$FASTQ_R1" "$FASTQ_R2" | \
        samtools view -@ $THREADS -bS - > "$BAM_FILE" || \
        log_error "Alignment failed"
    
    log_success "Alignment complete"
    
    log_info "Sorting BAM file..."
    samtools sort -@ $THREADS -o "$SORTED_BAM" "$BAM_FILE" || \
        log_error "Sorting failed"
    
    log_info "Indexing BAM file..."
    samtools index "$SORTED_BAM" || log_error "BAM indexing failed"
    
    # Clean up unsorted BAM to save space
    rm "$BAM_FILE"
    
    log_success "BAM file sorted and indexed"
else
    log_success "Sorted BAM file already exists"
fi

echo ""

################################################################################
# Extract APOE Region
################################################################################

APOE_BAM="${OUTPUT_DIR}/alignment/${SAMPLE_ID}.apoe.bam"

if [ ! -f "$APOE_BAM" ]; then
    log_info "Extracting APOE region (${APOE_REGION})..."
    samtools view -@ $THREADS -b "$SORTED_BAM" "$APOE_REGION" > "$APOE_BAM" || \
        log_error "Failed to extract APOE region"
    
    samtools index "$APOE_BAM"
    
    log_success "APOE region extracted"
    
    # Get coverage statistics
    log_info "Calculating coverage statistics..."
    COVERAGE=$(samtools depth "$APOE_BAM" | awk '{sum+=$3} END {print sum/NR}')
    log_info "Average coverage in APOE region: ${COVERAGE}x"
else
    log_success "APOE BAM already exists"
fi

echo ""

################################################################################
# Variant Calling
################################################################################

VCF_FILE="${OUTPUT_DIR}/variants/${SAMPLE_ID}.apoe.vcf.gz"

if [ ! -f "$VCF_FILE" ]; then
    log_info "Calling variants in APOE region..."
    
    bcftools mpileup -Ou -f "$REF_GENOME" -r "$APOE_REGION" "$APOE_BAM" | \
        bcftools call -mv -Oz -o "$VCF_FILE" || \
        log_error "Variant calling failed"
    
    log_info "Indexing VCF file..."
    bcftools index "$VCF_FILE" || log_error "VCF indexing failed"
    
    log_success "Variant calling complete"
else
    log_success "VCF file already exists"
fi

echo ""

################################################################################
# Extract Key APOE SNPs
################################################################################

log_info "Extracting key APOE SNPs (rs429358 and rs7412)..."

# Extract the two critical positions
bcftools view -r "${APOE_CHR}:${RS429358_POS}-${RS429358_POS}" "$VCF_FILE" > \
    "${OUTPUT_DIR}/variants/rs429358.vcf" || log_warning "Could not extract rs429358"

bcftools view -r "${APOE_CHR}:${RS7412_POS}-${RS7412_POS}" "$VCF_FILE" > \
    "${OUTPUT_DIR}/variants/rs7412.vcf" || log_warning "Could not extract rs7412"

log_success "SNPs extracted"

echo ""

################################################################################
# Genotype Interpretation
################################################################################

log_info "Interpreting APOE genotype..."

# Create Python script for genotype interpretation
cat > "${OUTPUT_DIR}/interpret_apoe.py" << 'PYTHON_SCRIPT'
#!/usr/bin/env python3
"""
APOE Genotype Interpretation Script

Determines APOE genotype (ε2, ε3, or ε4) from VCF files containing:
- rs429358 (chr19:45411941, C>T)
- rs7412 (chr19:45412079, C>T)

Genotype Logic:
ε2: rs429358=C/C, rs7412=T/T
ε3: rs429358=C/C, rs7412=C/C
ε4: rs429358=T/T, rs7412=C/C
"""

import sys
import re

def parse_vcf(vcf_file):
    """Extract genotype from VCF file"""
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
                
                # Parse genotype (GT field, typically first in FORMAT)
                format_fields = fields[8].split(':')
                sample_fields = fields[9].split(':')
                
                if 'GT' in format_fields:
                    gt_index = format_fields.index('GT')
                    genotype = sample_fields[gt_index]
                    
                    # Parse genotype (0/0, 0/1, 1/1, etc.)
                    gt_alleles = re.findall(r'\d+', genotype)
                    
                    if len(gt_alleles) >= 2:
                        allele1 = ref if gt_alleles[0] == '0' else alt.split(',')[int(gt_alleles[0])-1]
                        allele2 = ref if gt_alleles[1] == '0' else alt.split(',')[int(gt_alleles[1])-1]
                        
                        return (allele1, allele2), chrom, pos, ref, alt
                        
        return None, None, None, None, None
    except FileNotFoundError:
        return None, None, None, None, None

def determine_apoe_genotype(rs429358_gt, rs7412_gt):
    """
    Determine APOE alleles from genotypes
    
    APOE alleles are defined by two SNPs:
    - rs429358 (position 112)
    - rs7412 (position 158)
    
    Combinations:
    ε2 = Cys112, Cys158 (rs429358=C, rs7412=T)
    ε3 = Cys112, Arg158 (rs429358=C, rs7412=C)
    ε4 = Arg112, Arg158 (rs429358=T, rs7412=C)
    """
    
    def classify_allele(snp429358, snp7412):
        """Classify single APOE allele"""
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
    
    # Get both alleles
    allele1 = classify_allele(rs429358_gt[0], rs7412_gt[0])
    allele2 = classify_allele(rs429358_gt[1], rs7412_gt[1])
    
    # Sort to standardize output (e.g., ε3/ε4 instead of ε4/ε3)
    alleles = sorted([allele1, allele2])
    
    genotype = f"{alleles[0]}/{alleles[1]}"
    
    return allele1, allele2, genotype

def get_clinical_interpretation(genotype):
    """Provide clinical interpretation of APOE genotype"""
    
    interpretations = {
        'ε2/ε2': {
            'risk': 'REDUCED RISK',
            'description': 'Protective against Alzheimer\'s disease',
            'relative_risk': '~0.5x (about half the average risk)'
        },
        'ε2/ε3': {
            'risk': 'REDUCED RISK',
            'description': 'Somewhat protective against Alzheimer\'s disease',
            'relative_risk': '~0.6x (slightly below average risk)'
        },
        'ε2/ε4': {
            'risk': 'AVERAGE TO MODERATE RISK',
            'description': 'Mixed effect; ε2 is protective, ε4 increases risk',
            'relative_risk': '~2-3x (moderate increase)'
        },
        'ε3/ε3': {
            'risk': 'AVERAGE RISK',
            'description': 'Most common genotype; neutral effect',
            'relative_risk': '1x (population average)'
        },
        'ε3/ε4': {
            'risk': 'INCREASED RISK',
            'description': 'Moderately increased risk for Alzheimer\'s disease',
            'relative_risk': '~3x (three times average risk)'
        },
        'ε4/ε4': {
            'risk': 'SIGNIFICANTLY INCREASED RISK',
            'description': 'Substantially increased risk for Alzheimer\'s disease',
            'relative_risk': '~8-12x (eight to twelve times average risk)'
        }
    }
    
    return interpretations.get(genotype, {
        'risk': 'UNKNOWN',
        'description': 'Unable to determine clinical significance',
        'relative_risk': 'Unknown'
    })

def main():
    rs429358_file = sys.argv[1]
    rs7412_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Parse VCF files
    rs429358_gt, chr1, pos1, ref1, alt1 = parse_vcf(rs429358_file)
    rs7412_gt, chr2, pos2, ref2, alt2 = parse_vcf(rs7412_file)
    
    # Determine genotype
    allele1, allele2, genotype = determine_apoe_genotype(rs429358_gt, rs7412_gt)
    
    # Get clinical interpretation
    clinical = get_clinical_interpretation(genotype)
    
    # Generate report
    with open(output_file, 'w') as out:
        out.write("=" * 80 + "\n")
        out.write("APOE GENOTYPING REPORT\n")
        out.write("=" * 80 + "\n\n")
        
        out.write(f"Sample ID: {sys.argv[4] if len(sys.argv) > 4 else 'Unknown'}\n")
        out.write(f"Analysis Date: 2025-10-20\n")
        out.write(f"Reference: GRCh37/hg19\n\n")
        
        out.write("-" * 80 + "\n")
        out.write("RAW GENOTYPE DATA\n")
        out.write("-" * 80 + "\n\n")
        
        if rs429358_gt:
            out.write(f"rs429358 (chr19:45411941):\n")
            out.write(f"  Reference: {ref1}\n")
            out.write(f"  Alternate: {alt1}\n")
            out.write(f"  Genotype: {rs429358_gt[0]}/{rs429358_gt[1]}\n\n")
        else:
            out.write("rs429358: NOT FOUND\n\n")
        
        if rs7412_gt:
            out.write(f"rs7412 (chr19:45412079):\n")
            out.write(f"  Reference: {ref2}\n")
            out.write(f"  Alternate: {alt2}\n")
            out.write(f"  Genotype: {rs7412_gt[0]}/{rs7412_gt[1]}\n\n")
        else:
            out.write("rs7412: NOT FOUND\n\n")
        
        out.write("-" * 80 + "\n")
        out.write("APOE GENOTYPE\n")
        out.write("-" * 80 + "\n\n")
        
        if genotype and 'unknown' not in genotype.lower():
            out.write(f"APOE Genotype: {genotype}\n\n")
            
            out.write("-" * 80 + "\n")
            out.write("CLINICAL INTERPRETATION\n")
            out.write("-" * 80 + "\n\n")
            
            out.write(f"Risk Category: {clinical['risk']}\n")
            out.write(f"Description: {clinical['description']}\n")
            out.write(f"Relative Risk: {clinical['relative_risk']}\n\n")
            
            out.write("Important Notes:\n")
            out.write("- APOE genotype is just one of many risk factors for Alzheimer's disease\n")
            out.write("- Having ε4 alleles does NOT mean you will definitely develop Alzheimer's\n")
            out.write("- Many people with ε4 alleles never develop the disease\n")
            out.write("- Environmental and lifestyle factors also play important roles\n")
            out.write("- This information is for research purposes only\n")
            out.write("- Consult a genetic counselor or physician for clinical interpretation\n\n")
        else:
            out.write(f"Unable to determine APOE genotype\n")
            out.write(f"Reason: {genotype if genotype else 'Missing SNP data'}\n\n")
        
        out.write("=" * 80 + "\n")
        out.write("END OF REPORT\n")
        out.write("=" * 80 + "\n")
    
    # Print to console
    with open(output_file, 'r') as f:
        print(f.read())

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: python3 interpret_apoe.py <rs429358.vcf> <rs7412.vcf> <output.txt> [sample_id]")
        sys.exit(1)
    
    main()
PYTHON_SCRIPT

chmod +x "${OUTPUT_DIR}/interpret_apoe.py"

# Run interpretation
python3 "${OUTPUT_DIR}/interpret_apoe.py" \
    "${OUTPUT_DIR}/variants/rs429358.vcf" \
    "${OUTPUT_DIR}/variants/rs7412.vcf" \
    "${OUTPUT_DIR}/results/${SAMPLE_ID}_apoe_report.txt" \
    "${SAMPLE_ID}"

log_success "APOE genotype analysis complete!"

echo ""
echo "=================================="
log_success "PIPELINE COMPLETED SUCCESSFULLY!"
echo "=================================="
echo ""
log_info "Results saved to: ${OUTPUT_DIR}/results/${SAMPLE_ID}_apoe_report.txt"
log_info "Full VCF file: ${VCF_FILE}"
log_info "APOE BAM file: ${APOE_BAM}"
echo ""

# Display the report
cat "${OUTPUT_DIR}/results/${SAMPLE_ID}_apoe_report.txt"

