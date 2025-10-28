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
