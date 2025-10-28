# Quick Start Guide - APOE Analysis

## Your Files

**Location**: `E:\`
- `HRR024685_f1.fq.gz` - Forward reads (FASTQ format)
- `HRR024685_r2.fq.gz` - Reverse reads (FASTQ format)
- `HRR024685.sra` - Original SRA file
- `HRR024685.xml` - Sample metadata

**Type**: Paired-end whole genome sequencing (WGS) reads from CGGA/DDBJ

---

## Fastest Way to Get Started (Windows)

### 1. Install WSL (5 minutes)
```powershell
# Open PowerShell as Administrator
wsl --install
# Restart computer
```

### 2. Install Tools in WSL (5 minutes)
```bash
# Open Ubuntu from Start menu
sudo apt update
sudo apt install -y bwa samtools bcftools python3 wget
```

### 3. Copy and Run Pipeline (30-120 minutes)
```bash
# Copy pipeline to E:\
cp /mnt/d/APOE/apoe_pipeline.sh /mnt/e/

# Navigate to E:\
cd /mnt/e/

# Make executable
chmod +x apoe_pipeline.sh

# Run it!
./apoe_pipeline.sh
```

### 4. View Results
```bash
# Results will be in:
cat /mnt/e/apoe_analysis/results/HRR024685_apoe_report.txt
```

---

## What You'll Get

The pipeline will tell you:
- Your **APOE genotype**: ε2/ε2, ε2/ε3, ε2/ε4, ε3/ε3, ε3/ε4, or ε4/ε4
- Your **Alzheimer's disease risk** relative to average
- The specific SNP variants at rs429358 and rs7412

---

## Expected Timeline

| Step | Time | Disk Usage |
|------|------|------------|
| Setup WSL + tools | 10 min | ~500 MB |
| Download reference genome | 10-30 min | ~3 GB |
| Align reads (BWA) | 20-90 min | ~20-40 GB |
| Variant calling | 1-5 min | ~500 MB |
| Generate report | <1 min | ~10 KB |
| **TOTAL** | **30-120 min** | **~50 GB** |

*Time varies based on your CPU and internet speed*

---

## System Requirements

✅ **Required**:
- Windows 10/11 (version 2004 or higher for WSL)
- 8 GB RAM (16 GB recommended)
- 50-100 GB free space on E:\
- Internet connection (for reference genome download)

✅ **Recommended**:
- 4+ CPU cores
- SSD for E:\ drive (faster than HDD)

---

## Troubleshooting

### "wsl: command not found"
- You need Windows 10 version 2004+ or Windows 11
- Or manually install WSL: https://docs.microsoft.com/en-us/windows/wsl/install-manual

### "No space left on device"
- Check E:\ has 50+ GB free: `df -h /mnt/e/`
- Delete unnecessary files to free up space

### "Permission denied"
- Make script executable: `chmod +x apoe_pipeline.sh`
- Or run with: `bash apoe_pipeline.sh`

### Pipeline fails during alignment
- Reduce threads: Edit script, change `THREADS=4` to `THREADS=2`
- Ensure you have 8+ GB RAM available

### Reference genome download fails
- Download manually from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
- Place in: `E:\reference\human_g1k_v37.fasta`

---

## Alternative: Use Conda (No WSL)

If you prefer not to use WSL:

```bash
# Install Miniconda from: https://docs.conda.io/en/latest/miniconda.html

# Then:
conda create -n apoe python=3.9
conda activate apoe
conda install -c bioconda bwa samtools bcftools

# Copy script to E:\ and run
cd E:\
bash apoe_pipeline.sh
```

---

## Need Help?

1. Read `WINDOWS_SETUP.md` for detailed Windows instructions
2. Read `INSTALLATION_GUIDE.md` for tool installation details
3. Read `README.md` for complete pipeline documentation

---

## Understanding Your Results

Your report will show something like:

```
APOE Genotype: ε3/ε4

Risk Category: INCREASED RISK
Description: Moderately increased risk for Alzheimer's disease
Relative Risk: ~3x (three times average risk)
```

**Remember**:
- This is ONE of many risk factors
- Having ε4 doesn't mean you'll definitely get Alzheimer's
- This is for research purposes only
- Consult a genetic counselor for medical advice

---

## Files You'll Get

After running the pipeline:

```
E:\apoe_analysis\
├── results/
│   └── HRR024685_apoe_report.txt    ← YOUR MAIN RESULT
├── variants/
│   ├── HRR024685.apoe.vcf.gz        ← All variants in APOE region
│   ├── rs429358.vcf                  ← Key SNP 1
│   └── rs7412.vcf                    ← Key SNP 2
├── alignment/
│   ├── HRR024685.sorted.bam          ← Full genome alignment
│   └── HRR024685.apoe.bam            ← APOE region only
└── qc/
    └── *.html                        ← Quality control reports
```

---

## Ready to Start?

```bash
# 1. Open PowerShell as Admin
wsl --install

# 2. Restart computer

# 3. Open Ubuntu from Start menu

# 4. Install tools
sudo apt update && sudo apt install -y bwa samtools bcftools python3 wget

# 5. Run pipeline
cp /mnt/d/APOE/apoe_pipeline.sh /mnt/e/
cd /mnt/e/
chmod +x apoe_pipeline.sh
./apoe_pipeline.sh

# 6. Wait 30-120 minutes

# 7. View your APOE genotype!
cat apoe_analysis/results/HRR024685_apoe_report.txt
```

**Good luck with your analysis!** 🧬

