# APOE Genotyping Pipeline

This pipeline analyzes whole genome sequencing (WGS) data to determine APOE genotype, which is associated with Alzheimer's disease risk.

## Overview

The pipeline processes FASTQ files from your CGGA WGS data to identify the two key SNPs that define APOE alleles:
- **rs429358** (chr19:45411941, C>T)
- **rs7412** (chr19:45412079, C>T)

These SNPs determine whether you have ε2, ε3, or ε4 alleles.

## Targeted APOE patches (exon slices)

For rapid demonstrations or validation runs against existing BAM/CRAM files you
can use the new `apoe_analysis/apoe_patches.py` helper.  It wraps the PI
"exon 2" patch and the scientific validation track (exon 4 + rs429358/rs7412)
into reproducible commands:

```bash
python -m apoe_analysis.apoe_patches patch-a \
  --bam /path/to/sample.bam \
  --reference /path/to/GRCh38.fa \
  --gtf /path/to/genes.gtf \
  --output-dir patch_outputs

python -m apoe_analysis.apoe_patches patch-b \
  --bam /path/to/sample.bam \
  --reference /path/to/GRCh38.fa \
  --gtf /path/to/genes.gtf \
  --output-dir patch_outputs \
  --build hg38
```

* `patch-a` generates exon 2 read slices (`*_APOE_exon2.bam`), reference
  sequences, and depth summaries exactly as requested by the PI.
* `patch-b` adds exon 4 extraction plus targeted calls for rs429358/rs7412,
  producing BAM, FASTA, compressed VCF, and a ready-to-share TSV report.

Supply multiple `--bam` arguments (optionally paired with `--label`) to process
tumour/normal or cohort BAMs together.  Outputs are organised per-sample and a
JSON manifest summarises the generated artefacts.

## Your Data

Located in `E:\`:
- `HRR024685_f1.fq.gz` - Forward reads
- `HRR024685_r2.fq.gz` - Reverse reads

## Quick Start

### Step 1: Install Required Tools

See `INSTALLATION_GUIDE.md` for detailed instructions.

**Recommended (works on Windows):**
```bash
# Install Miniconda, then:
conda create -n apoe_analysis python=3.9
conda activate apoe_analysis
conda install -c bioconda bwa samtools bcftools fastqc
```

### Step 2: Run the Pipeline

**On Windows with WSL or Git Bash:**
```bash
cd /mnt/e/  # or E:/ in Git Bash
bash apoe_pipeline.sh
```

**On Windows with Conda:**
```bash
conda activate apoe_analysis
cd /d E:\
bash apoe_pipeline.sh
```

### Step 3: View Results

Results will be saved in `E:\apoe_analysis\results\HRR024685_apoe_report.txt`

## What the Pipeline Does

1. **Download Reference Genome** - GRCh37/hg19 human reference (~3 GB)
2. **Quality Control** - FastQC analysis of your FASTQ files
3. **Read Alignment** - BWA aligns reads to reference genome
4. **Extract APOE Region** - Focus on chromosome 19:45409039-45412650
5. **Variant Calling** - BCFtools identifies variants
6. **Genotype Determination** - Identifies your APOE alleles (ε2/ε3/ε4)
7. **Clinical Interpretation** - Provides Alzheimer's disease risk assessment

## APOE Genotypes and Risk

| Genotype | Alzheimer's Risk | Frequency |
|----------|------------------|-----------|
| ε2/ε2 | Reduced (~0.5x) | Rare (<1%) |
| ε2/ε3 | Reduced (~0.6x) | Uncommon (~10%) |
| ε2/ε4 | Average to Moderate (~2-3x) | Rare (~2%) |
| ε3/ε3 | Average (1x) | Most common (~60%) |
| ε3/ε4 | Increased (~3x) | Common (~20%) |
| ε4/ε4 | Significantly Increased (~8-12x) | Rare (~2%) |

## Important Notes

- APOE genotype is just ONE of many risk factors for Alzheimer's disease
- Having ε4 alleles does NOT guarantee you will develop Alzheimer's
- Many people with ε4/ε4 never develop the disease
- This is for research purposes only
- Consult a genetic counselor for clinical interpretation

## System Requirements

- **OS**: Windows 10/11 (with WSL or Conda), Linux, or macOS
- **RAM**: 8 GB minimum (16 GB recommended)
- **CPU**: 4+ cores recommended
- **Disk Space**: 50-100 GB free on E:\ drive
- **Runtime**: 30 minutes to 2 hours

## File Structure

```
E:\
├── HRR024685_f1.fq.gz              # Your input FASTQ files
├── HRR024685_r2.fq.gz
├── apoe_pipeline.sh                # Main pipeline script
├── reference/                       # Reference genome (auto-downloaded)
│   └── human_g1k_v37.fasta
└── apoe_analysis/                   # Output directory
    ├── qc/                          # FastQC reports
    ├── alignment/                   # BAM files
    ├── variants/                    # VCF files
    └── results/                     # Final APOE genotype report
```

## Troubleshooting

See `INSTALLATION_GUIDE.md` for common issues and solutions.

**Common Issues:**
- **Out of memory**: Reduce threads in script (change `THREADS=4` to `THREADS=2`)
- **Command not found**: Make sure conda environment is activated
- **Path errors on Windows**: Run from WSL where E:\ is `/mnt/e/`
- **Download fails**: Manually download reference genome (see installation guide)

## Technical Details

- **Reference**: GRCh37/hg19 (1000 Genomes Phase 2)
- **Aligner**: BWA-MEM
- **Variant Caller**: BCFtools mpileup + call
- **APOE Region**: chr19:45409039-45412650 (3,611 bp)
- **Gene**: APOE gene, exon 4

## Support

For issues or questions:
1. Check `INSTALLATION_GUIDE.md`
2. Ensure all tools are properly installed: `bwa`, `samtools`, `bcftools`, `python3`
3. Verify input files exist in E:\
4. Check you have sufficient disk space

## License

This pipeline is for research and educational purposes only.

## References

1. Genin et al. (2011) "APOE and Alzheimer's disease: a major gene with semi-dominant inheritance"
2. Farrer et al. (1997) "Effects of age, sex, and ethnicity on the association between apolipoprotein E genotype and Alzheimer disease"
3. Liu et al. (2013) "Apolipoprotein E and Alzheimer disease: risk, mechanisms and therapy"

