# APOE Analysis - Tool Installation Guide

## Required Bioinformatics Tools

### 1. **BWA** (Burrows-Wheeler Aligner)
Aligns sequencing reads to the reference genome.

**Installation:**
```bash
# Ubuntu/Debian
sudo apt-get install bwa

# macOS
brew install bwa

# Conda (recommended - works on all platforms including Windows)
conda install -c bioconda bwa
```

### 2. **SAMtools**
Manipulates BAM/SAM alignment files.

**Installation:**
```bash
# Ubuntu/Debian
sudo apt-get install samtools

# macOS
brew install samtools

# Conda
conda install -c bioconda samtools
```

### 3. **BCFtools**
Variant calling and VCF manipulation.

**Installation:**
```bash
# Ubuntu/Debian
sudo apt-get install bcftools

# macOS
brew install bcftools

# Conda
conda install -c bioconda bcftools
```

### 4. **FastQC** (Optional but recommended)
Quality control of sequencing data.

**Installation:**
```bash
# Ubuntu/Debian
sudo apt-get install fastqc

# macOS
brew install fastqc

# Conda
conda install -c bioconda fastqc
```

### 5. **Python 3**
For genotype interpretation script.

Usually pre-installed on Linux/macOS. Check with:
```bash
python3 --version
```

---

## Quick Installation with Conda (Recommended for Windows)

If you have Conda/Miniconda installed, you can install everything at once:

```bash
# Create a new environment
conda create -n apoe_analysis python=3.9

# Activate environment
conda activate apoe_analysis

# Install all tools
conda install -c bioconda bwa samtools bcftools fastqc

# Verify installations
bwa
samtools --version
bcftools --version
fastqc --version
```

---

## Installing Conda/Miniconda on Windows (if needed)

**Windows:**
1. Download Miniconda from: https://docs.conda.io/en/latest/miniconda.html
2. Run the installer: `Miniconda3-latest-Windows-x86_64.exe`
3. During installation, check "Add Miniconda to PATH"
4. Open a new Command Prompt or PowerShell and test: `conda --version`

**Alternative: WSL (Windows Subsystem for Linux):**
```powershell
# In PowerShell as Administrator
wsl --install

# After restart, in WSL terminal:
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

---

## Disk Space Requirements

- **Reference genome**: ~3 GB
- **Indexed reference**: ~5 GB total
- **BAM files**: Depends on coverage, typically 5-50 GB for WGS
- **Total recommended**: At least 50-100 GB free space on E:\ drive

---

## Running the Pipeline

Once tools are installed:

### On WSL or Linux:
```bash
cd /mnt/e/
chmod +x apoe_pipeline.sh
./apoe_pipeline.sh
```

### On Windows with Conda:
```bash
# Open Anaconda Prompt or Miniconda Prompt
conda activate apoe_analysis
cd /d E:\
bash apoe_pipeline.sh
```

The script will:
- Download the reference genome automatically
- Process your FASTQ files
- Call variants
- Determine APOE genotype

**Estimated runtime**: 30 minutes to 2 hours (depending on your computer and coverage depth)

---

## Troubleshooting

### "Command not found" errors
- Make sure tools are installed and in your PATH
- If using Conda, ensure the environment is activated: `conda activate apoe_analysis`

### Reference genome download fails
- Download manually from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/
- Place in `E:/reference/` directory

### Out of memory errors
- Reduce the number of threads in the script (change `THREADS=4` to `THREADS=2`)
- Ensure you have at least 8 GB RAM available

### Permission errors on Windows
- Run the script in WSL (Windows Subsystem for Linux)
- Or use Git Bash with admin privileges
- Or use Anaconda Prompt with admin privileges

### Path issues on Windows
The script uses Unix-style paths. If you get path errors:
- Run from WSL where E:\ is mounted at `/mnt/e/`
- Or edit the script to use Windows paths (E:\\ instead of E:/)

---

## System Requirements

- **OS**: Linux, macOS, or Windows 10/11 with WSL or Conda
- **RAM**: Minimum 8 GB (16 GB recommended)
- **CPU**: Multi-core processor (4+ cores recommended)
- **Disk**: 50-100 GB free space
- **Internet**: For downloading reference genome (~3 GB)

