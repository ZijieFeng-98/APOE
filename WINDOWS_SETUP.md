# Windows Setup Guide for APOE Pipeline

Since you're running Windows, you have three main options to run this bioinformatics pipeline:

## Option 1: Windows Subsystem for Linux (WSL) - RECOMMENDED

WSL provides a full Linux environment within Windows. This is the most compatible option.

### Install WSL:
```powershell
# Open PowerShell as Administrator
wsl --install

# Restart your computer when prompted
```

### After restart, set up the tools:
```bash
# Open Ubuntu from Start menu
# Update system
sudo apt update && sudo apt upgrade -y

# Install bioinformatics tools
sudo apt install -y bwa samtools bcftools fastqc wget python3

# Verify installations
bwa
samtools --version
bcftools --version
fastqc --version
```

### Run the pipeline:
```bash
# Your E:\ drive is accessible at /mnt/e/ in WSL
cd /mnt/e/

# Copy the pipeline script from D:\APOE to E:\
cp /mnt/d/APOE/apoe_pipeline.sh /mnt/e/

# Make it executable
chmod +x apoe_pipeline.sh

# Run it
./apoe_pipeline.sh
```

---

## Option 2: Conda/Miniconda (Works on Windows directly)

Conda works natively on Windows and can install most bioinformatics tools.

### Install Miniconda:
1. Download from: https://docs.conda.io/en/latest/miniconda.html
2. Run: `Miniconda3-latest-Windows-x86_64.exe`
3. During installation, check "Add Miniconda to PATH"
4. Open a new Command Prompt or PowerShell

### Install tools:
```bash
# Create environment
conda create -n apoe_analysis python=3.9 -y

# Activate environment
conda activate apoe_analysis

# Install bioinformatics tools
conda install -c bioconda bwa samtools bcftools fastqc -y

# Verify
bwa
samtools --version
bcftools --version
```

### Run the pipeline:
```bash
# Copy script to E:\ first (use File Explorer)
# Then:
cd E:\
bash apoe_pipeline.sh
```

**Note**: You may need to install Git for Windows to get the `bash` command:
https://git-scm.com/download/win

---

## Option 3: Docker (Advanced)

Docker provides containerized environments with all tools pre-installed.

### Install Docker Desktop for Windows:
1. Download from: https://www.docker.com/products/docker-desktop
2. Install and restart
3. Enable WSL 2 backend when prompted

### Use a biocontainers image:
```bash
# Pull a bioinformatics container
docker pull quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40

# Run the pipeline in the container
docker run -v E:/:/data -it quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40 bash

# Inside the container:
cd /data
bash apoe_pipeline.sh
```

---

## Comparison

| Option | Ease of Setup | Compatibility | Performance | Best For |
|--------|---------------|---------------|-------------|----------|
| **WSL** | Medium | Excellent | Excellent | Best overall |
| **Conda** | Easy | Good | Good | Quick start |
| **Docker** | Medium | Excellent | Excellent | Reproducibility |

---

## Recommended: WSL

For this pipeline, **WSL is recommended** because:
- Full Linux compatibility
- Best performance
- All tools work perfectly
- Easy access to Windows files
- Free and officially supported by Microsoft

---

## After Setup

Once you've chosen and set up your environment:

1. **Verify all tools are installed**:
```bash
bwa
samtools --version
bcftools --version
python3 --version
```

2. **Copy the pipeline script to E:\**:
```bash
# Copy from D:\APOE to E:\
cp D:\APOE\apoe_pipeline.sh E:\
# OR in WSL:
cp /mnt/d/APOE/apoe_pipeline.sh /mnt/e/
```

3. **Run the pipeline**:
```bash
cd E:\  # or /mnt/e/ in WSL
bash apoe_pipeline.sh
```

---

## Disk Space Reminder

Make sure E:\ has at least **50-100 GB free** for:
- Reference genome: ~5 GB
- BAM files: ~10-50 GB (depends on coverage)
- VCF files: ~1 GB
- Temporary files: ~10 GB

Check your space:
```bash
# Windows
dir E:\

# WSL
df -h /mnt/e/
```

---

## Getting Help

If you run into issues:
1. Check you have enough disk space
2. Make sure tools are in your PATH
3. Try running with WSL if Conda isn't working
4. Check file paths use correct format for your environment

---

## Quick Test

Test if your environment is working:

```bash
# Activate your environment (if using conda)
conda activate apoe_analysis

# Test each tool
bwa 2>&1 | head -5
samtools --version
bcftools --version
python3 --version

# Check if FASTQ files exist
ls -lh /mnt/e/HRR024685*.fq.gz  # WSL
dir E:\HRR024685*.fq.gz         # Windows CMD
```

All tools should respond without "command not found" errors.

