# APOE Analysis - Real-Time Status

## Current Status
**Pipeline:** ULTRA-FAST MODE (12 threads per sample)  
**Started:** Sunday, October 26, 2025 - 21:27:53 EDT  
**Target:** 24 samples (HRR024686-HRR024710)

## How to Check Progress

###  **From Windows:**
```powershell
# Check completed samples
wsl bash -c "ls /mnt/d/APOE/ultra_fast/results/*.txt 2>/dev/null | wc -l"

# See latest results
wsl bash -c "ls -t /mnt/d/APOE/ultra_fast/results/*.txt | head -3 | xargs -I {} bash -c 'echo; cat {}'"

# Check if pipeline is running
wsl bash -c "ps aux | grep ultra_fast_pipeline | grep -v grep"
```

### **From Ubuntu (WSL):**
```bash
cd /mnt/d/APOE

# Quick status
echo "Completed: $(ls ultra_fast/results/*.txt 2>/dev/null | wc -l) / 24"

# See all results  
for f in ultra_fast/results/*.txt; do
    sample=$(basename $f _result.txt)
    genotype=$(grep "APOE Genotype:" $f | cut -d: -f2)
    echo "$sample:$genotype"
done

# Watch pipeline live
tail -f ultra_fast.log
```

## Expected Timeline

- **Per sample:** ~15-20 minutes
- **Total time:** ~6-8 hours for all 24 samples
- **Completion:** ~3:00-5:00 AM EDT

## What's Happening

Each sample goes through 6 steps:
1. ⚡ **BWA Alignment** (10-15 min) - Aligns reads to reference genome
2. 🔄 **SAM to BAM** (2-3 min) - Converts format
3. 📊 **Sorting** (3-5 min) - Sorts alignments
4. 📇 **Indexing** (1 min) - Creates index
5. 🎯 **Extract APOE** (<1 min) - Gets chr19 region
6. 🧬 **Genotyping** (<1 min) - Calls variants

## Results Location

- **Individual reports:** `D:\APOE\ultra_fast\results\`
- **Logs:** `D:\APOE\ultra_fast\logs\`
- **Main log:** `D:\APOE\ultra_fast.log`

## If Something Goes Wrong

1. Check if processes are running:
   ```powershell
   wsl bash -c "ps aux | grep bwa"
   ```

2. Check last error:
   ```powershell
   wsl bash -c "grep -i error /mnt/d/APOE/ultra_fast.log | tail -5"
   ```

3. Restart from where it stopped (it will skip completed samples):
   ```powershell
   wsl bash -c "cd /mnt/d/APOE && bash ultra_fast_pipeline.sh 2>&1 | tee -a ultra_fast.log"
   ```


