# BATCH APOE ANALYSIS - Live Status

## Cohort Information

**Study:** CGGA WESeq  
**Total Patients:** 26  
**Sample IDs:** HRR024685 - HRR024710  
**Analysis Started:** 2025-10-26  

---

## Pipeline Workflow

### Stage 1: Data Processing ⏳
For each patient (26 samples):
1. Read alignment with BWA-MEM (~20-30 min per sample)
2. BAM sorting and indexing (~5-10 min per sample)
3. APOE region extraction (~1-2 min per sample)
4. Coverage calculation (~1 min per sample)

### Stage 2: Genotyping ⏳
For each patient:
1. Extract reads at rs429358 (chr19:45411941)
2. Extract reads at rs7412 (chr19:45412079)
3. Determine APOE genotype (ε2/ε3/ε4)
4. Calculate coverage at key positions

### Stage 3: Reporting ⏳
1. Generate individual patient reports
2. Create cohort summary with statistics
3. Analyze genotype distribution
4. Generate comprehensive report

---

## Estimated Timeline

| Task | Time per Sample | Total for 26 Samples |
|------|----------------|----------------------|
| Alignment | 20-30 min | 8-13 hours |
| Sorting | 5-10 min | 2-4 hours |
| APOE extraction | 1-2 min | 26-52 min |
| Genotyping | 1 min | 26 min |
| **TOTAL** | **~30-45 min** | **~11-18 hours** |

**Note:** Pipeline runs sequentially (one patient at a time)

---

## Real-Time Monitoring

### Check Progress:
```bash
# View live log
tail -f /mnt/d/APOE/batch_analysis.log

# Check which sample is processing
grep "Processing:" /mnt/d/APOE/batch_analysis.log | tail -1

# Check how many completed
ls /mnt/d/APOE/batch_analysis/results/*.txt | wc -l

# View latest sample result
ls -t /mnt/d/APOE/batch_analysis/results/*.txt | head -1 | xargs cat
```

### Check if running:
```bash
ps aux | grep batch_apoe_pipeline
```

---

## Output Files

### During Analysis:
```
/mnt/d/APOE/batch_analysis/
├── logs/
│   ├── HRR024685.log
│   ├── HRR024686.log
│   └── ... (one per sample)
├── results/
│   ├── HRR024685_summary.txt
│   ├── HRR024686_summary.txt
│   └── ... (one per sample)
└── HRR024685/
    └── alignment/
        ├── HRR024685.sorted.bam
        └── HRR024685.apoe.bam
```

### Final Report:
```
/mnt/d/APOE/batch_analysis/summary/
└── BATCH_SUMMARY_YYYYMMDD_HHMMSS.txt
```

---

## Expected Results

For each patient, you'll get:

### Individual Summary:
```
Sample: HRR024685
Status: SUCCESS

Coverage:
- APOE region: XX.Xx
- rs429358: XXx
- rs7412: XXx

Genotypes:
- rs429358: T/T (or C/T, or C/C)
- rs7412: C/C (or C/T, or T/T)

APOE Genotype: ε?/ε?
Risk Category: [REDUCED/AVERAGE/INCREASED]
Relative AD Risk: [X]x
```

### Cohort Summary:
```
Total: 26 patients
Genotype Distribution:
- ε2/ε2: X patients (X%)
- ε2/ε3: X patients (X%)
- ε2/ε4: X patients (X%)
- ε3/ε3: X patients (X%)
- ε3/ε4: X patients (X%)
- ε4/ε4: X patients (X%)
```

---

## Status Updates

**Current Status:** RUNNING ⏳

Progress will be updated here as samples complete...

---

## Quick Commands

### Stop the analysis (if needed):
```bash
pkill -f batch_apoe_pipeline
```

### Resume analysis:
```bash
cd /mnt/d/APOE
bash batch_apoe_pipeline.sh
```

### View summary so far:
```bash
cat /mnt/d/APOE/batch_analysis/summary/BATCH_SUMMARY_*.txt
```

---

**Last Updated:** Analysis started  
**Next Update:** When first sample completes (~30-45 minutes)


