# 📁 Project Structure (October 29, 2025)

The repository now separates **data**, **pipelines**, **results**, and
**documentation** so it is easy to navigate and maintain. This document mirrors
the final layout after the cleanup that followed the 25-patient analysis.

---

## 🗂️ Top-Level Layout

```
D:\APOE\
├── data/
│   ├── external/cgga/         # Original CGGA clinical spreadsheets & extras
│   ├── intermediate/forensic/ # Mapping evidence (headers, sex checks, etc.)
│   ├── processed/             # Final CSV/XLS tables used in analyses
│   ├── raw/                   # Large FASTQs (kept locally, not versioned)
│   └── reference/genome_build/# GRCh37 reference bundle (ignored by git)
│
├── pipelines/
│   ├── apoe_genotyping/       # Targeted APOE extraction helpers
│   ├── cohort_patch_validation/
│   └── patch_validation/
│
├── results/
│   ├── clinical/reclassified_cohort/  # ✅ Final 25-patient analysis package
│   ├── cox/therapy_stratified/        # Cox + KM outputs (all & GBM-only)
│   └── legacy/gbm_who2021_analysis/   # Archived WHO-2016 style results
│
├── docs/                   # Guides and cleanup notes (this file, etc.)
├── reports/                # Publication-ready narrative reports
├── scripts/
│   ├── analysis/           # Entry points to regenerate figures/tables
│   └── utilities/          # Checkers, installers, cleanup helpers
└── README.md               # Repository overview
```

---

## 🔑 Key Locations

- **Final clinical artefacts** → `results/clinical/reclassified_cohort/`
  - Start with `RECLASSIFIED_ANALYSIS/00_READ_ME_FIRST.md`
  - Figures: `apoe_grade_distribution_RECLASSIFIED.png`, `survival_IDH_WT_by_APOE.png`,
    `survival_IDH_WT_by_treatment.png`
  - Tables: `25_samples_RECLASSIFIED.csv`, `survival_cox_IDH_WT_ONLY.xlsx`

- **Therapy-stratified survival modelling** → `results/cox/therapy_stratified/`
  - `gbm_only/Fig_Treatment_by_APOE_Isoform.png`
  - `gbm_only/Fig3_Therapy_Stratified_KM_GBM_Only.png`

- **External source data** → `data/external/cgga/`
  - CGGA spreadsheets (`HRA000071.xlsx`, `HRA000073.xlsx`, `HRA000074.xlsx`)
  - Supplementary methylation files + MD5 manifests

- **Processed tables used by scripts** → `data/processed/`
  - `APOE_GENOTYPES.csv`
  - `FINAL_HRR_CGGA_MAPPING.csv`
  - `FINAL_HRR_CGGA_MAPPING_WITH_CLINICAL.csv`

- **Pipelines** → `pipelines/`
  - `apoe_genotyping/` holds the APOE patch helpers and validation artefacts
  - Validation slices (`cohort_patch_validation/`, `patch_validation/`)

- **Documentation & reports**
  - `docs/CLINICAL_ANALYSIS_GUIDE.md`
  - `docs/CLEANUP_COMPLETE.md`
  - `reports/FINAL_CLINICAL_ANALYSIS_REPORT.md`
  - `reports/COMPREHENSIVE_VALIDATION_REPORT.md`

---

## 🧮 Scripts Overview

- `scripts/analysis/run_clinical_analysis.py`
- `scripts/analysis/rerun_all_analyses_RECLASSIFIED.py`
- `scripts/analysis/survival_cox_IDH_WT_ONLY.py`
- `scripts/analysis/reclassify_IDH_WT_to_grade4.py`
- Utility helpers under `scripts/utilities/` (metadata check, installation, etc.)

All paths in these scripts default to the reorganised `data/` and `results/`
hierarchy.

---

## 🚀 Quick Navigation Cheatsheet

- **Presentations** → `results/clinical/reclassified_cohort/RECLASSIFIED_ANALYSIS/`
- **Cox modelling** → `results/cox/therapy_stratified/`
- **Legacy comparisons** → `results/legacy/gbm_who2021_analysis/`
- **Source data** → `data/external/cgga/`
- **Processed tables** → `data/processed/`
- **Genotyping validation** → `pipelines/apoe_genotyping/`
- **Reports** → `reports/`

---

## ✅ Status

- Cohort size: 25 patients (20 IDH-wildtype reclassified to WHO IV)
- ε4 carriers: 3 (all alive; `p = 0.0086` vs non-carriers)
- Figures and reports ready for Dr. Guo
- Repository cleaned and aligned with the current analysis workflow

For more context see the top-level `README.md` or `docs/CLINICAL_ANALYSIS_GUIDE.md`.

