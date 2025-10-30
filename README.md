# APOE Clinical Analysis Repository

End-to-end analysis of the 25 CGGA patients with validated APOE genotypes. The
repository now only tracks the clinical cohort deliverable (reclassified to
WHO-2021 rules) plus the supporting genotyping pipeline that produced the
APOE calls.

## Highlights

- ✅ **25 patients with confirmed APOE genotypes** (includes HRR024685 ε4/ε4)
- ✅ **IDH-wildtype reclassified to WHO Grade IV** for survival modelling
- ✅ **Kaplan–Meier and Cox regression** focussed on ε4 carrier status
- ✅ **Treatment-stratified figures and tables** ready for presentation
- ✅ **All raw/external data, scripts, results, and reports organised by role**

Key clinical result: **ε4 carriers show 100% survival in IDH-wildtype GBM**
(p = 0.0086), with the strongest effect in the RT+TMZ group (p = 0.0401).

## Repository Layout

```
data/
  external/cgga/          # Original CGGA clinical spreadsheets and manifests
  intermediate/forensic/  # Mapping evidence (BAM headers, sex checks, etc.)
  processed/              # Final tables (genotypes, ID mapping, metadata)
  raw/                    # Large FASTQ inputs (not tracked by git)
  reference/genome_build/ # GRCh37 reference bundle (ignored by git)

pipelines/
  apoe_genotyping/        # Targeted APOE extraction helpers
  cohort_patch_validation/
  patch_validation/

results/
  clinical/reclassified_cohort/   # Definitive 25-patient analysis package
  cox/therapy_stratified/         # Cox + KM outputs (all + GBM-only)
  legacy/gbm_who2021_analysis/    # Deprecated WHO-2016 era runs

docs/      # How-to guides and cleanup notes
reports/   # Narrative deliverables ready to share
scripts/
  analysis/   # End-to-end analysis entry points
  utilities/  # Helper scripts (metadata checks, install helpers)

pipelines/.../alignment, results, variants  # Supporting genotyping artefacts
```

## Data Inventory

- `data/external/cgga/` – Raw CGGA clinical spreadsheets, methylation files,
  and manifest MD5 checks.
- `data/processed/` – Harmonised tables used by the analysis scripts:
  `APOE_GENOTYPES.csv`, `FINAL_HRR_CGGA_MAPPING.csv`,
  `FINAL_HRR_CGGA_MAPPING_WITH_CLINICAL.csv`.
- `data/intermediate/forensic/` – Evidence produced while verifying HRR ↔ CGGA
  identities (BAM header matches, sex checks, preliminary matches).
- `data/raw/` – Sequencing FASTQs for HRR024685 (kept outside git, but stored
  locally for provenance).

## Running the Clinical Analysis

All orchestrator scripts live in `scripts/analysis/` and expect Python 3.9+.

1. Create / activate environment and install requirements (see
   `scripts/utilities/install_dependencies.sh` for exact commands).
2. Generate the merged dataset (optional – already in `data/processed/`):
   ```bash
   python scripts/analysis/reclassify_IDH_WT_to_grade4.py \
     --input data/processed/FINAL_HRR_CGGA_MAPPING_WITH_CLINICAL.csv \
     --output data/processed/25_samples_RECLASSIFIED.csv
   ```
3. Run the full analytical suite (recreates everything in `results/clinical/`):
   ```bash
   python scripts/analysis/rerun_all_analyses_RECLASSIFIED.py \
     --clinical-table data/processed/25_samples_RECLASSIFIED.csv \
     --output-dir results/clinical/reclassified_cohort
   ```
4. Therapy and Cox modelling (reproduces `results/cox/therapy_stratified/`):
   ```bash
   python scripts/analysis/survival_cox_IDH_WT_ONLY.py \
     --clinical-table data/processed/25_samples_RECLASSIFIED.csv \
     --output-dir results/cox/therapy_stratified
   ```

All figures, tables, and markdown summaries land in the corresponding results
directories and are ready for distribution.

## Reproducing Genotyping Artifacts

The `pipelines/apoe_genotyping/` module contains the targeted APOE patching
helpers used to validate HRR024685:

```bash
python -m pipelines.apoe_genotyping.apoe_patches patch-a \
  --bam pipelines/apoe_genotyping/alignment/HRR024685.sorted.bam \
  --reference data/reference/genome_build/human_g1k_v37.fasta \
  --gtf data/reference/genome_build/apoe_grch37_NCBI_CORRECT.gtf \
  --output-dir pipelines/apoe_genotyping/results
```

The validation bundles (`pipelines/cohort_patch_validation/` and
`pipelines/patch_validation/`) capture the exon-slice outputs and manifests.

## Key Outputs for Presentations

- `results/clinical/reclassified_cohort/RECLASSIFIED_ANALYSIS/`
  - `apoe_grade_distribution_RECLASSIFIED.png`
  - `survival_IDH_WT_by_APOE.png`
  - `survival_IDH_WT_by_treatment.png`
- `results/cox/therapy_stratified/gbm_only/`
  - `Fig_Treatment_by_APOE_Isoform.png`
  - `Fig3_Therapy_Stratified_KM_GBM_Only.png`
- `reports/FINAL_CLINICAL_ANALYSIS_REPORT.md`
  – Narratives, tables, and interpretation.

## Notes

- Raw sequencing inputs and reference genomes remain outside git history and
  live under `data/raw/` and `data/reference/` respectively.
- Pycache artefacts are ignored via `.gitignore`.
- Legacy WHO-2016 analyses are preserved under `results/legacy/` for traceability
  but are excluded from new deliverables.

For questions or further adjustments, start with `docs/CLINICAL_ANALYSIS_GUIDE.md`.

