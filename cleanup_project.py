#!/usr/bin/env python3
"""
Comprehensive project cleanup - remove temporary/old files, keep essential analyses
"""

import os
import shutil
from pathlib import Path

# Files and directories to KEEP (essential)
KEEP_ITEMS = {
    # Core data
    'DATA/',
    'reference/',
    'APOE_GENOTYPES.csv',
    'FINAL_HRR_CGGA_MAPPING.csv',
    'FINAL_HRR_CGGA_MAPPING_WITH_CLINICAL.csv',
    
    # Sequencing data (BAM files)
    'apoe_analysis/',
    
    # Final correct analyses
    'clinical_analysis_outputs/RECLASSIFIED_ANALYSIS/',
    'clinical_analysis_outputs/25_samples_metadata.xlsx',
    'clinical_analysis_outputs/25_samples_RECLASSIFIED_metadata.xlsx',
    'clinical_analysis_outputs/00_NAVIGATION_GUIDE.md',
    
    # Cox and GBM analyses
    'cox_therapy_analysis_outputs/',
    'gbm_who2021_analysis/',
    
    # Important validation reports
    'COHORT_PATCH_A_VALIDATION_REPORT.md',
    'COMPREHENSIVE_VALIDATION_REPORT.md',
    'REAL_EXON4_SEQUENCES_ALL_PATIENTS.md',
    'DATA_INTEGRITY_AUDIT_REPORT.md',
    'FORENSIC_MAPPING_REPORT.md',
    'MAPPING_SUCCESS_REPORT.md',
    
    # Key documentation
    'README.md',
    'INSTALLATION_GUIDE.md',
    'QUICK_START.md',
    
    # Reclassification files
    'reclassify_IDH_WT_to_grade4.py',
    'rerun_all_analyses_RECLASSIFIED.py',
    'survival_cox_IDH_WT_ONLY.py',
    'check_25_samples_metadata.py',
}

# Directories to REMOVE (temporary/old/superseded)
REMOVE_DIRS = [
    'fast_analysis',
    'ultra_fast',
    'working_analysis',
    'reanalysis',
    'final_analysis',
    'batch_analysis',
    'gbm_analysis_outputs',  # Superseded by gbm_who2021_analysis
]

# Files to REMOVE (temporary/old scripts and logs)
REMOVE_FILES = [
    # Old pipeline scripts
    'apoe_pipeline.sh',
    'batch_apoe_pipeline.sh',
    'bulletproof_pipeline.sh',
    'complete_pipeline.sh',
    'comprehensive_genotype_analysis.sh',
    'fast_pipeline.sh',
    'memory_efficient_pipeline.sh',
    'ultra_fast_pipeline.sh',
    'working_pipeline.sh',
    'setup_and_run.sh',
    'install_and_run.sh',
    
    # Monitoring scripts
    'auto_monitor.sh',
    'auto_monitor_and_report.sh',
    'check_progress.sh',
    'continuous_monitor.sh',
    'hourly_monitor.sh',
    'live_monitor.sh',
    'monitor_batch.sh',
    'monitor_fast.sh',
    'monitor_progress.sh',
    'simple_monitor.sh',
    'watchdog.sh',
    
    # Old logs
    'apoe_analysis.log',
    'batch_analysis.log',
    'bulletproof.log',
    'completion.log',
    'fast_analysis.log',
    'pipeline_output.log',
    'reanalysis.log',
    'ultra_fast.log',
    'working.log',
    
    # Old status/progress files
    'AUTO_STATUS.txt',
    'BATCH_ANALYSIS_STATUS.md',
    'BATCH_PROGRESS.md',
    'CHECK_STATUS.md',
    'COMMANDS.txt',
    'COMPLETION_LOG.txt',
    'LEAVE_SUMMARY.txt',
    'MONITORING_PLAN.md',
    'MONITORING_STATUS.md',
    'PROGRESS_REPORT.md',
    'PROGRESS_UPDATES.txt',
    'SLEEP_TIGHT.md',
    'START_HERE.md',
    'STATUS.txt',
    'WHEN_YOU_WAKE_UP.md',
    'WINDOWS_SETUP.md',
    
    # Old/superseded reports
    'COHORT_ANALYSIS_OVERVIEW.md',
    'CORRECTED_ANALYSIS_NOTES.md',
    '00_CORRECT_REPORTS_SUMMARY.md',
    'FINAL_CLEANUP_REPORT.md',
    'FIND_MAPPING_GUIDE.md',
    
    # Old analysis scripts
    'apoe_grade_distribution_analysis.py',
    'apoe_grade_distribution_REAL_DATA.py',
    'rerun_survival_cox_RECLASSIFIED.py',
    
    # Temporary mapping scripts
    'extract_mapping.py',
    'quick_check.py',
    'search_all_mapping.py',
    'find_missing_cgga_ids.py',
    'verify_clinical_coverage.py',
    
    # Excel mapping files (data already incorporated)
    'HRA000071.xlsx',
    'HRA000073.xlsx',
    'HRA000074.xlsx',
    'HRR_CGGA_ID_MAPPING.csv',  # Superseded by FINAL version
]

# Files in clinical_analysis_outputs to REMOVE (superseded)
REMOVE_CLINICAL_FILES = [
    'apoe_grade_distribution_REAL_DATA.png',
    'apoe_grade_REAL_DATA_tables.xlsx',
    'clinical_analysis_manifest.json',
    'clinical_analysis_results.xlsx',
    'clinical_overview_panel.png',
    'COMPREHENSIVE_CLINICAL_SUMMARY.md',
    'PANEL_CUSTOMIZATION_NOTES.md',
    'SURVIVAL_PLOT_UPDATE.md',
    'TREATED_COHORT_NOTES.md',
    'survival_by_Grade.png',
    'survival_by_apoe_genotype.png',
]

def cleanup():
    """Perform cleanup"""
    
    print("=" * 80)
    print("PROJECT CLEANUP - Removing Temporary/Old Files")
    print("=" * 80)
    print()
    
    removed_count = 0
    kept_count = 0
    
    # Remove directories
    print("Removing old/temporary directories...")
    for dirname in REMOVE_DIRS:
        if os.path.exists(dirname):
            try:
                shutil.rmtree(dirname)
                print(f"  ✓ Removed: {dirname}/")
                removed_count += 1
            except Exception as e:
                print(f"  ✗ Failed to remove {dirname}: {e}")
        else:
            print(f"  - Not found: {dirname}/")
    print()
    
    # Remove files from root
    print("Removing old scripts and logs from root...")
    for filename in REMOVE_FILES:
        if os.path.exists(filename):
            try:
                os.remove(filename)
                print(f"  ✓ Removed: {filename}")
                removed_count += 1
            except Exception as e:
                print(f"  ✗ Failed to remove {filename}: {e}")
    print()
    
    # Remove superseded clinical analysis files
    print("Removing superseded files from clinical_analysis_outputs...")
    for filename in REMOVE_CLINICAL_FILES:
        filepath = os.path.join('clinical_analysis_outputs', filename)
        if os.path.exists(filepath):
            try:
                os.remove(filepath)
                print(f"  ✓ Removed: {filepath}")
                removed_count += 1
            except Exception as e:
                print(f"  ✗ Failed to remove {filepath}: {e}")
    print()
    
    # Count kept items
    print("Verifying essential files are kept...")
    for item in KEEP_ITEMS:
        if os.path.exists(item):
            print(f"  ✓ Kept: {item}")
            kept_count += 1
        else:
            print(f"  ! Not found: {item}")
    print()
    
    print("=" * 80)
    print("CLEANUP SUMMARY")
    print("=" * 80)
    print(f"Files/directories removed: {removed_count}")
    print(f"Essential items verified: {kept_count}/{len(KEEP_ITEMS)}")
    print()
    print("✅ Cleanup complete!")
    print()

if __name__ == "__main__":
    cleanup()

