#!/usr/bin/env python3
"""
Run the complete clinical analysis WITHOUT risk category bias
Only analyze by individual APOE genotypes, not risk categories
"""
import sys
from pathlib import Path
import pandas as pd
import json

# Add apoe_analysis to path
sys.path.insert(0, str(Path(__file__).parent))

from apoe_analysis.clinical_analysis import (
    load_clinical_table,
    load_apoe_genotypes,
    _load_optional_mapping,
    merge_apoe_clinical,
    create_clinical_summary,
    stratified_analysis,
    perform_km_analysis,
    cox_regression_analysis,
    analyze_therapy_combinations,
    perform_logrank_tests,
    create_clinical_figure_panel,
)

def main():
    print("=" * 80)
    print("🔬 APOE CLINICAL ANALYSIS - NO RISK CATEGORY BIAS")
    print("=" * 80)
    print()
    print("Analysis will include:")
    print("  ✓ APOE genotype stratification (ε2/ε2, ε2/ε3, ε2/ε4, ε3/ε3, ε3/ε4, ε4/ε4)")
    print("  ✓ Cox regression")
    print("  ✓ Kaplan-Meier by genotype")
    print("  ✗ NO risk category (Protective/Baseline/Increased) - REMOVED")
    print()
    
    # Configuration
    clinical_path = Path("DATA/CGGA.WEseq_286_clinical.20200506.txt")
    output_dir = Path("clinical_analysis_outputs")
    apoe_manifest = Path("APOE_GENOTYPES.csv")
    id_mapping = Path("HRR_CGGA_ID_MAPPING.csv")
    time_column = "OS"
    event_column = "Censor (alive=0; dead=1)"
    group_column = "Grade"
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("📊 Loading data...")
    clinical_df = load_clinical_table(clinical_path)
    
    # Load APOE genotypes WITHOUT risk category assignment
    apoe_raw = pd.read_csv(apoe_manifest)
    # Add neutral Risk_Category (required by merge function, but we won't use it)
    apoe_raw['Risk_Category'] = 'Not_Analyzed'
    print(f"   ✓ Loaded {len(apoe_raw)} APOE genotypes")
    print(f"   Genotypes: {apoe_raw['APOE_genotype'].value_counts().to_dict()}")
    
    id_mapping_df = _load_optional_mapping(id_mapping)
    
    # Merge data
    print("\n🔗 Merging APOE genotypes with clinical data...")
    merged_df = merge_apoe_clinical(
        apoe_raw,
        clinical_df,
        id_mapping=id_mapping_df
    )
    print(f"   ✓ Merged cohort: {len(merged_df)} patients")
    
    # Overall summary - USE MERGED DATA (our 25 patients only!)
    print("\n📈 Generating cohort summaries...")
    overall_summary = create_clinical_summary(
        merged_df, time_column=time_column, event_column=event_column
    )
    
    # Stratified by Grade - USE MERGED DATA
    stratified_summary = stratified_analysis(
        merged_df,
        group_column,
        time_column=time_column,
        event_column=event_column,
    )
    
    # Stratified by APOE genotype
    apoe_summary = stratified_analysis(
        merged_df,
        "APOE_genotype",
        time_column=time_column,
        event_column=event_column,
    )
    
    # Kaplan-Meier plots
    print("\n📊 Generating Kaplan-Meier plots...")
    
    # By Grade - USE MERGED DATA (our 25 patients)
    import matplotlib.pyplot as plt
    km_axes = perform_km_analysis(
        merged_df,
        group_column,
        time_column=time_column,
        event_column=event_column,
    )
    km_figure = km_axes.get_figure()
    km_path = output_dir / f"survival_by_{group_column}.png"
    km_figure.savefig(km_path, dpi=300, bbox_inches="tight")
    print(f"   ✓ Saved: {km_path.name}")
    plt.close(km_figure)
    
    # By APOE genotype ONLY (no risk category)
    apoe_km_axes = perform_km_analysis(
        merged_df,
        "APOE_genotype",
        time_column=time_column,
        event_column=event_column,
    )
    apoe_km_figure = apoe_km_axes.get_figure()
    apoe_km_path = output_dir / "survival_by_apoe_genotype.png"
    apoe_km_figure.savefig(apoe_km_path, dpi=300, bbox_inches="tight")
    print(f"   ✓ Saved: {apoe_km_path.name}")
    plt.close(apoe_km_figure)
    
    # Cox regression - USE MERGED DATA
    print("\n📉 Running Cox regression...")
    cph = cox_regression_analysis(
        merged_df, time_column=time_column, event_column=event_column
    )
    print("   ✓ Cox model fitted")
    
    # Therapy analysis - USE MERGED DATA
    print("\n💊 Analyzing therapy combinations...")
    therapy_summary = analyze_therapy_combinations(
        merged_df, time_column=time_column, event_column=event_column
    )
    
    # Log-rank tests - USE MERGED DATA
    print("\n📊 Running log-rank tests...")
    logrank_grade = perform_logrank_tests(
        merged_df,
        group_column,
        time_column=time_column,
        event_column=event_column,
    )
    
    logrank_apoe = perform_logrank_tests(
        merged_df,
        "APOE_genotype",
        time_column=time_column,
        event_column=event_column,
    )
    
    # Clinical overview panel - USE MERGED DATA (our 25 patients, not all 286!)
    print("\n🎨 Creating overview panel...")
    panel_fig = create_clinical_figure_panel(merged_df, time_column=time_column)
    panel_path = output_dir / "clinical_overview_panel.png"
    panel_fig.savefig(panel_path, dpi=300, bbox_inches="tight")
    print(f"   ✓ Saved: {panel_path.name}")
    plt.close(panel_fig)
    
    # Save all results to Excel
    print("\n💾 Saving results to Excel...")
    results_book = output_dir / "clinical_analysis_results.xlsx"
    with pd.ExcelWriter(results_book) as writer:
        overall_summary.to_excel(writer, sheet_name="Overall_Summary", index=False)
        stratified_summary.to_excel(writer, sheet_name=f"By_{group_column}", index=False)
        apoe_summary.to_excel(writer, sheet_name="By_APOE_Genotype", index=False)
        therapy_summary.to_excel(writer, sheet_name="Therapy_Analysis", index=False)
        logrank_grade.to_excel(writer, sheet_name="LogRank_Grade", index=False)
        logrank_apoe.to_excel(writer, sheet_name="LogRank_APOE", index=False)
        merged_df.to_excel(writer, sheet_name="Merged_APOE_Clinical", index=False)
        cph.summary.to_excel(writer, sheet_name="Cox_Model")
    print(f"   ✓ Saved: {results_book.name}")
    
    # Create manifest
    manifest = {
        "note": "Analysis WITHOUT risk category bias - genotypes analyzed individually",
        "overall_summary": str(results_book),
        "survival_by_grade": str(km_path),
        "survival_by_apoe_genotype": str(apoe_km_path),
        "overview_panel": str(panel_path),
        "cox_model": str(results_book),
        "therapy_summary": str(results_book),
        "logrank_tests": str(results_book),
        "merged_cohort": str(results_book),
        "genotypes_included": apoe_raw['APOE_genotype'].value_counts().to_dict()
    }
    manifest_path = output_dir / "clinical_analysis_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    
    print()
    print("=" * 80)
    print("✅ ANALYSIS COMPLETE!")
    print("=" * 80)
    print()
    print(f"📁 Results saved to: {output_dir.absolute()}")
    print()
    print("Generated files:")
    print("  - clinical_analysis_results.xlsx (All tables)")
    print("  - survival_by_Grade.png (KM plot by WHO grade)")
    print("  - survival_by_apoe_genotype.png (KM plot by APOE genotype - NO RISK CATEGORIES)")
    print("  - clinical_overview_panel.png (Demographics overview)")
    print("  - clinical_analysis_manifest.json (File manifest)")
    print()
    print("📊 APOE Genotype Distribution:")
    for genotype, count in apoe_raw['APOE_genotype'].value_counts().items():
        print(f"   {genotype}: {count} patients")
    print()
    return 0

if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print()
        print("=" * 80)
        print("❌ ANALYSIS FAILED")
        print("=" * 80)
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

