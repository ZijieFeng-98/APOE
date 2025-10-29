#!/usr/bin/env python3
"""
Reclassify ALL IDH-wildtype patients as WHO Grade IV (GBM)
This follows biological reality: IDH-WT gliomas are aggressive regardless of histological grade
"""

import pandas as pd

# Load data
apoe = pd.read_csv('APOE_GENOTYPES.csv')
mapping = pd.read_csv('FINAL_HRR_CGGA_MAPPING.csv')
clinical = pd.read_csv('DATA/CGGA.WEseq_286_clinical.20200506.txt', sep='\t')

# Merge
df = apoe.merge(mapping, left_on='Patient_ID', right_on='HRR_ID', how='left')
df = df.merge(clinical, left_on='CGGA_ID', right_on='CGGA_ID', how='left')

print("=" * 100)
print("RECLASSIFICATION: ALL IDH-WILDTYPE → WHO GRADE IV")
print("=" * 100)
print()

# Show current distribution
print("BEFORE RECLASSIFICATION:")
print(f"  WHO II:  {(df['Grade'] == 'WHO II').sum()} patients")
print(f"  WHO III: {(df['Grade'] == 'WHO III').sum()} patients")
print(f"  WHO IV:  {(df['Grade'] == 'WHO IV').sum()} patients")
print()

# Identify IDH-WT patients that are NOT Grade IV
idh_wt_low_grade = df[(df['IDH_mut_status'] == 'Wildtype') & (df['Grade'] != 'WHO IV')]
print(f"IDH-WILDTYPE patients with Grade II/III: {len(idh_wt_low_grade)}")
if len(idh_wt_low_grade) > 0:
    print()
    print("Patients to be RECLASSIFIED to WHO IV:")
    for _, row in idh_wt_low_grade.iterrows():
        print(f"  {row['Patient_ID']} ({row['CGGA_ID']}) - {row['APOE_genotype']} - "
              f"Current: {row['Grade']} → New: WHO IV")
    print()

# RECLASSIFY: All IDH-wildtype becomes WHO IV
df['Grade_Original'] = df['Grade'].copy()
df.loc[df['IDH_mut_status'] == 'Wildtype', 'Grade'] = 'WHO IV'

print("AFTER RECLASSIFICATION:")
print(f"  WHO II:  {(df['Grade'] == 'WHO II').sum()} patients (all IDH-mutant)")
print(f"  WHO III: {(df['Grade'] == 'WHO III').sum()} patients (all IDH-mutant)")
print(f"  WHO IV:  {(df['Grade'] == 'WHO IV').sum()} patients (includes all IDH-WT)")
print()

# Verify: Grade IV should = all IDH-WT + IDH-mut Grade IV
idh_wt_count = (df['IDH_mut_status'] == 'Wildtype').sum()
print(f"✅ Total IDH-wildtype: {idh_wt_count}")
print(f"✅ Total WHO IV after reclassification: {(df['Grade'] == 'WHO IV').sum()}")
print()

# Show final distribution by IDH
print("=" * 100)
print("FINAL GRADE DISTRIBUTION BY IDH STATUS")
print("=" * 100)
print()
grade_idh = pd.crosstab(df['Grade'], df['IDH_mut_status'], margins=True)
print(grade_idh)
print()

# Show APOE by reclassified grade
print("=" * 100)
print("APOE GENOTYPE BY RECLASSIFIED GRADE")
print("=" * 100)
print()
apoe_grade = pd.crosstab(df['Grade'], df['APOE_genotype'], margins=True)
print(apoe_grade)
print()

# Save reclassified data
print("Saving reclassified data...")
df_save = df[['Patient_ID', 'CGGA_ID', 'APOE_genotype', 'Grade_Original', 'Grade', 
              'IDH_mut_status', 'Histology', 'Age', 'Gender', 'OS', 'Censor (alive=0; dead=1)']].copy()
df_save.columns = ['HRR_ID', 'CGGA_ID', 'APOE_Genotype', 'Grade_Original', 'Grade_Reclassified', 
                   'IDH_Status', 'Histology', 'Age', 'Gender', 'OS_days', 'Death']

with pd.ExcelWriter('clinical_analysis_outputs/25_samples_RECLASSIFIED_metadata.xlsx') as writer:
    df_save.to_excel(writer, sheet_name='Reclassified_Data', index=False)
    grade_idh.to_excel(writer, sheet_name='Grade_by_IDH')
    apoe_grade.to_excel(writer, sheet_name='APOE_by_Grade')

print("✅ Saved: clinical_analysis_outputs/25_samples_RECLASSIFIED_metadata.xlsx")
print()

# Also save as CSV for easy use
df_save.to_csv('clinical_analysis_outputs/25_samples_RECLASSIFIED.csv', index=False)
print("✅ Saved: clinical_analysis_outputs/25_samples_RECLASSIFIED.csv")
print()

print("=" * 100)
print("RECLASSIFICATION COMPLETE")
print("=" * 100)
print()
print("KEY POINTS:")
print("  • ALL IDH-wildtype patients now classified as WHO IV (GBM)")
print("  • WHO II and WHO III now contain ONLY IDH-mutant patients")
print(f"  • Total GBM (WHO IV): {(df['Grade'] == 'WHO IV').sum()} patients")
print(f"    - IDH-wildtype: {((df['Grade'] == 'WHO IV') & (df['IDH_mut_status'] == 'Wildtype')).sum()}")
print(f"    - IDH-mutant: {((df['Grade'] == 'WHO IV') & (df['IDH_mut_status'] == 'Mutant')).sum()}")
print()

