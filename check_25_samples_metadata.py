#!/usr/bin/env python3
"""
Check IDH status and WHO grade metadata for all 25 genotyped patients
"""

import pandas as pd

# Load data
apoe = pd.read_csv('APOE_GENOTYPES.csv')
mapping = pd.read_csv('FINAL_HRR_CGGA_MAPPING.csv')
clinical = pd.read_csv('DATA/CGGA.WEseq_286_clinical.20200506.txt', sep='\t')

# Merge
df = apoe.merge(mapping, left_on='Patient_ID', right_on='HRR_ID', how='left')
df = df.merge(clinical, left_on='CGGA_ID', right_on='CGGA_ID', how='left')

# Select relevant columns
metadata = df[['Patient_ID', 'CGGA_ID', 'APOE_genotype', 'Grade', 'IDH_mut_status', 
               'Histology', 'Age', 'Gender', 'OS', 'Censor (alive=0; dead=1)']].copy()
metadata.columns = ['HRR_ID', 'CGGA_ID', 'APOE_Genotype', 'WHO_Grade', 'IDH_Status', 
                    'Histology', 'Age', 'Gender', 'OS_days', 'Death']

# Display all
print('=' * 100)
print('METADATA FOR ALL 25 GENOTYPED PATIENTS')
print('=' * 100)
print()
print(metadata.to_string(index=False))
print()

# Summary statistics
print('=' * 100)
print('SUMMARY STATISTICS')
print('=' * 100)
print()
print('WHO GRADE DISTRIBUTION:')
print(metadata['WHO_Grade'].value_counts().sort_index())
print()
print('IDH STATUS DISTRIBUTION:')
print(metadata['IDH_Status'].value_counts())
print()
print('HISTOLOGY DISTRIBUTION:')
print(metadata['Histology'].value_counts())
print()

# Cross-tabulation
print('=' * 100)
print('IDH STATUS BY WHO GRADE')
print('=' * 100)
print()
idh_grade = pd.crosstab(metadata['WHO_Grade'], metadata['IDH_Status'], margins=True)
print(idh_grade)
print()

# APOE by IDH
print('=' * 100)
print('APOE GENOTYPE BY IDH STATUS')
print('=' * 100)
print()
apoe_idh = pd.crosstab(metadata['APOE_Genotype'], metadata['IDH_Status'], margins=True)
print(apoe_idh)
print()

# WHO 2021 GBM criteria check
print('=' * 100)
print('WHO 2021 GBM CRITERIA CHECK')
print('=' * 100)
print()
print('TRUE GBM (IDH-wildtype Grade IV):')
true_gbm = metadata[(metadata['WHO_Grade'] == 'WHO IV') & (metadata['IDH_Status'] == 'Wildtype')]
print(f"  Count: {len(true_gbm)}")
print()
print('IDH-mutant Grade IV (NOT GBM under WHO 2021):')
idh_mut_iv = metadata[(metadata['WHO_Grade'] == 'WHO IV') & (metadata['IDH_Status'] == 'Mutant')]
print(f"  Count: {len(idh_mut_iv)}")
if len(idh_mut_iv) > 0:
    print('  Patients:')
    for _, row in idh_mut_iv.iterrows():
        print(f"    {row['HRR_ID']} ({row['CGGA_ID']}) - {row['APOE_Genotype']}")
print()

# Save to file
print('Saving detailed metadata to: clinical_analysis_outputs/25_samples_metadata.xlsx')
with pd.ExcelWriter('clinical_analysis_outputs/25_samples_metadata.xlsx') as writer:
    metadata.to_excel(writer, sheet_name='Patient_Metadata', index=False)
    idh_grade.to_excel(writer, sheet_name='IDH_by_Grade')
    apoe_idh.to_excel(writer, sheet_name='APOE_by_IDH')
    
    # Add WHO 2021 GBM sheet
    who2021_summary = pd.DataFrame({
        'Category': ['True GBM (IDH-wt Grade IV)', 'IDH-mut Grade IV (not GBM)', 
                    'Grade II', 'Grade III', 'Total'],
        'Count': [
            len(metadata[(metadata['WHO_Grade'] == 'WHO IV') & (metadata['IDH_Status'] == 'Wildtype')]),
            len(metadata[(metadata['WHO_Grade'] == 'WHO IV') & (metadata['IDH_Status'] == 'Mutant')]),
            len(metadata[metadata['WHO_Grade'] == 'WHO II']),
            len(metadata[metadata['WHO_Grade'] == 'WHO III']),
            len(metadata)
        ]
    })
    who2021_summary.to_excel(writer, sheet_name='WHO2021_Classification', index=False)

print('✅ Complete!')

