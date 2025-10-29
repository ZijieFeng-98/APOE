#!/usr/bin/env python3
"""Recreate main figure: Treatment | ε4 Carrier Status"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
import warnings
warnings.filterwarnings('ignore')

plt.style.use('seaborn-v0_8-darkgrid')

df = pd.read_csv('clinical_analysis_outputs/RECLASSIFIED_ANALYSIS/25_samples_RECLASSIFIED.csv')
df = df[df['IDH_Status'] == 'Wildtype'].copy()
df = df[df['OS_days'].notna() & (df['OS_days'] > 0)].copy()
df['OS_months'] = df['OS_days'] / 30.44

clinical = pd.read_csv('DATA/CGGA.WEseq_286_clinical.20200506.txt', sep='\t')
df = df.merge(clinical[['CGGA_ID', 'Radio_status (treated=1;un-treated=0)', 
                        'Chemo_status (TMZ treated=1;un-treated=0)']], on='CGGA_ID', how='left')

df['Radio'] = pd.to_numeric(df['Radio_status (treated=1;un-treated=0)'], errors='coerce')
df['Chemo'] = pd.to_numeric(df['Chemo_status (TMZ treated=1;un-treated=0)'], errors='coerce')

df['Treatment_Group'] = 'Unknown'
df.loc[(df['Radio'] == 1) & (df['Chemo'] == 1), 'Treatment_Group'] = 'Radio+Chemo'
df.loc[(df['Radio'] == 1) & (df['Chemo'] == 0), 'Treatment_Group'] = 'Radio only'
df.loc[(df['Radio'] == 0) & (df['Chemo'] == 1), 'Treatment_Group'] = 'Chemo only'

df['e4_carrier'] = df['APOE_Genotype'].str.contains('ε4', na=False).astype(int)
df['Carrier_Status'] = df['e4_carrier'].map({1: 'ε4 carrier', 0: 'ε4 non-carrier'})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
kmf = KaplanMeierFitter()

# Panel A: Treatment
colors_treat = {'Radio+Chemo': '#2E7D32', 'Radio only': '#1976D2', 'Chemo only': '#F57C00'}
for treatment in ['Radio+Chemo', 'Radio only', 'Chemo only']:
    treat_df = df[df['Treatment_Group'] == treatment]
    if len(treat_df) > 0:
        kmf.fit(treat_df['OS_months'], treat_df['Death'], label=f"{treatment} (n={len(treat_df)})")
        kmf.plot_survival_function(ax=ax1, ci_show=True, color=colors_treat[treatment], linewidth=2.5)

ax1.set_xlabel('Time (months)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Survival Probability', fontsize=12, fontweight='bold')
ax1.set_title('A. Survival by Treatment Modality\n(GBM Cohort, n=20, Reclassified)', fontsize=12, fontweight='bold', pad=10)
ax1.legend(loc='best', fontsize=10, frameon=True, shadow=True)
ax1.grid(True, alpha=0.3)
ax1.set_ylim([0, 1.05])

# Panel B: ε4 Carrier
colors_carrier = {'ε4 non-carrier': '#1976D2', 'ε4 carrier': '#C62828'}
for status in ['ε4 non-carrier', 'ε4 carrier']:
    status_df = df[df['Carrier_Status'] == status]
    if len(status_df) > 0:
        kmf.fit(status_df['OS_months'], status_df['Death'], label=f"{status} (n={len(status_df)})")
        kmf.plot_survival_function(ax=ax2, ci_show=True, color=colors_carrier[status], linewidth=3.0)

ax2.set_xlabel('Time (months)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Survival Probability', fontsize=12, fontweight='bold')
ax2.set_title('B. Survival by APOE ε4 Carrier Status\n(GBM Cohort, n=20, Reclassified)', fontsize=12, fontweight='bold', pad=10)
ax2.legend(loc='best', fontsize=10, frameon=True, shadow=True)
ax2.grid(True, alpha=0.3)
ax2.set_ylim([0, 1.05])

results_apoe = multivariate_logrank_test(df['OS_months'], df['Carrier_Status'], df['Death'])
ax2.text(0.98, 0.05, f'Log-rank p={results_apoe.p_value:.4f}',
        transform=ax2.transAxes, fontsize=11, fontweight='bold',
        verticalalignment='bottom', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))

plt.suptitle('Therapy-Stratified Kaplan-Meier Survival Curves (GBM Only)\n' +
            'IDH-Wildtype Patients - Grade Reclassified to WHO IV',
            fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('cox_therapy_analysis_outputs/gbm_only/Fig3_Therapy_Stratified_KM_GBM_Only.png', dpi=300, bbox_inches='tight')
print("✅ Recreated: Fig3_Therapy_Stratified_KM_GBM_Only.png")

