#!/usr/bin/env python3
"""
Survival, Cox, and Treatment Analysis - IDH-WILDTYPE PATIENTS ONLY
20 IDH-wildtype patients (true GBM cohort)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("Set2")

def load_idh_wt_data():
    """Load IDH-wildtype patients only"""
    df = pd.read_csv('clinical_analysis_outputs/RECLASSIFIED_ANALYSIS/25_samples_RECLASSIFIED.csv')
    
    # Filter to IDH-wildtype only
    df_wt = df[df['IDH_Status'] == 'Wildtype'].copy()
    
    # Filter to patients with survival data
    df_wt = df_wt[df_wt['OS_days'].notna() & (df_wt['OS_days'] > 0)].copy()
    
    # Convert survival to months
    df_wt['OS_months'] = df_wt['OS_days'] / 30.44
    
    print(f"Loaded {len(df_wt)} IDH-wildtype (true GBM) patients")
    print(f"Deaths: {df_wt['Death'].sum():.0f}")
    print(f"Censored: {(df_wt['Death'] == 0).sum()}")
    print()
    
    return df_wt

def perform_km_by_apoe(df):
    """Kaplan-Meier by APOE genotype - IDH-WT only"""
    
    print("=" * 80)
    print("KAPLAN-MEIER: IDH-WILDTYPE PATIENTS BY APOE GENOTYPE")
    print("=" * 80)
    print()
    
    kmf = KaplanMeierFitter()
    fig, ax = plt.subplots(figsize=(12, 8))
    
    genotype_styles = {
        'ε2/ε2': {'color': '#2E7D32', 'linestyle': '-', 'marker': 'o'},
        'ε2/ε3': {'color': '#388E3C', 'linestyle': '-', 'marker': 's'},
        'ε2/ε4': {'color': '#FBC02D', 'linestyle': '--', 'marker': '^'},
        'ε3/ε3': {'color': '#1976D2', 'linestyle': '-', 'marker': 'D'},
        'ε3/ε4': {'color': '#F57C00', 'linestyle': '--', 'marker': 'v'},
        'ε4/ε4': {'color': '#C62828', 'linestyle': '-', 'marker': 'p'}
    }
    
    summary_data = []
    
    for genotype in sorted(df['APOE_Genotype'].unique()):
        genotype_df = df[df['APOE_Genotype'] == genotype]
        
        if len(genotype_df) > 0:
            kmf.fit(
                genotype_df['OS_months'], 
                genotype_df['Death'],
                label=f"{genotype} (n={len(genotype_df)})"
            )
            
            style = genotype_styles.get(genotype, {'color': 'gray', 'linestyle': '-', 'marker': 'o'})
            kmf.plot_survival_function(
                ax=ax,
                ci_show=True,
                color=style['color'],
                linestyle=style['linestyle'],
                linewidth=2.5,
                marker=style['marker'],
                markevery=5,
                markersize=8
            )
            
            median_survival = kmf.median_survival_time_
            deaths = genotype_df['Death'].sum()
            
            summary_data.append({
                'Genotype': genotype,
                'N': len(genotype_df),
                'Deaths': int(deaths),
                'Death_Rate': f"{deaths}/{len(genotype_df)} ({deaths/len(genotype_df)*100:.1f}%)",
                'Median_Survival_Months': f"{median_survival:.1f}" if not np.isnan(median_survival) else "Not reached"
            })
            
            print(f"{genotype}: n={len(genotype_df)}, Deaths={deaths:.0f}, Median={median_survival:.1f} months" if not np.isnan(median_survival) else f"{genotype}: n={len(genotype_df)}, Deaths={deaths:.0f}, Median=Not reached")
    
    ax.set_xlabel('Time (months)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Overall Survival Probability', fontsize=12, fontweight='bold')
    ax.set_title('Kaplan-Meier Survival: IDH-Wildtype GBM by APOE Genotype\n(n=20 IDH-WT patients)', 
                 fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=10, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    plt.savefig('clinical_analysis_outputs/RECLASSIFIED_ANALYSIS/survival_IDH_WT_by_APOE.png', 
                dpi=300, bbox_inches='tight')
    print("\n✅ Saved: survival_IDH_WT_by_APOE.png")
    
    # Log-rank test
    if len(df['APOE_Genotype'].unique()) > 1:
        results = multivariate_logrank_test(
            df['OS_months'],
            df['APOE_Genotype'],
            df['Death']
        )
        print(f"\nLog-rank test: Chi²={results.test_statistic:.2f}, p={results.p_value:.4f}")
    
    return fig, pd.DataFrame(summary_data)

def perform_km_by_treatment(df):
    """Kaplan-Meier by treatment - IDH-WT only"""
    
    print("\n" + "=" * 80)
    print("KAPLAN-MEIER: IDH-WILDTYPE PATIENTS BY TREATMENT")
    print("=" * 80)
    print()
    
    # Load original clinical data to get treatment info
    clinical = pd.read_csv('DATA/CGGA.WEseq_286_clinical.20200506.txt', sep='\t')
    df_merged = df.merge(clinical[['CGGA_ID', 'Radio_status (treated=1;un-treated=0)', 
                                    'Chemo_status (TMZ treated=1;un-treated=0)']], 
                         on='CGGA_ID', how='left')
    
    # Create treatment groups
    df_merged['Radio'] = pd.to_numeric(df_merged['Radio_status (treated=1;un-treated=0)'], errors='coerce')
    df_merged['Chemo'] = pd.to_numeric(df_merged['Chemo_status (TMZ treated=1;un-treated=0)'], errors='coerce')
    
    df_merged['Treatment_Group'] = 'Unknown'
    df_merged.loc[(df_merged['Radio'] == 1) & (df_merged['Chemo'] == 1), 'Treatment_Group'] = 'Radio+Chemo'
    df_merged.loc[(df_merged['Radio'] == 1) & (df_merged['Chemo'] == 0), 'Treatment_Group'] = 'Radio only'
    df_merged.loc[(df_merged['Radio'] == 0) & (df_merged['Chemo'] == 1), 'Treatment_Group'] = 'Chemo only'
    df_merged.loc[(df_merged['Radio'] == 0) & (df_merged['Chemo'] == 0), 'Treatment_Group'] = 'No treatment'
    
    print("Treatment distribution:")
    print(df_merged['Treatment_Group'].value_counts())
    print()
    
    kmf = KaplanMeierFitter()
    fig, ax = plt.subplots(figsize=(12, 8))
    
    colors = {'Radio+Chemo': '#2E7D32', 'Radio only': '#1976D2', 
              'Chemo only': '#F57C00', 'No treatment': '#C62828', 'Unknown': 'gray'}
    
    summary_data = []
    
    for treatment in ['Radio+Chemo', 'Radio only', 'Chemo only', 'No treatment', 'Unknown']:
        treat_df = df_merged[df_merged['Treatment_Group'] == treatment]
        
        if len(treat_df) > 0:
            kmf.fit(
                treat_df['OS_months'], 
                treat_df['Death'],
                label=f"{treatment} (n={len(treat_df)})"
            )
            
            kmf.plot_survival_function(
                ax=ax,
                ci_show=True,
                color=colors[treatment],
                linewidth=3
            )
            
            median_survival = kmf.median_survival_time_
            deaths = treat_df['Death'].sum()
            
            summary_data.append({
                'Treatment': treatment,
                'N': len(treat_df),
                'Deaths': int(deaths),
                'Death_Rate': f"{deaths}/{len(treat_df)} ({deaths/len(treat_df)*100:.1f}%)",
                'Median_Survival_Months': f"{median_survival:.1f}" if not np.isnan(median_survival) else "Not reached"
            })
            
            print(f"{treatment}: n={len(treat_df)}, Deaths={deaths:.0f}, Median={median_survival:.1f} months" if not np.isnan(median_survival) else f"{treatment}: n={len(treat_df)}, Deaths={deaths:.0f}, Median=Not reached")
    
    ax.set_xlabel('Time (months)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Overall Survival Probability', fontsize=12, fontweight='bold')
    ax.set_title('Kaplan-Meier Survival: IDH-Wildtype GBM by Treatment\n(n=20 IDH-WT patients)', 
                 fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=10, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    plt.savefig('clinical_analysis_outputs/RECLASSIFIED_ANALYSIS/survival_IDH_WT_by_treatment.png', 
                dpi=300, bbox_inches='tight')
    print("\n✅ Saved: survival_IDH_WT_by_treatment.png")
    
    return fig, pd.DataFrame(summary_data), df_merged

def perform_cox_idh_wt(df):
    """Cox regression - IDH-WT only"""
    
    print("\n" + "=" * 80)
    print("COX REGRESSION: IDH-WILDTYPE PATIENTS ONLY")
    print("=" * 80)
    print()
    
    # Prepare data
    cox_df = df.copy()
    cox_df['e4_carrier'] = cox_df['APOE_Genotype'].str.contains('ε4', na=False).astype(int)
    cox_df['Age'] = pd.to_numeric(cox_df['Age'], errors='coerce')
    cox_df['Male'] = (cox_df['Gender'] == 'Male').astype(int)
    
    # Load treatment data
    clinical = pd.read_csv('DATA/CGGA.WEseq_286_clinical.20200506.txt', sep='\t')
    cox_df = cox_df.merge(clinical[['CGGA_ID', 'Radio_status (treated=1;un-treated=0)', 
                                     'Chemo_status (TMZ treated=1;un-treated=0)']], 
                          on='CGGA_ID', how='left')
    
    cox_df['Radio'] = pd.to_numeric(cox_df['Radio_status (treated=1;un-treated=0)'], errors='coerce')
    cox_df['Chemo'] = pd.to_numeric(cox_df['Chemo_status (TMZ treated=1;un-treated=0)'], errors='coerce')
    cox_df['Radio_Chemo'] = ((cox_df['Radio'] == 1) & (cox_df['Chemo'] == 1)).astype(int)
    
    # Clean data
    cox_df_clean = cox_df.dropna(subset=['OS_months', 'Death', 'Age'])
    
    print(f"Cox analysis on {len(cox_df_clean)} IDH-WT patients\n")
    
    all_results = []
    
    # Model 1: APOE ε4 carrier only
    print("MODEL 1: APOE ε4 Carrier Status")
    cph1 = CoxPHFitter()
    try:
        cph1.fit(
            cox_df_clean[['OS_months', 'Death', 'e4_carrier']],
            duration_col='OS_months',
            event_col='Death'
        )
        print(cph1.summary[['coef', 'exp(coef)', 'se(coef)', 'p']])
        model1 = cph1.summary[['coef', 'exp(coef)', 'se(coef)', 'p']].copy()
        model1['Model'] = 'ε4_carrier_only'
        all_results.append(model1)
    except Exception as e:
        print(f"Model 1 failed: {e}")
    
    # Model 2: Clinical factors (Age, Gender)
    print("\nMODEL 2: Age + Gender")
    cph2 = CoxPHFitter()
    try:
        cph2.fit(
            cox_df_clean[['OS_months', 'Death', 'Age', 'Male']],
            duration_col='OS_months',
            event_col='Death'
        )
        print(cph2.summary[['coef', 'exp(coef)', 'se(coef)', 'p']])
        model2 = cph2.summary[['coef', 'exp(coef)', 'se(coef)', 'p']].copy()
        model2['Model'] = 'Age_Gender'
        all_results.append(model2)
    except Exception as e:
        print(f"Model 2 failed: {e}")
    
    # Model 3: Treatment only
    print("\nMODEL 3: Treatment (Radio+Chemo combined)")
    cox_df_treat = cox_df_clean.dropna(subset=['Radio_Chemo'])
    if len(cox_df_treat) > 0:
        cph3 = CoxPHFitter()
        try:
            cph3.fit(
                cox_df_treat[['OS_months', 'Death', 'Radio_Chemo']],
                duration_col='OS_months',
                event_col='Death'
            )
            print(cph3.summary[['coef', 'exp(coef)', 'se(coef)', 'p']])
            model3 = cph3.summary[['coef', 'exp(coef)', 'se(coef)', 'p']].copy()
            model3['Model'] = 'Treatment'
            all_results.append(model3)
        except Exception as e:
            print(f"Model 3 failed: {e}")
    
    # Model 4: Multivariable (Age, Gender, ε4, Treatment)
    print("\nMODEL 4: Multivariable (Age + Gender + ε4 + Treatment)")
    cox_df_multi = cox_df_clean.dropna(subset=['Radio_Chemo'])
    if len(cox_df_multi) > 5:
        cph4 = CoxPHFitter()
        try:
            cph4.fit(
                cox_df_multi[['OS_months', 'Death', 'Age', 'Male', 'e4_carrier', 'Radio_Chemo']],
                duration_col='OS_months',
                event_col='Death'
            )
            print(cph4.summary[['coef', 'exp(coef)', 'se(coef)', 'p']])
            model4 = cph4.summary[['coef', 'exp(coef)', 'se(coef)', 'p']].copy()
            model4['Model'] = 'Multivariable'
            all_results.append(model4)
        except Exception as e:
            print(f"Model 4 failed: {e}")
    
    if all_results:
        combined = pd.concat(all_results)
        return combined
    else:
        return None

def main():
    """Execute IDH-WT only analysis"""
    
    print("=" * 80)
    print("SURVIVAL & COX ANALYSIS - IDH-WILDTYPE PATIENTS ONLY")
    print("True GBM Cohort (n=20)")
    print("=" * 80)
    print()
    
    # Load IDH-WT patients
    df = load_idh_wt_data()
    
    # KM by APOE
    fig1, summary_apoe = perform_km_by_apoe(df)
    
    # KM by treatment
    fig2, summary_treat, df_treat = perform_km_by_treatment(df)
    
    # Cox regression
    cox_results = perform_cox_idh_wt(df)
    
    # Save results
    print("\n" + "=" * 80)
    print("SAVING RESULTS")
    print("=" * 80)
    
    with pd.ExcelWriter('clinical_analysis_outputs/RECLASSIFIED_ANALYSIS/survival_cox_IDH_WT_ONLY.xlsx') as writer:
        summary_apoe.to_excel(writer, sheet_name='KM_by_APOE', index=False)
        summary_treat.to_excel(writer, sheet_name='KM_by_Treatment', index=False)
        if cox_results is not None:
            cox_results.to_excel(writer, sheet_name='Cox_Regression')
        
        # Cohort summary
        summary = pd.DataFrame({
            'Metric': [
                'Total IDH-WT Patients',
                'Deaths',
                'Censored',
                'Median Follow-up (months)',
                'ε4 Carriers',
                'Radio+Chemo treated'
            ],
            'Value': [
                len(df),
                int(df['Death'].sum()),
                int((df['Death'] == 0).sum()),
                f"{df['OS_months'].median():.1f}",
                df['APOE_Genotype'].str.contains('ε4', na=False).sum(),
                df_treat['Radio_Chemo'].sum() if 'Radio_Chemo' in df_treat.columns else 'N/A'
            ]
        })
        summary.to_excel(writer, sheet_name='Cohort_Summary', index=False)
    
    print("\n✅ Saved: survival_cox_IDH_WT_ONLY.xlsx")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE - IDH-WILDTYPE ONLY")
    print("=" * 80)
    print("\nGenerated files:")
    print("  1. survival_IDH_WT_by_APOE.png")
    print("  2. survival_IDH_WT_by_treatment.png")
    print("  3. survival_cox_IDH_WT_ONLY.xlsx")
    print()

if __name__ == "__main__":
    main()

