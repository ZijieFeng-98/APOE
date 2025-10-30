#!/usr/bin/env python3
"""
Rerun ALL analyses using RECLASSIFIED grades (IDH-WT → WHO IV)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

def load_reclassified_data():
    """Load the reclassified data"""
    df = pd.read_csv('clinical_analysis_outputs/25_samples_RECLASSIFIED.csv')
    
    # Extract alleles
    def extract_alleles(genotype):
        if pd.isna(genotype):
            return None, None
        alleles = genotype.split('/')
        return alleles[0], alleles[1]
    
    df[['Allele1', 'Allele2']] = df['APOE_Genotype'].apply(
        lambda x: pd.Series(extract_alleles(x))
    )
    
    print(f"Loaded {len(df)} patients with reclassified grades")
    print(f"Grade distribution:")
    print(df['Grade_Reclassified'].value_counts().sort_index())
    print()
    
    return df

def analyze_genotype_by_grade(df):
    """Analyze APOE genotype distribution across reclassified WHO grades"""
    
    # Filter to valid grades
    df_valid = df[df['Grade_Reclassified'].isin(['WHO II', 'WHO III', 'WHO IV'])].copy()
    
    if len(df_valid) == 0:
        return None, None, None, None
    
    # Create contingency table
    contingency = pd.crosstab(
        df_valid['Grade_Reclassified'], 
        df_valid['APOE_Genotype'],
        margins=True,
        margins_name='Total'
    )
    
    # Calculate percentages
    contingency_pct = pd.crosstab(
        df_valid['Grade_Reclassified'], 
        df_valid['APOE_Genotype'],
        normalize='index'
    ) * 100
    
    # Chi-square test
    try:
        chi2, p_value, dof, expected = stats.chi2_contingency(
            contingency.iloc[:-1, :-1]
        )
    except:
        chi2, p_value = None, None
    
    return contingency, contingency_pct, chi2, p_value

def analyze_allele_frequency_by_grade(df):
    """Calculate allele frequencies for reclassified grades"""
    
    df_valid = df[df['Grade_Reclassified'].isin(['WHO II', 'WHO III', 'WHO IV'])].copy()
    
    if len(df_valid) == 0:
        return pd.DataFrame()
    
    allele_freq = {}
    
    for grade in ['WHO II', 'WHO III', 'WHO IV']:
        grade_df = df_valid[df_valid['Grade_Reclassified'] == grade]
        
        if len(grade_df) == 0:
            continue
            
        all_alleles = list(grade_df['Allele1'].dropna()) + list(grade_df['Allele2'].dropna())
        total_alleles = len(all_alleles)
        
        if total_alleles > 0:
            freq = {
                'ε2': all_alleles.count('ε2') / total_alleles * 100,
                'ε3': all_alleles.count('ε3') / total_alleles * 100,
                'ε4': all_alleles.count('ε4') / total_alleles * 100,
                'Total_alleles': total_alleles,
                'N_patients': len(grade_df)
            }
        else:
            freq = {'ε2': 0, 'ε3': 0, 'ε4': 0, 'Total_alleles': 0, 'N_patients': 0}
        
        allele_freq[grade] = freq
    
    return pd.DataFrame(allele_freq).T

def create_visualization_panel(df, contingency_pct, allele_freq):
    """Create comprehensive visualization panel"""
    
    fig = plt.figure(figsize=(16, 10))
    
    # 1. Genotype distribution by grade (stacked bar)
    ax1 = plt.subplot(2, 3, 1)
    if contingency_pct is not None and len(contingency_pct) > 0:
        contingency_pct.plot(kind='bar', stacked=True, ax=ax1)
        ax1.set_title('APOE Genotype Distribution by WHO Grade\n(IDH-WT Reclassified to Grade IV)', 
                      fontsize=11, fontweight='bold')
        ax1.set_xlabel('WHO Grade')
        ax1.set_ylabel('Percentage (%)')
        ax1.legend(title='APOE Genotype', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=0)
    
    # 2. Allele frequency by grade
    ax2 = plt.subplot(2, 3, 2)
    if allele_freq is not None and len(allele_freq) > 0:
        allele_data = allele_freq[['ε2', 'ε3', 'ε4']]
        allele_data.plot(kind='bar', ax=ax2)
        ax2.set_title('APOE Allele Frequencies by WHO Grade\n(Reclassified)', 
                      fontsize=11, fontweight='bold')
        ax2.set_xlabel('WHO Grade')
        ax2.set_ylabel('Allele Frequency (%)')
        ax2.legend(title='Allele', fontsize=8)
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=0)
    
    # 3. Sample size by grade
    ax3 = plt.subplot(2, 3, 3)
    df_valid = df[df['Grade_Reclassified'].isin(['WHO II', 'WHO III', 'WHO IV'])]
    if len(df_valid) > 0:
        grade_counts = df_valid['Grade_Reclassified'].value_counts().sort_index()
        colors = ['green' if g == 'WHO II' else 'orange' if g == 'WHO III' else 'red' 
                  for g in grade_counts.index]
        grade_counts.plot(kind='bar', ax=ax3, color=colors)
        ax3.set_title('Sample Distribution\n(All IDH-WT → WHO IV)', 
                      fontsize=11, fontweight='bold')
        ax3.set_xlabel('WHO Grade (Reclassified)')
        ax3.set_ylabel('Number of Patients')
        ax3.set_xticklabels(ax3.get_xticklabels(), rotation=0)
        
        for i, v in enumerate(grade_counts.values):
            ax3.text(i, v + 0.5, str(v), ha='center', va='bottom', fontweight='bold', fontsize=12)
    
    # 4. IDH status by grade
    ax4 = plt.subplot(2, 3, 4)
    if len(df_valid) > 0:
        idh_by_grade = pd.crosstab(df_valid['Grade_Reclassified'], df_valid['IDH_Status'])
        idh_by_grade.plot(kind='bar', stacked=True, ax=ax4, color=['lightcoral', 'lightblue'])
        ax4.set_title('IDH Status by WHO Grade\n(After Reclassification)', 
                      fontsize=11, fontweight='bold')
        ax4.set_xlabel('WHO Grade')
        ax4.set_ylabel('Number of Patients')
        ax4.legend(title='IDH Status', fontsize=8)
        ax4.set_xticklabels(ax4.get_xticklabels(), rotation=0)
    
    # 5. ε4 carrier status by grade
    ax5 = plt.subplot(2, 3, 5)
    if len(df_valid) > 0:
        df_valid_copy = df_valid.copy()
        df_valid_copy['ε4_carrier'] = df_valid_copy['APOE_Genotype'].str.contains('ε4', na=False)
        e4_by_grade = pd.crosstab(df_valid_copy['Grade_Reclassified'], df_valid_copy['ε4_carrier'])
        
        e4_by_grade.plot(kind='bar', ax=ax5, color=['lightblue', 'darkred'], stacked=False)
        ax5.set_title('APOE ε4 Carrier Status by WHO Grade\n(Reclassified)', 
                      fontsize=11, fontweight='bold')
        ax5.set_xlabel('WHO Grade')
        ax5.set_ylabel('Number of Patients')
        ax5.legend(title='ε4 Carrier', labels=['Non-carrier', 'ε4 Carrier'], fontsize=8)
        ax5.set_xticklabels(ax5.get_xticklabels(), rotation=0)
    
    # 6. Heatmap of genotype counts
    ax6 = plt.subplot(2, 3, 6)
    if len(df_valid) > 0:
        contingency_for_heatmap = pd.crosstab(
            df_valid['Grade_Reclassified'], 
            df_valid['APOE_Genotype']
        )
        if len(contingency_for_heatmap) > 0:
            sns.heatmap(contingency_for_heatmap, annot=True, fmt='d', cmap='YlOrRd', 
                       ax=ax6, cbar_kws={'label': 'Patient Count'})
            ax6.set_title('APOE Genotype Counts Heatmap\n(Reclassified)', 
                         fontsize=11, fontweight='bold')
            ax6.set_xlabel('APOE Genotype')
            ax6.set_ylabel('WHO Grade (Reclassified)')
    
    plt.suptitle('APOE Analysis with Reclassified Grades: ALL IDH-Wildtype → WHO Grade IV\n25 Genotyped Patients', 
                 fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    return fig

def generate_summary_report(df, contingency, contingency_pct, allele_freq, chi2, p_value):
    """Generate text summary report"""
    
    report = []
    report.append("=" * 100)
    report.append("APOE GENOTYPE ANALYSIS - RECLASSIFIED GRADES")
    report.append("ALL IDH-WILDTYPE PATIENTS → WHO GRADE IV (GBM)")
    report.append("=" * 100)
    report.append("")
    
    report.append("RECLASSIFICATION RATIONALE:")
    report.append("  IDH-wildtype gliomas are biologically aggressive regardless of histological grade.")
    report.append("  All IDH-WT patients reclassified to WHO Grade IV for clinical accuracy.")
    report.append("")
    
    # Overall summary
    report.append("COHORT OVERVIEW:")
    report.append(f"Total Genotyped Patients: {len(df)}")
    df_valid = df[df['Grade_Reclassified'].isin(['WHO II', 'WHO III', 'WHO IV'])]
    report.append(f"Patients with Valid WHO Grade: {len(df_valid)}")
    report.append("")
    
    # IDH status
    report.append("IDH STATUS DISTRIBUTION:")
    for status in ['Wildtype', 'Mutant']:
        count = (df['IDH_Status'] == status).sum()
        pct = count / len(df) * 100
        report.append(f"  {status}: {count} ({pct:.1f}%)")
    report.append("")
    
    # Overall genotype distribution
    report.append("OVERALL GENOTYPE DISTRIBUTION (All 25 patients):")
    overall_geno = df['APOE_Genotype'].value_counts().sort_index()
    for geno, count in overall_geno.items():
        pct = count / len(df) * 100
        report.append(f"  {geno}: {count} ({pct:.1f}%)")
    report.append("")
    
    if len(df_valid) > 0:
        # Grade distribution
        report.append("RECLASSIFIED GRADE DISTRIBUTION:")
        for grade in ['WHO II', 'WHO III', 'WHO IV']:
            count = (df_valid['Grade_Reclassified'] == grade).sum()
            if len(df_valid) > 0:
                pct = count / len(df_valid) * 100
                # Show IDH breakdown
                idh_wt = ((df_valid['Grade_Reclassified'] == grade) & (df_valid['IDH_Status'] == 'Wildtype')).sum()
                idh_mut = ((df_valid['Grade_Reclassified'] == grade) & (df_valid['IDH_Status'] == 'Mutant')).sum()
                report.append(f"  {grade}: {count} ({pct:.1f}%) - IDH-WT: {idh_wt}, IDH-mut: {idh_mut}")
        report.append("")
        
        # Genotype counts by grade
        if contingency is not None:
            report.append("GENOTYPE DISTRIBUTION BY RECLASSIFIED GRADE:")
            report.append("")
            report.append(contingency.to_string())
            report.append("")
        
        # Genotype percentages by grade
        if contingency_pct is not None and len(contingency_pct) > 0:
            report.append("GENOTYPE PERCENTAGES BY RECLASSIFIED GRADE:")
            report.append("")
            contingency_pct_formatted = contingency_pct.round(1).astype(str) + '%'
            report.append(contingency_pct_formatted.to_string())
            report.append("")
        
        # Allele frequencies
        if allele_freq is not None and len(allele_freq) > 0:
            report.append("ALLELE FREQUENCIES BY RECLASSIFIED GRADE:")
            report.append("")
            allele_freq_formatted = allele_freq[['ε2', 'ε3', 'ε4']].round(1).astype(str) + '%'
            allele_freq_formatted['N_patients'] = allele_freq['N_patients'].astype(int)
            report.append(allele_freq_formatted.to_string())
            report.append("")
        
        # Statistical test
        if chi2 is not None and p_value is not None:
            report.append("STATISTICAL ANALYSIS:")
            report.append(f"Chi-square test for independence:")
            report.append(f"  Chi-square statistic: {chi2:.2f}")
            report.append(f"  P-value: {p_value:.4f}")
            if p_value < 0.05:
                report.append("  Result: Significant association between APOE genotype and WHO grade")
            else:
                report.append("  Result: No significant association detected")
            report.append("")
        
        # Key findings
        report.append("KEY FINDINGS:")
        if contingency_pct is not None and len(contingency_pct) > 0:
            for grade in ['WHO II', 'WHO III', 'WHO IV']:
                if grade in contingency_pct.index:
                    grade_data = contingency_pct.loc[grade]
                    if len(grade_data) > 0:
                        most_common = grade_data.idxmax()
                        freq = grade_data.max()
                        n_patients = (df_valid['Grade_Reclassified'] == grade).sum()
                        report.append(f"  {grade} (n={n_patients}): Most common is {most_common} ({freq:.1f}%)")
        
        report.append("")
        
        # ε4 carrier analysis
        report.append("APOE ε4 CARRIER ANALYSIS BY RECLASSIFIED GRADE:")
        for grade in ['WHO II', 'WHO III', 'WHO IV']:
            grade_df = df_valid[df_valid['Grade_Reclassified'] == grade]
            e4_carriers = grade_df['APOE_Genotype'].str.contains('ε4', na=False).sum()
            total = len(grade_df)
            if total > 0:
                pct = e4_carriers / total * 100
                report.append(f"  {grade}: {e4_carriers}/{total} ({pct:.1f}%) are ε4 carriers")
        
        report.append("")
    
    report.append("=" * 100)
    report.append("✅ This analysis uses RECLASSIFIED data")
    report.append("✅ ALL IDH-wildtype patients classified as WHO Grade IV (GBM)")
    report.append("✅ Reflects biological behavior, not just histology")
    report.append("=" * 100)
    
    return "\n".join(report)

def main():
    """Execute complete analysis with reclassified grades"""
    
    print("=" * 100)
    print("RERUNNING ALL ANALYSES WITH RECLASSIFIED GRADES")
    print("=" * 100)
    print()
    
    # Load reclassified data
    df = load_reclassified_data()
    
    # Analyze genotype by grade
    print("Analyzing genotype distribution by reclassified WHO grade...")
    contingency, contingency_pct, chi2, p_value = analyze_genotype_by_grade(df)
    
    # Analyze allele frequencies
    print("Calculating allele frequencies by reclassified grade...")
    allele_freq = analyze_allele_frequency_by_grade(df)
    
    # Create visualizations
    print("Generating visualization panel...")
    fig = create_visualization_panel(df, contingency_pct, allele_freq)
    plt.savefig('clinical_analysis_outputs/apoe_grade_distribution_RECLASSIFIED.png', 
                dpi=300, bbox_inches='tight')
    print("✅ Saved: clinical_analysis_outputs/apoe_grade_distribution_RECLASSIFIED.png")
    print()
    
    # Generate text report
    print("Generating summary report...")
    report = generate_summary_report(df, contingency, contingency_pct, allele_freq, chi2, p_value)
    
    with open('clinical_analysis_outputs/apoe_grade_distribution_RECLASSIFIED_report.txt', 'w') as f:
        f.write(report)
    print("✅ Saved: clinical_analysis_outputs/apoe_grade_distribution_RECLASSIFIED_report.txt")
    print()
    
    # Save detailed data tables
    print("Saving detailed data tables...")
    with pd.ExcelWriter('clinical_analysis_outputs/apoe_grade_RECLASSIFIED_tables.xlsx') as writer:
        df.to_excel(writer, sheet_name='Full_Patient_Data', index=False)
        
        if contingency is not None:
            contingency.to_excel(writer, sheet_name='Genotype_Counts')
        if contingency_pct is not None and len(contingency_pct) > 0:
            contingency_pct.to_excel(writer, sheet_name='Genotype_Percentages')
        if allele_freq is not None and len(allele_freq) > 0:
            allele_freq.to_excel(writer, sheet_name='Allele_Frequencies')
        
        summary_df = pd.DataFrame({
            'Metric': ['Total Patients', 'WHO II', 'WHO III', 'WHO IV (GBM)', 
                      'IDH-Wildtype', 'IDH-Mutant', 'Chi-square', 'P-value'],
            'Value': [
                len(df),
                (df['Grade_Reclassified'] == 'WHO II').sum(),
                (df['Grade_Reclassified'] == 'WHO III').sum(),
                (df['Grade_Reclassified'] == 'WHO IV').sum(),
                (df['IDH_Status'] == 'Wildtype').sum(),
                (df['IDH_Status'] == 'Mutant').sum(),
                chi2 if chi2 is not None else 'N/A',
                p_value if p_value is not None else 'N/A'
            ]
        })
        summary_df.to_excel(writer, sheet_name='Summary', index=False)
    
    print("✅ Saved: clinical_analysis_outputs/apoe_grade_RECLASSIFIED_tables.xlsx")
    print()
    
    # Print summary
    print("=" * 100)
    print("ANALYSIS COMPLETE - RECLASSIFIED DATA")
    print("=" * 100)
    print(f"Total patients: {len(df)}")
    print()
    print("Reclassified Grade Distribution:")
    print(f"  WHO II:  {(df['Grade_Reclassified'] == 'WHO II').sum()} (all IDH-mutant)")
    print(f"  WHO III: {(df['Grade_Reclassified'] == 'WHO III').sum()} (all IDH-mutant)")
    print(f"  WHO IV:  {(df['Grade_Reclassified'] == 'WHO IV').sum()} (20 IDH-WT + 3 IDH-mutant)")
    print()
    if chi2 is not None and p_value is not None:
        print(f"Chi-square p-value: {p_value:.4f}")
    print()
    print("✅ All analyses regenerated with reclassified grades!")
    print("=" * 100)
    
    return df, contingency, allele_freq

if __name__ == "__main__":
    df, contingency, allele_freq = main()

