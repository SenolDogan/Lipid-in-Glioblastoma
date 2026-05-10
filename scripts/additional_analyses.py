#!/usr/bin/env python3
"""
Additional Advanced Analyses for Lipidomics-Metabolomics Integration
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score
import warnings
import os

warnings.filterwarnings('ignore')

output_dir = "/Users/senol/Desktop/GBM/Lipidomics_Metabolomics_Integration"

# Figure typography: Times New Roman and larger, manuscript-ready
plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "font.size": 16,
        "axes.titlesize": 22,
        "axes.labelsize": 18,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 16,
        "legend.title_fontsize": 16,
    }
)
sns.set_theme(style="whitegrid", font="Times New Roman")

print("="*80)
print("ADDITIONAL ADVANCED ANALYSES")
print("="*80)

# Load previous results
try:
    metabolome = pd.read_csv("/Users/senol/Desktop/GBM/metabolome_quantification.txt", sep="\t")
    lipidome_pos = pd.read_csv("/Users/senol/Desktop/GBM/lipidome_positive_quantification.txt", sep="\t")
    lipidome_neg = pd.read_csv("/Users/senol/Desktop/GBM/lipidome_negative_quantification.txt", sep="\t")
    clinical = pd.read_csv("/Users/senol/Desktop/GBM/clinical_patient.txt", sep="\t", skiprows=4)
    cluster_assignments = pd.read_csv(f"{output_dir}/results/cluster_assignments.csv")
    
    print("\n[1/6] Lipid Class Analysis...")
    
    # Extract lipid classes from names
    def categorize_lipid(lipid_name):
        """Categorize lipid by class"""
        name_lower = str(lipid_name).lower()
        
        if 'cer' in name_lower or 'ceramide' in name_lower:
            return 'Ceramide'
        elif 'ce' in name_lower and 'carnitine' not in name_lower:
            return 'Cholesterol Ester'
        elif 'carnitine' in name_lower:
            return 'Carnitine'
        elif 'pc' in name_lower or 'phosphatidylcholine' in name_lower:
            return 'Phosphatidylcholine'
        elif 'pe' in name_lower or 'phosphatidylethanolamine' in name_lower:
            return 'Phosphatidylethanolamine'
        elif 'tg' in name_lower or 'triacylglycerol' in name_lower:
            return 'Triacylglycerol'
        elif 'dg' in name_lower or 'diacylglycerol' in name_lower:
            return 'Diacylglycerol'
        elif 'sm' in name_lower or 'sphingomyelin' in name_lower:
            return 'Sphingomyelin'
        elif 'anandamide' in name_lower:
            return 'Anandamide'
        else:
            return 'Other'
    
    # Get sample columns and ensure proper Patient ID matching
    lipid_samples = [col for col in lipidome_pos.columns if col.startswith('C3')]
    metabolome_samples = [col for col in metabolome.columns if col.startswith('C3')]
    clinical_patient_ids = clinical['PATIENT_ID'].dropna().unique().tolist()
    
    # Find common samples across all datasets
    common_samples = list(set(lipid_samples) & set(metabolome_samples) & 
                         set(clinical_patient_ids) & set(cluster_assignments['sample_id']))
    common_samples = sorted(common_samples)  # Sort for consistency
    print(f"  Using {len(common_samples)} common patients across all datasets")
    
    # Prepare lipid data
    lipidome_pos_clean = lipidome_pos.set_index('ENTITY_STABLE_ID')
    lipidome_neg_clean = lipidome_neg.set_index('ENTITY_STABLE_ID')
    
    lipid_data = pd.concat([
        lipidome_pos_clean[common_samples].T,
        lipidome_neg_clean[common_samples].T
    ], axis=1)
    lipid_data = lipid_data.apply(pd.to_numeric, errors='coerce').fillna(lipid_data.median())
    
    # Categorize lipids
    lipid_classes = {}
    for col in lipid_data.columns:
        # Try to get name
        if col in lipidome_pos_clean.index:
            name = lipidome_pos.loc[lipidome_pos['ENTITY_STABLE_ID'] == col, 'NAME'].values[0] if col in lipidome_pos['ENTITY_STABLE_ID'].values else col
        elif col in lipidome_neg_clean.index:
            name = lipidome_neg.loc[lipidome_neg['ENTITY_STABLE_ID'] == col, 'NAME'].values[0] if col in lipidome_neg['ENTITY_STABLE_ID'].values else col
        else:
            name = col
        
        lipid_class = categorize_lipid(name)
        if lipid_class not in lipid_classes:
            lipid_classes[lipid_class] = []
        lipid_classes[lipid_class].append(col)
    
    # Calculate class abundances
    lipid_class_abundance = pd.DataFrame(index=common_samples)
    for lipid_class, lipids in lipid_classes.items():
        if len(lipids) > 0:
            available_lipids = [l for l in lipids if l in lipid_data.columns]
            if len(available_lipids) > 0:
                lipid_class_abundance[lipid_class] = lipid_data[available_lipids].mean(axis=1)
    
    # Merge with cluster assignments
    cluster_df = cluster_assignments.set_index('sample_id')
    lipid_class_abundance = lipid_class_abundance.merge(
        cluster_df[['nmf_subtype']], left_index=True, right_index=True, how='inner'
    )
    # Use manuscript subtype labels (A/B) instead of 0/1
    subtype_label_map = {0: "A", 1: "B", "0": "A", "1": "B"}
    lipid_class_abundance["Subtype"] = lipid_class_abundance["nmf_subtype"].map(subtype_label_map).fillna(
        lipid_class_abundance["nmf_subtype"].astype(str)
    )
    
    # Plot lipid class distribution by subtype
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Box plot
    lipid_class_melted = lipid_class_abundance.melt(
        id_vars=['Subtype'],
        value_vars=[c for c in lipid_class_abundance.columns if c not in ('nmf_subtype', 'Subtype')],
        var_name='Lipid Class',
        value_name='Abundance'
    )
    
    sns.boxplot(data=lipid_class_melted, x='Lipid Class', y='Abundance', 
               hue='Subtype', ax=axes[0])
    axes[0].set_title('Lipid Class Abundance by Subtype', fontweight='bold')
    axes[0].set_xlabel('Lipid Class')
    axes[0].set_ylabel('Abundance')
    axes[0].tick_params(axis='x', rotation=45)
    axes[0].legend(title='Subtype')
    
    # Stacked bar plot
    lipid_class_avg = lipid_class_abundance.groupby('Subtype')[
        [c for c in lipid_class_abundance.columns if c not in ('nmf_subtype', 'Subtype')]
    ].mean()
    
    lipid_class_avg.plot(kind='bar', stacked=True, ax=axes[1], colormap='tab20')
    axes[1].set_title('Average Lipid Class Composition by Subtype', fontweight='bold')
    axes[1].set_xlabel('Subtype')
    axes[1].set_ylabel('Relative Abundance')
    axes[1].legend(title='Lipid Class', bbox_to_anchor=(1.05, 1), loc='upper left')
    axes[1].tick_params(axis='x', rotation=0)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/plots/lipid_class_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save results
    lipid_class_abundance.to_csv(f"{output_dir}/results/lipid_class_abundance.csv")
    
    print(f"  Analyzed {len(lipid_classes)} lipid classes")
    
    print("\n[2/6] Metabolite Category Analysis...")
    
    # Categorize metabolites
    def categorize_metabolite(met_name):
        """Categorize metabolite by class"""
        name_lower = str(met_name).lower()
        
        if any(kw in name_lower for kw in ['amino', 'acid', 'alanine', 'glutamine', 'glycine', 'serine', 'threonine', 'valine', 'leucine', 'isoleucine', 'phenylalanine', 'tyrosine', 'tryptophan', 'histidine', 'lysine', 'arginine', 'aspartate', 'asparagine', 'glutamate', 'proline', 'methionine', 'cysteine']):
            return 'Amino Acid'
        elif any(kw in name_lower for kw in ['glucose', 'fructose', 'galactose', 'sugar', 'glycol', 'pentose', 'ribose', 'lactate', 'pyruvate']):
            return 'Carbohydrate/Glycolysis'
        elif any(kw in name_lower for kw in ['fatty', 'acyl', 'butyric', 'palmitic', 'oleic', 'linoleic']):
            return 'Fatty Acid'
        elif any(kw in name_lower for kw in ['nucleotide', 'purine', 'pyrimidine', 'adenine', 'uracil', 'guanine', 'cytosine', 'thymine']):
            return 'Nucleotide'
        elif any(kw in name_lower for kw in ['atp', 'adp', 'amp', 'nad', 'fad', 'coenzyme', 'phosphate']):
            return 'Energy/Coenzyme'
        elif any(kw in name_lower for kw in ['cholesterol', 'sterol', 'steroid']):
            return 'Sterol'
        else:
            return 'Other'
    
    # Prepare metabolite data
    metabolome_clean = metabolome.set_index('ENTITY_STABLE_ID')
    met_samples = [col for col in metabolome.columns if col.startswith('C3')]
    common_met_samples = list(set(met_samples) & set(cluster_assignments['sample_id']))
    
    met_data = metabolome_clean[common_met_samples].T
    met_data = met_data.apply(pd.to_numeric, errors='coerce').fillna(met_data.median())
    
    # Categorize metabolites
    met_categories = {}
    for col in met_data.columns:
        if col in metabolome['ENTITY_STABLE_ID'].values:
            name = metabolome.loc[metabolome['ENTITY_STABLE_ID'] == col, 'NAME'].values[0] if len(metabolome.loc[metabolome['ENTITY_STABLE_ID'] == col]) > 0 else col
        else:
            name = col
        
        met_category = categorize_metabolite(name)
        if met_category not in met_categories:
            met_categories[met_category] = []
        met_categories[met_category].append(col)
    
    # Calculate category abundances
    met_category_abundance = pd.DataFrame(index=common_met_samples)
    for category, mets in met_categories.items():
        if len(mets) > 0:
            available_mets = [m for m in mets if m in met_data.columns]
            if len(available_mets) > 0:
                met_category_abundance[category] = met_data[available_mets].mean(axis=1)
    
    # Merge with clusters
    met_category_abundance = met_category_abundance.merge(
        cluster_df[['nmf_subtype']], left_index=True, right_index=True, how='inner'
    )
    # Use manuscript subtype labels (A/B) instead of 0/1
    met_category_abundance["Subtype"] = met_category_abundance["nmf_subtype"].map(subtype_label_map).fillna(
        met_category_abundance["nmf_subtype"].astype(str)
    )
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    met_category_melted = met_category_abundance.melt(
        id_vars=['Subtype'],
        value_vars=[c for c in met_category_abundance.columns if c not in ('nmf_subtype', 'Subtype')],
        var_name='Metabolite Category',
        value_name='Abundance'
    )
    
    sns.boxplot(data=met_category_melted, x='Metabolite Category', y='Abundance',
               hue='Subtype', ax=axes[0])
    axes[0].set_title('Metabolite Category Abundance by Subtype', fontweight='bold')
    axes[0].set_xlabel('Metabolite Category')
    axes[0].set_ylabel('Abundance')
    axes[0].tick_params(axis='x', rotation=45)
    axes[0].legend(title='Subtype')
    
    met_category_avg = met_category_abundance.groupby('Subtype')[
        [c for c in met_category_abundance.columns if c not in ('nmf_subtype', 'Subtype')]
    ].mean()
    
    met_category_avg.plot(kind='bar', stacked=True, ax=axes[1], colormap='Set3')
    axes[1].set_title('Average Metabolite Category Composition by Subtype', fontweight='bold')
    axes[1].set_xlabel('Subtype')
    axes[1].set_ylabel('Relative Abundance')
    axes[1].legend(title='Category', bbox_to_anchor=(1.05, 1), loc='upper left')
    axes[1].tick_params(axis='x', rotation=0)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/plots/metabolite_category_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    met_category_abundance.to_csv(f"{output_dir}/results/metabolite_category_abundance.csv")
    
    print(f"  Analyzed {len(met_categories)} metabolite categories")
    
    print("\n[3/6] Principal Component Analysis...")
    
    # Combine all data
    combined_all = pd.concat([met_data, lipid_data], axis=1)
    combined_all = combined_all.loc[common_samples]
    combined_all = combined_all.fillna(combined_all.median())
    
    # Standardize
    scaler = StandardScaler()
    combined_scaled = pd.DataFrame(
        scaler.fit_transform(combined_all),
        index=combined_all.index,
        columns=combined_all.columns
    )
    
    # PCA
    pca = PCA(n_components=10)
    pca_results = pca.fit_transform(combined_scaled)
    
    # Merge with clusters
    pca_df = pd.DataFrame(pca_results, index=combined_scaled.index, 
                        columns=[f'PC{i+1}' for i in range(10)])
    pca_df = pca_df.merge(cluster_df[['nmf_subtype']], left_index=True, right_index=True, how='inner')
    
    # Plot PCA
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # PC1 vs PC2
    for subtype in sorted(pca_df['nmf_subtype'].unique()):
        subtype_data = pca_df[pca_df['nmf_subtype'] == subtype]
        subtype_label = subtype_label_map.get(subtype, str(subtype))
        axes[0].scatter(subtype_data['PC1'], subtype_data['PC2'], 
                       label=f'Subtype {subtype_label}', s=110, alpha=0.75)
    
    axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontweight='bold')
    axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontweight='bold')
    axes[0].set_title('PCA: PC1 vs PC2 by Subtype', fontweight='bold')
    axes[0].legend(title="Subtype")
    axes[0].grid(True, alpha=0.3)
    
    # Explained variance
    axes[1].bar(range(1, 11), pca.explained_variance_ratio_[:10], color='steelblue', alpha=0.7)
    axes[1].plot(range(1, 11), np.cumsum(pca.explained_variance_ratio_[:10]), 
                'ro-', linewidth=2, markersize=8, label='Cumulative')
    axes[1].set_xlabel('Principal Component', fontweight='bold')
    axes[1].set_ylabel('Explained Variance Ratio', fontweight='bold')
    axes[1].set_title('PCA: Explained Variance', fontweight='bold')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/plots/pca_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save PCA loadings
    pca_loadings = pd.DataFrame(
        pca.components_[:5].T,
        index=combined_scaled.columns,
        columns=[f'PC{i+1}' for i in range(5)]
    )
    pca_loadings.to_csv(f"{output_dir}/results/pca_loadings.csv")
    
    print(f"  First 5 PCs explain {np.cumsum(pca.explained_variance_ratio_[:5])[-1]*100:.1f}% of variance")
    
    print("\n[4/6] Correlation Matrix Analysis...")
    
    # Calculate correlation matrix between lipid classes and metabolite categories
    lipid_class_avg_clean = lipid_class_abundance[
        [c for c in lipid_class_abundance.columns if c not in ('nmf_subtype', 'Subtype')]
    ]
    met_category_avg_clean = met_category_abundance[
        [c for c in met_category_abundance.columns if c not in ('nmf_subtype', 'Subtype')]
    ]
    
    # Align samples
    common_corr_samples = list(set(lipid_class_avg_clean.index) & set(met_category_avg_clean.index))
    lipid_class_avg_clean = lipid_class_avg_clean.loc[common_corr_samples]
    met_category_avg_clean = met_category_avg_clean.loc[common_corr_samples]
    
    # Correlation matrix
    corr_matrix = pd.DataFrame(index=lipid_class_avg_clean.columns, 
                               columns=met_category_avg_clean.columns)
    
    for lipid_class in lipid_class_avg_clean.columns:
        for met_category in met_category_avg_clean.columns:
            corr, pval = stats.pearsonr(
                lipid_class_avg_clean[lipid_class],
                met_category_avg_clean[met_category]
            )
            corr_matrix.loc[lipid_class, met_category] = corr
    
    corr_matrix = corr_matrix.astype(float)
    
    # Plot correlation heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
               square=True, linewidths=0.5, cbar_kws={"shrink": 0.8})
    plt.title('Lipid Class - Metabolite Category Correlation Matrix', 
             fontsize=14, fontweight='bold')
    plt.xlabel('Metabolite Category', fontsize=12)
    plt.ylabel('Lipid Class', fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/plots/lipid_metabolite_category_correlation.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    corr_matrix.to_csv(f"{output_dir}/results/lipid_metabolite_category_correlation.csv")
    
    print("  Correlation matrix computed and saved")
    
    print("\n[5/6] Survival Prediction Models...")
    
    # Load survival data
    clinical_clean = clinical[clinical['PATIENT_ID'].notna()].copy()
    clinical_clean = clinical_clean.set_index('PATIENT_ID')
    
    # Prepare survival data
    survival_data = []
    for sample_id in common_samples:
        if sample_id in clinical_clean.index:
            row = clinical_clean.loc[sample_id]
            if pd.notna(row.get('PATH_DIAG_TO_DEATH_DAYS')):
                survival_data.append({
                    'sample_id': sample_id,
                    'time': row['PATH_DIAG_TO_DEATH_DAYS'],
                    'event': 1
                })
            elif pd.notna(row.get('PATH_DIAG_TO_LAST_CONTACT_DAYS')):
                survival_data.append({
                    'sample_id': sample_id,
                    'time': row['PATH_DIAG_TO_LAST_CONTACT_DAYS'],
                    'event': 0
                })
    
    survival_df = pd.DataFrame(survival_data)
    
    if len(survival_df) > 10:
        # Prepare features
        ml_features = pd.concat([
            lipid_class_avg_clean.loc[survival_df['sample_id']],
            met_category_avg_clean.loc[survival_df['sample_id']]
        ], axis=1)
        ml_features = ml_features.fillna(ml_features.median())
        
        # Train models
        X = ml_features.values
        y_time = survival_df['time'].values
        y_event = survival_df['event'].values
        
        # Random Forest for survival time prediction
        rf = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1)
        rf_scores = cross_val_score(rf, X, y_time, cv=5, scoring='neg_mean_squared_error')
        print(f"  Random Forest RMSE: {np.sqrt(-rf_scores.mean()):.2f} days")
        
        # Feature importance
        rf.fit(X, y_time)
        feature_importance = pd.DataFrame({
            'feature': ml_features.columns,
            'importance': rf.feature_importances_
        }).sort_values('importance', ascending=False)
        
        feature_importance.to_csv(f"{output_dir}/results/survival_feature_importance.csv", index=False)
        
        # Plot feature importance
        plt.figure(figsize=(10, 8))
        top_features = feature_importance.head(20)
        plt.barh(range(len(top_features)), top_features['importance'].values)
        plt.yticks(range(len(top_features)), top_features['feature'].values)
        plt.xlabel('Importance', fontsize=12)
        plt.title('Top 20 Features for Survival Time Prediction', fontsize=14, fontweight='bold')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(f"{output_dir}/plots/survival_feature_importance.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    print("\n[6/6] Summary Statistics...")
    
    # Create summary report
    summary_stats = {
        'Total Samples': len(common_samples),
        'Total Metabolites': len(met_data.columns),
        'Total Lipids': len(lipid_data.columns),
        'Lipid Classes': len(lipid_classes),
        'Metabolite Categories': len(met_categories),
        'Subtypes Identified': len(set(cluster_assignments['nmf_subtype'])),
        'Samples with Survival Data': len(survival_df) if 'survival_df' in locals() else 0
    }
    
    summary_df = pd.Series(summary_stats)
    summary_df.to_csv(f"{output_dir}/results/summary_statistics.csv", header=['Value'])
    
    print("\n" + "="*80)
    print("ADDITIONAL ANALYSES COMPLETE!")
    print("="*80)
    for key, value in summary_stats.items():
        print(f"  {key}: {value}")
    print("="*80)
    
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()

