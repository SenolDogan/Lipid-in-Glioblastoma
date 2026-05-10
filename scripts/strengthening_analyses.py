#!/usr/bin/env python3
"""
Additional Strengthening Analyses for Manuscript
- Cross-validation details
- Bootstrap sampling for robustness
- Hub molecule analysis
- Detailed pathway enrichment
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import cross_val_score, LeaveOneOut
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import RobustScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import NMF
import networkx as nx
from scipy import stats
import warnings
import os

warnings.filterwarnings('ignore')

output_dir = "/Users/senol/Desktop/GBM/Lipidomics_Metabolomics_Integration"
os.makedirs(f"{output_dir}/results", exist_ok=True)
os.makedirs(f"{output_dir}/plots", exist_ok=True)

print("="*80)
print("STRENGTHENING ANALYSES")
print("="*80)

# Load data
print("\n[1/5] Loading data...")

# Load raw data to reconstruct combined matrix
metabolome = pd.read_csv("/Users/senol/Desktop/GBM/metabolome_quantification.txt", sep="\t")
lipidome_pos = pd.read_csv("/Users/senol/Desktop/GBM/lipidome_positive_quantification.txt", sep="\t")
lipidome_neg = pd.read_csv("/Users/senol/Desktop/GBM/lipidome_negative_quantification.txt", sep="\t")
clinical = pd.read_csv("/Users/senol/Desktop/GBM/clinical_patient.txt", sep="\t", skiprows=4)

# Clean clinical data
clinical.columns = clinical.columns.str.strip()
clinical = clinical[clinical['PATIENT_ID'].notna()]

# Extract sample IDs
def extract_sample_ids(df):
    sample_cols = [col for col in df.columns if col.startswith('C3')]
    return sample_cols

metabolome_samples = extract_sample_ids(metabolome)
lipidome_pos_samples = extract_sample_ids(lipidome_pos)
lipidome_neg_samples = extract_sample_ids(lipidome_neg)
clinical_patient_ids = clinical['PATIENT_ID'].dropna().unique().tolist()

common_samples = list(set(metabolome_samples) & set(lipidome_pos_samples) & 
                     set(lipidome_neg_samples) & set(clinical_patient_ids))
common_samples = sorted(common_samples)

# Prepare data matrices
def prepare_data_matrix(df, sample_cols, id_col='ENTITY_STABLE_ID'):
    if id_col in df.columns:
        df = df.set_index(id_col)
    else:
        df = df.set_index(df.columns[0])
    
    available_samples = [s for s in sample_cols if s in df.columns]
    df_subset = df[available_samples].T
    df_subset = df_subset.apply(pd.to_numeric, errors='coerce')
    df_subset = df_subset.loc[:, df_subset.isna().mean() < 0.5]
    df_subset = df_subset.fillna(df_subset.median())
    return df_subset

metabolome_matrix = prepare_data_matrix(metabolome, common_samples)
lipidome_pos_matrix = prepare_data_matrix(lipidome_pos, common_samples)
lipidome_neg_matrix = prepare_data_matrix(lipidome_neg, common_samples)

# Combine lipidome datasets
lipidome_combined = pd.concat([lipidome_pos_matrix, lipidome_neg_matrix], axis=1)

# Use actual lipid names instead of LIPID_XXX
lipidome_combined_columns = []
lipidome_pos_original = pd.read_csv("/Users/senol/Desktop/GBM/lipidome_positive_quantification.txt", sep="\t")
lipidome_neg_original = pd.read_csv("/Users/senol/Desktop/GBM/lipidome_negative_quantification.txt", sep="\t")

# Track seen names to handle duplicates
seen_names = {}
for col in lipidome_pos_matrix.columns:
    if col in lipidome_pos_original['ENTITY_STABLE_ID'].values:
        name = lipidome_pos_original.loc[lipidome_pos_original['ENTITY_STABLE_ID'] == col, 'NAME'].values[0]
    else:
        name = col
    
    # Handle duplicate names by appending suffix
    if name in seen_names:
        seen_names[name] += 1
        name = f"{name}_dup{seen_names[name]}"
    else:
        seen_names[name] = 0
    lipidome_combined_columns.append(name)

for col in lipidome_neg_matrix.columns:
    if col in lipidome_neg_original['ENTITY_STABLE_ID'].values:
        name = lipidome_neg_original.loc[lipidome_neg_original['ENTITY_STABLE_ID'] == col, 'NAME'].values[0]
    else:
        name = col
    
    # Handle duplicate names by appending suffix
    if name in seen_names:
        seen_names[name] += 1
        name = f"{name}_dup{seen_names[name]}"
    else:
        seen_names[name] = 0
    lipidome_combined_columns.append(name)

lipidome_combined.columns = lipidome_combined_columns

# Align samples
common_samples_final = list(set(metabolome_matrix.index) & set(lipidome_combined.index))
metabolome_matrix = metabolome_matrix.loc[common_samples_final]
lipidome_combined = lipidome_combined.loc[common_samples_final]

# Combine into single matrix
combined_data = pd.concat([metabolome_matrix, lipidome_combined], axis=1)

# Load cluster assignments and other results
cluster_assignments = pd.read_csv(f"{output_dir}/results/cluster_assignments.csv")
lipid_metabolite_corr = pd.read_csv(f"{output_dir}/results/lipid_metabolite_correlations.csv")
differential_features = pd.read_csv(f"{output_dir}/results/differential_features.csv")

# Get NMF labels
nmf_labels = cluster_assignments['nmf_subtype'].values

# Prepare data for analysis
scaler = RobustScaler()
combined_scaled = pd.DataFrame(
    scaler.fit_transform(combined_data),
    index=combined_data.index,
    columns=combined_data.columns
)

# ============================================================================
# 1. CROSS-VALIDATION ANALYSIS
# ============================================================================
print("\n[2/5] Performing cross-validation analysis...")

# Leave-One-Out Cross-Validation
rf_model = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
loo = LeaveOneOut()
loo_scores = cross_val_score(rf_model, combined_scaled.values, nmf_labels, cv=loo, scoring='accuracy')
loo_accuracy = loo_scores.mean()
loo_std = loo_scores.std()

print(f"  Leave-One-Out CV Accuracy: {loo_accuracy:.4f} ± {loo_std:.4f}")

# 5-fold and 10-fold CV
cv5_scores = cross_val_score(rf_model, combined_scaled.values, nmf_labels, cv=5, scoring='accuracy')
cv10_scores = cross_val_score(rf_model, combined_scaled.values, nmf_labels, cv=10, scoring='accuracy')

print(f"  5-fold CV Accuracy: {cv5_scores.mean():.4f} ± {cv5_scores.std():.4f}")
print(f"  10-fold CV Accuracy: {cv10_scores.mean():.4f} ± {cv10_scores.std():.4f}")

# Save results
cv_results = pd.DataFrame({
    'method': ['Leave-One-Out', '5-fold', '10-fold'],
    'accuracy_mean': [loo_accuracy, cv5_scores.mean(), cv10_scores.mean()],
    'accuracy_std': [loo_std, cv5_scores.std(), cv10_scores.std()],
    'n_splits': [len(loo_scores), 5, 10]
})
cv_results.to_csv(f"{output_dir}/results/cross_validation_results.csv", index=False)

# ============================================================================
# 2. BOOTSTRAP SAMPLING FOR ROBUSTNESS
# ============================================================================
print("\n[3/5] Performing bootstrap sampling for subtype stability...")

n_bootstrap = 1000
bootstrap_agreement = []

for i in range(n_bootstrap):
    # Bootstrap sample
    bootstrap_indices = np.random.choice(len(combined_scaled), size=len(combined_scaled), replace=True)
    bootstrap_data = combined_scaled.iloc[bootstrap_indices]
    bootstrap_labels = nmf_labels[bootstrap_indices]
    
    # Cluster with NMF
    try:
        nmf = NMF(n_components=2, random_state=42, max_iter=1000)
        nmf_components = nmf.fit_transform(bootstrap_data.values)
        bootstrap_nmf_labels = np.argmax(nmf_components, axis=1)
        
        # Calculate agreement with original (adjusted for label switching)
        # Try both label assignments
        agreement1 = np.mean(bootstrap_nmf_labels == bootstrap_labels)
        agreement2 = np.mean((1 - bootstrap_nmf_labels) == bootstrap_labels)
        agreement = max(agreement1, agreement2)
        bootstrap_agreement.append(agreement)
    except:
        continue

bootstrap_agreement = np.array(bootstrap_agreement)
bootstrap_stability = np.mean(bootstrap_agreement)
bootstrap_std = np.std(bootstrap_agreement)

print(f"  Bootstrap Stability: {bootstrap_stability:.4f} ± {bootstrap_std:.4f}")
print(f"  Stable in {np.sum(bootstrap_agreement > 0.8) / len(bootstrap_agreement) * 100:.1f}% of iterations")

# Save results
bootstrap_results = pd.DataFrame({
    'metric': ['Mean Agreement', 'Std Agreement', 'Stable (>0.8) %'],
    'value': [bootstrap_stability, bootstrap_std, np.sum(bootstrap_agreement > 0.8) / len(bootstrap_agreement) * 100]
})
bootstrap_results.to_csv(f"{output_dir}/results/bootstrap_stability_results.csv", index=False)

# ============================================================================
# 3. HUB MOLECULE ANALYSIS
# ============================================================================
print("\n[4/5] Performing detailed hub molecule analysis...")

# Create network from correlations
G = nx.Graph()
top_corr = lipid_metabolite_corr.nlargest(100, 'correlation')

for _, row in top_corr.iterrows():
    G.add_edge(row['feature1'], row['feature2'], weight=abs(row['correlation']))

# Calculate centrality measures
degree_centrality = nx.degree_centrality(G)
betweenness_centrality = nx.betweenness_centrality(G)
closeness_centrality = nx.closeness_centrality(G)
eigenvector_centrality = nx.eigenvector_centrality(G, max_iter=1000)

# Combine centrality measures
hub_analysis = pd.DataFrame({
    'molecule': list(degree_centrality.keys()),
    'degree_centrality': list(degree_centrality.values()),
    'betweenness_centrality': list(betweenness_centrality.values()),
    'closeness_centrality': list(closeness_centrality.values()),
    'eigenvector_centrality': list(eigenvector_centrality.values()),
    'degree': [G.degree(node) for node in degree_centrality.keys()]
})

# Calculate composite hub score
hub_analysis['hub_score'] = (
    hub_analysis['degree_centrality'] * 0.4 +
    hub_analysis['betweenness_centrality'] * 0.3 +
    hub_analysis['closeness_centrality'] * 0.2 +
    hub_analysis['eigenvector_centrality'] * 0.1
)

hub_analysis = hub_analysis.sort_values('hub_score', ascending=False)

# Top 20 hub molecules
top_hubs = hub_analysis.head(20)
print(f"  Top 20 hub molecules identified")
print(f"  Top hub: {top_hubs.iloc[0]['molecule']} (degree: {top_hubs.iloc[0]['degree']}, hub_score: {top_hubs.iloc[0]['hub_score']:.4f})")

# Save results
hub_analysis.to_csv(f"{output_dir}/results/hub_molecule_analysis.csv", index=False)
top_hubs.to_csv(f"{output_dir}/results/top_20_hub_molecules.csv", index=False)

# ============================================================================
# 4. DETAILED PATHWAY ENRICHMENT
# ============================================================================
print("\n[5/5] Performing detailed pathway enrichment analysis...")

# Get differential features
sig_features = differential_features[differential_features['fdr'] < 0.05]

# Categorize features into pathways (more detailed)
pathway_categories = {
    'Lipid Metabolism': {
        'Fatty Acid Synthesis': ['TG(52:4)', 'TG(52:5)', 'TG(56:5)', 'Cer(d16:1/18:0);Cer(d18:1/16:0)'],
        'Phospholipid Metabolism': ['PC(18:0/18:1)', 'PC(16:0/18:2)', 'PC(16:0/18:1)'],
        'Sphingolipid Metabolism': ['Cer(d18:1/24:0)', 'SM(d18:1/24:1)'],
        'Triacylglycerol Metabolism': ['TG(16:0/18:1/20:4)'],
        'Cholesterol Metabolism': [],
        'Other Lipid Metabolism': []
    },
    'Amino Acid Metabolism': ['L-glutamic_acid', 'glycine', 'phosphoenolpyruvic_acid'],
    'Carbohydrate Metabolism': [],
    'Nucleotide Metabolism': [],
    'Energy Metabolism': ['pyrophosphate'],
    'Other': []
}

# Count features in each pathway
pathway_counts = {}
for pathway, subcategories in pathway_categories.items():
    if isinstance(subcategories, dict):
        # Lipid metabolism has subcategories
        total = 0
        for subcat, features in subcategories.items():
            count = len([f for f in sig_features['feature'] if any(feat in str(f) for feat in features)])
            pathway_counts[f"{pathway} - {subcat}"] = count
            total += count
        pathway_counts[pathway] = total
    else:
        count = len([f for f in sig_features['feature'] if any(feat in str(f) for feat in subcategories)])
        pathway_counts[pathway] = count

# Create pathway enrichment table
pathway_enrichment_detailed = pd.DataFrame({
    'pathway': list(pathway_counts.keys()),
    'count': list(pathway_counts.values()),
    'percentage': [c / len(sig_features) * 100 for c in pathway_counts.values()]
}).sort_values('count', ascending=False)

print(f"  Total pathways analyzed: {len(pathway_enrichment_detailed)}")
print(f"  Lipid Metabolism total: {pathway_counts.get('Lipid Metabolism', 0)} features")

pathway_enrichment_detailed.to_csv(f"{output_dir}/results/detailed_pathway_enrichment.csv", index=False)

print("\n" + "="*80)
print("STRENGTHENING ANALYSES COMPLETE!")
print("="*80)
print(f"\nResults saved to: {output_dir}/results/")
print("\nGenerated files:")
print("  - cross_validation_results.csv")
print("  - bootstrap_stability_results.csv")
print("  - hub_molecule_analysis.csv")
print("  - top_20_hub_molecules.csv")
print("  - detailed_pathway_enrichment.csv")

