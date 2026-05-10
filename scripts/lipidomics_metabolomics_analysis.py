#!/usr/bin/env python3
"""
Comprehensive Lipidomics-Metabolomics Integration Analysis for GBM
This script performs all statistical and machine learning analyses
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.decomposition import NMF, PCA
from sklearn.cluster import KMeans, AgglomerativeClustering, DBSCAN
from sklearn.manifold import TSNE
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import LassoCV, ElasticNetCV
from sklearn.svm import SVC
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.metrics import pairwise_distances
from sklearn.metrics import roc_curve, auc, confusion_matrix, classification_report
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import warnings
import os
import re
from itertools import combinations
import networkx as nx
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')

# Set style
try:
    plt.style.use('seaborn-v0_8-darkgrid')
except:
    plt.style.use('seaborn-darkgrid')
sns.set_palette("husl")

# Create output directories
output_dir = "/Users/senol/Desktop/GBM/Lipidomics_Metabolomics_Integration"
os.makedirs(f"{output_dir}/plots", exist_ok=True)
os.makedirs(f"{output_dir}/results", exist_ok=True)
os.makedirs(f"{output_dir}/networks", exist_ok=True)

print("="*80)
print("GBM LIPIDOMICS-METABOLOMICS INTEGRATION ANALYSIS")
print("="*80)

# ============================================================================
# 1. DATA LOADING AND PREPROCESSING
# ============================================================================
print("\n[1/8] Loading and preprocessing data...")

# Load datasets
metabolome = pd.read_csv("/Users/senol/Desktop/GBM/metabolome_quantification.txt", sep="\t")
lipidome_pos = pd.read_csv("/Users/senol/Desktop/GBM/lipidome_positive_quantification.txt", sep="\t")
lipidome_neg = pd.read_csv("/Users/senol/Desktop/GBM/lipidome_negative_quantification.txt", sep="\t")
clinical = pd.read_csv("/Users/senol/Desktop/GBM/clinical_patient.txt", sep="\t", skiprows=4)

# Clean clinical data
clinical.columns = clinical.columns.str.strip()
clinical = clinical[clinical['PATIENT_ID'].notna()]

# Extract sample IDs from column names
def extract_sample_ids(df):
    """Extract sample IDs from dataframe columns"""
    # First column is usually ID/name, rest are samples
    sample_cols = [col for col in df.columns if col.startswith('C3')]
    return sample_cols

# Get common samples across all datasets
metabolome_samples = extract_sample_ids(metabolome)
lipidome_pos_samples = extract_sample_ids(lipidome_pos)
lipidome_neg_samples = extract_sample_ids(lipidome_neg)

# Get clinical patient IDs
clinical_patient_ids = clinical['PATIENT_ID'].dropna().unique().tolist()

# Find common samples across ALL datasets including clinical
common_samples = list(set(metabolome_samples) & set(lipidome_pos_samples) & 
                     set(lipidome_neg_samples) & set(clinical_patient_ids))
common_samples = sorted(common_samples)  # Sort for consistency

print(f"Total patients in clinical data: {len(clinical_patient_ids)}")
print(f"Total patients in metabolome: {len(metabolome_samples)}")
print(f"Total patients in lipidome_pos: {len(lipidome_pos_samples)}")
print(f"Total patients in lipidome_neg: {len(lipidome_neg_samples)}")
print(f"Found {len(common_samples)} common patients across ALL datasets (metabolome + lipidome + clinical)")
print(f"Using {len(common_samples)} patients for analysis")

# Prepare data matrices
def prepare_data_matrix(df, sample_cols, id_col='ENTITY_STABLE_ID', name_col='NAME'):
    """Prepare data matrix with samples as columns and features as rows"""
    if id_col in df.columns:
        df = df.set_index(id_col)
    elif name_col in df.columns:
        df = df.set_index(name_col)
    else:
        df = df.set_index(df.columns[0])
    
    # Select only common samples
    available_samples = [s for s in sample_cols if s in df.columns]
    df_subset = df[available_samples].T
    
    # Convert to numeric, handling NA
    df_subset = df_subset.apply(pd.to_numeric, errors='coerce')
    
    # Remove features with >50% missing values
    missing_threshold = 0.5
    df_subset = df_subset.loc[:, df_subset.isna().mean() < missing_threshold]
    
    # Fill remaining missing values with median
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

print(f"Final dataset: {len(common_samples_final)} samples")
print(f"  - Metabolites: {metabolome_matrix.shape[1]}")
print(f"  - Lipids: {lipidome_combined.shape[1]}")

# Prepare clinical data
clinical_clean = clinical[clinical['PATIENT_ID'].isin(common_samples_final)].copy()
clinical_clean = clinical_clean.set_index('PATIENT_ID')

# Create survival data
def prepare_survival_data(clinical_df, sample_ids):
    """Prepare survival data"""
    survival_data = pd.DataFrame(index=sample_ids)
    
    for sample_id in sample_ids:
        if sample_id in clinical_df.index:
            row = clinical_df.loc[sample_id]
            
            # Time to event (days)
            if pd.notna(row.get('PATH_DIAG_TO_DEATH_DAYS')):
                survival_data.loc[sample_id, 'time'] = row['PATH_DIAG_TO_DEATH_DAYS']
                survival_data.loc[sample_id, 'event'] = 1
            elif pd.notna(row.get('PATH_DIAG_TO_LAST_CONTACT_DAYS')):
                survival_data.loc[sample_id, 'time'] = row['PATH_DIAG_TO_LAST_CONTACT_DAYS']
                survival_data.loc[sample_id, 'event'] = 0
            else:
                survival_data.loc[sample_id, 'time'] = np.nan
                survival_data.loc[sample_id, 'event'] = 0
            
            # Recurrence
            if pd.notna(row.get('TUMOR_REOCCUR_AFTER_TREATMENT')):
                survival_data.loc[sample_id, 'recurrence'] = 1 if str(row['TUMOR_REOCCUR_AFTER_TREATMENT']).upper() == 'TRUE' else 0
            else:
                survival_data.loc[sample_id, 'recurrence'] = 0
            
            # Recurrence time
            if pd.notna(row.get('PATH_DIAG_TO_TUMOR_REOCCUR_DAYS')):
                survival_data.loc[sample_id, 'recurrence_time'] = row['PATH_DIAG_TO_TUMOR_REOCCUR_DAYS']
            else:
                survival_data.loc[sample_id, 'recurrence_time'] = np.nan
        else:
            survival_data.loc[sample_id, 'time'] = np.nan
            survival_data.loc[sample_id, 'event'] = 0
            survival_data.loc[sample_id, 'recurrence'] = 0
            survival_data.loc[sample_id, 'recurrence_time'] = np.nan
    
    return survival_data.dropna(subset=['time'])

survival_data = prepare_survival_data(clinical_clean, common_samples_final)
print(f"Survival data: {len(survival_data)} samples with survival information")

# ============================================================================
# 2. LIPID-METABOLITE CORRELATION NETWORKS
# ============================================================================
print("\n[2/8] Computing lipid-metabolite correlation networks...")

def compute_correlation_network(df1, df2, method='spearman', threshold=0.3, fdr_threshold=0.05):
    """Compute correlation network between two datasets"""
    # Align samples
    common_idx = df1.index.intersection(df2.index)
    df1_aligned = df1.loc[common_idx]
    df2_aligned = df2.loc[common_idx]
    
    # Compute correlations
    correlations = []
    pvalues = []
    
    for col1 in df1_aligned.columns:
        for col2 in df2_aligned.columns:
            # Remove NaN pairs
            mask = ~(df1_aligned[col1].isna() | df2_aligned[col2].isna())
            if mask.sum() < 10:  # Need at least 10 samples
                continue
            
            if method == 'spearman':
                corr, pval = stats.spearmanr(df1_aligned[col1][mask], df2_aligned[col2][mask])
            else:
                corr, pval = stats.pearsonr(df1_aligned[col1][mask], df2_aligned[col2][mask])
            
            correlations.append({
                'feature1': col1,
                'feature2': col2,
                'correlation': corr,
                'pvalue': pval
            })
            pvalues.append(pval)
    
    corr_df = pd.DataFrame(correlations)
    
    # FDR correction
    if len(pvalues) > 0:
        _, corr_df['fdr'], _, _ = multipletests(pvalues, method='fdr_bh')
        corr_df = corr_df[(abs(corr_df['correlation']) >= threshold) & (corr_df['fdr'] < fdr_threshold)]
    
    return corr_df

# Compute correlations
lipid_metabolite_corr = compute_correlation_network(
    lipidome_combined, metabolome_matrix, 
    method='spearman', threshold=0.3, fdr_threshold=0.05
)

print(f"Found {len(lipid_metabolite_corr)} significant correlations (|r|>=0.3, FDR<0.05)")

# Save correlation results
lipid_metabolite_corr.to_csv(f"{output_dir}/results/lipid_metabolite_correlations.csv", index=False)

def _safe_read_feature_name_mapping(mapping_xlsx_path: str):
    """Best-effort feature name mapping for prettier plots.

    Returns dict[str, str]. If file/engine missing, returns empty dict.
    """
    try:
        if not os.path.exists(mapping_xlsx_path):
            return {}
        df_map = pd.read_excel(mapping_xlsx_path)
        # try common column names
        cols = {c.lower(): c for c in df_map.columns}
        key_col = cols.get('feature') or cols.get('feature_id') or cols.get('id') or df_map.columns[0]
        val_col = cols.get('name') or cols.get('feature_name') or cols.get('full_name') or df_map.columns[-1]
        m = (
            df_map[[key_col, val_col]]
            .dropna()
            .astype(str)
            .set_index(key_col)[val_col]
            .to_dict()
        )
        return {k: v for k, v in m.items() if k and v}
    except Exception:
        return {}


def _normalize_lipid_like_label(label: str) -> str:
    """Standardize lipid-like labels for consistent plotting.

    - Removes trailing ' A'/' B' suffixes (e.g., 'PC(34:1) A' -> 'PC(34:1)').
    - Standardizes TG with chain-level composition to total composition:
      'TG(14:0/20:4/20:4)' -> 'TG(54:8)'.
    - Standardizes spacing and bracket usage.
    """
    if label is None:
        return label
    s = str(label).strip()

    # If there are multiple alternative annotations (e.g., "X;Y" or "X:Y"),
    # collapse to a single shorthand label when both alternatives map to the same sum composition.
    if ";" in s or (":" in s and s.count("(") >= 2 and s.count(")") >= 2):
        sep = ";" if ";" in s else ":"
        parts = [p.strip() for p in s.split(sep) if p.strip()]
        if len(parts) == 2:
            a = _normalize_lipid_like_label(parts[0])
            b = _normalize_lipid_like_label(parts[1])
            if a == b:
                return a

    # Remove common trailing A/B flags (keep internal letters untouched)
    s = re.sub(r"\s+([AB])\s*$", "", s)

    # Normalize "TG 46:2" -> "TG(46:2)"
    s = re.sub(r"^TG\s+(\d{2}:\d{1,2})$", r"TG(\1)", s)

    # Collapse ceramide detailed composition to sum composition:
    # Cer(d16:1/18:0) -> Cer(d34:1)
    m = re.match(r"^(Cer)\(d(\d+):(\d+)/(\d+):(\d+)\)$", s)
    if m:
        head = m.group(1)
        base_c = int(m.group(2))
        base_d = int(m.group(3))
        acyl_c = int(m.group(4))
        acyl_d = int(m.group(5))
        return f"{head}(d{base_c+acyl_c}:{base_d+acyl_d})"

    # Convert TG(chain/chain/chain) -> TG(total:total)
    m = re.match(r"^(TG)\(\s*([0-9]{1,2}:[0-9]{1,2})\s*/\s*([0-9]{1,2}:[0-9]{1,2})\s*/\s*([0-9]{1,2}:[0-9]{1,2})\s*\)\s*$", s)
    if m:
        head = m.group(1)
        parts = [m.group(2), m.group(3), m.group(4)]
        carbons = 0
        doubles = 0
        ok = True
        for p in parts:
            try:
                c, d = p.split(":")
                carbons += int(c)
                doubles += int(d)
            except Exception:
                ok = False
                break
        if ok:
            s = f"{head}({carbons}:{doubles})"

    # Convert phospholipids with chain-level composition to sum composition:
    # PC(18:1/24:1) -> PC(42:2)
    m = re.match(r"^([A-Z]{2,4})\(\s*([0-9]{1,2}):([0-9]{1,2})\s*/\s*([0-9]{1,2}):([0-9]{1,2})\s*\)\s*$", s)
    if m:
        head = m.group(1)
        if head in {"PC", "PE", "PS", "PI", "PG", "PA"}:
            c1, d1, c2, d2 = int(m.group(2)), int(m.group(3)), int(m.group(4)), int(m.group(5))
            s = f"{head}({c1+c2}:{d1+d2})"

    # Convert sphingomyelin detailed chain composition to sum composition:
    # SM(d18:2/25:1) -> SM(d43:3)
    m = re.match(r"^(SM)\(d(\d+):(\d+)/(\d+):(\d+)\)\s*$", s)
    if m:
        head = m.group(1)
        base_c, base_d = int(m.group(2)), int(m.group(3))
        acyl_c, acyl_d = int(m.group(4)), int(m.group(5))
        s = f"{head}(d{base_c+acyl_c}:{base_d+acyl_d})"

    # Standardize class formatting like "PC 20:4/36:5" -> "PC(20:4/36:5)"
    s = re.sub(r"^([A-Z]{2,4})\s+([0-9]{1,2}:[0-9]{1,2}/[0-9]{1,2}:[0-9]{1,2})$", r"\1(\2)", s)

    # Remove redundant spaces inside parentheses
    s = re.sub(r"\(\s+", "(", s)
    s = re.sub(r"\s+\)", ")", s)
    return s


def _pretty_feature_label(raw: str, mapping: dict) -> str:
    if raw is None:
        return raw
    s = str(raw)
    # Apply mapping first (Unknown_### -> identified names if present)
    s = mapping.get(s, s)
    # Remove metabolomics annotation artifacts from display labels
    # Examples: "cholesterol_(spectral_match)" / "aminomalonic_acid_(spectralmatch)" / "(spectral match)"
    s = re.sub(r"\s*\(\s*spectral\s*match\s*\)\s*", "", s, flags=re.IGNORECASE)
    s = re.sub(r"_\(\s*spectral[_\s]*match\s*\)\s*", "", s, flags=re.IGNORECASE)
    s = re.sub(r"_\(\s*spectralmatch\s*\)\s*", "", s, flags=re.IGNORECASE)
    # Remove internal standard and duplicate suffixes commonly appended during preprocessing
    s = re.sub(r"(_IS_dup\d+|_IS_dup\d*|_IS)$", "", s)
    # Remove trailing batch/isomer flags like "_B" used in some exported matrices
    # Handle both end-of-string and delimiter cases (e.g., "... )_A;...").
    s = re.sub(r"(?<=\))_A(?=\b|$|[;:])", "", s)
    s = re.sub(r"(?<=\))_B(?=\b|$|[;:])", "", s)
    s = re.sub(r"_A$", "", s)
    s = re.sub(r"_B$", "", s)
    s = _normalize_lipid_like_label(s)
    return s


# Create correlation heatmaps
if len(lipid_metabolite_corr) > 0:
    feature_name_map = _safe_read_feature_name_mapping(os.path.join(output_dir, "Feature_Names_Mapping.xlsx"))

    # (A) "Top associations" heatmap (top 100 by absolute correlation)
    top_corr = lipid_metabolite_corr.assign(abs_corr=lambda d: d["correlation"].abs()).nlargest(100, "abs_corr")
    pivot_data = top_corr.pivot_table(values="correlation", index="feature1", columns="feature2")

    pivot_data = pivot_data.rename(
        index=lambda x: _pretty_feature_label(x, feature_name_map),
        columns=lambda x: _pretty_feature_label(x, feature_name_map),
    )

    fig, ax = plt.subplots(figsize=(22, 16))
    hm = sns.heatmap(
        pivot_data,
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        square=False,
        linewidths=0.4,
        cbar_kws={"shrink": 0.8},
        ax=ax,
    )
    # Rasterize the heavy artist for better readability/export performance in vector workflows
    try:
        hm.collections[0].set_rasterized(True)
    except Exception:
        pass

    plt.xlabel("Metabolites", fontsize=12)
    plt.ylabel("Lipids", fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/plots/lipid_metabolite_correlation_heatmap.png", dpi=450, bbox_inches="tight")
    plt.close()

    # (B) "Figure 3-style" heatmap: top 20 lipids (lowest FDR) and their significant metabolite correlations
    # This matches the caption logic used in the Word draft comments.
    lipid_top20 = (
        lipid_metabolite_corr.groupby("feature1")["fdr"].min().sort_values().head(20).index.tolist()
    )
    fig3_df = lipid_metabolite_corr[lipid_metabolite_corr["feature1"].isin(lipid_top20)].copy()

    # Remove specific unidentified metabolite features from the figure for clarity
    _excluded_unknowns = {"Unknown_003", "Unknown_012", "Unknown_018"}
    fig3_df = fig3_df[~fig3_df["feature2"].isin(_excluded_unknowns)].copy()

    # Keep metabolites that have at least one significant link to the top lipids; order by connectivity then by best FDR
    metab_order = (
        fig3_df.groupby("feature2")
        .agg(n_links=("feature1", "count"), best_fdr=("fdr", "min"), best_abs_corr=("correlation", lambda s: s.abs().max()))
        .sort_values(["n_links", "best_fdr", "best_abs_corr"], ascending=[False, True, False])
        .index.tolist()
    )

    fig3_pivot = fig3_df.pivot_table(values="correlation", index="feature1", columns="feature2")
    fig3_pivot = fig3_pivot.reindex(index=lipid_top20, columns=metab_order)

    fig3_pivot = fig3_pivot.rename(
        index=lambda x: _pretty_feature_label(x, feature_name_map),
        columns=lambda x: _pretty_feature_label(x, feature_name_map),
    )

    fig, ax = plt.subplots(figsize=(30, 16))
    hm = sns.heatmap(
        fig3_pivot,
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        square=False,
        linewidths=0.6,
        linecolor="white",
        cbar_kws={"shrink": 0.85},
        ax=ax,
    )
    try:
        hm.collections[0].set_rasterized(True)
    except Exception:
        pass

    # Improve readability: larger fonts and stronger contrast
    ax.set_xlabel("Metabolites", fontsize=20, fontweight="bold", labelpad=10)
    ax.set_ylabel("Top 20 lipids (lowest FDR)", fontsize=20, fontweight="bold", labelpad=10)
    ax.tick_params(axis="both", which="major", labelsize=18, colors="black")
    for t in ax.get_xticklabels():
        t.set_rotation(55)
        t.set_horizontalalignment("right")
    for t in ax.get_yticklabels():
        t.set_rotation(0)

    # Colorbar label and ticks
    try:
        cbar = hm.collections[0].colorbar
        cbar.ax.tick_params(labelsize=16, colors="black")
        cbar.set_label("Spearman ρ", fontsize=18, fontweight="bold", labelpad=12)
    except Exception:
        pass

    plt.tight_layout()
    plt.savefig(f"{output_dir}/plots/figure3_top20lipids_correlation_heatmap.png", dpi=450, bbox_inches="tight")
    plt.close()

# Create network visualization
def create_correlation_network(corr_df, top_n=50):
    """Create network graph from correlations"""
    G = nx.Graph()
    
    # Add top correlations
    top_corr = corr_df.nlargest(top_n, 'correlation')
    
    for _, row in top_corr.iterrows():
        G.add_edge(
            row['feature1'],
            row['feature2'],
            weight=abs(row['correlation']),
            correlation=row['correlation'],
        )
    
    return G

if len(lipid_metabolite_corr) > 0:
    G = create_correlation_network(lipid_metabolite_corr, top_n=100)
    
    plt.figure(figsize=(20, 20))
    # Layout: keep it as ONE figure (no separate panels), but reduce whitespace by packing components
    # so disconnected parts don't drift far apart.
    components = [list(c) for c in nx.connected_components(G)]
    components.sort(key=len, reverse=True)

    pos = {}
    # deterministic offsets to keep components close to center
    offsets = [(0.0, 0.0), (0.9, 0.0), (-0.9, 0.0), (0.0, 0.9), (0.0, -0.9),
               (0.7, 0.7), (-0.7, 0.7), (0.7, -0.7), (-0.7, -0.7)]

    for i, nodes in enumerate(components):
        sub = G.subgraph(nodes)
        # smaller k -> tighter cluster; higher iterations -> stable layout
        sub_pos = nx.spring_layout(sub, k=0.22, iterations=300, seed=42)
        ox, oy = offsets[i] if i < len(offsets) else (0.0, 0.0)
        for n, (x, y) in sub_pos.items():
            pos[n] = (x + ox, y + oy)

    # Final rescale to use canvas area better (fills right/left/top/bottom without changing figure size)
    pos = nx.rescale_layout_dict(pos, scale=3.2)
    
    # Determine node types (lipid vs metabolite) for coloring
    lipid_nodes = set(lipid_metabolite_corr["feature1"].astype(str).unique().tolist())
    metabolite_nodes = set(lipid_metabolite_corr["feature2"].astype(str).unique().tolist())
    main_hubs = {"Unknown_048", "creatinine", "L-mimosine"}
    node_colors = []
    node_sizes = []
    for n in G.nodes():
        # Base size scaled by degree (fills space without changing figure size)
        deg = G.degree(n)
        base_size = 1100 + 55 * deg

        if str(n) in main_hubs:
            node_colors.append("#006400")  # dark green for main hubs
            node_sizes.append(max(base_size, 3200))
            continue

        if n in lipid_nodes and n not in metabolite_nodes:
            node_colors.append("#ff7f0e")  # lipid: orange
        elif n in metabolite_nodes and n not in lipid_nodes:
            node_colors.append("#1f77b4")  # metabolite: blue
        else:
            node_colors.append("#7f7f7f")  # fallback/overlap

        node_sizes.append(base_size)

    # Edge color by sign (correlation), and width by magnitude (bold for high |corr|)
    edges = list(G.edges())
    edge_corr = np.array([G[u][v]["correlation"] for u, v in edges], dtype=float)
    edge_abs = np.abs(edge_corr)
    if len(edge_abs) > 0 and float(edge_abs.max()) != float(edge_abs.min()):
        w_min, w_max = 0.6, 4.0
        edge_widths = w_min + (edge_abs - edge_abs.min()) / (edge_abs.max() - edge_abs.min()) * (w_max - w_min)
    else:
        edge_widths = np.full_like(edge_abs, 2.0, dtype=float)

    nx.draw_networkx_edges(
        G,
        pos,
        alpha=0.55,
        width=edge_widths.tolist(),
        edge_color=edge_corr.tolist(),
        edge_cmap=plt.cm.RdBu_r,
        edge_vmin=-1,
        edge_vmax=1,
    )
    nx.draw_networkx_nodes(
        G,
        pos,
        node_size=node_sizes,
        node_color=node_colors,
        alpha=0.9,
        edgecolors="white",
        linewidths=1.0,
    )
    # Clean display labels (remove spectral match / _IS / _B artifacts) for readability
    try:
        label_map = {n: _pretty_feature_label(n, feature_name_map) for n in G.nodes()}
    except Exception:
        label_map = {n: str(n) for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=label_map, font_size=13, font_weight='bold')
    
    # Title removed - caption will be in manuscript
    plt.axis('off')
    plt.margins(0.01)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/plots/lipid_metabolite_network.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save network statistics
    network_stats = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'density': nx.density(G),
        'avg_clustering': nx.average_clustering(G)
    }
    pd.Series(network_stats).to_csv(f"{output_dir}/results/network_statistics.csv")

# ============================================================================
# 3. METABOLIC-LIPIDOMIC SUBTYPE CLASSIFICATION
# ============================================================================
print("\n[3/8] Performing metabolic-lipidomic subtype classification...")

# Combine metabolome and lipidome for clustering
combined_data = pd.concat([metabolome_matrix, lipidome_combined], axis=1)
combined_data = combined_data.fillna(combined_data.median())

# Standardize
scaler = StandardScaler()
combined_scaled = pd.DataFrame(
    scaler.fit_transform(combined_data),
    index=combined_data.index,
    columns=combined_data.columns
)

# Determine optimal number of clusters using multiple methods
def find_optimal_clusters(data, max_k=10):
    """Find optimal number of clusters using multiple metrics"""
    results = []
    
    for k in range(2, max_k + 1):
        # KMeans
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        kmeans_labels = kmeans.fit_predict(data)
        kmeans_silhouette = silhouette_score(data, kmeans_labels)
        
        # Hierarchical
        hierarchical = AgglomerativeClustering(n_clusters=k)
        hier_labels = hierarchical.fit_predict(data)
        hier_silhouette = silhouette_score(data, hier_labels)
        
        results.append({
            'k': k,
            'kmeans_silhouette': kmeans_silhouette,
            'hierarchical_silhouette': hier_silhouette
        })
    
    return pd.DataFrame(results)

cluster_metrics = find_optimal_clusters(combined_scaled, max_k=8)

# Plot cluster metrics
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
axes[0].plot(cluster_metrics['k'], cluster_metrics['kmeans_silhouette'], 
             marker='o', linewidth=2, label='KMeans')
axes[0].plot(cluster_metrics['k'], cluster_metrics['hierarchical_silhouette'], 
             marker='s', linewidth=2, label='Hierarchical')
axes[0].set_xlabel('Number of Clusters (k)', fontsize=12)
axes[0].set_ylabel('Silhouette Score', fontsize=12)
axes[0].set_title('Cluster Quality Metrics', fontsize=14, fontweight='bold')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Optimal k (highest silhouette)
optimal_k = cluster_metrics.loc[cluster_metrics['kmeans_silhouette'].idxmax(), 'k']
axes[0].axvline(x=optimal_k, color='r', linestyle='--', alpha=0.5, label=f'Optimal k={int(optimal_k)}')
axes[0].legend()

# Elbow method (inertia)
inertias = []
for k in range(2, 9):
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    kmeans.fit(combined_scaled)
    inertias.append(kmeans.inertia_)

axes[1].plot(range(2, 9), inertias, marker='o', linewidth=2)
axes[1].set_xlabel('Number of Clusters (k)', fontsize=12)
axes[1].set_ylabel('Inertia', fontsize=12)
axes[1].set_title('Elbow Method', fontsize=14, fontweight='bold')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{output_dir}/plots/cluster_optimization.png", dpi=300, bbox_inches='tight')
plt.close()

# Perform clustering with optimal k
optimal_k = int(optimal_k)
print(f"Using k={optimal_k} clusters")

# KMeans clustering
kmeans_final = KMeans(n_clusters=optimal_k, random_state=42, n_init=20)
kmeans_labels = kmeans_final.fit_predict(combined_scaled)
combined_data['kmeans_subtype'] = kmeans_labels

# Hierarchical clustering
hierarchical_final = AgglomerativeClustering(n_clusters=optimal_k)
hier_labels = hierarchical_final.fit_predict(combined_scaled)
combined_data['hierarchical_subtype'] = hier_labels

# NMF-based clustering
def nmf_clustering(data, n_components=optimal_k, n_runs=10):
    """Perform NMF clustering"""
    best_score = -1
    best_labels = None
    
    # NMF requires non-negative values, so we need to shift the data
    data_min = data.min().min()
    if data_min < 0:
        data_shifted = data - data_min + 1e-10
    else:
        data_shifted = data + 1e-10
    
    for _ in range(n_runs):
        try:
            nmf = NMF(n_components=n_components, random_state=None, max_iter=500, init='random')
            W = nmf.fit_transform(data_shifted)
            H = nmf.components_
            
            # Assign to cluster with highest weight
            labels = np.argmax(W, axis=1)
            score = silhouette_score(data, labels)  # Use original data for silhouette
            
            if score > best_score:
                best_score = score
                best_labels = labels
        except:
            continue
    
    return best_labels, best_score

nmf_labels, nmf_score = nmf_clustering(combined_scaled, n_components=optimal_k)
combined_data['nmf_subtype'] = nmf_labels

print(f"NMF clustering silhouette score: {nmf_score:.3f}")

# Visualize clusters using t-SNE
print("Computing t-SNE embedding...")
tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(combined_scaled)-1), max_iter=1000)
tsne_results = tsne.fit_transform(combined_scaled)

# Plot t-SNE with different clustering methods
fig, axes = plt.subplots(1, 3, figsize=(21, 6))

def _subtype_label(cluster_value):
    # Map numeric cluster labels to A/B for readability in figures
    try:
        c = int(cluster_value)
    except Exception:
        return f"Subtype {cluster_value}"
    if c == 0:
        return "Subtype A"
    if c == 1:
        return "Subtype B"
    return f"Subtype {c}"

# Make subgroup colors clearly separable (A=red)
subtype_colors = {
    0: "#d62728",  # red
    1: "#1f77b4",  # blue
}

def _permanova_pvalue(X, groups, n_perm=999, seed=42):
    """
    PERMANOVA for multivariate separation (Euclidean distance), returning (pseudo-F, p).
    p is permutation-based; with 999 permutations, minimum reportable p is 0.001.
    """
    rng = np.random.default_rng(seed)
    g = np.asarray(groups).astype(int)
    D = pairwise_distances(X, metric="euclidean")
    N = D.shape[0]
    SS_total = (D ** 2).sum() / N

    def ss_within(gv):
        ss = 0.0
        for lev in np.unique(gv):
            idx = np.where(gv == lev)[0]
            if len(idx) <= 1:
                continue
            Dij = D[np.ix_(idx, idx)]
            ss += (Dij ** 2).sum() / len(idx)
        return ss

    SS_w = ss_within(g)
    SS_b = SS_total - SS_w
    k = len(np.unique(g))
    df_b = k - 1
    df_w = N - k
    F = (SS_b / df_b) / (SS_w / df_w)

    cnt = 0
    for _ in range(n_perm):
        gp = rng.permutation(g)
        SS_wp = ss_within(gp)
        SS_bp = SS_total - SS_wp
        Fp = (SS_bp / df_b) / (SS_wp / df_w)
        if Fp >= F:
            cnt += 1
    p = (cnt + 1) / (n_perm + 1)
    return float(F), float(p)

for idx, (method, labels) in enumerate([
    ('KMeans', kmeans_labels),
    ('Hierarchical', hier_labels),
    ('NMF', nmf_labels)
]):
    unique_labels = sorted(set(int(x) for x in labels))

    for lab in unique_labels:
        mask = (labels == lab)
        color = subtype_colors.get(int(lab), "#7f7f7f")
        axes[idx].scatter(
            tsne_results[mask, 0],
            tsne_results[mask, 1],
            s=95,
            alpha=0.85,
            color=color,
            edgecolors="white",
            linewidths=0.6,
            label=_subtype_label(lab),
        )

    axes[idx].set_title(f'{method} Clustering (k={optimal_k})', fontsize=20, fontweight='bold')
    axes[idx].set_xlabel('t-SNE 1', fontsize=18)
    axes[idx].set_ylabel('t-SNE 2', fontsize=18)
    axes[idx].tick_params(axis='both', which='major', labelsize=14)
    axes[idx].legend(frameon=False, fontsize=14, loc="best", handletextpad=0.4)

    # Add PERMANOVA p-value in bottom-right corner
    try:
        _, p_perm = _permanova_pvalue(combined_scaled, labels, n_perm=999, seed=42)
        p_txt = "<0.001" if p_perm <= 0.001 else f"{p_perm:.3g}"
        axes[idx].text(
            0.98,
            0.02,
            f"PERMANOVA p={p_txt}",
            transform=axes[idx].transAxes,
            ha="right",
            va="bottom",
            fontsize=14,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="none", alpha=0.75),
        )
    except Exception:
        pass

plt.tight_layout()
plt.savefig(f"{output_dir}/plots/clustering_tsne.png", dpi=300, bbox_inches='tight')
plt.close()

# Save cluster assignments
cluster_assignments = pd.DataFrame({
    'sample_id': combined_data.index,
    'kmeans_subtype': kmeans_labels,
    'hierarchical_subtype': hier_labels,
    'nmf_subtype': nmf_labels
})
cluster_assignments.to_csv(f"{output_dir}/results/cluster_assignments.csv", index=False)

# ============================================================================
# 4. SURVIVAL ANALYSIS BY SUBTYPES
# ============================================================================
print("\n[4/8] Performing survival analysis by subtypes...")

# Merge survival data with cluster assignments
survival_clusters = survival_data.merge(
    cluster_assignments.set_index('sample_id'),
    left_index=True, right_index=True,
    how='inner'
)

# Perform survival analysis for each clustering method
def plot_survival_curves(survival_df, cluster_col, title, save_path):
    """Plot Kaplan-Meier survival curves by cluster"""
    kmf = KaplanMeierFitter()
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    clusters = sorted(survival_df[cluster_col].unique())
    colors = plt.cm.tab10(np.linspace(0, 1, len(clusters)))
    
    def _subtype_label(cluster_value):
        # Map numeric cluster labels to A/B for readability in figures
        try:
            c = int(cluster_value)
        except Exception:
            return f"Subtype {cluster_value}"
        if c == 0:
            return "Subtype A"
        if c == 1:
            return "Subtype B"
        return f"Subtype {c}"

    for cluster, color in zip(clusters, colors):
        cluster_data = survival_df[survival_df[cluster_col] == cluster]
        cluster_data_clean = cluster_data.dropna(subset=['time', 'event'])
        
        if len(cluster_data_clean) > 0:
            kmf.fit(cluster_data_clean['time'], cluster_data_clean['event'], 
                   label=f"{_subtype_label(cluster)} (n={len(cluster_data_clean)})")
            kmf.plot_survival_function(ax=ax, color=color, linewidth=2)
    
    # Log-rank test
    p_value = None
    if len(clusters) > 1:
        groups = []
        times = []
        events = []
        
        for cluster in clusters:
            cluster_data = survival_df[survival_df[cluster_col] == cluster].dropna(subset=['time', 'event'])
            if len(cluster_data) > 0:
                groups.extend([cluster] * len(cluster_data))
                times.extend(cluster_data['time'].tolist())
                events.extend(cluster_data['event'].tolist())
        
        if len(set(groups)) > 1:
            results = logrank_test(
                [t for g, t in zip(groups, times) if g == clusters[0]],
                [t for g, t in zip(groups, times) if g == clusters[1]],
                event_observed_A=[e for g, e in zip(groups, events) if g == clusters[0]],
                event_observed_B=[e for g, e in zip(groups, events) if g == clusters[1]]
            )
            p_value = results.p_value
    
    # Add p-value to plot if significant
    if p_value is not None:
        if p_value < 0.05:
            ax.text(0.05, 0.05, f'Log-rank p-value = {p_value:.4f}***', 
                   transform=ax.transAxes, fontsize=12, fontweight='bold',
                   bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
        else:
            ax.text(0.05, 0.05, f'Log-rank p-value = {p_value:.4f}', 
                   transform=ax.transAxes, fontsize=12,
                   bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.7))
    
    ax.set_xlabel('Time (days)', fontsize=12)
    ax.set_ylabel('Survival Probability', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return p_value

# Plot survival curves for each clustering method
for method in ['kmeans_subtype', 'hierarchical_subtype', 'nmf_subtype']:
    pval = plot_survival_curves(
        survival_clusters, 
        method,
        f'Survival Analysis: {method.replace("_", " ").title()}',
        f"{output_dir}/plots/survival_{method}.png"
    )
    if pval:
        print(f"  {method}: Log-rank p-value = {pval:.4f}")

# ============================================================================
# 5. DIFFERENTIAL ANALYSIS BETWEEN SUBTYPES
# ============================================================================
print("\n[5/8] Performing differential analysis between subtypes...")

def differential_analysis(data, labels, method='kmeans'):
    """Perform differential analysis between clusters"""
    results = []
    
    clusters = sorted(set(labels))
    
    for cluster in clusters:
        cluster_samples = data[labels == cluster]
        other_samples = data[labels != cluster]
        
        for feature in data.columns:
            cluster_values = cluster_samples[feature].dropna()
            other_values = other_samples[feature].dropna()
            
            if len(cluster_values) < 3 or len(other_values) < 3:
                continue
            
            # Mann-Whitney U test
            stat, pval = stats.mannwhitneyu(cluster_values, other_values, alternative='two-sided')
            
            # Effect size (Cohen's d)
            pooled_std = np.sqrt(((len(cluster_values) - 1) * cluster_values.std()**2 + 
                                (len(other_values) - 1) * other_values.std()**2) / 
                               (len(cluster_values) + len(other_values) - 2))
            cohens_d = (cluster_values.mean() - other_values.mean()) / pooled_std if pooled_std > 0 else 0
            
            results.append({
                'cluster': cluster,
                'feature': feature,
                'mean_cluster': cluster_values.mean(),
                'mean_other': other_values.mean(),
                'fold_change': cluster_values.mean() / other_values.mean() if other_values.mean() != 0 else np.inf,
                'pvalue': pval,
                'cohens_d': cohens_d
            })
    
    results_df = pd.DataFrame(results)
    
    # FDR correction
    if len(results_df) > 0:
        _, results_df['fdr'], _, _ = multipletests(results_df['pvalue'], method='fdr_bh')
    
    return results_df

# Perform differential analysis for NMF subtypes (most novel)
diff_results = differential_analysis(combined_data.drop(columns=['kmeans_subtype', 'hierarchical_subtype', 'nmf_subtype']), 
                                     nmf_labels, method='nmf')

# Filter significant results
significant_features = diff_results[
    (diff_results['fdr'] < 0.05) & (abs(diff_results['cohens_d']) > 0.5)
].sort_values('fdr')

print(f"Found {len(significant_features)} significantly different features (FDR<0.05, |d|>0.5)")

# Save results
significant_features.to_csv(f"{output_dir}/results/differential_features.csv", index=False)

# Create heatmap of top differential features
if len(significant_features) > 0:
    top_features = significant_features.head(50)['feature'].unique()
    top_features = [f for f in top_features if f in combined_data.columns]
    
    if len(top_features) > 0:
        # Load mapping for nicer display labels (and shorthand normalization)
        feature_name_map = _safe_read_feature_name_mapping(os.path.join(output_dir, "Feature_Names_Mapping.xlsx"))

        # Create heatmap data
        heatmap_data = combined_data[top_features].copy()
        heatmap_data['subtype'] = nmf_labels
        
        # Average by subtype
        heatmap_avg = heatmap_data.groupby('subtype')[top_features].mean()
        
        # Normalize for visualization
        heatmap_avg_norm = (heatmap_avg - heatmap_avg.mean()) / heatmap_avg.std()
        
        # Shorten feature names for readability in the figure (sum composition)
        y_labels = [_pretty_feature_label(f, feature_name_map) for f in heatmap_avg_norm.columns]

        # Improve readability + Times New Roman for publication
        plt.rcParams['font.family'] = 'Times New Roman'

        fig, ax = plt.subplots(figsize=(14, 11))
        sns.heatmap(
            heatmap_avg_norm.T,
            cmap='RdBu_r',
            center=0,
            annot=True,
            fmt='.2f',
            cbar_kws={"shrink": 0.8},
            xticklabels=[
                ('Subtype A' if int(i) == 0 else 'Subtype B' if int(i) == 1 else f'Subtype {int(i)}')
                for i in heatmap_avg_norm.index
            ],
            yticklabels=y_labels,
            ax=ax,
        )
        
        plt.title('Top 50 Differential Features by Subtype (Z-score normalized)', 
                 fontsize=24, fontweight='bold')
        plt.xlabel('Subtype', fontsize=20)
        plt.ylabel('Feature', fontsize=20)
        ax.tick_params(axis='both', which='major', labelsize=18)

        # Make the annotation numbers larger for readability
        for t in ax.texts:
            t.set_fontsize(15)

        try:
            cbar = ax.collections[0].colorbar
            cbar.ax.tick_params(labelsize=16)
        except Exception:
            pass
        plt.tight_layout()
        plt.savefig(f"{output_dir}/plots/differential_features_heatmap.png", dpi=300, bbox_inches='tight')
        plt.close()

# ============================================================================
# 6. PATHWAY ENRICHMENT ANALYSIS
# ============================================================================
print("\n[6/8] Performing pathway enrichment analysis...")

# This is a simplified version - in real analysis, you'd use KEGG/Reactome APIs
def simple_pathway_analysis(feature_names):
    """Simple pathway categorization based on feature names"""
    pathways = {
        'Lipid Metabolism': [],
        'Amino Acid Metabolism': [],
        'Carbohydrate Metabolism': [],
        'Nucleotide Metabolism': [],
        'Energy Metabolism': [],
        'Other': []
    }
    
    lipid_keywords = ['lipid', 'cer', 'ce', 'carnitine', 'cholesterol', 'fatty', 'acyl']
    aa_keywords = ['amino', 'acid', 'alanine', 'glutamine', 'glycine', 'serine', 'threonine']
    carb_keywords = ['glucose', 'fructose', 'galactose', 'sugar', 'glycol', 'pentose']
    nucleotide_keywords = ['nucleotide', 'purine', 'pyrimidine', 'adenine', 'uracil']
    energy_keywords = ['atp', 'adp', 'nad', 'fad', 'coenzyme', 'phosphate']
    
    for feature in feature_names:
        feature_lower = str(feature).lower()
        categorized = False
        
        if any(kw in feature_lower for kw in lipid_keywords):
            pathways['Lipid Metabolism'].append(feature)
            categorized = True
        elif any(kw in feature_lower for kw in aa_keywords):
            pathways['Amino Acid Metabolism'].append(feature)
            categorized = True
        elif any(kw in feature_lower for kw in carb_keywords):
            pathways['Carbohydrate Metabolism'].append(feature)
            categorized = True
        elif any(kw in feature_lower for kw in nucleotide_keywords):
            pathways['Nucleotide Metabolism'].append(feature)
            categorized = True
        elif any(kw in feature_lower for kw in energy_keywords):
            pathways['Energy Metabolism'].append(feature)
            categorized = True
        
        if not categorized:
            pathways['Other'].append(feature)
    
    return pathways

# Analyze pathways for differential features
if len(significant_features) > 0:
    pathway_analysis = simple_pathway_analysis(significant_features['feature'].unique())
    
    # Create pathway summary
    pathway_summary = pd.DataFrame([
        {'pathway': pathway, 'count': len(features), 'features': ', '.join(features[:10])}
        for pathway, features in pathway_analysis.items()
    ])
    pathway_summary.to_csv(f"{output_dir}/results/pathway_enrichment.csv", index=False)
    
    # Plot pathway distribution
    plt.rcParams['font.family'] = 'Times New Roman'
    fig, ax = plt.subplots(figsize=(10, 6.5))
    pathway_counts = {k: len(v) for k, v in pathway_analysis.items()}
    ax.barh(list(pathway_counts.keys()), list(pathway_counts.values()), color='steelblue')
    ax.set_xlabel('Number of Features', fontsize=20)
    ax.set_title('Pathway Distribution of Differential Features', fontsize=24, fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=18)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/plots/pathway_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()

# ============================================================================
# 7. MACHINE LEARNING MODELS
# ============================================================================
print("\n[7/8] Building machine learning models...")

# Prepare data for ML
X = combined_scaled.values
y_survival = survival_clusters['event'].values if 'event' in survival_clusters.columns else None
y_subtype = nmf_labels

# Subtype prediction model (if we have external validation)
# For now, we'll do cross-validation
if len(set(y_subtype)) > 1:
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_subtype, test_size=0.3, random_state=42, stratify=y_subtype
    )
    
    # Random Forest
    rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
    rf.fit(X_train, y_train)
    rf_score = rf.score(X_test, y_test)
    print(f"  Random Forest accuracy: {rf_score:.3f}")
    
    # Feature importance
    feature_importance = pd.DataFrame({
        'feature': combined_scaled.columns,
        'importance': rf.feature_importances_
    }).sort_values('importance', ascending=False)
    
    feature_importance.head(50).to_csv(f"{output_dir}/results/feature_importance_rf.csv", index=False)
    
    # Plot top features
    plt.figure(figsize=(10, 8))
    top_features = feature_importance.head(30)
    plt.barh(range(len(top_features)), top_features['importance'].values)
    plt.yticks(range(len(top_features)), top_features['feature'].values)
    plt.xlabel('Importance', fontsize=12)
    plt.title('Top 30 Features by Random Forest Importance', fontsize=14, fontweight='bold')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(f"{output_dir}/plots/feature_importance.png", dpi=300, bbox_inches='tight')
    plt.close()

# Survival prediction model
if y_survival is not None and len(set(y_survival)) > 1:
    # Use Cox regression for survival
    survival_ml_data = survival_clusters[['time', 'event']].copy()
    survival_ml_data = survival_ml_data.dropna()
    
    if len(survival_ml_data) > 10:
        # Select top features for Cox model
        if len(significant_features) > 0:
            top_survival_features = significant_features.nsmallest(20, 'fdr')['feature'].unique()
            top_survival_features = [f for f in top_survival_features if f in combined_scaled.columns]
            
            if len(top_survival_features) > 0:
                cox_data = combined_scaled[top_survival_features].copy()
                cox_data = cox_data.merge(survival_ml_data, left_index=True, right_index=True, how='inner')
                
                if len(cox_data) > 10:
                    try:
                        cph = CoxPHFitter()
                        cph.fit(cox_data, duration_col='time', event_col='event')
                        
                        # Save Cox model results
                        cph.summary.to_csv(f"{output_dir}/results/cox_regression_results.csv")
                        
                        # Plot hazard ratios
                        plt.figure(figsize=(10, max(6, len(top_survival_features)//3)))
                        hazard_ratios = cph.hazard_ratios_
                        ci_lower = cph.confidence_intervals_.iloc[:, 0]
                        ci_upper = cph.confidence_intervals_.iloc[:, 1]
                        
                        y_pos = np.arange(len(hazard_ratios))
                        err_lower = np.maximum(0, hazard_ratios.values - ci_lower.values)
                        err_upper = np.maximum(0, ci_upper.values - hazard_ratios.values)
                        plt.errorbar(hazard_ratios.values, y_pos, 
                                   xerr=[err_lower, err_upper],
                                   fmt='o', capsize=5, capthick=2)
                        plt.axvline(x=1, color='r', linestyle='--', alpha=0.5)
                        plt.yticks(y_pos, hazard_ratios.index)
                        plt.xlabel('Hazard Ratio', fontsize=12)
                        plt.title('Cox Regression: Hazard Ratios for Top Features', 
                                fontsize=14, fontweight='bold')
                        plt.grid(True, alpha=0.3)
                        plt.tight_layout()
                        plt.savefig(f"{output_dir}/plots/cox_hazard_ratios.png", dpi=300, bbox_inches='tight')
                        plt.close()
                    except Exception as e:
                        print(f"  Cox regression failed: {e}")

# ============================================================================
# 8. RECURRENCE ANALYSIS
# ============================================================================
print("\n[8/8] Analyzing recurrence patterns...")

if 'recurrence' in survival_clusters.columns:
    # Compare metabolites/lipids between recurrence and no recurrence
    recurrence_data = survival_clusters[survival_clusters['recurrence'].notna()].copy()
    
    if len(recurrence_data) > 0:
        # Get metabolite and lipid data for recurrence analysis
        recurrence_metabolites = metabolome_matrix.loc[recurrence_data.index]
        recurrence_lipids = lipidome_combined.loc[recurrence_data.index]
        
        # Differential analysis: recurrence vs no recurrence
        recurrence_results = []
        
        for feature in list(recurrence_metabolites.columns) + list(recurrence_lipids.columns):
            if feature in recurrence_metabolites.columns:
                values = recurrence_metabolites[feature]
            else:
                values = recurrence_lipids[feature]
            
            recur_values = values[recurrence_data['recurrence'] == 1].dropna()
            no_recur_values = values[recurrence_data['recurrence'] == 0].dropna()
            
            if len(recur_values) >= 3 and len(no_recur_values) >= 3:
                stat, pval = stats.mannwhitneyu(recur_values, no_recur_values, alternative='two-sided')
                
                pooled_std = np.sqrt(((len(recur_values) - 1) * recur_values.std()**2 + 
                                    (len(no_recur_values) - 1) * no_recur_values.std()**2) / 
                                   (len(recur_values) + len(no_recur_values) - 2))
                cohens_d = (recur_values.mean() - no_recur_values.mean()) / pooled_std if pooled_std > 0 else 0
                
                recurrence_results.append({
                    'feature': feature,
                    'mean_recurrence': recur_values.mean(),
                    'mean_no_recurrence': no_recur_values.mean(),
                    'fold_change': recur_values.mean() / no_recur_values.mean() if no_recur_values.mean() != 0 else np.inf,
                    'pvalue': pval,
                    'cohens_d': cohens_d
                })
        
        if len(recurrence_results) > 0:
            recurrence_df = pd.DataFrame(recurrence_results)
            _, recurrence_df['fdr'], _, _ = multipletests(recurrence_df['pvalue'], method='fdr_bh')
            
            recurrence_sig = recurrence_df[
                (recurrence_df['fdr'] < 0.05) & (abs(recurrence_df['cohens_d']) > 0.5)
            ].sort_values('fdr')
            
            print(f"Found {len(recurrence_sig)} features associated with recurrence (FDR<0.05)")
            
            recurrence_sig.to_csv(f"{output_dir}/results/recurrence_associated_features.csv", index=False)
            
            # Plot top recurrence-associated features
            if len(recurrence_sig) > 0:
                top_recur = recurrence_sig.head(20)
                
                fig, ax = plt.subplots(figsize=(10, 8))
                y_pos = np.arange(len(top_recur))
                colors = ['red' if fc > 1 else 'blue' for fc in top_recur['fold_change']]
                ax.barh(y_pos, top_recur['fold_change'].values, color=colors, alpha=0.7)
                ax.axvline(x=1, color='black', linestyle='--', linewidth=1)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(top_recur['feature'].values)
                ax.set_xlabel('Fold Change (Recurrence / No Recurrence)', fontsize=12)
                ax.set_title('Top 20 Recurrence-Associated Features', fontsize=14, fontweight='bold')
                ax.grid(True, alpha=0.3, axis='x')
                plt.tight_layout()
                plt.savefig(f"{output_dir}/plots/recurrence_features.png", dpi=300, bbox_inches='tight')
                plt.close()

# ============================================================================
# SUMMARY REPORT
# ============================================================================
print("\n" + "="*80)
print("ANALYSIS COMPLETE!")
print("="*80)
print(f"\nResults saved to: {output_dir}")
print("\nGenerated files:")
print("  - plots/: All visualization figures")
print("  - results/: All statistical results and tables")
print("  - networks/: Network analysis results")
print("\nKey findings:")
print(f"  - {len(lipid_metabolite_corr)} significant lipid-metabolite correlations")
print(f"  - {optimal_k} metabolic-lipidomic subtypes identified")
print(f"  - {len(significant_features)} differential features between subtypes")
if 'recurrence_sig' in locals():
    print(f"  - {len(recurrence_sig)} recurrence-associated features")
print("\n" + "="*80)

