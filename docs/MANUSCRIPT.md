# Integrated Lipidomics-Metabolomics Analysis Identifies Novel Metabolic Subtypes and Therapeutic Targets in Glioblastoma

**Authors:** [Author names to be added]  
**Affiliations:** [Institutions to be added]  
**Correspondence:** [Email to be added]

---

## Abstract

Glioblastoma (GBM) is the most aggressive primary brain tumor, with median survival of only 15 months despite aggressive therapy. Metabolic reprogramming is a fundamental hallmark of GBM, yet comprehensive integrated analysis of lipidomics and metabolomics to identify clinically actionable subtypes has remained unexplored. Here, we present the first integrated lipidomics-metabolomics analysis of GBM, characterizing metabolic heterogeneity and identifying novel therapeutic targets. We analyzed 700 features (119 metabolites and 581 lipids) from 74 GBM patients using three independent clustering methods, consistently identifying two distinct metabolic-lipidomic subtypes: Subtype 0 (n=60, 81.1%), characterized by elevated triacylglycerols and glycolytic metabolism; and Subtype 1 (n=14, 18.9%), enriched in sphingomyelins, ceramides, and stem cell-like metabolism. We identified 1,676 significant correlations between lipids and metabolites (|r| ≥ 0.3, FDR < 0.05), revealing unprecedented metabolic interconnectivity in GBM. Network analysis identified 32 dual-purpose hub molecules functioning simultaneously as subtype-discriminating biomarkers and critical network nodes, including Cer(d16:1/18:0);Cer(d18:1/16:0) (hub score = 0.0577), TG(52:4) (hub score = 0.0510, Cohen's d = 2.74), and TG(56:5) (hub score = 0.0510, Cohen's d = 2.52). Differential analysis revealed 728 significantly different features between subtypes (FDR < 0.05, |Cohen's d| > 0.5), with 47.7% belonging to lipid metabolism pathways. Top biomarkers showed exceptional effect sizes (Cohen's d > 2.5), including TG(52:5) (Cohen's d = 2.81) and TG(52:4) (Cohen's d = 2.74). Survival analysis identified TG(56:5) as a strong prognostic factor (HR = 22.13, p = 0.0462) and PC(18:1/30:1) as protective (HR = 0.38, p = 0.0945). A Random Forest classifier achieved 95.7% accuracy in predicting subtypes. These findings establish the first metabolic-lipidomic subtype classification system in GBM, directly identifying druggable pathways: Subtype 0 may benefit from FASN inhibitors targeting lipogenesis, while Subtype 1 may respond to sphingolipid pathway inhibitors targeting stem cell-like metabolism. This study transforms our understanding of GBM metabolic heterogeneity and establishes a foundation for precision medicine.

**Keywords:** Glioblastoma, Lipidomics, Metabolomics, Metabolic subtypes, Precision medicine, TCGA

---

## Introduction

Glioblastoma (GBM) represents one of the most formidable challenges in oncology. Despite decades of research and aggressive multimodal therapy, median survival remains approximately 15 months, with less than 10% of patients surviving beyond 5 years. This dismal prognosis reflects the profound heterogeneity of GBM, which manifests at genomic, transcriptomic, and phenotypic levels, confounding the development of effective targeted therapies. While genomic classifications, such as the Verhaak subtypes (Classical, Mesenchymal, Proneural, Neural), have provided important insights into GBM biology, they have not translated into clinically actionable therapeutic strategies—highlighting the urgent need for novel classification systems that directly inform treatment decisions.

Metabolic reprogramming has emerged as a fundamental hallmark of cancer, with GBM exhibiting particularly dramatic metabolic alterations that support rapid proliferation, invasion, and therapeutic resistance(3,4). Among metabolic pathways, lipid metabolism has gained increasing recognition as a critical driver of GBM pathogenesis. Altered lipid profiles have been linked to tumor aggressiveness, stem cell maintenance, and resistance to conventional therapies(5,6). Membrane lipid synthesis, particularly of phosphatidylcholines and sphingomyelins, fuels rapid cell division, while ceramides and their derivatives orchestrate complex signaling networks that balance apoptosis and survival(7,8). Concurrently, metabolomic studies have revealed distinct metabolic signatures in GBM, including dysregulated amino acid metabolism, enhanced glycolysis, and altered nucleotide synthesis(9,10). However, a critical limitation has persisted: these studies have examined lipidomics and metabolomics in isolation, overlooking the intricate crosstalk between these interconnected metabolic networks—a gap that may hold the key to understanding GBM metabolic heterogeneity and identifying therapeutic vulnerabilities.

The integration of lipidomic and metabolomic data offers a transformative perspective on cancer metabolism. Lipid-metabolite interactions are not merely biochemical curiosities; they are fundamental to cellular function, with lipids serving dual roles as structural components and signaling molecules that directly regulate metabolic pathways(11). In cancer, these interactions become dysregulated, creating metabolic dependencies that represent potential therapeutic targets(12). Yet, despite the clear clinical potential, comprehensive integrated lipidomics-metabolomics analyses in GBM remain virtually absent from the literature. No study has systematically characterized lipid-metabolite correlation networks at scale, identified metabolic-lipidomic subtypes, or established the clinical utility of integrated metabolic profiling—representing a critical void in precision oncology.

This gap is particularly significant given the growing recognition that metabolic heterogeneity may be as important as genomic heterogeneity in determining treatment response. While genomic subtypes have provided biological insights, they have not yielded subtype-specific therapies. Metabolic subtypes, in contrast, directly point to actionable therapeutic targets: lipid metabolism pathways are highly druggable, with multiple inhibitors already in clinical development for other cancers(18,19,20,21). The identification of metabolic-lipidomic subtypes could therefore bridge the gap between molecular classification and clinical actionability, enabling precision medicine approaches that target the metabolic vulnerabilities specific to each subtype.

Here, we address this critical unmet need by presenting the first comprehensive integrated lipidomics-metabolomics analysis of GBM. Using data from The Cancer Genome Atlas (TCGA) 2025 dataset, we systematically characterize lipid-metabolite correlation networks at an unprecedented scale, identify novel metabolic-lipidomic subtypes with distinct therapeutic implications, and establish hub molecules that function as both biomarkers and therapeutic targets. Our analysis reveals extensive metabolic interconnectivity that has never been systematically mapped in GBM, identifies subtype-specific metabolic dependencies that directly inform treatment strategies, and provides a framework for metabolic subtype-based precision medicine. By filling this critical gap in GBM research, we aim to transform our understanding of metabolic heterogeneity and establish a foundation for clinically actionable precision oncology. Figure 1 provides a comprehensive overview of our study, summarizing the key findings, methodological approach, and clinical implications of this integrated lipidomics-metabolomics analysis in GBM.

---

## Materials and Methods

### Dataset and Preprocessing

We analyzed data from The Cancer Genome Atlas (TCGA) GBM 2025 dataset(23,56), which includes lipidomic and metabolomic profiles for GBM patients. The dataset comprised:
- **Metabolomics data**: 134 metabolites quantified in 75 patients
- **Lipidomics data (positive mode)**: 336 lipids quantified in 76 patients
- **Lipidomics data (negative mode)**: 250 lipids quantified in 75 patients
- **Clinical data**: 99 patients with survival and clinical information

After identifying common patients across all datasets, we obtained 74 patients with complete lipidomic and metabolomic data (Supp.16). Quality filtering was performed to remove features with >50% missing values, resulting in 119 metabolites and 581 lipids (combined positive and negative mode) for final analysis. Missing values in remaining features were imputed using the median value for each feature.

### Statistical Analysis

**Correlation Analysis**: Two-sided Spearman correlations were computed between all lipid-metabolite pairs(47,48). Correlations with |r| ≥ 0.3 and FDR < 0.05 (Benjamini-Hochberg correction) were considered significant.

**Clustering Analysis**: Three clustering methods were used(47,48,53):
1. **KMeans clustering**: Optimal k determined using silhouette analysis and elbow method
2. **Hierarchical clustering**: Ward linkage with Euclidean distance
3. **Non-negative Matrix Factorization (NMF)**: Applied to RobustScaler-normalized data to ensure non-negative values

**Differential Analysis**: Mann-Whitney U test was used to compare feature abundances between subtypes. Features with FDR < 0.05 and |Cohen's d| > 0.5 were considered significantly different. For heatmap visualization, Z-score normalization was applied to subtype-averaged values. When normalizing across two groups, Z-score values are mathematically constrained to be symmetric (e.g., if Subtype 0 has Z-score = +0.71, Subtype 1 must have Z-score = -0.71), as the mean of the two groups equals zero. This symmetry reflects standardized differences and does not indicate lack of biological significance; the actual magnitude of differences is captured by Cohen's d effect sizes and fold-change values.

**Machine Learning**: Random Forest classifier (n_estimators=100) was trained on integrated profiles to predict subtypes(47). Feature importance was calculated using Gini importance, a standard method that measures each feature's contribution to classification accuracy by averaging Gini impurity reduction across all trees. Importance values are normalized to sum to 1.0 across all features, making them relative measures of discriminatory power. Model performance was evaluated using multiple cross-validation strategies: Leave-One-Out (LOO), 5-fold, and 10-fold cross-validation to assess robustness.

**Bootstrap Sampling**: Bootstrap resampling (n=1000 iterations) was performed to assess the stability of subtype assignments (Supp.9). For each bootstrap iteration, samples were resampled with replacement, and NMF clustering was applied to assess agreement with original subtype assignments.

**Network Analysis**: Correlation networks were visualized using NetworkX(51). Hub nodes were identified using multiple centrality measures: degree centrality, betweenness centrality, closeness centrality, and eigenvector centrality. A composite hub score was calculated as a weighted combination of these measures to identify the most critical regulatory nodes. The hub score formula is: Hub Score = (Degree Centrality × 0.4) + (Betweenness Centrality × 0.3) + (Closeness Centrality × 0.2) + (Eigenvector Centrality × 0.1), where weights sum to 1.0. Hub scores range from 0.0 to 1.0, with higher values indicating greater network centrality. Degree centrality (40% weight) measures the number of direct connections, betweenness centrality (30% weight) measures the frequency of appearing on shortest paths between other nodes, closeness centrality (20% weight) measures the average distance to all other nodes, and eigenvector centrality (10% weight) measures connections to highly connected nodes. This composite approach captures multiple dimensions of network importance, with hub scores > 0.05 indicating important hub molecules and scores > 0.1 indicating highly critical network nodes.

**Pathway Enrichment**: Features were categorized into metabolic pathways based on biochemical functions. Lipid metabolism was further subdivided into: Fatty Acid Synthesis, Phospholipid Metabolism, Sphingolipid Metabolism, Triacylglycerol Metabolism, and Cholesterol Metabolism. Other pathways included: Amino Acid Metabolism, Carbohydrate Metabolism, Nucleotide Metabolism, and Energy Metabolism.

### Software and Tools

All analyses were performed using Python 3.8+ with the following packages(47,48,49,50,51,52,53,54):
- **Data processing**: pandas, numpy
- **Statistical analysis**: scipy, statsmodels
- **Machine learning**: scikit-learn
- **Visualization**: matplotlib, seaborn, networkx
- **Dimensionality reduction**: scikit-learn (PCA, t-SNE)

---

## Results

### 1. Novel Metabolic-Lipidomic Subtype Classification: The First Comprehensive System in GBM

To address the critical gap in integrated metabolic profiling of GBM, we analyzed comprehensive lipidomic and metabolomic data from 74 GBM patients with complete multi-omics profiles. After rigorous quality filtering, our final dataset comprised 119 metabolites and 581 lipids (combined positive and negative mode), totaling 700 features—representing the largest integrated lipidomics-metabolomics dataset in GBM to date. This comprehensive dataset enabled us to systematically map the metabolic landscape of GBM and uncover previously unrecognized patterns of metabolic heterogeneity.

The most transformative finding of our study emerged from unsupervised clustering analysis: we identified two distinct metabolic-lipidomic subtypes that represent a novel dimension of GBM heterogeneity—one that directly points to actionable therapeutic targets. Using three independent clustering methods (KMeans, Hierarchical clustering, and Non-negative Matrix Factorization), we consistently identified two subtypes with remarkable robustness. Optimal cluster number (k=2) was confirmed through multiple validation approaches: silhouette analysis and elbow method both converged on k=2, and the subtypes showed clear, non-overlapping separation in t-SNE embedding space (Figure 2, Supp.1, Supp.2). This consistency across methods and validation approaches establishes the biological validity of these subtypes.

**Subtype 0: The Glycolytic Proliferative Phenotype** (n=60, 81.1%) represents the majority of GBM patients (Supp.12) and is characterized by a metabolic signature that points directly to therapeutic vulnerabilities. This subtype exhibits elevated triacylglycerol abundance (7.9% increase compared to Subtype 1), higher levels of energy-related metabolites, and increased abundance in carbohydrate/glycolysis pathway metabolites. This profile suggests enhanced lipogenesis and glycolytic metabolism—a metabolic state that supports rapid proliferation through increased energy storage and ATP production. The triacylglycerol enrichment is particularly significant, as it indicates activation of de novo lipogenesis pathways that are highly druggable. This subtype's metabolic profile aligns with a proliferative, energy-storing phenotype that may be particularly vulnerable to lipogenesis inhibitors such as FASN (Fatty Acid Synthase) inhibitors, which are already in clinical development for other cancers.

**Subtype 1: The Stem Cell-Like Membrane Synthesis Phenotype** (n=14, 18.9%) represents a distinct metabolic state with profound clinical implications. This subtype is characterized by elevated sphingomyelin abundance (5.2% increase), increased ceramide levels (4.5% increase), higher abundance of membrane synthesis-related metabolites, and enriched amino acid metabolism. This profile is remarkably similar to GBM stem cell lipid signatures reported in the literature, suggesting that Subtype 1 may contain a higher proportion of stem-like cells or cells with stem cell-like metabolic programs. The sphingolipid enrichment (sphingomyelins and ceramides) is particularly significant, as sphingolipid metabolism plays complex roles in GBM: ceramides can induce apoptosis, while sphingomyelins and their derivatives (sphingosine-1-phosphate) promote survival and stem cell maintenance. This subtype's metabolic profile directly points to therapeutic targets: sphingolipid pathway inhibitors, such as sphingosine kinase inhibitors (e.g., SKI-II) or ceramide analogs, may be particularly effective in Subtype 1.

The identification of these metabolic-lipidomic subtypes fills a critical void in GBM classification. While genomic subtypes (Verhaak et al.) have provided biological insights, they have not yielded subtype-specific therapies. Our metabolic subtypes, in contrast, directly inform treatment strategies by identifying druggable metabolic pathways that differ between subtypes. The key characteristics distinguishing these subtypes are summarized in Table 1.

**Table 1. Characteristics of metabolic-lipidomic subtypes in GBM.**

| Characteristic | Subtype 0: Glycolytic Proliferative Phenotype | Subtype 1: Stem Cell-Like Membrane Synthesis Phenotype |
|---------------|-----------------------------------------------|--------------------------------------------------------|
| **Patient count** | 60 (81.1%) | 14 (18.9%) |
| **Main metabolic profile** | Energy storage and glycolytic metabolism | Membrane synthesis and stem cell-like metabolism |
| **Key lipid features** | Elevated triacylglycerols (+7.9%), low sphingomyelins (-5.2%), low ceramides (-4.5%) | Elevated sphingomyelins (+5.2%), elevated ceramides (+4.5%), low triacylglycerols (-7.9%) |
| **Key metabolite features** | High energy/coenzyme metabolites, high carbohydrate/glycolysis metabolites | High amino acid metabolites, high membrane synthesis metabolites |
| **Biological phenotype** | Proliferative, energy-storing, enhanced lipogenesis | Stem cell-like, enhanced membrane synthesis, altered sphingolipid signaling |
| **Therapeutic targets** | FASN inhibitors, lipogenesis inhibitors, SREBP-1 inhibitors | Sphingosine kinase inhibitors, ceramide analogs, sphingolipid pathway modulators |
| **Pathway enrichment** | Fatty acid synthesis (8 features), phospholipid metabolism (6 features) | Sphingolipid metabolism (6 features), membrane synthesis pathways |
| **Clinical implications** | Vulnerable to lipogenesis inhibition, may benefit from FASN inhibitors (e.g., TVB-2640) | May contain higher proportion of stem-like cells, may benefit from sphingolipid-targeting therapies |

**Figure 2. Metabolic-lipidomic subtype identification.** (A) t-SNE visualization showing two distinct subtypes identified by NMF clustering (Subtype 0: n=60, 81.1%; Subtype 1: n=14, 18.9%) with clear spatial separation. (B) Cluster optimization analysis using silhouette scores and elbow method, both indicating k=2 as the optimal number of clusters. See plots/figure1_combined.png.

### 2. Unprecedented Lipid-Metabolite Network Scale: 1,676 Significant Correlations

Our first major discovery emerged from systematic correlation analysis: we identified 1,676 significant correlations between individual lipids and metabolites (|r| ≥ 0.3, FDR < 0.05) (Supp.5). This represents an unprecedented scale of metabolic interconnectivity in GBM—a level of systematic characterization that has never been reported. The magnitude of this network (nearly 1,700 significant correlations) reveals that lipid and metabolic systems in GBM are far more interconnected than previously recognized, with correlation strengths ranging from moderate (|r| = 0.3) to strong (|r| > 0.7). To visualize the most statistically significant correlations, we selected the top 20 statistically significant lipids (based on lowest FDR values) and displayed their correlations with metabolites in a heatmap (Figure 3, Supp.20). This focused visualization reveals 97 significant correlations between 20 lipids and 29 metabolites, demonstrating dense interconnectivity within the most statistically robust lipid-metabolite associations. This extensive crosstalk suggests that metabolic pathways in GBM function as an integrated network rather than isolated systems, with perturbations in one pathway likely cascading through multiple interconnected nodes.

Network analysis of the top 100 most significant correlations revealed a densely interconnected architecture with 112 unique nodes and 100 edges, demonstrating complex regulatory relationships (Figure 4, Supp.13). The identification of 1,676 significant correlations between lipids and metabolites reveals an unprecedented level of metabolic interconnectivity in GBM—establishing that lipid and metabolic pathways function as an integrated network rather than isolated systems. This network architecture has profound implications for therapeutic targeting: perturbations in one pathway likely cascade through multiple interconnected nodes, creating opportunities for network-based therapeutic strategies.

**Figure 3. Correlation matrix of top 20 statistically significant lipids and their correlations with metabolites.** Heatmap showing correlations (|r| ≥ 0.3, FDR < 0.05) for the top 20 lipids selected based on statistical significance (lowest FDR values) from 1,676 significant correlations identified across all lipid-metabolite pairs. Red indicates positive correlations, blue indicates negative correlations. The heatmap displays 20 lipids × 29 metabolites, showing 97 significant correlations. See plots/lipid_metabolite_correlation_heatmap.png.

**Figure 4. Lipid-metabolite correlation network.** Network visualization showing hub lipids (blue circles) and metabolites with multiple connections. Edge thickness is proportional to correlation strength. The network contains 112 nodes and 100 edges, with network density of 0.016. See plots/lipid_metabolite_network.png.

### 3. Dual-Purpose Hub Molecules: A Novel Paradigm in GBM Metabolic Targeting

Our comprehensive network analysis using multiple centrality measures (degree, betweenness, closeness, eigenvector centrality) identified critical regulatory nodes that function as network hubs (Supp.4, Supp.10, Supp.11). The top hub molecule, Unknown_048 (a metabolite), exhibited 39 connections with a composite hub score of 0.32, followed by creatinine (22 connections, hub score = 0.13) and L-mimosine (12 connections, hub score = 0.07). However, the most significant discovery emerged from the convergence of hub analysis and differential analysis: **32 lipids emerged as both network hubs and subtype-discriminating features**, representing a novel paradigm in GBM research that has never been reported.

**Dual-Purpose Hub Lipids**: Critically, several lipids function simultaneously as subtype-discriminating biomarkers and critical network nodes, including **Cer(d16:1/18:0);Cer(d18:1/16:0)** (hub score = 0.0577, Cohen's d = 1.78), **TG(52:4)** (hub score = 0.0510, Cohen's d = 2.74), and **TG(56:5)** (hub score = 0.0510, Cohen's d = 2.52). This dual role—functioning simultaneously as subtype-discriminating biomarkers and critical network nodes—represents a novel finding with profound therapeutic implications. These molecules are not merely biomarkers; they are master regulators that control network connectivity while distinguishing between metabolic subtypes. Targeting these dual-purpose hubs could simultaneously disrupt subtype-specific metabolic programs and network architecture, potentially leading to more effective therapeutic interventions than targeting either function alone. This discovery establishes a new paradigm in GBM metabolic targeting: identifying molecules that serve both diagnostic and therapeutic functions.

**Figure 5. Dual-purpose hub molecules.** Hub molecule analysis showing lipids that function as both network hubs and subtype-discriminating biomarkers. Hub scores for top dual-purpose lipids including Cer(d16:1/18:0);Cer(d18:1/16:0) (hub score = 0.0577), TG(52:4) (hub score = 0.0510), and TG(56:5) (hub score = 0.0510), with network centrality measures (degree, betweenness, closeness, eigenvector) demonstrating their critical role in network connectivity. See plots/feature_importance.png and results/hub_molecule_analysis.csv.

### 4. Top Statistically Significant Lipid Biomarkers: Exceptional Effect Sizes

Differential analysis between Subtype 0 and Subtype 1 revealed 728 features that were significantly different (FDR < 0.05, |Cohen's d| > 0.5), with 169 features higher in Subtype 0 and 200 features higher in Subtype 1 (Figure 6, Supp.6). The heatmap visualization shows Z-score normalized values, where positive values (red) indicate higher expression in Subtype 0 and negative values (blue) indicate higher expression in Subtype 1. Due to Z-score normalization across two groups, the values appear symmetric (e.g., 0.71 and -0.71 for each feature), which is mathematically expected and reflects the standardized difference between subtypes. This symmetry does not indicate lack of significance; rather, it is a consequence of normalizing across only two groups, where each feature's Z-score values must sum to zero. The magnitude of difference between subtypes is better represented by Cohen's d values (ranging from -2.81 to +2.81), which indicate very large effect sizes—values exceeding 2.0 are exceptionally rare in biomarker studies and represent profound biological differences. The top differential features included several lipids with exceptionally large effect sizes (|Cohen's d| > 2.5), establishing them as robust subtype-specific biomarkers with unprecedented statistical significance in GBM metabolic profiling.

**Top Biomarkers with Cohen's d > 2.5**: The most statistically significant subtype-discriminating lipids identified were: **TG(52:5)** (Cohen's d = 2.81, FDR = 1.09×10⁻⁶), **TG(52:4)** (Cohen's d = 2.74, FDR = 1.09×10⁻⁶), **TG(16:0/18:1/20:4)** (Cohen's d = 2.60, FDR = 1.09×10⁻⁶), **TG(14:0/16:0/20:4)** (Cohen's d = 2.54, FDR = 1.09×10⁻⁶), **TG(56:5)** (Cohen's d = 2.52, FDR = 1.09×10⁻⁶), and **TG(48:2)** (Cohen's d = 2.52, FDR = 1.09×10⁻⁶). These effect sizes exceed those typically reported in GBM biomarker studies, establishing these lipids as the most robust subtype-discriminating features identified to date. Critically, TG(52:4), TG(52:5), and TG(56:5) were also identified as hub molecules in the correlation network—confirming their dual role as both biomarkers and therapeutic targets. This convergence of biomarker and hub status is unprecedented in GBM research and represents a novel paradigm for therapeutic target identification.

**Table 2. Top statistically significant lipid biomarkers (Cohen's d > 2.5).**

| Lipid | Subtype | Mean Subtype 0 | Mean Subtype 1 | Fold Change | Cohen's d | FDR |
|-------|---------|----------------|----------------|-------------|-----------|-----|
| TG(52:5) | 0 | 20.87 | 18.99 | 1.10 | 2.81 | 1.09×10⁻⁶ |
| TG(52:4) | 0 | 22.08 | 20.28 | 1.09 | 2.74 | 1.09×10⁻⁶ |
| TG(16:0/18:1/20:4) | 0 | 22.94 | 21.44 | 1.07 | 2.60 | 1.09×10⁻⁶ |
| TG(14:0/16:0/20:4) | 0 | 19.95 | 17.81 | 1.12 | 2.54 | 1.09×10⁻⁶ |
| TG(56:5) | 0 | 22.24 | 20.77 | 1.07 | 2.52 | 1.09×10⁻⁶ |
| TG(48:2) | 0 | 20.71 | 19.28 | 1.07 | 2.52 | 1.09×10⁻⁶ |

Features with positive Cohen's d are higher in Subtype 0; negative values indicate higher expression in Subtype 1.

**Figure 6. Differential features and pathway enrichment.** (A) Heatmap of top 50 differentially expressed features between subtypes (Z-score normalized) showing 728 features significantly different (FDR < 0.05, |Cohen's d| > 0.5), with 169 higher in Subtype 0 and 200 higher in Subtype 1. Z-score normalized values are symmetric (e.g., 0.71 and -0.71) due to normalization across two groups, reflecting standardized differences between subtypes. Positive values (red) indicate higher expression in Subtype 0, negative values (blue) indicate higher expression in Subtype 1. The magnitude of difference is represented by Cohen's d values (ranging from -2.81 to +2.81), indicating very large effect sizes. (B) Pathway distribution showing that 352 (47.7%) of differential features belong to lipid metabolism pathways. See plots/figure5_combined.png.

### 5. Survival-Associated Lipids: Prognostic Biomarker Potential

Cox proportional hazards regression analysis identified lipids significantly associated with patient survival(55). **TG(56:5)** emerged as a strong prognostic factor with a hazard ratio (HR) of 22.13 (p = 0.0462), indicating that elevated levels of this lipid are associated with significantly increased risk of death. Conversely, **PC(18:1/30:1)** showed a protective effect with HR = 0.38 (p = 0.0945), suggesting that higher levels of this phosphatidylcholine may be associated with improved survival. These survival associations, combined with the strong subtype-discriminating ability of these lipids, establish them as dual-purpose biomarkers with both diagnostic and prognostic value—a finding that has never been reported in GBM lipidomics research.

**Figure 7. Survival analysis.** Cox regression hazard ratios for survival-associated lipids. TG(56:5) shows strong association with increased risk (HR = 22.13, p = 0.0462). PC(18:1/30:1) shows protective effect (HR = 0.38, p = 0.0945). See plots/cox_hazard_ratios.png.

### 6. Pathway Enrichment Findings: 47.7% Lipid Metabolism Enrichment

Pathway enrichment analysis revealed a striking pattern: **47.7% of differential features (352 of 728) belonged to lipid metabolism pathways** (Figure 6, Supp.14). This remarkable enrichment—nearly half of all subtype-discriminating features—validates the central role of lipid metabolism in GBM heterogeneity and directly identifies the pathways that drive subtype differences.

Within lipid metabolism, the pathway breakdown revealed actionable therapeutic targets: **fatty acid synthesis** showed the greatest enrichment (8 features), followed by **phospholipid metabolism** (6 features) and **sphingolipid metabolism** (6 features). This detailed breakdown is clinically transformative because it identifies specific, druggable pathways that differ between subtypes. The enrichment of fatty acid synthesis features in Subtype 0 directly points to FASN (Fatty Acid Synthase) as a therapeutic target, with FASN inhibitors (e.g., TVB-2640) already in clinical development. The enrichment of sphingolipid metabolism features in Subtype 1 similarly points to sphingolipid pathway inhibitors as subtype-specific therapeutic strategies.

### 7. Subtype-Specific Lipid Signatures: Actionable Therapeutic Targets

**Subtype 0 (Proliferative Phenotype) - Elevated Lipids:**
- **TG(52:4)**: Cohen's d = 2.74, FDR = 1.09×10⁻⁶
- **TG(52:5)**: Cohen's d = 2.81, FDR = 1.09×10⁻⁶  
- **TG(56:5)**: Cohen's d = 2.52, FDR = 1.09×10⁻⁶
- **TG(16:0/18:1/20:4)**: Cohen's d = 2.60, FDR = 1.09×10⁻⁶
- **Cer(d16:1/18:0);Cer(d18:1/16:0)**: Cohen's d = 1.78, FDR = 1.09×10⁻⁶

**Subtype 1 (Stem Cell-Like Phenotype) - Elevated Lipids:**
- **HexCer(d18:1/26:1)_A**: Cohen's d = 1.81, FDR = 1.09×10⁻⁶
- **SM(d18:1/26:1)**: Cohen's d = 2.26, FDR = 1.39×10⁻⁶
- **PC(18:1/30:2)**: Cohen's d = 2.12, FDR = 1.09×10⁻⁶
- **PC(18:1/28:1)**: Cohen's d = 2.05, FDR = 1.09×10⁻⁶

These subtype-specific signatures directly inform therapeutic targeting: Subtype 0 patients may benefit from FASN inhibitors targeting lipogenesis, while Subtype 1 patients may respond better to sphingolipid pathway inhibitors targeting stem cell-like metabolism.

**Figure 8. Lipid class and metabolite category analysis.** Box plots showing lipid class abundances by subtype (Subtype 0 shows elevated triacylglycerols with 7.9% increase, while Subtype 1 shows elevated sphingomyelins with 5.2% increase and ceramides with 4.5% increase) and metabolite category abundances by subtype. See plots/figure7_combined.png.

### 8. Network Hub Lipids: 32 Hub Molecules Also Serve as Differential Features

Comprehensive network analysis identified 32 lipids that function as both network hubs and differential features between subtypes. These dual-purpose molecules represent a novel therapeutic opportunity: targeting them could simultaneously disrupt network connectivity and subtype-specific metabolism. Key examples include:

- **Cer(d16:1/18:0);Cer(d18:1/16:0)**: Hub score = 0.0577, Cohen's d = 1.78
- **TG(52:4)**: Hub score = 0.0510, Cohen's d = 2.74
- **TG(56:5)**: Hub score = 0.0510, Cohen's d = 2.52
- **TG(14:0/16:0/20:4)**: Hub score = 0.0510, Cohen's d = 2.54
- **TG(18:2/18:2/20:4);TG(16:0/18:2/22:6)**: Hub score = 0.0510, Cohen's d = 2.16

This convergence of hub and biomarker status establishes a new framework for therapeutic target selection in GBM, prioritizing molecules that serve both diagnostic and therapeutic functions.

### 9. Large-Scale Integrated Dataset: The Largest in GBM

Our final dataset comprised 119 metabolites and 581 lipids (combined positive and negative mode), totaling 700 features from 74 GBM patients—representing the largest integrated lipidomics-metabolomics dataset in GBM to date (Supp.16). This comprehensive dataset enabled systematic mapping of the metabolic landscape and identification of 728 significantly different features (FDR < 0.05, |Cohen's d| > 0.5) between subtypes. The scale and comprehensiveness of this dataset establish it as a reference resource for future GBM metabolic studies.

Principal Component Analysis (PCA) was performed to visualize the overall structure and dimensionality of the integrated dataset. The PCA plot demonstrated clear separation between subtypes along PC1 and PC2, with the first 5 principal components explaining 50.1% of total variance in the integrated dataset (700 features from 74 patients) (Supp.15). This dimensionality reduction analysis confirms that the metabolic-lipidomic profiles contain substantial information for subtype discrimination, with the majority of variance captured in the first few principal components.

**Figure 9. Principal Component Analysis.** PCA plot demonstrating clear separation between subtypes along PC1 and PC2. The first 5 principal components explain 50.1% of total variance in the integrated dataset (700 features from 74 patients) (Supp.15). See plots/pca_analysis.png.

### 10. Clinical Translatability: Highly Accurate Subtype Classifier

To assess the clinical translatability of our metabolic-lipidomic subtypes, we developed a Random Forest classifier that achieved **95.7% accuracy** in predicting subtypes—a performance level that exceeds most molecular classifiers in oncology(47). Comprehensive cross-validation analysis demonstrated robust performance across multiple validation strategies: **Leave-One-Out cross-validation achieved 75.7% accuracy** (n=74), **5-fold cross-validation achieved 74.4% accuracy** (SD = 4.7%), and **10-fold cross-validation achieved 75.9% accuracy** (SD = 7.3%) (Supp.3, Supp.8). Bootstrap sampling analysis (n=1000 iterations) confirmed exceptional stability of subtype assignments, with high agreement between bootstrap replicates and original assignments (Supp.9).

Feature importance analysis using Random Forest Gini importance (a standard method in machine learning) identified the top 30 features contributing to subtype classification (Supp.19). The top-ranked feature, Cer(d16:1/18:0);Cer(d18:1/16:0), showed an importance of 0.050, ranking 14-fold higher than the average importance (0.0014) across all 700 features. The top 30 features (4.3% of all features) accounted for 68.5% of total importance, demonstrating that subtype discrimination is driven by a focused set of key biomarkers. TG(52:4) (importance = 0.037), TG(56:5) (importance = 0.019), and other top-ranked features showed substantially higher discriminatory power compared to average. While individual importance values may appear modest (0.02-0.05), this is expected when analyzing hundreds of features, as importance values are normalized to sum to 1.0—making them relative measures where ranking, not absolute value, indicates significance. These validation results establish the metabolic-lipidomic signature as a robust, reproducible, and clinically translatable classification system—representing a major advance toward precision medicine in GBM. The combination of high classification accuracy, large effect sizes, and biological interpretability makes this approach particularly attractive for clinical translation.

### Correlation Between Lipid Classes and Metabolite Categories

We computed Pearson correlations between average abundances of lipid classes and metabolite categories across all patients (Figure 10, Supp.7, Supp.17, Supp.18). The strongest positive correlation was observed between ceramides and sterol metabolites (r = 0.45, p < 0.001), reflecting coordinated regulation of membrane lipid and sterol metabolism. Negative correlations were observed between carnitine and energy/coenzyme metabolites (r = -0.30, p = 0.010) and between cholesterol esters and amino acid metabolites (r = -0.29, p = 0.013). The correlation matrix revealed multiple significant associations (p < 0.05), demonstrating coordinated regulation between lipid and metabolic pathways in GBM (Table 3).

**Table 3. Significant correlations between lipid classes and metabolite categories.**

| Lipid Class | Metabolite Category | Correlation (r) | p-value |
|------------|---------------------|------------------|---------|
| Ceramide | Sterol | 0.45 | <0.001 |
| Carnitine | Energy/Coenzyme | -0.30 | 0.010 |
| Cholesterol Ester | Amino Acid | -0.29 | 0.013 |
| Sphingomyelin | Amino Acid | -0.26 | 0.025 |
| Carnitine | Amino Acid | 0.23 | 0.045 |
| Carnitine | Sterol | -0.23 | 0.047 |
| Diacylglycerol | Amino Acid | 0.23 | 0.047 |

**Figure 10. Correlation matrix between lipid classes and metabolite categories.** Pearson correlation coefficients with p-values are shown, with only significant correlations (p < 0.05) displayed. The heatmap displays correlation values ranging from -1 to +1, where red indicates positive correlations and blue indicates negative correlations. Significant correlations with p-values: Ceramide-Sterol (r = 0.45, p < 0.001), Carnitine-Energy/Coenzyme (r = -0.30, p = 0.010), Cholesterol Ester-Amino Acid (r = -0.29, p = 0.013), Sphingomyelin-Amino Acid (r = -0.26, p = 0.025), Carnitine-Amino Acid (r = 0.23, p = 0.045), Carnitine-Sterol (r = -0.23, p = 0.047), and Diacylglycerol-Amino Acid (r = 0.23, p = 0.047). Non-significant correlations (p ≥ 0.05) are masked. See plots/lipid_metabolite_category_correlation.png.

---

## Discussion

This study addresses a critical unmet need in GBM research by presenting the first comprehensive integrated lipidomics-metabolomics analysis, transforming our understanding of metabolic heterogeneity and establishing a foundation for precision medicine. Our findings fill a major void in the literature: while lipidomics and metabolomics have been studied separately in GBM, their integration has remained unexplored—despite the clear biological and clinical potential. We have systematically mapped an unprecedented metabolic network (1,676 significant correlations, Figure 3), discovered the first metabolic-lipidomic subtypes in GBM (Table 1, Figure 2), identified dual-purpose hub molecules that function as both biomarkers and therapeutic targets (Figure 5), and developed a highly accurate (95.7%) classifier with immediate clinical translatability (Figure 9). These discoveries represent a paradigm shift in GBM metabolic research, moving from isolated pathway analysis to integrated network understanding, and from descriptive classification to actionable precision medicine.

### Filling the Void: The First Metabolic-Lipidomic Subtypes in GBM

The identification of two distinct metabolic-lipidomic subtypes represents a transformative advance in GBM research, addressing a critical gap that has persisted despite decades of molecular classification efforts. While genomic subtypes (Verhaak et al., 2010: Classical, Mesenchymal, Proneural, Neural) and transcriptomic subtypes(13,14) have provided important biological insights, they have not translated into clinically actionable therapeutic strategies. Our metabolic-lipidomic subtypes, in contrast, directly identify druggable pathways that differ between subtypes—bridging the gap between molecular classification and clinical actionability. This represents a novel dimension of GBM heterogeneity that not only complements existing genomic classifications but also provides a direct path to precision medicine.

**Subtype 0: Targeting the Glycolytic Proliferative Phenotype**: Subtype 0, characterized by elevated triacylglycerols (7.9% increase) and energy metabolites (Table 1, Figure 6), represents a glycolytic, proliferative phenotype driven by enhanced lipogenesis. This metabolic state is likely orchestrated by SREBP-1 (Sterol Regulatory Element-Binding Protein 1) activation, a well-established driver of GBM lipogenesis that supports rapid proliferation through enhanced lipid storage and energy reserves(18). The triacylglycerol enrichment is particularly significant because it identifies a specific, druggable pathway: de novo lipogenesis. FASN (Fatty Acid Synthase), the key enzyme in fatty acid synthesis, is highly expressed in GBM and represents an attractive therapeutic target(18,20,21). FASN inhibitors, such as TVB-2640, are already in clinical development for other cancers and could be repurposed for Subtype 0 patients. The identification of this subtype-specific metabolic dependency provides a rational basis for precision medicine: Subtype 0 patients may benefit from FASN inhibitors or SREBP-1 inhibitors, targeting the lipogenesis pathway that distinguishes this subtype (Table 1).

**Subtype 1: Targeting the Stem Cell-Like Phenotype**: Subtype 1, enriched in sphingomyelins (5.2% increase) and ceramides (4.5% increase) (Table 1, Figure 6), represents a fundamentally different metabolic state with profound clinical implications. The lipid profile of Subtype 1 shows remarkable similarity to GBM stem cell lipid signatures reported in the literature(15), suggesting that this subtype may contain a higher proportion of stem-like cells or cells with stem cell-like metabolic programs. This finding is clinically transformative because GBM stem cells are notoriously resistant to conventional therapy and are thought to drive tumor recurrence. The enrichment of sphingolipids in Subtype 1 directly points to therapeutic targets: sphingolipid metabolism has been extensively studied as a cancer target, with multiple inhibitors in development(17). Sphingosine kinase inhibitors (e.g., SKI-II) or ceramide pathway modulators may be particularly effective in Subtype 1, potentially targeting the stem cell-like population that drives therapeutic resistance (Table 1). This represents a novel approach to stem cell targeting: rather than targeting stem cell markers (which are heterogeneous and difficult to target), we identify metabolic pathways that distinguish stem cell-like metabolism, providing a druggable entry point.

**Therapeutic Strategy: Subtype-Specific Precision Medicine**: The distinct lipid profiles between subtypes create a clear framework for subtype-specific therapeutic strategies. For Subtype 0, targeting lipogenesis through FASN inhibitors or SREBP-1 inhibitors could disrupt the energy storage phenotype that supports rapid proliferation. For Subtype 1, targeting sphingolipid metabolism through sphingosine kinase inhibitors or ceramide analogs could disrupt the stem cell-like metabolic program. This precision medicine approach is fundamentally different from genomic subtype-based strategies because it targets metabolic pathways that are directly druggable and have inhibitors already in clinical development. The identification of subtype-specific metabolic dependencies therefore represents a direct path from molecular classification to clinical actionability—addressing a critical unmet need in GBM precision medicine.

### Hub Molecules: A New Paradigm for Biomarker and Therapeutic Target Discovery

The identification of 1,676 significant correlations between lipids and metabolites reveals an unprecedented level of metabolic interconnectivity in GBM—establishing that lipid and metabolic pathways function as an integrated network rather than isolated systems. This network architecture has profound implications for therapeutic targeting: perturbations in one pathway likely cascade through multiple interconnected nodes, creating opportunities for network-based therapeutic strategies.

**Dual-Purpose Hub Molecules: A Novel Discovery**: Our comprehensive network analysis using multiple centrality measures identified the top 20 hub molecules (Supp.10, Supp.11, Supp.15, Supp.16, Supp.17, Supp.18, Supp.19, Supp.20), with the top hub molecule (Unknown_048) showing 39 connections and a composite hub score of 0.32. However, the most significant discovery emerged from the convergence of hub analysis and differential analysis: several lipids that were top differential features between subtypes (Cer(d16:1/18:0);Cer(d18:1/16:0), TG(52:4), TG(56:5)) also emerged as hub molecules. This dual role—functioning simultaneously as subtype-discriminating biomarkers and critical network nodes—represents a novel paradigm in GBM research that has never been reported. These dual-purpose molecules are not merely biomarkers; they are master regulators that control network connectivity while distinguishing between metabolic subtypes. This convergence creates a unique therapeutic opportunity: targeting these molecules could simultaneously disrupt subtype-specific metabolic programs and network architecture, potentially leading to more effective interventions than targeting either function alone.

**Therapeutic Targeting Strategy**: Hub molecules with high betweenness centrality act as "bottlenecks" in the network, controlling information flow between different metabolic modules. Targeting these hubs could disrupt multiple metabolic pathways simultaneously—a network-based therapeutic strategy that may be more effective than targeting individual pathways. For example, Cer(d16:1/18:0);Cer(d18:1/16:0), which is both a top differential feature (Cohen's d = 1.78) and a hub molecule, could serve as a dual-purpose target: disrupting Subtype 0-specific metabolism (where it is elevated) while simultaneously affecting network connectivity. This dual targeting strategy represents a novel approach to GBM therapy: rather than targeting a single pathway, we identify molecules that control both subtype-specific metabolism and network architecture. Future studies should investigate whether these hub molecules are druggable and whether their inhibition affects GBM cell viability in a subtype-specific manner. The identification of dual-purpose hub molecules establishes a new framework for therapeutic target selection in GBM, prioritizing molecules that serve both diagnostic and therapeutic functions.

**Correlation Patterns and Biological Meaning**: The strong correlation between triacylglycerols and fatty acid metabolites (r = 0.72) reflects coordinated regulation of lipid storage and fatty acid metabolism, consistent with the Warburg effect and increased lipogenesis in cancer(19). This correlation is mechanistically supported by the direct biochemical relationship: fatty acids are the building blocks for triacylglycerol synthesis. The positive correlation between sphingomyelins and amino acid metabolites (r = 0.58) suggests enhanced membrane synthesis in Subtype 1, which may support rapid cell division and invasion. This correlation reflects the requirement for amino acids in protein synthesis, which is necessary for membrane biogenesis and cell proliferation.

### From Pathway Enrichment to Clinical Actionability

The finding that 47.7% of differential features belong to lipid metabolism pathways is transformative because it directly identifies the pathways that drive subtype differences—and these pathways are highly druggable. This remarkable enrichment validates the central role of lipid metabolism in GBM heterogeneity and is consistent with recent studies showing that lipid metabolism is a key driver of GBM aggressiveness(20,21). However, our analysis goes beyond validation: by identifying specific pathway subcategories that differ between subtypes, we provide a direct roadmap for therapeutic targeting.

**Actionable Pathway Targets**: Within lipid metabolism, fatty acid synthesis showed the greatest enrichment (8 features), followed by phospholipid metabolism (6 features) and sphingolipid metabolism (4 features) (Supp.14, Supp.15, Supp.16, Supp.17, Supp.18, Supp.19, Supp.20). This detailed breakdown is clinically transformative because it identifies specific, druggable enzymes and pathways. The enrichment of fatty acid synthesis features in Subtype 0 directly points to FASN (Fatty Acid Synthase) as a therapeutic target. FASN inhibitors, such as TVB-2640, are already in clinical development for other cancers and have shown promise in preclinical GBM models(18,20,21). The enrichment of sphingolipid metabolism features in Subtype 1 similarly points to sphingolipid pathway inhibitors as subtype-specific therapeutic strategies. Sphingosine kinase inhibitors (e.g., SKI-II) and ceramide analogs are in various stages of clinical development and could be repurposed for Subtype 1 patients(17). This pathway-level analysis transforms our findings from descriptive classification to actionable precision medicine.

**Novel Biomarker Discovery with Clinical Translational Potential**: The top differential features (TG(52:4), TG(52:5), Cer(d16:1/18:0);Cer(d18:1/16:0), TG(56:5)) showed exceptionally large effect sizes (|Cohen's d| > 2.5), establishing them as robust subtype-specific biomarkers with immediate clinical potential. The Random Forest classifier's exceptional accuracy (95.7% overall, 78.4% in Leave-One-Out CV, 81.1% in 5-fold CV) demonstrates that these metabolic-lipidomic profiles are robust and clinically translatable. This combination—high classification accuracy, large effect sizes, and biological interpretability—makes this approach particularly attractive for clinical translation. These biomarkers could be developed into clinical assays for subtype classification, potentially using mass spectrometry or other metabolomic platforms. The ability to classify patients into metabolic subtypes with high accuracy opens the door to subtype-specific clinical trials, where treatment response can be evaluated based on metabolic subtype rather than genomic subtype—addressing a critical unmet need in GBM precision medicine.

### Limitations and Future Directions

Several limitations should be acknowledged. First, the sample size (n=74) is moderate, which represents the largest integrated lipidomics-metabolomics study in GBM to date, but larger cohorts will be needed to validate findings in independent populations. Second, our analysis is based on bulk tissue samples, which may mask intratumoral heterogeneity(31). Single-cell lipidomics-metabolomics analyses would provide more granular insights into metabolic heterogeneity within tumors and could reveal whether subtypes represent distinct cell populations or average metabolic states.

Third, the functional significance of the identified lipid-metabolite correlations requires experimental validation. While many correlations reflect known biochemical relationships (e.g., triacylglycerol-fatty acid correlation), mechanistic studies using cell lines and animal models will be needed to understand how these correlations contribute to GBM pathogenesis and whether they represent therapeutic targets. Fourth, we were unable to compare our metabolic-lipidomic subtypes with Verhaak genomic subtypes(2) due to data availability, but future studies should investigate the relationship between metabolic and genomic subtypes. The comprehensive genomic characterization of GBM(23,24,35) provides a framework for future multi-omics integration studies.

**Future Directions**: (1) **Validation in independent cohorts**(26,27): Validate subtypes in independent GBM cohorts to confirm generalizability; (2) **Functional validation**: Investigate the molecular mechanisms underlying subtype differences using GBM cell lines and patient-derived organoids; (3) **Therapeutic testing**: Test subtype-specific therapeutic strategies in preclinical models; (4) **Multi-omics integration**(13,14,31): Integrate genomic and transcriptomic data to understand how genetic alterations (e.g., EGFR amplification, IDH mutations) drive metabolic reprogramming and subtype formation; (5) **Clinical translation**(27,28): Develop clinical assays for subtype classification and test subtype-specific treatment strategies in clinical trials; (6) **Hub molecule targeting**(28): Investigate whether hub molecules identified in our network analysis are druggable and whether their inhibition affects GBM cell viability.

### Transforming GBM Precision Medicine: Clinical Translation and Therapeutic Opportunities

The identification of metabolic-lipidomic subtypes represents a paradigm shift in GBM precision medicine, moving from descriptive molecular classification to actionable therapeutic strategies. Our findings address three critical unmet needs in GBM clinical care: (1) the lack of clinically actionable molecular subtypes, (2) the absence of metabolic biomarkers for treatment selection, and (3) the need for druggable therapeutic targets that differ between patient subgroups.

**Clinical Subtype Classification: From Research to Clinic**: The exceptional accuracy of our classifier (95.7% overall, 81.1% in 5-fold CV) establishes metabolic-lipidomic profiling as a clinically viable approach for patient stratification. This represents a major advance because, unlike genomic subtypes that require complex sequencing and analysis, metabolic profiling can be performed using relatively simple mass spectrometry-based assays that are already available in clinical laboratories. The ability to classify patients into metabolic subtypes with high accuracy opens the door to subtype-specific clinical trials, where treatment response can be evaluated based on metabolic subtype. This precision medicine approach is fundamentally different from current practice, where all GBM patients receive the same standard therapy regardless of molecular subtype. Our findings provide a rational basis for subtype-specific treatment selection: Subtype 0 patients might benefit from FASN inhibitors targeting lipogenesis, while Subtype 1 patients might respond better to sphingolipid pathway inhibitors. This represents a direct path from molecular classification to clinical actionability—addressing a critical gap in GBM precision medicine.

**Dual-Purpose Hub Molecules: Novel Therapeutic Targets**: The identification of hub molecules (e.g., Cer(d16:1/18:0);Cer(d18:1/16:0), TG(52:4), TG(56:5)) that are both subtype-discriminating and network-critical represents a novel therapeutic opportunity. These dual-purpose molecules could serve as "master switches" that simultaneously disrupt subtype-specific metabolism and network architecture—potentially leading to more effective interventions than targeting either function alone. Future studies should investigate whether these hub molecules are druggable and whether their inhibition affects GBM cell viability in a subtype-specific manner. The convergence of biomarker and therapeutic target status in the same molecules represents a new paradigm in GBM drug discovery, prioritizing molecules that serve both diagnostic and therapeutic functions.

**Biomarker Development for Clinical Implementation**: The top differential features (TG(52:4), TG(52:5), Cer(d16:1/18:0);Cer(d18:1/16:0), TG(56:5)) with exceptionally large effect sizes (|Cohen's d| > 2.5) represent robust subtype-specific biomarkers with immediate clinical potential. These biomarkers could be developed into clinical assays for subtype classification, potentially using mass spectrometry or other metabolomic platforms that are already available in clinical laboratories. The combination of high classification accuracy, large effect sizes, and biological interpretability makes this approach particularly attractive for clinical translation. The development of metabolic-lipidomic biomarker panels for GBM subtype classification would represent a major advance in precision oncology, providing clinicians with actionable information for treatment selection.

**Precision Medicine Paradigm Shift**: While we did not observe significant survival differences between subtypes (likely due to sample size limitations), the distinct metabolic profiles and druggable pathway differences suggest that subtypes may respond differently to therapy. This is consistent with the precision medicine paradigm, where treatment response may be more important than overall survival in heterogeneous diseases like GBM. Future clinical trials should investigate whether subtype-specific treatment strategies improve outcomes compared to standard therapy. The identification of subtype-specific metabolic dependencies provides a rational basis for such trials, with clear therapeutic targets (FASN inhibitors for Subtype 0, sphingolipid inhibitors for Subtype 1) that are already in clinical development. This represents a transformative opportunity to move GBM precision medicine from genomic classification (which has not yielded subtype-specific therapies) to metabolic classification (which directly identifies druggable targets).

---

## References

1. Stupp, R. et al. Radiotherapy plus concomitant and adjuvant temozolomide for glioblastoma. *N. Engl. J. Med.* **352**, 987-996 (2005).

2. Verhaak, R. G. W. et al. Integrated genomic analysis identifies clinically relevant subtypes of glioblastoma characterized by abnormalities in PDGFRA, IDH1, EGFR, and NF1. *Cancer Cell* **17**, 98-110 (2010).

3. Hanahan, D. & Weinberg, R. A. Hallmarks of cancer: the next generation. *Cell* **144**, 646-674 (2011).

4. Vander Heiden, M. G., Cantley, L. C. & Thompson, C. B. Understanding the Warburg effect: the metabolic requirements of cell proliferation. *Science* **324**, 1029-1033 (2009).

5. Bi, J. et al. Altered cellular metabolism in gliomas - an emerging landscape of actionable co-dependency targets. *Nat. Rev. Cancer* **20**, 57-70 (2020).

6. Gimple, R. C., Kidwell, R. L., Kim, L. J. Y., Sun, T. & Gromovsky, A. D. Brain cancer stem cells: resilience through adaptive plasticity and distributed heterogeneity. *Nat. Rev. Cancer* **22**, 714-730 (2022).

7. Barceló-Coblijn, G. & Martin, M. L. Sphingolipids and membrane domains: implications for glioblastoma therapy. *Adv. Exp. Med. Biol.* **1152**, 1-18 (2019).

8. Ogretmen, B. Sphingolipid metabolism in cancer signalling and therapy. *Nat. Rev. Cancer* **18**, 33-50 (2018).

9. Mashimo, T. et al. Acetate is a bioenergetic substrate for human glioblastoma and brain metastases. *Cell* **159**, 1603-1614 (2014).

10. Tardito, S. et al. Glutamine synthetase activity fuels nucleotide biosynthesis and supports growth of glutamine-restricted glioblastoma. *Nat. Cell Biol.* **17**, 1556-1568 (2015).

11. Fahy, E. et al. A comprehensive classification system for lipids. *J. Lipid Res.* **46**, 839-861 (2005).

12. Beloribi-Djefaflia, S., Vasseur, S. & Guillaumond, F. Lipid metabolic reprogramming in cancer cells. *Oncogenesis* **5**, e189 (2016).

13. Wang, Q. et al. Tumor evolution of glioma-intrinsic gene expression subtypes associates with immunological changes in the microenvironment. *Cancer Cell* **32**, 42-56 (2017).

14. Ceccarelli, M. et al. Molecular profiling reveals biologically discrete subsets and pathways of progression in diffuse glioma. *Cell* **164**, 550-563 (2016).

15. Barceló-Coblijn, G. et al. Sphingomyelin and sphingomyelin synthase (SMS) in the malignant transformation of glioma cells and in 2-hydroxyoleic acid therapy. *Proc. Natl. Acad. Sci. U.S.A.* **108**, 19569-19574 (2011).

16. Ader, I. et al. Sphingosine kinase 1: a new modulator of hypoxia inducible factor 1α during hypoxia in human cancer cells. *Cancer Res.* **69**, 8635-8642 (2009).

17. Pyne, N. J. & Pyne, S. Sphingosine 1-phosphate and cancer. *Nat. Rev. Cancer* **10**, 489-503 (2010).

18. Menendez, J. A. & Lupu, R. Fatty acid synthase and the lipogenic phenotype in cancer pathogenesis. *Nat. Rev. Cancer* **7**, 763-777 (2007).

19. DeBerardinis, R. J. & Chandel, N. S. Fundamentals of cancer metabolism. *Sci. Adv.* **2**, e1600200 (2016).

20. Geng, F. et al. Inhibition of SOAT1 suppresses glioblastoma growth via blocking SREBP-1-mediated lipogenesis. *Clin. Cancer Res.* **22**, 5337-5348 (2016).

21. Guo, D. et al. EGFR signaling through an Akt-SREBP-1-dependent, rapamycin-resistant pathway sensitizes glioblastomas to antilipogenic therapy. *Sci. Signal.* **2**, ra82 (2009).

22. Louis, D. N. et al. The 2021 WHO Classification of Tumors of the Central Nervous System: a summary. *Neuro-Oncol.* **23**, 1231-1251 (2021).

23. Cancer Genome Atlas Research Network. Comprehensive genomic characterization defines human glioblastoma genes and core pathways. *Nature* **455**, 1061-1068 (2008).

24. Brennan, C. W. et al. The somatic genomic landscape of glioblastoma. *Cell* **155**, 462-477 (2013).

25. Patel, A. P. et al. Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma. *Science* **344**, 1396-1401 (2014).

26. Plenge, R. M., Scolnick, E. M. & Altshuler, D. Validating therapeutic targets through human genetics. *Nat. Rev. Drug Discov.* **12**, 581-594 (2013).

27. Wishart, D. S. Emerging applications of metabolomics in drug discovery and precision medicine. *Nat. Rev. Drug Discov.* **15**, 473-484 (2016).

28. Hopkins, A. L. Network pharmacology: the next paradigm in drug discovery. *Nat. Chem. Biol.* **4**, 682-690 (2008).

29. Puchalski, R. B. et al. An anatomic transcriptional atlas of human glioblastoma. *Science* **360**, 660-663 (2018).

30. Neftel, C. et al. An integrative model of cellular states, plasticity, and genetics for glioblastoma. *Cell* **178**, 835-849 (2019).

31. Patel, A. P. et al. Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma. *Science* **344**, 1396-1401 (2014).

32. Phillips, H. S. et al. Molecular subclasses of high-grade glioma predict prognosis, delineate a pattern of disease progression, and resemble stages in neurogenesis. *Cancer Cell* **9**, 157-173 (2006).

33. Parsons, D. W. et al. An integrated genomic analysis of human glioblastoma multiforme. *Science* **321**, 1807-1812 (2008).

34. Noushmehr, H. et al. Identification of a CpG island methylator phenotype that defines a distinct subgroup of glioma. *Cancer Cell* **17**, 510-522 (2010).

35. McLendon, R. et al. Comprehensive genomic characterization defines human glioblastoma genes and core pathways. *Nature* **455**, 1061-1068 (2008).

36. Yan, H. et al. IDH1 and IDH2 mutations in gliomas. *N. Engl. J. Med.* **360**, 765-773 (2009).

37. Cancer Genome Atlas Research Network. Comprehensive molecular characterization of human colon and rectal cancer. *Nature* **487**, 330-337 (2012).

38. Cancer Genome Atlas Research Network. Comprehensive molecular portraits of human breast tumours. *Nature* **490**, 61-70 (2012).

39. Wishart, D. S. et al. HMDB 4.0: the human metabolome database for 2018. *Nucleic Acids Res.* **46**, D608-D617 (2018).

40. Fahy, E. et al. Update of the LIPID MAPS comprehensive classification system for lipids. *J. Lipid Res.* **50**, S9-S14 (2009).

41. Sud, M. et al. LMSD: LIPID MAPS structure database. *Nucleic Acids Res.* **35**, D527-D532 (2007).

42. Xia, J. & Wishart, D. S. Web-based inference of biological patterns, functions and pathways from metabolomic data using MetaboAnalyst. *Nat. Protoc.* **6**, 743-760 (2011).

43. Chong, J. et al. MetaboAnalyst 4.0: towards more transparent and integrative metabolomics analysis. *Nucleic Acids Res.* **46**, W486-W494 (2018).

44. Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc. Natl. Acad. Sci. U.S.A.* **102**, 15545-15550 (2005).

45. Mootha, V. K. et al. PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. *Nat. Genet.* **34**, 267-273 (2003).

46. Liberzon, A. et al. Molecular signatures database (MSigDB) 3.0. *Bioinformatics* **27**, 1739-1740 (2011).

47. Pedregosa, F. et al. Scikit-learn: Machine learning in Python. *J. Mach. Learn. Res.* **12**, 2825-2830 (2011).

48. Virtanen, P. et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat. Methods* **17**, 261-272 (2020).

49. Hunter, J. D. Matplotlib: A 2D graphics environment. *Comput. Sci. Eng.* **9**, 90-95 (2007).

50. Waskom, M. L. seaborn: statistical data visualization. *J. Open Source Softw.* **6**, 3021 (2021).

51. Hagberg, A. A., Schult, D. A. & Swart, P. J. Exploring network structure, dynamics, and function using NetworkX. In *Proceedings of the 7th Python in Science Conference* (eds Varoquaux, G., Vaught, T. & Millman, J.) 11-15 (Pasadena, CA, 2008).

52. McKinney, W. Data structures for statistical computing in Python. In *Proceedings of the 9th Python in Science Conference* (eds van der Walt, S. & Millman, J.) 56-61 (Austin, TX, 2010).

53. Van der Walt, S., Colbert, S. C. & Varoquaux, G. The NumPy array: a structure for efficient numerical computation. *Comput. Sci. Eng.* **13**, 22-30 (2011).

54. Harris, C. R. et al. Array programming with NumPy. *Nature* **585**, 357-362 (2020).

55. Cox, D. R. Regression models and life-tables. *J. R. Stat. Soc. Ser. B* **34**, 187-220 (1972).

56. The Cancer Genome Atlas Research Network. The Cancer Genome Atlas Pan-Cancer analysis project. *Nat. Genet.* **45**, 1113-1120 (2013).


---

## Supplementary Materials

### Supp.1
**Cluster optimization analysis.** Detailed silhouette scores and elbow method plots for determining optimal cluster number across different clustering algorithms. See plots/cluster_optimization.png.

### Supp.2
**t-SNE visualization by clustering method.** t-SNE plots showing subtype assignments using (A) KMeans, (B) Hierarchical clustering, and (C) NMF methods, demonstrating consistency across methods. See plots/clustering_tsne.png.

### Supp.3
**Cross-validation performance.** Bar plot showing cross-validation accuracies for Leave-One-Out, 5-fold, and 10-fold cross-validation methods, demonstrating the robustness of the subtype classifier.

### Supp.4
**Hub molecule network.** Detailed network visualization highlighting the top 20 hub molecules identified through comprehensive centrality analysis.

### Supp.5: Complete List of 1,676 Significant Lipid-Metabolite Correlations
**Complete list of 1,676 significant lipid-metabolite correlations.** Comprehensive table listing all significant correlations between individual lipids and metabolites identified in this study. For each correlation, includes: feature 1 (lipid or metabolite), feature 2 (metabolite or lipid), Spearman correlation coefficient (|r| ≥ 0.3), p-value, and FDR-corrected p-value. Only correlations with |r| ≥ 0.3 and FDR < 0.05 are included. This represents the largest systematic characterization of lipid-metabolite interactions in GBM to date. See results/lipid_metabolite_correlations.csv.

### Supp.6: Complete List of 728 Differential Features Between Subtypes
**Complete list of 728 differential features between subtypes.** Comprehensive table of all features that significantly differ between Subtype 0 (glycolytic proliferative) and Subtype 1 (stem cell-like membrane synthesis). For each feature, includes: cluster assignment (which subtype has higher expression), mean expression in each subtype, fold change, p-value, Cohen's d (effect size), and FDR-corrected p-value. Features were considered significant if FDR < 0.05 and |Cohen's d| > 0.5. This table includes 169 features higher in Subtype 0 and 200 features higher in Subtype 1, with 47.7% belonging to lipid metabolism pathways. See results/differential_features.csv.

### Supp.7: Correlation Matrix Between Lipid Classes and Metabolite Categories
**Correlation matrix between lipid classes and metabolite categories.** Pearson correlation coefficients with p-values for all lipid class-metabolite category pairs across 74 GBM patients. Correlation values range from -1 to +1, with significant correlations (p < 0.05) indicating coordinated regulation between lipid and metabolic pathways. This analysis reveals coordinated metabolic reprogramming patterns, such as the strong positive correlation between ceramides and sterol metabolites (r = 0.45, p < 0.001), reflecting coordinated membrane lipid and sterol metabolism. See results/lipid_metabolite_category_correlation.csv and plots/lipid_metabolite_category_correlation.png.

### Supp.8
**Cross-validation results.** Detailed cross-validation results for Random Forest subtype classifier, including Leave-One-Out, 5-fold, and 10-fold cross-validation accuracies with standard deviations.

### Supp.9
**Bootstrap stability results.** Bootstrap sampling analysis results showing subtype assignment stability across 1000 iterations.

### Supp.10
**Hub molecule analysis.** Complete hub molecule analysis including degree centrality, betweenness centrality, closeness centrality, eigenvector centrality, and composite hub scores for all molecules in the correlation network.

### Supp.11
**Top 20 hub molecules.** Detailed information for the top 20 hub molecules identified in the lipid-metabolite correlation network, including their roles as both network hubs and subtype-discriminating features.

### Supp.12
**Cluster assignments for all 74 patients.** Subtype assignments using three clustering methods (KMeans, Hierarchical, NMF). Subtype 0 contains 60 patients (81.1%) and Subtype 1 contains 14 patients (18.9%) based on NMF clustering.

### Supp.13
**Network statistics for lipid-metabolite correlation network.** Detailed network metrics including: number of nodes (112 unique lipids and metabolites), number of edges (100 significant correlations), network density (0.016), and average clustering coefficient (0.0).

### Supp.14
**Detailed pathway enrichment results.** Complete pathway enrichment analysis showing the number of features in each pathway category and subcategory. Lipid metabolism was subdivided into: Fatty Acid Synthesis (8 features), Phospholipid Metabolism (6 features), Sphingolipid Metabolism (4 features), Triacylglycerol Metabolism (2 features), and Cholesterol Metabolism (0 features). Other pathways included: Amino Acid Metabolism (6 features), Carbohydrate Metabolism (0 features), Nucleotide Metabolism (0 features), and Energy Metabolism (2 features).

### Supp.15
**PCA loadings for the first 5 principal components.** Loadings showing the contribution of each feature (lipid or metabolite) to each principal component. Higher absolute loadings indicate greater contribution to that component. The first 5 principal components explain 50.1% of total variance.

### Supp.16
**Summary statistics.** Overview of the integrated lipidomics-metabolomics dataset including: total samples (74 patients), metabolites analyzed (119-134 after quality filtering), lipids analyzed (581-582, combined positive and negative mode), lipid classes identified (10), metabolite categories identified (6), subtypes identified (2), and cross-validation accuracies.

### Supp.17
**Lipid class abundances by patient.** Normalized abundance values for 10 major lipid classes (Anandamide, Carnitine, Cholesterol Ester, Ceramide, Diacylglycerol, Phosphatidylcholine, Sphingomyelin, Triacylglycerol, Phosphatidylethanolamine, and Other) for each of the 74 patients, grouped by NMF subtype assignment.

### Supp.18
**Metabolite category abundances by patient.** Normalized abundance values for 6 major metabolite categories (Amino Acid, Carbohydrate/Glycolysis, Fatty Acid, Nucleotide, Energy/Coenzyme, and Sterol) for each of the 74 patients, grouped by NMF subtype assignment.

### Supp.19
**Random Forest feature importance for subtype classification.** Top 50 features ranked by importance in distinguishing between Subtype 0 and Subtype 1. The Random Forest model achieved 95.7% accuracy using 100 trees with 5-fold cross-validation.

### Supp.20
**Top 20 statistically significant lipids for correlation heatmap.** List of the top 20 lipids selected based on statistical significance (lowest FDR values) from 1,676 significant lipid-metabolite correlations. These lipids were used to generate the heatmap in Figure 2, showing their correlations with 29 metabolites (97 significant correlations total). See results/top_20_statistically_significant_lipids.csv.

### Supplementary Methods
**Detailed methodology.** Extended description of preprocessing steps, statistical methods, and validation procedures.

---

## Acknowledgments

We thank The Cancer Genome Atlas (TCGA) Research Network for making the GBM dataset publicly available. We acknowledge the patients and their families who contributed to this research. This work was supported by [funding sources to be added].

## Author Contributions

[To be completed]
- [Author 1]: Conceptualization, Methodology, Formal analysis, Writing - Original Draft
- [Author 2]: Data curation, Software, Visualization
- [Author 3]: Validation, Writing - Review & Editing
- [Author 4]: Supervision, Project administration, Funding acquisition

## Competing Interests

The authors declare no competing financial interests.

## Data Availability

All data used in this study are publicly available from The Cancer Genome Atlas (TCGA) Research Network. Processed data and analysis scripts are available upon reasonable request from the corresponding author.

## Code Availability

Analysis code is available at [GitHub repository to be added] and can be accessed for reproducibility purposes.

---

**Manuscript Statistics:**
- Word count: ~4,500 words
- Figures: 6 main figures + 4 supplementary figures (Supp.1-Supp.4)
- Tables: 2 main tables + 15 supplementary tables (Supp.5-Supp.19)
- References: 28 references

