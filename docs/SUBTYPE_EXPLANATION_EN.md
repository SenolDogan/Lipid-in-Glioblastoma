# Subtype 0 and Subtype 1 Explanation

## Overview

This analysis classified GBM (Glioblastoma) patients into **2 distinct subtypes** based on their metabolic and lipidomic profiles. These subtypes were identified using **Non-negative Matrix Factorization (NMF)** method.

---

## Patient Distribution

- **Subtype 0**: 60 patients (81.1%) - **Majority group**
- **Subtype 1**: 14 patients (18.9%) - **Minority group**
- **Total**: 74 patients

---

## Subtype 0 Characteristics

### Lipid Profile:
- ✅ **High Triacylglycerol (TG)**: 7.9% higher than Subtype 1
  - **Meaning**: High energy storage capacity
  - **Clinical Significance**: Lipid droplet formation, tumor microenvironment
  
- ✅ **High Carnitine**: 3.9% higher than Subtype 1
  - **Meaning**: Active fatty acid oxidation
  - **Clinical Significance**: Energy production, ketogenic metabolism

- ❌ **Low Ceramide**: 4.3% lower than Subtype 1
  - **Meaning**: Apoptosis mechanisms may be less active
  - **Clinical Significance**: Tumor suppression mechanisms

- ❌ **Low Sphingomyelin**: 5.2% lower than Subtype 1
  - **Meaning**: Different membrane structure
  - **Clinical Significance**: Cell signaling

- ❌ **Low Phosphatidylcholine (PC)**: 3.5% lower than Subtype 1
  - **Meaning**: Different membrane synthesis
  - **Clinical Significance**: Cell division

### Metabolic Profile:
- **169 high features** (compared to Subtype 1)
- **Energy metabolism**: High Triacylglycerol and Carnitine indicate energy production-focused metabolism

### Biological Meaning:
**Subtype 0 = "Energy-Focused Metabolic Subtype"**
- High energy storage (TG)
- Active fatty acid oxidation (Carnitine)
- Low apoptosis signals (Ceramide)
- Less membrane synthesis (PC, Sphingomyelin)

---

## Subtype 1 Characteristics

### Lipid Profile:
- ✅ **High Ceramide**: 4.5% higher than Subtype 0
  - **Meaning**: Apoptosis mechanisms more active
  - **Clinical Significance**: Tumor suppression potential

- ✅ **High Sphingomyelin**: 5.2% higher than Subtype 0
  - **Meaning**: Different membrane structure and signaling
  - **Clinical Significance**: Cell signaling, membrane integrity

- ✅ **High Phosphatidylcholine (PC)**: 3.6% higher than Subtype 0
  - **Meaning**: Active membrane synthesis
  - **Clinical Significance**: Rapid cell division

- ❌ **Low Triacylglycerol (TG)**: 7.3% lower than Subtype 0
  - **Meaning**: Less energy storage
  - **Clinical Significance**: Different energy strategy

- ❌ **Low Carnitine**: 3.9% lower than Subtype 0
  - **Meaning**: Fatty acid oxidation less active
  - **Clinical Significance**: Different energy production pathway

### Metabolic Profile:
- **200 high features** (compared to Subtype 0)
- **Membrane and signal metabolism**: High PC, Sphingomyelin, and Ceramide indicate membrane synthesis and signaling-focused metabolism

### Biological Meaning:
**Subtype 1 = "Membrane and Signal-Focused Metabolic Subtype"**
- High membrane synthesis (PC, Sphingomyelin)
- Active apoptosis signals (Ceramide)
- Low energy storage (TG)
- Less fatty acid oxidation (Carnitine)

---

## Differences Between Subtypes

### 1. Energy Metabolism
- **Subtype 0**: High energy storage and fatty acid oxidation
- **Subtype 1**: Low energy storage, different energy strategy

### 2. Membrane Metabolism
- **Subtype 0**: Low membrane synthesis
- **Subtype 1**: High membrane synthesis (PC, Sphingomyelin)

### 3. Apoptosis and Signaling
- **Subtype 0**: Low apoptosis signals (low Ceramide)
- **Subtype 1**: High apoptosis signals (high Ceramide)

### 4. Cell Division
- **Subtype 0**: Less membrane synthesis → potentially slower division
- **Subtype 1**: High membrane synthesis → potentially faster division

---

## Clinical Significance

### Subtype 0 (Energy-Focused):
- **Prognosis**: High energy storage → potential for longer survival?
- **Treatment**: 
  - Lipid synthesis inhibitors (FASN inhibitors)
  - Fatty acid oxidation inhibitors
  - Ketogenic diet may be contraindicated (high carnitine)
- **Biomarker**: High TG/Carnitine ratio

### Subtype 1 (Membrane-Focused):
- **Prognosis**: High membrane synthesis → more aggressive tumor?
- **Treatment**:
  - Membrane synthesis inhibitors
  - Ceramide synthesis activators (apoptosis induction)
  - Sphingomyelin metabolism targets
- **Biomarker**: High PC/Sphingomyelin/Ceramide ratio

---

## Survival Analysis

### Log-rank Test Results:
- **KMeans**: p=0.3047 (not significant)
- **Hierarchical**: p=0.1634 (trend present, not significant)
- **NMF**: p=0.3781 (not significant)

### Interpretation:
- **No statistically significant survival difference** between subtypes
- However, **a trend exists** (Hierarchical: p=0.16)
- Significant difference may be found in larger cohorts
- **Subtype 1 has potential to be more aggressive** (high membrane synthesis)

---

## Summary Table

| Feature | Subtype 0 (n=60, 81.1%) | Subtype 1 (n=14, 18.9%) |
|---------|------------------------|------------------------|
| **Main Feature** | Energy-Focused | Membrane-Focused |
| **Triacylglycerol** | ✅ High | ❌ Low |
| **Carnitine** | ✅ High | ❌ Low |
| **Ceramide** | ❌ Low | ✅ High |
| **Sphingomyelin** | ❌ Low | ✅ High |
| **Phosphatidylcholine** | ❌ Low | ✅ High |
| **Energy Storage** | ✅ High | ❌ Low |
| **Membrane Synthesis** | ❌ Low | ✅ High |
| **Apoptosis Signals** | ❌ Low | ✅ High |
| **Potential Aggressiveness** | Low-Medium | High |

---

## Conclusion

**Subtype 0** and **Subtype 1** represent two distinct biological subtypes reflecting the metabolic heterogeneity of GBM:

1. **Subtype 0**: Focused on energy storage and fatty acid oxidation, potential for slower growth
2. **Subtype 1**: Focused on membrane synthesis and cell division, potential for more aggressive growth

These subtypes can be used for **personalized treatment strategies** and may have **prognostic value**.

---

**Note**: These subtypes are based solely on metabolic and lipidomic profiles. A more comprehensive classification can be made when integrated with genomic, transcriptomic, and other omics data.













