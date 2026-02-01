![Fig1](https://github.com/user-attachments/assets/da5bacf5-4f92-4760-9e91-94f9ff4cbe35)# Proteomics_Migraine

**Plasma proteomics identifies proteins and pathways associated with incident migraine in 50,668 adults**

This repository contains analysis code supporting a large-scale prospective plasma proteomics study of incident migraine based on the UK Biobank.

The project integrates **longitudinal proteomics**, **survival analysis**, **trajectory modeling**, **single-cell transcriptomics**, **genetic analyses (PRS, MR, colocalization)**, **brain MRI**, and **machine learning–based risk prediction** to characterize pre-diagnostic molecular signatures of migraine.

---

## 📌 Study Overview

- **Cohort**: UK Biobank (n = 50,668, migraine-free at baseline)
- **Proteomics**: Olink Explore (≈2,900 plasma proteins)
- **Outcome**: Incident migraine (ICD-10 G43)
- **Follow-up**: Mean ~13 years
- **Key methods**:
  - Cox proportional hazards models
  - Pre-diagnostic protein trajectory analysis
  - Functional enrichment & PPI networks
  - Single-cell RNA-seq integration
  - Brain MRI association analysis
  - Polygenic risk scores (PRS)
  - Mendelian randomization & colocalization
  - Machine learning prediction (LightGBM)

---

## 📁 Repository Structure

## 🧬 Key Analyses

### 1. Plasma proteome–wide association
- Cox regression for incident migraine
- Sex- and age-stratified analyses
- FDR correction for multiple testing

### 2. Pre-diagnostic protein trajectories
- Incidence-density matching
- Residualization of protein levels
- LOESS-based temporal profiling
- Unsupervised trajectory clustering

### 3. Functional & network analyses
- GO / KEGG enrichment
- Protein–protein interaction networks
- Hub protein identification
- Transcription factor target analysis

### 4. Single-cell transcriptomics
- scRNA-seq integration (GSE269117)
- Immune cell annotation
- Cell-type–specific expression of migraine-associated genes

### 5. Brain MRI integration
- Global brain structure phenotypes (WMH, CGV, TGV, SGV)
- Protein–brain structure associations
- Adjustment for TIV and assessment center

### 6. Genetic analyses
- Migraine PRS–protein associations
- Mendelian randomization
- Bayesian colocalization

### 7. Risk prediction
- LightGBM classifiers
- Cross-validated feature importance
- Stepwise AUC evaluation
- ROC curve visualization

---

## 🔒 Data Availability

- **UK Biobank data** are available via application to UK Biobank  
- **Individual-level data are not shared** due to access restrictions
- scRNA-seq data: **GEO GSE269117**

---

## 🧪 Software & Environment

### R
- R ≥ 4.2.0  
- Key packages:
  - `survival`, `MatchIt`, `limma`, `clusterProfiler`
  - `Seurat`, `harmony`, `SingleR`
  - `TwoSampleMR`, `coloc`

### Python
- Python ≥ 3.8  
- Key packages:
  - `pandas`, `numpy`, `scikit-learn`
  - `lightgbm`, `matplotlib`

---

## ⭐ Notes

- This repository focuses on **analysis pipelines**, not raw data
- Scripts are modular and can be adapted to other proteomic or longitudinal disease studies
- Contributions and issues are welcome
