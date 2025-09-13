# About: Protein Differential Expression (Protein DEX) Shiny App

The **Protein DEX Shiny App** is an interactive tool for analyzing and visualizing proteomics differential expression data.  
It is designed to streamline the process of **data preprocessing, normalization, statistical testing, and visualization** for label-based proteomics datasets.

This project was developed by **Rithwik Nambiar** (2025) as part of personal and academic explorations in bioinformatics.

### ✨ Features
- **Data Input**: Supports CSV and Excel proteomics data.  
- **Preprocessing**:  
  - Impute missing values.  
  - Remove contaminants.  
  - Filter by detection thresholds.  
- **Normalization**:  
  - Median normalization.  
  - Quantile normalization.  
- **Statistical Analysis**:  
  - Fold change calculation.  
  - t-test based p-values.  
  - Significance classification (Up/Down/Not Sig).  
- **Interactive Visualizations**:  
  - PCA plots.  
  - Volcano plots with customizable highlights.  
  - Heatmaps of top proteins.  
  - Violin + boxplots (before vs. after normalization).  
- **Downloads**:  
  - Normalized data.  
  - Differential expression results.  
  - Upregulated and downregulated proteins.  
  - Publication-ready plots (PCA, volcano, boxplot, heatmap).

### 🛠️ Tech Stack
- **R** with **Shiny**  
- Supporting packages: `tidyr`, `dplyr`, `readxl`, `writexl`, `preprocessCore`, `ggplot2`, `DT`, `ggrepel`, `reshape2`, `pheatmap`.

### 🎯 Target Users
This tool is intended for:
- Bioinformatics researchers working with **proteomics datasets**.  
- Students and scientists who want a **user-friendly, reproducible pipeline** for protein differential expression analysis.  

---
