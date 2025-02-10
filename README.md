# Integrating Single-cell Sequencing and Clinical Insights to Explore Malignant Transformation in Odontogenic Keratocyst

This repository contains the code for the paper titled **"Integrating Single-cell Sequencing and Clinical Insights to Explore Malignant Transformation in Odontogenic Keratocyst."**

## File Structure

```
scripts/
    cellchat.R                     # Cell-cell communication analysis using CellChat
    Dotplot.R                       # Visualization of gene expression through dot plot
    GSVA and survival analysis.R    # GSVA analysis and survival analysis
    Heatmap.R                       # Heatmap visualization of gene expression data
    Processing.R                    # Data preprocessing for single-cell analysis
    slingshot.R                     # Trajectory inference using Slingshot
README.md                         # This file
```

## Requirements

### R Version
- R (version 4.2 or higher)

### Required R packages:
- **Seurat** - For single-cell RNA sequencing analysis.
- **ggplot2** - For data visualization.
- **dplyr** - For data manipulation.
- **Harmony** - For integrating multiple single-cell datasets.
- **monocle** - For single-cell trajectory analysis.
- **GSVA** - For Gene Set Variation Analysis.
- **ComplexHeatmap** - For generating heatmaps with annotations.
- **circlize** - For circular visualizations and complex data visualization.

## Installation

To install the required packages, run the following R code:

```R
install.packages(c("ggplot2", "dplyr", "circlize", "ComplexHeatmap"))
install.packages("Seurat")
devtools::install_github("jokergoo/ComplexHeatmap")
install.packages("GSVA")
install.packages("monocle")
install.packages("harmony")
install.packages("CellChat")
```

## Description

This repository includes scripts for processing single-cell RNA sequencing data, performing various analyses such as clustering, differential expression analysis, and gene set enrichment analysis, as well as visualizing the results through heatmaps, dot plots, and survival analysis.

The main analysis includes integrating single-cell data with clinical insights to explore the malignant transformation of odontogenic keratocyst (OKC) and its progression to cancerous forms (COKC).

### Key Components:
- **cellchat.R**: Analyzes cell-cell communication using the CellChat package, providing insights into the interactions between different cell types in OKC and COKC.
- **Dotplot.R**: Generates dot plots for visualizing the expression of selected genes across clusters or conditions.
- **GSVA and survival analysis.R**: Performs GSVA to analyze gene set activity in single-cell data and survival analysis to correlate gene signatures with patient survival.
- **Heatmap.R**: Creates heatmaps to display the Z-scores of gene expression levels and their clustering across different samples.
- **Processing.R**: Prepares the data, normalizes gene expression, and filters low-quality cells before analysis.
- **slingshot.R**: Uses Slingshot for trajectory inference, modeling the differentiation paths of various cell types.

## Usage

After ensuring that you have all the required packages installed, you can use these scripts to process your single-cell RNA sequencing data.

1. **Preprocess the data**:
    - Run `Processing.R` to prepare your data by normalizing gene expression and filtering out unwanted cells.
    
2. **Perform analyses**:
    - Run `GSVA and survival analysis.R` to perform GSVA analysis and assess the correlation between gene signatures and survival outcomes.
    - Use `cellchat.R` to analyze cell-cell communication and `slingshot.R` for trajectory analysis.

3. **Visualize results**:
    - Generate heatmaps and plots by running `Dotplot.R` and `Heatmap.R`.
