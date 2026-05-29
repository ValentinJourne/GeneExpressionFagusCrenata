# Gene Expression Analysis of Mast Flowering in *Fagus crenata*

This repository contains the code used to analyse seasonal transcriptomic variation and identify candidate genes associated with mast flowering in *Fagus crenata*.

The workflow combines gene expression data, flowering records, climate variables, and nutrient measurements to investigate the molecular mechanisms underlying mast flowering.

## Repository structure

### Main analysis

**Full_analysis_gene_expression_v2026.R**

Main script used to reproduce all analyses presented in the manuscript, including:

* Processing and formatting transcriptomic datasets
* Calculation of gene synchrony and inter-annual variability metrics
* Logistic regression analyses linking gene expression to flowering
* Machine-learning feature selection using Boruta
* Gene ontology (GO) enrichment analyses
* Figure generation

### Functions

**functions_gene_exp_analysis_viz.R**

Collection of custom functions used throughout the analysis, including:

* Data formatting and cleaning
* Gene synchrony and variability calculations
* Climate window extraction
* Regression model fitting and summarisation
* GO enrichment analyses
* Plotting utilities

### Package management

**packages_loading.R**

Script used to load and install all required R packages prior to running the analyses.

## Data requirements

The analyses require:

* Gene expression data from leaf tissues collected between 2014 and 2023
* Flowering intensity records for individual trees
* Daily climate data (temperature, precipitation, sunshine duration)
* Nitrogen, sulfur, amino acid measurements

Raw sequencing data and associated metadata are available through the repositories listed in the manuscript.

## Reproducibility

Analyses were performed using R (version ≥ 4.3). To reproduce the analyses:

```r
source("packages_loading.R")
source("functions_gene_exp_analysis_viz.R")
source("Full_analysis_gene_expression_v2026.R")
```

## Citation

If you use this code, please cite:

Journé V., Satake A., et al. *Sulfur–nitrogen imbalance drives mass flowering in forest tree via FLC–FT pathway activation*.
