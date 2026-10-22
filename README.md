# Lipid Atlas Demo (R Shiny)
**Live demo:** https://qiliu-chn.shinyapps.io/lipid-atlas-demo/
An interactive prototype for exploring quantitative lipidomics data: you can search lipid
species, compare concentrations between groups, run univariate and multivariate statistics,
and test a cross-validated classifier.

> **All data in this repository are simulated** (`R/simulate_data.R`). They do not contain real patient data.

## Features

| Tab | What it does |
|---|---|
| Explorer | Search for a lipid species by class or name. Shows a box plot by group (optionally split by sex), a summary table, and a CSV download. |
| Volcano | Welch t-tests on log2 data with fold change (Cancer/Control) and a choice of Bonferroni or Benjamini–Hochberg correction. Shows an interactive volcano plot and a results table. |
| PCA | Log2 transform, mean-centring and Pareto scaling, followed by a PCA scores plot and the top loadings. |
| Classifier | Logistic regression with k-fold cross-validation. Lipid selection happens *inside* each training fold to avoid information leakage. Shows a cross-validated ROC/AUC, sensitivity, specificity and accuracy, and a label-permutation test of the AUC. |

## Run locally

```r
install.packages(c("shiny", "ggplot2", "DT"))
shiny::runApp()          # run from the project folder
```

To regenerate the simulated data, run `source("R/simulate_data.R")`.

## Project structure

```
app.R                 # the Shiny app (UI + server)
R/simulate_data.R     # creates data/demo_lipids.csv
data/demo_lipids.csv  # simulated data: 300 samples x 37 lipid species
```

## Built with

R, shiny, ggplot2, DT
