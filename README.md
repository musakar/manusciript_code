# Comparative Transcriptomic Analysis of Pesticide Stress in Arabidopsis

This repository contains the R scripts for the meta-analysis of plant transcriptomic responses to nine different pesticides.

## Description

This study performs an integrated analysis of nine publicly available microarray datasets to understand how pesticides with diverse Modes of Action (MoAs) affect the *Arabidopsis thaliana* transcriptome. We use a robust bioinformatics pipeline to correct for inter-experimental batch effects and identify both conserved core responses and compound-specific transcriptional signatures.

## Workflow

1.  **Data Acquisition**: Raw `.CEL` files for the 9 GSE series must be downloaded from the NCBI GEO database and placed inside the `data/` folder.
2.  **Environment Setup**: Run the `01_Environment_Setup.R` script in R to install all required packages.
3.  **Analysis**: Run the `02_Analysis_and_Figures.R` script to perform the entire analysis pipeline and generate all outputs in the `results/` folder.
