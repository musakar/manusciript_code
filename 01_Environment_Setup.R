# -------------------------------------------------- #
# Script 01: Environment Setup
# Purpose: Install all required R packages for the analysis.
# This script should be run once on any machine before the main analysis.
# -------------------------------------------------- #

# Increase the timeout limit for downloads
options(timeout = 600)

# Install BiocManager if it's not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Define the list of all required packages
required_packages <- c(
  "affy",             # For reading .CEL files
  "Biobase",          # Core bioinformatics data structures
  "sva",              # For ComBat batch effect correction
  "limma",            # For differential gene expression analysis
  "tidyverse",        # For data manipulation and ggplot2
  "ggplot2",          # For plotting
  "ggrepel",          # For non-overlapping plot labels
  "patchwork",        # For combining ggplot plots
  "UpSetR",           # For generating UpSet plots
  "clusterProfiler",  # For GO and KEGG enrichment analysis
  "AnnotationDbi",    # Base for annotation packages
  "org.At.tair.db",   # General annotation database for Arabidopsis
  "ath1121501.db",    # Chip-specific annotation for ID mapping
  "ath1121501cdf",    # Chip definition file for RMA normalization
  "knitr"             # For creating formatted tables
)

# Install all packages
BiocManager::install(required_packages)

cat("\nEnvironment setup complete. All required packages are installed or updated.\n")
cat("You can now run the '02_Analysis_and_Figures.R' script.\n")
