############################################################
# Microbiome Analysis in R: Building a Phyloseq Object
############################################################

# This script prepares microbiome sequencing data and constructs
# a phyloseq object for downstream microbiome analyses.
#
# Required input files:
#   1. feature_table.tsv : Feature/ASV abundance table
#   2. taxonomy.tsv      : Taxonomic annotation table
#   3. metadata          : Sample metadata
#   4. tree.nwk          : Rooted phylogenetic tree
#
# Note:
# The input tables should contain consistent sample and feature
# identifiers to ensure that the corresponding data can be merged
# correctly when constructing the phyloseq object.


# ------------------------------
# 1. Install Required Packages
# ------------------------------

# Install BiocManager if it is not already installed.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install rhdf5 from Bioconductor.
# rhdf5 is required by some Bioconductor packages and dependencies.
BiocManager::install("rhdf5")

# Install phyloseq from Bioconductor.
# phyloseq provides the data structure used to integrate
# feature abundance, taxonomy, sample metadata, and phylogenetic
# information for microbiome analyses.
BiocManager::install("phyloseq")


# Install required packages from CRAN.
# These packages are used for data manipulation, file import,
# phylogenetic tree handling, and statistical analyses.
install.packages(c("ape", "dplyr", "readxl", "tibble", "rio", "car"))
