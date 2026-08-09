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


# ------------------------------
# Load Required Packages
# ------------------------------

# Load packages required for data import, manipulation,
# phylogenetic tree handling, and downstream microbiome analyses.
library(phyloseq)
library(ape)
library(readxl)
library(dplyr)
library(tibble)
library(rio)
library(car)




# ------------------------------
# 2. Import Input Files
# ------------------------------

# Set the working directory to the folder containing the raw
# microbiome sequencing data and associated input files.
setwd("~/Documents/HillTribe_NGS/1.Raw_data")


# Define the path to the rooted phylogenetic tree.
tree_file      <- "tree.nwk"              # Rooted phylogenetic tree

# Import the feature abundance table.
# The first row is skipped because the input table contains an
# additional header/comment line. check.names = FALSE preserves
# the original feature or sample names.
feature_table   <- read.delim("feature-table.tsv",
                              skip = 1,
                              header = TRUE,
                              check.names = FALSE
                              )

# Import the raw taxonomy table.
# comment.char = "#" ignores lines beginning with "#" in the input file.
tax_table_raw   <- read.delim("~/Documents/HillTribe_NGS/1.Raw_data/taxonomy.tsv", comment.char="#")

# Import the sample metadata table.
metadata        <- read.delim("~/Documents/HillTribe_NGS/1.Raw_data/metadata.txt")

# Import the rooted phylogenetic tree.
phy_tree        <- read_tree(tree_file)




# ------------------------------
# Initial Data Inspection
# ------------------------------

# Inspect the first few rows of each imported dataset to verify
# that the files were imported correctly and that their structure
# is consistent with the expected input format.
head(feature_table)
head(tax_table_raw)
head(metadata)



# ----------------------------------------------------------
# 3. Prepare Feature Table for Phyloseq
# ----------------------------------------------------------

# Load required package for data-frame and tibble operations.
library(tibble)



# ----------------------------------------------------------
# Clean Feature Table
# ----------------------------------------------------------

# Rename the OTU identifier column to a standardized name.
# This identifier will later be used as the row name of the
# feature table when constructing the phyloseq object.
feature_table <- feature_table %>%
  dplyr::rename(ASV_ID = `#OTU ID`)

# Convert the imported feature table to a standard data frame
# for subsequent manipulation and phyloseq preparation.
feature_table_mat <- feature_table %>% as.data.frame()

# Inspect column names to verify the structure of the feature table.
colnames(feature_table_mat)

# Set ASV_ID as row names, as required for constructing the
# phyloseq OTU abundance table.
rownames(feature_table_mat) <- feature_table_mat$ASV_ID


# Calculate the total number of sequencing reads (library size)
# for each sample and display samples in descending order.
# This provides an initial assessment of sequencing depth
# across samples.
feature_table_mat %>%
  summarise(across(where(is.numeric), sum)) %>%
  tidyr::pivot_longer(everything(),
                      names_to = "SampleID",
                      values_to = "LibrarySize") %>%
  arrange(desc(LibrarySize)) %>%
  print(n = 130) # <- you may adjust the number (n) according to your sample size.



# ----------------------------------------------------------
# Extract Sample Order from Feature Table
# ----------------------------------------------------------

# Extract sample IDs from the feature table in their original
# column order. ASV_ID and taxonomy columns are excluded.
# The resulting order is used to align the sample metadata
# with the feature table.
sample_order <- colnames(feature_table_mat[2:127])
sample_order



# ----------------------------------------------------------
# Verify Sample Order Consistency
# ----------------------------------------------------------

# Check whether the sample order in the feature table is
# identical to the sample order currently present in the metadata.
# FALSE indicates that the two datasets are not in the same order.
identical(sample_order, metadata$Smaple_ID)
# [FALSE]



# ----------------------------------------------------------
# Reorder Metadata to Match Feature Table
# ----------------------------------------------------------

# Reorder metadata rows according to the sample order in the
# feature table. Matching is performed using the Sample_ID column.
# This ensures that each metadata row corresponds to the correct
# sample in the feature table.
metadata <- metadata %>%
  slice(match(sample_order, Sample_ID)) 

# Inspect the reordered metadata to verify the result.
head(metadata)



# ----------------------------------------------------------
# Assign Sample IDs as Row Names
# ----------------------------------------------------------

# Set Sample_ID as the row names of the metadata table.
# This allows sample identifiers to be used to match metadata
# with the corresponding samples in the phyloseq object.
rownames(metadata) <- metadata$Sample_ID

head(metadata)



# ----------------------------------------------------------
# Final Validation
# ----------------------------------------------------------

# Confirm that the sample order in the feature table now
# exactly matches the sample order in the metadata.
# TRUE confirms that the two datasets are aligned.
identical(sample_order, metadata$Sample_ID)
# [TRUE]




# ----------------------------------------------------------
# 4. Prepare Taxonomy Table for Phyloseq
# ----------------------------------------------------------

# ----------------------------------------------------------
# Select and Format Taxonomy Columns
# ----------------------------------------------------------

# Retain only the feature identifier and the corresponding
# taxonomy annotation from the raw taxonomy table.
taxonomy <- tax_table_raw[, c("Feature.ID", "Taxon")]

# Convert the taxonomy annotation to character format to ensure
# consistent handling during taxonomy parsing.
taxonomy$Taxon <- as.character(taxonomy$Taxon)

# Inspect the formatted taxonomy table.
head(taxonomy)



# ----------------------------------------------------------
# Split Taxonomy String into Ranks
# ----------------------------------------------------------

# Split each taxonomy string into individual taxonomic ranks
# using the semicolon (;) as the delimiter.
taxonomy_split <- strsplit(taxonomy$Taxon, ";")

# Inspect the resulting list of taxonomic assignments.
head(taxonomy_split)



# ----------------------------------------------------------
# Construct Taxonomy Matrix
# ----------------------------------------------------------

# Convert the list of taxonomic assignments into a matrix with
# six taxonomic ranks (Kingdom, phylum, class, order, family, and genus). 
# Missing ranks are retained as NA to ensure
# that all features have the same number of taxonomic levels.
taxonomy_matrix <- do.call(rbind, lapply(taxonomy_split, function(x) {
  length(x) <- 6
  return(x)
}))

# Assign feature IDs as row names so that taxonomy assignments
# can be matched to the corresponding features in the abundance table.
rownames(taxonomy_matrix) <- taxonomy$Feature.ID

# Inspect the resulting taxonomy matrix.
head(taxonomy_matrix)



# ----------------------------------------------------------
# Validate Feature–Taxonomy Alignment
# ----------------------------------------------------------

# Confirm that the feature identifiers and their order are
# identical between the taxonomy matrix and feature table.
# TRUE confirms that each feature has a corresponding taxonomy
# assignment in the same order.
identical(rownames(taxonomy_matrix), rownames(feature_table_mat))
# [TRUE]




# ------------------------------
# 5. Build the Phyloseq Object
# ------------------------------

# Combine the feature abundance table, taxonomy matrix,
# sample metadata, and rooted phylogenetic tree into a single
# phyloseq object for downstream microbiome analyses.
ps <- merge_phyloseq(otu_table(feature_table_mat[-1],  taxa_are_rows = TRUE), 
                     tax_table(taxonomy_matrix),
                     sample_data(metadata),
                     phy_tree(phy_tree))

# Inspect the resulting phyloseq object to confirm that all
# components have been successfully integrated.
ps




# Assign standard taxonomic rank names to the columns of the
# taxonomy table, from Kingdom through Genus.
colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", 
                             "Order", "Family", "Genus")

# Confirm the taxonomic rank names assigned to the phyloseq object.
rank_names(ps)


# Summarize the number of features per sample to assess sequencing
# depth (library size) across samples.
summary(sample_sums(ps))




# ----------------------------------------------------------
# Remove Unassigned Taxa and Potential Contaminants
# ----------------------------------------------------------

# Remove features that are unassigned at the Kingdom level,
# classified as Archaea, lack a Phylum or Order assignment,
# or are annotated as chloroplast or mitochondrial sequences.
#
# These filtering criteria are applied to reduce non-target
# sequences and retain taxa relevant to the bacterial microbiome
# analysis.
ps_filtered <- ps %>%
  subset_taxa(
    Kingdom != "Unassigned" & 
      Kingdom != "k__Archaea" &
      Phylum != "" &
      Order != " o__Chloroplast" &
      Family != " f__Mitochondria"
  )

# Inspect the filtered phyloseq object.
ps_filtered



# ------------------------------
# 6. Rarefy Sequencing Depth
# ------------------------------

# Calculate the sequencing depth (library size) for each sample
# before rarefaction.
library_sizes <- sample_sums(ps_filtered)
summary(library_sizes)

# Rarefy all samples to the minimum observed sequencing depth.
# This randomly subsamples reads without replacement so that all
# samples have an equal number of reads for downstream analyses
# that require equal sequencing depth.
ps_rare <- rarefy_even_depth(
  ps_filtered,
  sample.size = min(sample_sums(ps_filtered)),  # Minimum sequencing depth across samples
  rngseed = 123,                                # Set random seed for reproducibility
  replace = FALSE,                              # Sample reads without replacement
  verbose = TRUE
)

# Confirm that all samples have the same sequencing depth after
# rarefaction.
summary(sample_sums(ps_rare))




# ------------------------------
# 7. Explore Phyloseq Object
# ------------------------------

# Inspect the basic components and structure of the rarefied
# phyloseq object.
sample_names(ps_rare)[1:5]
tax_table(ps_rare)[1:5, ]
otu_table(ps_rare)[1:5, ]
rank_names(ps_rare)
sample_variables(ps_rare)
summary(sample_sums(ps_rare))


# Check for missing values in the taxonomy table.
# This provides an overview of the completeness of taxonomic
# assignments across all taxonomic ranks.
table(is.na(tax_table(ps_rare)))




# ------------------------------
# 8. Summarize Taxonomy Table
# ------------------------------

# Convert the phyloseq taxonomy table to a standard data frame
# for summarization and export.
taxa_df <- as.data.frame(tax_table(ps_rare))

# Count the number of unique assigned taxa at each taxonomic rank.
# Missing taxonomic assignments are excluded from the count.
summary_table <- data.frame(
  Rank = colnames(taxa_df),
  Unique_Taxa = sapply(taxa_df, function(x) length(unique(na.omit(x))))
)

# Display the number of unique taxa identified at each rank.
print(summary_table)





# ------------------------------
# 9. Export Taxonomy Tables
# ------------------------------

# Add feature IDs to the taxonomy table for easier identification
# of individual features in the exported file.
taxa_df_asv <- data.frame(ASV_ID = rownames(taxa_df), taxa_df)

# Combine feature IDs, taxonomy assignments, and the corresponding
# feature abundance table into a single data frame.
taxa_rank_asv <- data.frame(
  ASV_ID = rownames(taxa_df),
  taxa_df,
  as.data.frame(otu_table(ps_rare))
)

# Export the taxonomy-only and taxonomy-plus-abundance tables
# as Excel files for downstream inspection or reporting.
Export(taxa_df_asv, "taxa_df_asv.xlsx")
Export(taxa_rank_asv, "taxa_rank_asv.xlsx")




# ------------------------------
# 10. Rarefaction Curve
# ------------------------------

# Load vegan for rarefaction curve analysis.
library(vegan)

# Extract the feature abundance table from the filtered phyloseq object
# and convert it to a standard matrix.
otu <- as(otu_table(ps_filtered), "matrix")

# vegan::rarecurve() expects samples as rows and taxa/features as columns.
# Transpose the matrix when the phyloseq object stores taxa as rows.
if (taxa_are_rows(ps_filtered)) {
  
  otu <- t(otu)
  
}

# Check the dimensions of the resulting sample-by-feature matrix.
dim(otu)



# Define colors for each study group.
group_colors <- c(
  "Akha-CM" = "gold1",
  "Lahu-CM" = "lightblue2",
  "Khuen-CM" = "darkorange",
  
  "Akha-CR" = "dodgerblue4",
  "Lahu-CR" = "coral1",
  "Lisu-CR" = "grey40"
)

# Create a combined group variable from ethnicity and study area.
# The resulting labels are used to assign colors to samples in
# the rarefaction curve.
metadata$Group <- paste(metadata$Ethnicity, metadata$Area, sep = "-")


# Match metadata-derived group labels to the sample order in the
# abundance matrix. This ensures that each sample receives the
# correct group color in the rarefaction plot.
group <- metadata$Group[
  match(rownames(otu), metadata$Sample_ID)
  ]

# Confirm that the group labels are correctly aligned with the
# sample order in the abundance matrix.
identical(group, metadata$Group)


# Generate rarefaction curves for all samples.
# Curves are evaluated at intervals of 100 sequencing reads and
# extended to the minimum library size across samples.
# Each curve is colored according to the corresponding study group.
rarecurve(
  otu,
  step = 100,
  label = FALSE,
  sample = min(rowSums(otu)),
  xlab = "Sequencing depth",
  ylab = "Observed ASVs",
  col = group_colors[group],
  grid = FALSE
)


# Add a legend identifying the study groups represented by each color.
legend(
  "bottomright",
  legend = names(group_colors),
  col = group_colors,
  lwd = 3,
  bty = "n",
  cex = 0.7
)
