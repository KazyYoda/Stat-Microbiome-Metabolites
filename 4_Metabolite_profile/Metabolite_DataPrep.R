############################################################
# Metabolite Analysis
############################################################
#
# Purpose:
# Prepare the metabolomics dataset for downstream statistical
# analysis by integrating metabolite abundance data with
# compound identification and chemical classification
# annotations.
#
# Main workflow:
#   1. Load required packages and input data
#   2. Clean and prepare metabolite abundance data
#   3. Integrate PubChem and ClassyFire annotations
#   4. Remove known analytical contaminants
#   5. Align metabolite annotations with abundance data
#   6. Compare compound composition between datasets
#   7. Import and align sample metadata
#   8. Convert the metabolite matrix to sample-wise format
#   9. Verify sample and compound ordering
#  10. Export metabolite annotation information
############################################################


# ============================
# 1. Setup
# ============================

# Set the working directory for metabolomics analyses.
setwd("~/Documents/HillTribe_NGS/5.Metabolites")


# Load required packages.
library(dplyr)
library(tidyr)
library(car)
library(readxl)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggh4x)


#========================================================
# 2. Load and Prepare Metabolomics Data
#========================================================

#--------------------------------------------------------
# 2.1 Import Compound Mapping and Annotation Tables
#--------------------------------------------------------
#
# Import compound identification results from PubChem mapping
# and chemical classification results from ClassyFire for the
# 0.05 dataset.

PubMed_05 <- read_excel("PubMed_0.5.xlsx")
ClassyFireBatch_05 <- read_excel("ClassyFireBatch_0.5.xlsx")


#--------------------------------------------------------
# 2.2 Import Raw Metabolite Abundance Tables
#--------------------------------------------------------
#
# Import the metabolite abundance table used for downstream
# statistical analysis.

Metabolites_05 <- read_excel("Metabolites_0.05.xlsx")


#--------------------------------------------------------
# 2.3 Clean Column Names and Remove Unnecessary Fields
#--------------------------------------------------------
#
# Sample identifiers contain the "-AllHits" suffix, which is
# removed to match the sample identifiers used in the metadata.
#
# Identification-related columns are removed because only the
# metabolite abundance measurements are required for downstream
# statistical analysis.

# Remove "-AllHits" suffix from sample IDs
colnames(Metabolites_05) <- gsub("-AllHits", "", colnames(Metabolites_05))

# Remove identification-related columns and retain the
# metabolite abundance matrix.
Metabolites_05 <- Metabolites_05 %>% 
  select(-Mass:-`METLIN ID`)



#--------------------------------------------------------
# 2.4 Define and Remove Contaminants
#--------------------------------------------------------
#
# Define known derivatization products and background
# contaminants identified during the metabolite-processing
# workflow. These compounds are excluded from downstream
# analyses.

contaminants <- c(
  "Bis(heptamethylcyclotetrasiloxy)siloxane",
  "Cyclotrisiloxane, hexamethyl-",
  "2-Vinylphenol",
  "Dimethyl phthalate",
  "Naphthalene, 2-methyl-",
  "Acetonitrile, trifluoro-"
)



#--------------------------------------------------------
# 2.5 Merge Annotation Sources
#--------------------------------------------------------
#
# Combine PubChem compound-identification information with
# ClassyFire chemical classification using the cleaned
# compound name and InChIKey.
#
# Missing annotation values are explicitly labeled as
# "Not found" to distinguish missing annotations from
# successfully identified compounds.

Pub_Classy_05 <- PubMed_05 %>%
  left_join(ClassyFireBatch_05,
            by = c("Compound_Clean", "InChIKey")) %>%
  mutate(across(where(is.character),
                ~ replace_na(.x, "Not found")))



#--------------------------------------------------------
# 2.6 Clean Metabolite Abundance Table (0.05)
#--------------------------------------------------------
#
# Remove known contaminants from the abundance table.
# Assign a unique internal metabolite code (Met1, Met2, ...)
# to each remaining compound for convenient reference in
# downstream analyses.

Metabolites_05_clean <- Metabolites_05 %>%
  filter(!Compound %in% contaminants) %>%
  mutate(Code = paste0("Met", row_number())) %>%
  relocate(Code, .before = 2)


#--------------------------------------------------------
# 2.7 Align Annotation Order to Abundance Table
#--------------------------------------------------------
#
# Remove the same contaminants from the annotation table and
# reorder the annotation records to exactly match the compound
# order in the cleaned abundance table.
#
# The same internal metabolite codes are then assigned to the
# annotation table.

Pub_Classy_05_clean <- Pub_Classy_05 %>%
  filter(!Compound %in% contaminants) %>%
  slice(match(Metabolites_05_clean$Compound, Compound)) %>%
  mutate(Code = paste0("Met", row_number())) %>%
  relocate(Code, .before = 3)


# Verify that compounds occur in the same order in the
# abundance and annotation tables.
identical(Metabolites_05_clean$Compound, Pub_Classy_05_clean$Compound)




#========================================================
# 3. Load Metadata and Create Group Variable
#========================================================
#
# Import sample metadata and create a combined group variable
# representing ethnicity and geographic area.

sample_metadata <- read_excel("sample_metadata.xlsx")

# Create combined group label (Ethnicity-Area)
sample_metadata <- sample_metadata %>%
  mutate(Group = paste(Ethnicity, Area, sep = "-")) %>%
  relocate(Group, .before = 2)



#========================================================
# 3.1 Identify Missing Samples (Metabolite_05 Dataset)
#========================================================
#
# Identify samples listed in the metadata that are not present
# in the cleaned metabolite abundance table. These samples are
# excluded from the metadata used for downstream metabolomics
# analyses.

# Extract Sample_IDs present in the metabolite table
sample_met05 <- names(Metabolites_05_clean[, -c(1:2)])

# Identify samples absent from the metabolite dataset
# due to insufficient material
drop_sampleID <- sample_metadata %>%
  filter(!Sample_ID %in% sample_met05) %>%
  pull(Sample_ID)

# Remove missing samples from metadata
sample_metadata_drop <- sample_metadata %>%
  filter(!Sample_ID %in% drop_sampleID) 

# Check the number of samples remaining in each group
sample_metadata_drop %>% 
  group_by(Group) %>% 
  count()


#========================================================
# 3.2 Align Metadata Order to Metabolite Table
#========================================================
#
# Reorder the metadata so that sample order exactly matches
# the sample order in the metabolite abundance matrix.
#
# Matching sample order is essential because subsequent
# analyses assume that each row of the metabolite matrix
# corresponds to the same sample in the metadata.

Met05_match <- match_metadata_order(sample_metadata_drop, sample_met05)

# Verify identical Sample_ID order
identical(sample_met05, Met05_match$Sample_ID)
# [TRUE]


#----------
# Notes:
#----------
match_metadata_order() # can be found in HelperFn_Metabolite.R



#========================================================
# 3.3 Transpose Metabolite Matrix
#========================================================
#
# Convert the metabolite abundance table from a compound-wise
# format to a sample-wise format:
#
#   Before: compounds = rows, samples = columns
#   After:  samples   = rows, compounds = columns
#
# The internal metabolite codes are used as column names.

# Transpose metabolite abundance matrix (samples as rows)
metabo_05_t <- data.frame(t(Metabolites_05_clean[, -c(1:2)]))

# Assign metabolite codes as column names
colnames(metabo_05_t) <- Metabolites_05_clean$Code

# Confirm row order matches metadata
all(rownames(metabo_05_t) == Met05_match$Sample_ID)



#========================================================
# 3.4 Export Metabolite Annotation Table
#========================================================
#
# Extract the primary metabolite annotation fields and export
# them as a separate table for reference during downstream
# analyses and interpretation.

# Extract metabolite annotation columns
metabo_05_descp <- data.frame(
  Pub_Classy_05_clean[1:3]
)

# Export metabolite description table
Export(metabo_05_descp, "metabo_05_descp.txt")

