############################################################
# Microbiome Analysis in R: Compositional Profiles
#
# Purpose:
# Generate taxonomic composition tables from the rarefied
# phyloseq object at the phylum, class, order, family, and
# genus levels. Both rarefied count tables and relative
# abundance tables are exported for downstream analyses.
############################################################


# ------------------------------
# 1. Set Working Directory
# ------------------------------
# Set the directory for compositional-profile analysis
# outputs.
setwd("~/Documents/HillTribe_NGS/3.Compositional_profile")


# -----------------------------------
# 2. Load Phyloseq Object & Packages
# -----------------------------------

# Load the phyloseq object generated during the
# phyloseq-object construction step.
load("~/Documents/HillTribe_NGS/1.Raw_data/Building_phyloseq.RData")

# Load packages required for data manipulation,
# visualization, statistical analysis, and file export.
library(phyloseq)
library(dplyr)
library(ggplot2)
library(rio)
library(car)  
library(readxl)




# -------------------------------------------
# 3. Export ASV Table with Taxonomic Information
# -------------------------------------------

# Extract the ASV count table, taxonomy table, and
# sample metadata from the phyloseq object.
asv_count_table <- as.data.frame(otu_table(ps_rare))
taxonomy_table <- as.data.frame(tax_table(ps_rare))
sample_metadata <- as(sample_data(ps_rare), "data.frame")


# Calculate total counts for each ASV across all samples
# and sort ASVs from highest to lowest total abundance.
asv_count_table$TotalCounts <- rowSums(asv_count_table)
asv_taxonomy_table_sorted <- asv_count_table[order(-asv_count_table$TotalCounts), ]
#rownames(asv_taxonomy_table_sorted) <- NULL
asv_taxonomy_table_sorted <- data.frame(ASV_ID = rownames(asv_taxonomy_table_sorted),
                                        asv_taxonomy_table_sorted)

# Inspect the sorted ASV table.
head(asv_taxonomy_table_sorted)

# Add the corresponding taxonomic annotations to the
# ASV count table.
asv_taxonomy_table_sorted <- asv_taxonomy_table_sorted %>% left_join(taxonomy,
                                        by = c("ASV_ID" = "Feature.ID"))

# Export the ASV count table with taxonomic annotations.
Export(asv_taxonomy_table_sorted, "asv_taxonomy_table_sorted.xlsx")


# Define the standard taxonomic prefixes used in the
# taxonomy annotations for each taxonomic rank.
prefix_map <- c(
  Phylum = "p__",
  Class = "c__",
  Order = "o__",
  Family = "f__",
  Genus = "g__"
)






# ----------------------------------------------------------
# Helper Function: Summarize and Export Taxonomic Profiles
# ----------------------------------------------------------
# Purpose:
# Aggregate ASV-level counts at a specified taxonomic rank,
# merge the resulting table with sample metadata, and export
# both rarefied counts and sample-level relative abundances.
process_taxonomic_level <- function(level, prefix, 
                                    asv_count_table, 
                                    taxonomy_table, 
                                    sample_metadata, 
                                    prefix_map) {
  message("\nProcessing: ", level)
  
  # Assign the selected taxonomic rank to each ASV based on
  # the corresponding feature ID.
  asv_count_table[[level]] <- taxonomy_table[[level]][match(rownames(asv_count_table), rownames(taxonomy_table))]
  
  # Aggregate ASV counts by the selected taxonomic rank.
  # Counts of ASVs assigned to the same taxon are summed
  # within each sample.
  tax_table_level <- asv_count_table %>%
    group_by(.data[[level]]) %>%
    summarise(across(where(is.numeric), ~sum(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(TotalCounts = rowSums(across(where(is.numeric)))) %>%
    arrange(desc(TotalCounts)) %>%
    select(-TotalCounts)
  
  # Transpose the taxonomic count table so that samples
  # are represented as rows and taxa as columns.
  tax_table_t <- as.data.frame(t(tax_table_level[,-1]))
  colnames(tax_table_t) <- tax_table_level[[level]]
  tax_table_t$Sample_ID <- rownames(tax_table_t)
  
  # Replace missing taxonomic assignments with "unclassified".
  colnames(tax_table_t)[is.na(colnames(tax_table_t))] <- "unclassified"
  
  
  # Merge the taxonomic abundance table with sample metadata.
  tax_with_metadata <- sample_metadata %>%
    left_join(tax_table_t, by = "Sample_ID")
  
  # Remove the standard taxonomic prefix from taxon names
  # for cleaner column names.
  prefix_to_remove <- prefix_map[[level]]
  colnames(tax_with_metadata) <- gsub(paste0("^", prefix_to_remove), "", colnames(tax_with_metadata))
  
  # Replace missing or empty column names with "unclassified".
  colnames(tax_with_metadata) <- ifelse(
    is.na(colnames(tax_with_metadata)) | colnames(tax_with_metadata) == "",
    "unclassified",
    colnames(tax_with_metadata)
  )
  
  # Ensure all taxon column names are unique, including
  # multiple unclassified taxa.
  colnames(tax_with_metadata) <- make.unique(colnames(tax_with_metadata), sep = "_")
  
  
  # Export the aggregated rarefied count table.
  Export(tax_with_metadata, paste0(prefix, "_Counts_", level, ".xlsx"), row.names = FALSE)
  
  # Calculate relative abundance for each taxon within
  # each sample.
  tax_with_metadata_rel <- tax_with_metadata %>%
    rowwise() %>%
    mutate(across(
      .cols = -(Sample_ID:DBP), # Metadata columns
      .fns = ~ .x / sum(c_across(-(Sample_ID:DBP))),
      .names = "{.col}"
    )) %>%
    ungroup()
  
  # Export the sample-level relative abundance table.
  Export(tax_with_metadata_rel, paste0(prefix, "_Rel_", level, ".xlsx"), row.names = FALSE)
}


# Define the taxonomic ranks to be processed and their
# corresponding output prefixes.
levels <- c("Phylum", "Class", "Order", "Family", "Genus")
prefixes <- c("1", "2", "3", "4", "5")

# Process and export each taxonomic rank using the helper
# function defined above.
for (i in seq_along(levels)) {
  process_taxonomic_level(levels[i], 
                          prefixes[i], 
                          asv_count_table = asv_count_table,
                          taxonomy_table = taxonomy_table,
                          sample_metadata = sample_metadata,
                          prefix_map = prefix_map)
}
