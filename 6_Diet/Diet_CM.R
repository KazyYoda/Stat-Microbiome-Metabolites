############################################################
# Diet Analysis
#
# Purpose:
# This script evaluates differences in dietary variables among
# study groups within the Chiang Mai (CM) area.
#
# Statistical workflow:
# 1. Import and prepare dietary and metadata datasets
# 2. Perform Kruskal-Wallis tests for each dietary variable
# 3. Apply Dunn's post-hoc test to variables with significant
#    Kruskal-Wallis results
# 4. Summarize group-level medians to support interpretation
#    of the direction of pairwise differences
# 5. Visualize significant pairwise differences using a bubble plot
############################################################


# ============================================================
# 1. Environment Setup
# ============================================================

# Set the working directory containing the dietary analysis files.
setwd("~/Documents/HillTribe_NGS/8.Diet")

# Load packages required for data import, manipulation,
# statistical analysis, and visualization.
library(dplyr)
library(tidyr)
library(car)
library(readxl)
library(readr)
library(ggplot2)


# ============================================================
# 2. Import and Prepare Data
# ============================================================

# Import the study metadata and dietary data.
# The metadata file contains sample-level information, while
# Diet_CM.xlsx contains dietary variables measured for samples
# from the Chiang Mai (CM) area.
metadata <- read_delim("metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Diet_CM <- read_excel("Diet_CM.xlsx")

# Restrict the metadata to samples from the Chiang Mai area.
metadata_CM <- metadata %>% filter(Area == "CM")

# Check whether any CM samples in the metadata are missing
# from the dietary dataset.
#
# An empty result indicates that all CM samples represented in
# the metadata are available in the dietary dataset.
setdiff(metadata_CM$Sample_ID, Diet_CM$Sample_ID)

# Retain dietary observations corresponding to samples present
# in the CM metadata.
#
# Fresh_Meat is excluded from the analysis, and the resulting
# object is converted to a data.frame for subsequent analyses.
dataDiet_CM <- Diet_CM %>%
  filter(Sample_ID %in% metadata_CM$Sample_ID) %>%
  select(-Fresh_Meat) %>%
  as.data.frame()

# Use Sample_ID as row names and convert grouping variables
# to factors for statistical analysis.
rownames(dataDiet_CM) <- dataDiet_CM$Sample_ID
dataDiet_CM$Group <- factor(dataDiet_CM$Group)
dataDiet_CM$Ethnicity <- factor(dataDiet_CM$Ethnicity)


# ============================================================
# 3. Statistical Analysis
# ============================================================

# Load packages required for the statistical testing workflow.
library(dplyr)
library(tidyr)
library(ggplot2)
library(FSA)


# ------------------------------------------------------------
# 3.1 Kruskal-Wallis Test for Each Dietary Variable
# ------------------------------------------------------------

# Perform an independent Kruskal-Wallis test for each dietary
# variable.
#
# The Kruskal-Wallis test is used to evaluate whether the
# distribution of a dietary variable differs among the CM groups.
#
# Variables are selected from columns 6 through the final column
# of dataDiet_CM.
krus_results <- lapply(names(dataDiet_CM[6:ncol(dataDiet_CM)]), function(var) {
  
  test <- kruskal.test(dataDiet_CM[[var]], dataDiet_CM$Group)
  
  data.frame(
    Variable = var,
    P_value = test$p.value
  )
})

# Combine the individual test results into a single data frame.
krus_df <- do.call(rbind, krus_results)

# Retain dietary variables showing a statistically significant
# overall difference among groups at the nominal significance
# threshold of P < 0.05.
#
# These variables are subsequently evaluated using pairwise
# post-hoc comparisons.
krus_sig <- krus_df %>%
  filter(P_value < 0.05)

krus_sig


# ------------------------------------------------------------
# 3.2 Dunn's Post-hoc Test
# ------------------------------------------------------------

# Perform pairwise Dunn's tests for dietary variables that showed
# a significant Kruskal-Wallis result.
#
# The Benjamini-Hochberg (BH) procedure is applied within each
# Dunn's test to adjust the pairwise P-values for multiple
# comparisons.
dunn_results <- lapply(krus_sig$Variable, function(var) {
  
  dunn_res <- dunnTest(
    dataDiet_CM[[var]] ~ Group,
    data = dataDiet_CM,
    method = "bh"
  )
  
  dunn_df <- as.data.frame(dunn_res$res)
  dunn_df$Variable <- var
  
  return(dunn_df)
})

# Combine the Dunn's test results for all significant dietary
# variables into a single data frame.
dunn_df <- do.call(rbind, dunn_results)

# Retain pairwise comparisons with BH-adjusted P-values below 0.05.
dunn_sig <- dunn_df %>%
  filter(P.adj < 0.05)

dunn_sig


# ------------------------------------------------------------
# 3.3 Calculate Group Medians
# ------------------------------------------------------------

# Reshape dietary data to long format so that each row represents
# one sample, dietary variable, and corresponding value.
#
# This format facilitates calculation of group-level summaries
# and can also be used for downstream visualization.
diet_long <- dataDiet_CM %>%
  pivot_longer(
    cols = Egg:Sausage,
    names_to = "Variable",
    values_to = "Value"
  )

# Calculate the median value of each dietary variable within
# each study group.
#
# Medians are used because the primary statistical analysis is
# non-parametric and is based on rank distributions.
group_medians <- diet_long %>%
  group_by(Group, Variable) %>%
  summarise(Median = median(Value, na.rm = TRUE), .groups = "drop")


# ============================================================
# 4. Bubble Plot of Significant Pairwise Comparisons
# ============================================================

# Define the order in which pairwise group comparisons are
# displayed in the plot.
dunn_sig$Comparison <- factor(dunn_sig$Comparison,
                              levels = c("Akha-CM - Lahu-CM",
                                         "Khuen-CM - Lahu-CM",
                                         "Akha-CM - Khuen-CM"))

# Generate a bubble plot showing significant pairwise
# differences identified by Dunn's test.
#
# Bubble size represents the strength of statistical evidence,
# expressed as -log10 of the BH-adjusted P-value.
#
# Bubble fill indicates the direction of the Dunn test statistic:
#   Z > 0  = higher value in the first group of the comparison
#   Z < 0  = lower value in the first group of the comparison
#
# Only pairwise comparisons meeting the adjusted significance
# threshold (P.adj < 0.05) are displayed.
ggplot(dunn_sig, aes(
  x = Comparison,
  y = Variable
)) +
  
  geom_point(aes(
    size = -log10(P.adj),
    fill = Z > 0
  ),
  shape = 21,
  color = "grey30",
  stroke = 0.3
  ) +
  
  scale_fill_manual(
    values = c(
      "TRUE" = "#FC766A",
      "FALSE" = "steelblue3"
    ),
    labels = c(
      "FALSE" = "Lower in first group",
      "TRUE" = "Higher in first group"
    ),
    name = ""
  ) +
  
  scale_size_continuous(
    name = expression(-log[10](q-value)),
    range = c(1, 5)
  ) +
  
  theme_bw() +
  
  labs(
    x = "Group comparison",
    y = "Diet variable",
    title = "Dietary pattern differences across CM groups"
  ) +
  
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(
      size = 6,
      angle = 45,
      hjust = 1,   # helps align rotated labels
      vjust = 1
    ),
    axis.title = element_text(size = 6),
    plot.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    panel.grid = element_blank()
  )







############################################################
# CM Location
# Genus-Level Feature Processing, Diet–Microbiota Correlation,
# CLR Transformation, and HAllA Input Preparation
#
# Workflow:
#   1. Load genus abundance and diet data.
#   2. Match samples across genus, diet, and metadata tables.
#   3. Perform group-specific Spearman correlations between
#      dietary variables and genus-level features.
#   4. Prepare genus-level data for CLR transformation.
#   5. Generate CLR-transformed genus tables for each CM group.
#   6. Prepare diet and microbiota tables for HAllA analysis.
############################################################


# ==========================================================
# 1. Load Genus Abundance Data
# ==========================================================

# Import genus-level count and relative abundance tables.
Counts_Genus <- read_excel("~/Documents/HillTribe_NGS/3.Compositional_profile/5_Counts_Genus.xlsx")
Rel_Genus <- read_excel("~/Documents/HillTribe_NGS/3.Compositional_profile/5_Rel_Genus.xlsx")

# Retain only samples that are present in the CM diet dataset.
# Rel_Genus is used as the input for the downstream genus-level
# analyses in this section.
Counts_Genus_CM <- Rel_Genus %>% # Counts_Genus or Rel_Genus data
  filter(Sample_ID %in% dataDiet_CM$Sample_ID)

# Convert the filtered genus abundance table to a data frame.
Genus_abs_CM <- Counts_Genus_CM %>%
  as.data.frame()

# Assign Sample_ID as row names for sample-level matching.
rownames(Genus_abs_CM) <- Counts_Genus_CM$Sample_ID


# ==========================================================
# 2. Align Sample Order Across Data Tables
# ==========================================================

# Store the sample order from the genus abundance table.
# All subsequent datasets are reordered to this same sequence.
sample_order <- Genus_abs_CM$Sample_ID

# Reorder the diet dataset to match the genus abundance table.
# Matching samples by Sample_ID ensures that each row represents
# the same biological sample across datasets.
dataDiet_CM <- dataDiet_CM %>%
  slice(match(sample_order, Sample_ID)) 


# Confirm that sample order is identical between the
# genus abundance and diet datasets.
identical(rownames(Genus_abs_CM), dataDiet_CM$Sample_ID)

# Reorder metadata to the same sample order and construct
# the combined group identifier using Ethnicity and Area.
meta_CM <- metadata_CM %>%
  slice(match(sample_order, Sample_ID)) %>%
  mutate(Group = paste(Ethnicity, Area, sep = "-"))

# Confirm that sample order is consistent across the
# genus abundance, metadata, and diet datasets.
identical(rownames(Genus_abs_CM), meta_CM$Sample_ID)
identical(dataDiet_CM$Sample_ID, dataDiet_CM$Sample_ID)


# ==========================================================
# 3. Group-Specific Spearman Correlation Analysis
# ==========================================================

# Spearman correlation is used to evaluate monotonic associations
# between dietary variables and genus-level abundance.
library(dplyr)
library(tidyr)


# -----------------------------------
# 3.1 Define Diet and Genus Features
# -----------------------------------

# Identify dietary variables from the diet dataset.
diet_vars <- colnames(dataDiet_CM[6:ncol(dataDiet_CM)])

# Identify genus-level features from the genus abundance table.
# Update the column range if the structure of the input table changes.
genus_vars <- colnames(Genus_abs_CM)[13:100]  
# adjust if genus starts at another column


# -----------------------------------
# 3.2 Merge Diet and Genus Data
# -----------------------------------

# Combine dietary variables and genus-level features using Sample_ID.
# The join preserves the correspondence between dietary measurements
# and microbiota profiles from the same samples.
corr_data <- dataDiet_CM %>%
  select(Sample_ID, Group, all_of(diet_vars)) %>%
  left_join(
    Genus_abs_CM %>%
      select(Sample_ID, all_of(genus_vars)),
    by = "Sample_ID"
  )


# -----------------------------------
# 3.3 Calculate Spearman Correlations
#     Separately Within Each Group
# -----------------------------------

# Initialize a list for storing correlation results from each group.
spearman_results <- list()

# Perform all diet–genus correlations independently within each group.
for (grp in unique(corr_data$Group)) {
  
  # Subset samples belonging to the current group.
  group_data <- corr_data %>%
    filter(Group == grp)
  
  # Generate all possible combinations of diet variables
  # and genus-level features.
  group_results <- expand.grid(
    Diet = diet_vars,
    Genus = genus_vars,
    stringsAsFactors = FALSE
  ) %>%
    rowwise() %>%
    mutate(
      
      # Calculate Spearman's rank correlation coefficient (rho).
      # Complete observations are used for each feature pair.
      rho = suppressWarnings(
        cor(
          group_data[[Diet]],
          group_data[[Genus]],
          method = "spearman",
          use = "complete.obs"
        )
      ),
      
      # Calculate the corresponding statistical significance.
      # exact = FALSE is used for the Spearman test.
      p_value = suppressWarnings(
        cor.test(
          group_data[[Diet]],
          group_data[[Genus]],
          method = "spearman",
          exact = FALSE
        )$p.value
      ),
      
      # Record the group associated with each correlation.
      Group = grp
    ) %>%
    ungroup()
  
  # Apply Benjamini–Hochberg (BH) correction for multiple testing.
  # The adjustment is performed independently within each group.
  group_results$p_adj <- p.adjust(
    group_results$p_value,
    method = "BH"
  )
  
  # Store the results for the current group.
  spearman_results[[grp]] <- group_results
}


# -----------------------------------
# 3.4 Combine Results Across Groups
# -----------------------------------

# Combine the group-specific correlation results into one table.
spearman_df <- bind_rows(spearman_results)


# -----------------------------------
# 3.5 Extract Significant Associations
# -----------------------------------

# Retain diet–genus associations with BH-adjusted p-values < 0.05.
sig_spearman <- spearman_df %>%
  filter(p_adj < 0.05)

# View significant diet–genus associations.
sig_spearman


# ==========================================================
# 4. Format Genus Data for CLR Processing
# ==========================================================

# Transpose the genus abundance table so that:
#   - rows = genus-level features
#   - columns = samples
#
# The first 12 columns contain sample metadata and are excluded.
Genus_abs_CM_t <- as.data.frame(t(Genus_abs_CM[-c(1:12)]))  


# ==========================================================
# 5. Define Analysis Inputs
# ==========================================================

# Specify the taxonomic resolution used for the analysis.
rank_level <- c("Genus")

# Store the genus table in a list so that the same processing
# framework can be extended to additional taxonomic levels.
tables <- list(Genus = Genus_abs_CM_t)


# ==========================================================
# 6. Define Target CM Subgroups
# ==========================================================

# Define the CM groups to be processed independently.
# Each element contains the group label used in the metadata.
unique_group <- list(
  Akha_CM  = c("Akha-CM"),
  Lahu_CM  = c("Lahu-CM"),
  Khuen_CM = c("Khuen-CM")
)


# ==========================================================
# 7. Initialize CLR Processing Output
# ==========================================================

# Initialize an empty list for storing CLR-transformed
# genus-level tables generated for each subgroup.
results_list <- list()


# ==========================================================
# 8. CLR Transformation by CM Group
# ==========================================================

# Apply the CLR preparation function separately to each
# predefined CM subgroup.
#
# The prepare_clr_table() function is assumed to be defined
# in the previously loaded analysis workspace/script.
for (data_type in rank_level) {
  for (comp_name in names(unique_group)) {
    
    # Extract the group label for the current analysis.
    group_pair <- unique_group[[comp_name]]
    
    # Generate the CLR-transformed genus feature table.
    res <- prepare_clr_table(
      feature_table = tables[[rank_level]],
      sample_metadata = meta_CM,
      data_type = rank_level,
      group_col = "Group",
      group = group_pair,
      factor = "Eth"
    )
    
    # Store each result using a structured name indicating
    # the taxonomic level and subgroup.
    results_list[[paste(data_type, comp_name, sep = "_")]] <- res
    
    # Display progress during processing.
    message(paste("Completed:", data_type, comp_name))
  }
}


####------ Helper Fn -----#######
prepare_clr_table <- function(feature_table, 
                              sample_metadata, 
                              data_type,
                              group,
                              group_col = "Group", 
                              id_col = "Sample_ID",
                              factor = "Geo") {
  library(dplyr)
  library(compositions)
  
  message(paste("Processing:", data_type, "|", group))
  
  # Select samples for the group
  samples <- sample_metadata %>%
    filter(.data[[group_col]] %in% c(group)) %>%
    pull(.data[[id_col]])
  
  # Subset count table
  feature <- feature_table %>%
    dplyr::select(all_of(samples)) %>%
    t() %>%
    as.data.frame()
  
  feature_matrix <- as.matrix(feature)
  
  if (any(feature_matrix == 0, na.rm = TRUE)) {
    warning("Zero values detected; applying global pseudocount.")
  }
  
  # Global pseudocount for consistency
  feature_matrix <- feature_matrix + 1
  
  # CLR transformation
  clr_matrix <- clr(feature_matrix)
  clr_tbl <- t(clr_matrix)
  clr_tbl <- as.data.frame(clr_tbl)
  
  
  # Optional output file writing
  file_name <- paste0("feature_", data_type, "_", group, "_", factor, ".txt")
  
  write.table(clr_tbl,
              file = file_name,
              sep = "\t",
              row.names = TRUE,
              col.names = NA,
              quote = FALSE)
  
  return(clr_tbl)
}





# ==========================================================
# 9. Prepare Diet Data for HAllA
# ==========================================================

# Prepare the complete CM diet table for HAllA.
# The table is transposed so that:
#   - rows = diet variables
#   - columns = samples
#
# Sample identifiers and metadata columns through Age are excluded.
diet_HAllA <- dataDiet_CM %>% 
  select(-Sample_ID:-Age) %>%
  t() %>%
  as.data.frame()

# Prepare the diet table specifically for Akha-CM samples.
diet_Akha <- dataDiet_CM %>% 
  filter(Group == "Akha-CM") %>%
  select(-Sample_ID:-Age) %>%
  t() %>%
  as.data.frame()

# Prepare the diet table specifically for Khuen-CM samples.
diet_Khuen <- dataDiet_CM %>% 
  filter(Group == "Khuen-CM") %>%
  select(-Sample_ID:-Age) %>%
  t() %>%
  as.data.frame()

# Prepare the diet table specifically for Lahu-CM samples.
diet_Lahu <- dataDiet_CM %>% 
  filter(Group == "Lahu-CM") %>%
  select(-Sample_ID:-Age) %>%
  t() %>%
  as.data.frame()


# ==========================================================
# 10. Export Diet Tables for HAllA
# ==========================================================

# Export each group-specific diet table as a tab-delimited file.
# These files are used as the X input for HAllA.
write.table(
  diet_Akha,
  "diet_Akha.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

write.table(
  diet_Khuen,
  "diet_Khuen.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

write.table(
  diet_Lahu,
  "diet_Lahu.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)


# ==========================================================
# 11. HAllA Analysis
# ==========================================================

# HAllA is used to identify associations between dietary
# variables and CLR-transformed genus-level microbiota features.
#
# Method:
#   - Association measure: Spearman correlation
#   - Multiple-testing threshold: FDR alpha = 0.05
#
# The following commands are run externally from the terminal
# or command line using the HAllA software.

# HAllA Diet CM
# halla -x diet_Akha.tsv -y feature_Genus_Akha-CM_Eth.txt -o Diet_AkhaCM -m spearman --fdr_alpha 0.05 
# halla -x diet_Khuen.tsv -y feature_Genus_Khuen-CM_Eth.txt -o Diet_KhuenCM -m spearman --fdr_alpha 0.05 
# halla -x diet_Lahu.tsv -y feature_Genus_Lahu-CM_Eth.txt -o Diet_LahuCM -m spearman --fdr_alpha 0.05 









###################################################################
# CM Ethnic Group
# Diet–Metabolite Association Analysis
# Project: Hill Tribe NGS – Gut Microbiota
#
# Purpose:
#   Prepare metabolite abundance data for each ethnic group in the
#   Chiang Mai (CM) area and generate input files for HAllA analysis
#   of diet–metabolite associations.
#
# Workflow:
#   1. Load previously processed metabolite data.
#   2. Align metabolite, diet, and metadata sample IDs.
#   3. Transpose metabolite data into HAllA-compatible format.
#   4. Process metabolite data using log2 transformation and scaling.
#   5. Generate subgroup-specific diet tables.
#   6. Provide HAllA commands for diet–metabolite association analysis.
#
# Note:
#   The functions used below (e.g., prepare_log2_scaled_metabo) are
#   custom functions defined elsewhere in the project and must be
#   available in the R environment before running this script.
###################################################################

load("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Metabolite_Genus/HAllA_Metabolite_Genus.RData")


# ==========================================================
# 1. Define Functional Resolution Level
# ==========================================================

# Specify the feature type to be processed in this analysis.
feature_level <- c("Metabolite")


# ==========================================================
# 2. Load and Align Input Data
# ==========================================================

# Identify samples present in both the metabolite table and diet dataset.
# This ensures that all downstream analyses use the same set and order
# of samples across metabolite, diet, and metadata tables.
sample_order <- intersect(rownames(met_05_CM), dataDiet_CM$Sample_ID)


# Reorder the diet data to match the metabolite sample order.
dataDiet_CM_filtered <- dataDiet_CM %>%
  slice(match(sample_order, Sample_ID)) 


# Reorder the metabolite data to match the same sample order.
met_05_CM_filtered <- met_05_CM %>%
  slice(match(sample_order, rownames(met_05_CM))) 


# Confirm that sample IDs are in identical order in both datasets.
identical(rownames(met_05_CM_filtered), dataDiet_CM_filtered$Sample_ID)



# Reconstruct Sample_ID in the metadata from the row names and
# reorder metadata to match the metabolite sample order.
meta_CM <- metadata_CM %>%
  mutate(Sample_ID = rownames(metadata_CM)) %>%
  slice(match(sample_order, Sample_ID)) 


# Confirm consistent sample IDs and ordering across all three datasets.
identical(rownames(met_05_CM_filtered), meta_CM$Sample_ID)
identical(dataDiet_CM_filtered$Sample_ID, meta_CM$Sample_ID)



# Transpose the metabolite table so that:
#   - rows = metabolites
#   - columns = samples
#
# This format is required for the subsequent HAllA workflow.
met_05_CM_t <- as.data.frame(t(met_05_CM_filtered))  


# Store the metabolite table in a named list to allow the same
# processing workflow to be applied iteratively to different
# functional feature levels.
tables <- list(Metabolite  = met_05_CM_t)




# ==========================================================
# 3. Define Analysis Subgroups
# ==========================================================

# Define the CM ethnic groups to be analyzed independently.
# The values must correspond exactly to the levels in the
# "Group" column of the metadata.
unique_group <- list(
  Akha_CM  = c("Akha-CM"),
  Lahu_CM  = c("Lahu-CM"),
  Khuen_CM = c("Khuen-CM")
)



# ==========================================================
# 4. Initialize Output Container
# ==========================================================

# Store processed metabolite tables for each subgroup.
# Results are named according to feature level and subgroup.
results_list <- list()



# ==========================================================
# 5. Process Metabolite Data by Subgroup
# ==========================================================

# Apply the predefined metabolite preprocessing function separately
# to each CM ethnic group.
#
# The preprocessing function performs the required log2 transformation
# and scaling of metabolite abundances before downstream association
# analysis.
for (data_type in feature_level) {
  for (comp_name in names(unique_group)) {
    group_pair <- unique_group[[comp_name]]
    
    res <- prepare_log2_scaled_metabo(
      data_table = tables[[data_type]],
      sample_metadata = meta_CM,
      data_type = data_type,
      group_col = "Group",
      group = group_pair,
      factor = "Eth"
    )
    
    results_list[[paste(data_type, comp_name, sep = "_")]] <- res
    message(paste("Completed:", data_type, comp_name))
  }
}



####------ Helper Fn -----#######
prepare_log2_scaled_metabo <- function(data_table, 
                                       sample_metadata, 
                                       data_type,
                                       group,
                                       group_col = "Group", 
                                       id_col = "Sample_ID",
                                       factor = "Geo") {
  
  library(dplyr)
  
  message(paste("Processing:", data_type, "|", group, "| transform: log2 + scale"))
  
  # Select samples for the group
  samples <- sample_metadata %>%
    filter(.data[[group_col]] %in% group) %>%
    pull(.data[[id_col]])
  
  # Subset table (columns = samples)
  feature <- data_table %>%
    dplyr::select(all_of(samples))%>%
    t() %>%
    as.data.frame()
  
  feature_matrix <- as.matrix(feature)
  
  if (!any(feature_matrix == 0, na.rm = TRUE)) {
    message("No zeros detected; no pseudocount applied.")
  }
  
  
  # Log2 transform
  log_matrix <- log2(feature_matrix)
  
  
  # Remove zero variance metabolite to avoid NA after scaling
  log_matrix <- log_matrix[, apply(log_matrix, 2, sd, na.rm = TRUE) > 0]
  
  
  # Scale per feature (col-wise standardization)
  # Each metabolite (column) is standardized across all samples (rows)
  scaled_matrix <- t(scale(log_matrix))
  
  transformed_tbl <- as.data.frame(scaled_matrix)
  
  # Output file
  file_name <- paste0("feature_", data_type, "_", group,  "_", factor, ".txt")
  
  write.table(transformed_tbl,
              file = file_name,
              sep = "\t",
              row.names = TRUE,
              col.names = NA,
              quote = FALSE)
  
  return(transformed_tbl)
}





# ==========================================================
# 6. Prepare Diet Data for HAllA
# ==========================================================

# Prepare diet variables for HAllA analysis against metabolite features.
# HAllA expects variables as rows and samples as columns.
diet_HAllA2 <- dataDiet_CM_filtered %>% 
  filter(Sample_ID %in% meta_CM$Sample_ID) 


# Prepare diet data separately for each CM ethnic group.
diet_Akha2 <- dataDiet_CM_filtered %>% 
  filter(Group == "Akha-CM") %>%
  select(-Sample_ID:-Age) %>%
  t() %>%
  as.data.frame()

diet_Khuen2 <- dataDiet_CM_filtered %>% 
  filter(Group == "Khuen-CM") %>%
  select(-Sample_ID:-Age) %>%
  t() %>%
  as.data.frame()

diet_Lahu2 <- dataDiet_CM_filtered %>% 
  filter(Group == "Lahu-CM") %>%
   select(-Sample_ID:-Age) %>%
  t() %>%
  as.data.frame()



# ==========================================================
# 7. Export HAllA Input Files
# ==========================================================

# Export subgroup-specific diet tables as tab-delimited files.
# These files can be supplied to HAllA as the X (diet) matrix.
write.table(
  diet_Akha2,
  "diet_Akha2.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

write.table(
  diet_Khuen2,
  "diet_Khuen2.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

write.table(
  diet_Lahu2,
  "diet_Lahu2.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)



# ==========================================================
# 8. HAllA Commands
# ==========================================================

# Run HAllA separately for each ethnic group.
#
#   -x: diet input matrix
#   -y: processed metabolite feature matrix
#   -o: output directory/name
#   -m: association method (Spearman correlation)
#   --fdr_alpha: false discovery rate threshold
#
# The metabolite feature files (feature_Metabolite_*.txt) are generated
# by the metabolite preprocessing workflow above.

# halla -x diet_Akha2.tsv -y feature_Metabolite_Akha-CM_Eth.txt -o Diet_AkhaCM2 -m spearman --fdr_alpha 0.05 
# halla -x diet_Khuen2.tsv -y feature_Metabolite_Khuen-CM_Eth.txt -o Diet_KhuenCM2 -m spearman --fdr_alpha 0.05 
# halla -x diet_Lahu2.tsv -y feature_Metabolite_Lahu-CM_Eth.txt -o Diet_LahuCM2 -m spearman --fdr_alpha 0.05








# ==========================================================
# PCA of Dietary Patterns
# ==========================================================

# Purpose:
#   Summarize overall variation in dietary patterns among participants
#   from the three ethnic groups in the Chiang Mai (CM) area.
#
#   PCA is performed on standardized dietary variables. Samples are
#   visualized according to ethnic group, with confidence ellipses
#   used to illustrate the distribution of each group in PCA space.


# Load required packages
library(dplyr)
library(FactoMineR)
library(factoextra)


# ==========================================================
# 1. Select Dietary Variables
# ==========================================================

# Select the dietary variables from Egg through Sausage for PCA.
# These variables represent the dietary features included in the
# analysis.
diet_pca_data <- dataDiet_CM %>%
  select(
    Egg:Sausage
  )


# ==========================================================
# 2. Perform Principal Component Analysis (PCA)
# ==========================================================

# Perform PCA using standardized dietary variables.
#
# scale.unit = TRUE:
#   Each dietary variable is centered and scaled before PCA so that
#   variables with different measurement scales contribute comparably.
#
# graph = FALSE:
#   Suppress the default PCA plot because the visualization is created
#   separately below using factoextra.
diet_pca <- PCA(
  diet_pca_data,
  scale.unit = TRUE,   # standardize variables
  graph = FALSE
)


# ==========================================================
# 3. Define Group Colors
# ==========================================================

# Define colors for the three ethnic groups in the Chiang Mai area.
# The names must match the levels of dataDiet_CM$Group.
group_colors <- c(
  "Akha-CM" = "gold1",
  "Lahu-CM" = "lightblue2",
  "Khuen-CM" = "darkorange"
)


# ==========================================================
# 4. Visualize Individual Samples in PCA Space
# ==========================================================

# Plot individual samples using the first two principal components.
#
# Samples are colored according to ethnic group.
# Confidence ellipses are added to visualize the distribution of
# samples within each group.
fviz_pca_ind(
  diet_pca,
  geom = "point",
  habillage = dataDiet_CM$Group,   # color by group
  addEllipses = TRUE,
  ellipse.type = "confidence",
  palette = group_colors,
  pointsize = 2,
  repel = TRUE
) +
  
  theme_bw() +
  
  labs(
    title = "PCA of dietary patterns across ethnic groups (CM)"
  ) +
  
  theme(
    plot.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    panel.grid = element_blank()
  )
