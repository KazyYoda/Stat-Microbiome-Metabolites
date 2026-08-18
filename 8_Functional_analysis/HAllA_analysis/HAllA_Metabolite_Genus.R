############################################################
# Microbiome Analysis in R: HAllA Metabolite+Genus
############################################################


#--------------------------------------------
# Data Preparation for HAllA: Genus+KO
#--------------------------------------------

setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Metabolite_Genus")
load("~/Documents/HillTribe_NGS/7.RDA/RDA.RData")
load("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Metabolite_Genus/HAllA_Metabolite_Genus.RData")

# Load required packages
library(dplyr)
library(tidyr)
library(readxl)
library(car)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(ggh4x)
library(gridExtra)



#----------- Helper Function with Conditional Pseudocount ---------------
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



############################################################
# Akha Ethnic Group
# Genus-Level Feature Processing
# Generate CLR-transformed genus count by group
############################################################

# ==========================================================
# 1. Load Genus Abundance Data
# ==========================================================

# Metadata 
drop_Akha <- sample_metadata_drop %>%
  filter(Ethnicity == "Akha") %>%
  mutate(Group = paste(Ethnicity, Area, sep = "-")) %>%
  relocate(Group, .before = 2) %>%
  as.data.frame()

# Data preparation
Counts_Genus_Akha <- Counts_Genus %>%
  filter(Sample_ID %in% drop_Akha$Sample_ID)

# Extract numeric count columns (remove metadata)
Genus_abs_Akha <- Counts_Genus_Akha %>%
  as.data.frame()

# Assign Sample IDs as row names
rownames(Genus_abs_Akha) <- Counts_Genus_Akha$Sample_ID



# ==========================================================
# 2. Prepare metabolite data
# ==========================================================
met_05_Akha <- metabo_05_t[drop_Akha$Sample_ID, ]



# Confirm sample order matches metabolite table structure
identical(rownames(Genus_abs_Akha), drop_Akha$Sample_ID)
identical(rownames(Genus_abs_Akha), rownames(met_05_Akha))


# ==========================================================
# 3. Format for CLR Processing
# ==========================================================

# Transpose table: samples → columns, taxa → rows
Genus_abs_Akha_t <- as.data.frame(t(Genus_abs_Akha[-c(1:11)]))  



# ==========================================================
# 4. Define Analysis Inputs
# ==========================================================

# Taxonomic resolution level
rank_level <- c("Genus")

# Store table in iterable structure
tables <- list(Genus = Genus_abs_t)


# ==========================================================
# 5. Define Target Subgroups
# ==========================================================

# Each subgroup processed independently
unique_group <- list(
  Akha_CM  = c("Akha-CM"),
  Akha_CR  = c("Akha-CR")
)



# ==========================================================
# 6. Initialize Output Container
# ==========================================================

results_list <- list()



# ==========================================================
# 7. CLR Processing Loop (Genus Level)
# ==========================================================

for (data_type in rank_level) {
  for (comp_name in names(unique_group)) {
    
    # Extract subgroup label
    group_pair <- unique_group[[comp_name]]
    
    # Generate CLR-transformed genus table
    res <- prepare_clr_table(
      feature_table = tables[[rank_level]],
      sample_metadata = drop_Akha,
      data_type = rank_level,
      group_col = "Group",
      group = group_pair
    )
    
    # Store result using structured naming
    results_list[[paste(data_type, comp_name, sep = "_")]] <- res
    
    # Progress indicator
    message(paste("Completed:", data_type, comp_name))
  }
}




###################################################################
# Akha Ethnic Group
# Functional Profile Processing Loop (Metabolites)
# Generate CLR-transformed abundance tables by subgroup
##################################################################


# ==========================================================
# 1. Define Functional Resolution Levels
# ==========================================================

# Functional hierarchy levels to process
feature_level <- c("Metabolite")



# ==========================================================
# 2. Input Data Sources
# ==========================================================

# Confirm sample order matches metabolite metadata
identical(rownames(met_05_Akha), drop_Akha$Sample_ID)

# Transpose table: samples → columns, metabolite → rows
met_05_Akha_t <- as.data.frame(t(met_05_Akha))  

# Store count tables in a named list for iteration
tables <- list(Metabolite  = met_05_Akha_t)




# ==========================================================
# 3. Define Analysis Subgroups
# ==========================================================

# Each subgroup is processed independently
# Values must match the "Group" column in metadata
unique_group <- list(
  Akha_CM  = c("Akha-CM"),
  Akha_CR  = c("Akha-CR")
)



# ==========================================================
# 4. Initialize Output Container
# ==========================================================

#----------- Helper Function Log2 Transformation and scaling for Metabolites ---------------
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



# Stores processed CLR tables with structured naming
results_list <- list()



# ==========================================================
# 5. Processing Loop
# ==========================================================

for (data_type in feature_level) {
  for (comp_name in names(unique_group)) {
    group_pair <- unique_group[[comp_name]]
    
    res <- prepare_log2_scaled_metabo(
      data_table = tables[[data_type]],
      sample_metadata = drop_Akha,
      data_type = data_type,
      group_col = "Group",
      group = group_pair
    )
    
    results_list[[paste(data_type, comp_name, sep = "_")]] <- res
    message(paste("Completed:", data_type, comp_name))
  }
}











############################################################
# Lahu Ethnic Group
# Genus-Level Feature Processing
# Generate CLR-transformed genus count by group
############################################################

# ==========================================================
# 1. Load Genus Abundance Data
# ==========================================================

# Metadata 
drop_Lahu <- sample_metadata_drop %>%
  filter(Ethnicity == "Lahu") %>%
  mutate(Group = paste(Ethnicity, Area, sep = "-")) %>%
  relocate(Group, .before = 2) %>%
  as.data.frame()

# Data preparation
Counts_Genus_Lahu <- Counts_Genus %>%
  filter(Sample_ID %in% drop_Lahu$Sample_ID)

# Extract numeric count columns (remove metadata)
Genus_abs_Lahu <- Counts_Genus_Lahu %>%
  as.data.frame()

# Assign Sample IDs as row names
rownames(Genus_abs_Lahu) <- Counts_Genus_Lahu$Sample_ID



# ==========================================================
# 2. Prepare metabolite data
# ==========================================================
met_05_Lahu <- metabo_05_t[drop_Lahu$Sample_ID, ]



# Confirm sample order matches metabolite table structure
identical(rownames(Genus_abs_Lahu), drop_Lahu$Sample_ID)
identical(rownames(Genus_abs_Lahu), rownames(met_05_Lahu))


# ==========================================================
# 3. Format for CLR Processing
# ==========================================================

# Transpose table: samples → columns, taxa → rows
Genus_abs_Lahu_t <- as.data.frame(t(Genus_abs_Lahu[-c(1:11)]))  



# ==========================================================
# 4. Define Analysis Inputs
# ==========================================================

# Taxonomic resolution level
rank_level <- c("Genus")

# Store table in iterable structure
tables <- list(Genus = Genus_abs_t)


# ==========================================================
# 5. Define Target Subgroups
# ==========================================================

# Each subgroup processed independently
unique_group <- list(
  Lahu_CM  = c("Lahu-CM"),
  Lahu_CR  = c("Lahu-CR")
)



# ==========================================================
# 6. Initialize Output Container
# ==========================================================

results_list <- list()



# ==========================================================
# 7. CLR Processing Loop (Genus Level)
# ==========================================================

for (data_type in rank_level) {
  for (comp_name in names(unique_group)) {
    
    # Extract subgroup label
    group_pair <- unique_group[[comp_name]]
    
    # Generate CLR-transformed genus table
    res <- prepare_clr_table(
      feature_table = tables[[rank_level]],
      sample_metadata = drop_Lahu,
      data_type = rank_level,
      group_col = "Group",
      group = group_pair
    )
    
    # Store result using structured naming
    results_list[[paste(data_type, comp_name, sep = "_")]] <- res
    
    # Progress indicator
    message(paste("Completed:", data_type, comp_name))
  }
}




###################################################################
# Lahu Ethnic Group
# Functional Profile Processing Loop (Metabolites)
# Generate CLR-transformed abundance tables by subgroup
##################################################################


# ==========================================================
# 1. Define Functional Resolution Levels
# ==========================================================

# Functional hierarchy levels to process
feature_level <- c("Metabolite")



# ==========================================================
# 2. Input Data Sources
# ==========================================================

# Confirm sample order matches metabolite metadata
identical(rownames(met_05_Lahu), drop_Lahu$Sample_ID)

# Transpose table: samples → columns, metabolite → rows
met_05_Lahu_t <- as.data.frame(t(met_05_Lahu))  

# Store count tables in a named list for iteration
tables <- list(Metabolite  = met_05_Lahu_t)




# ==========================================================
# 3. Define Analysis Subgroups
# ==========================================================

# Each subgroup is processed independently
# Values must match the "Group" column in metadata
unique_group <- list(
  Lahu_CM  = c("Lahu-CM"),
  Lahu_CR  = c("Lahu-CR")
)



# ==========================================================
# 4. Initialize Output Container
# ==========================================================

# Stores processed CLR tables with structured naming
results_list <- list()



# ==========================================================
# 5. Processing Loop
# ==========================================================

for (data_type in feature_level) {
  for (comp_name in names(unique_group)) {
    group_pair <- unique_group[[comp_name]]
    
    res <- prepare_log2_scaled_metabo(
      data_table = tables[[data_type]],
      sample_metadata = drop_Lahu,
      data_type = data_type,
      group_col = "Group",
      group = group_pair
    )
    
    results_list[[paste(data_type, comp_name, sep = "_")]] <- res
    message(paste("Completed:", data_type, comp_name))
  }
}














############################################################
# CM location
# Genus-Level Feature Processing
# Generate CLR-transformed genus count by group
############################################################

# ==========================================================
# 1. Load Genus Abundance Data
# ==========================================================

# Metadata 
drop_CM <- sample_metadata_drop %>%
  filter(Area == "CM") %>%
  mutate(Group = paste(Ethnicity, Area, sep = "-")) %>%
  relocate(Group, .before = 2) %>%
  as.data.frame()

# Data preparation
Counts_Genus_CM <- Counts_Genus %>%
  filter(Sample_ID %in% drop_CM$Sample_ID)

# Extract numeric count columns (remove metadata)
Genus_abs_CM <- Counts_Genus_CM %>%
  as.data.frame()

# Assign Sample IDs as row names
rownames(Genus_abs_CM) <- Counts_Genus_CM$Sample_ID



# ==========================================================
# 2. Prepare metabolite data
# ==========================================================
met_05_CM <- metabo_05_t[drop_CM$Sample_ID, ]



# Confirm sample order matches metabolite table structure
identical(rownames(Genus_abs_CM), drop_CM$Sample_ID)
identical(rownames(Genus_abs_CM), rownames(met_05_CM))


# ==========================================================
# 3. Format for CLR Processing
# ==========================================================

# Transpose table: samples → columns, taxa → rows
Genus_abs_CM_t <- as.data.frame(t(Genus_abs_CM[-c(1:11)]))  



# ==========================================================
# 4. Define Analysis Inputs
# ==========================================================

# Taxonomic resolution level
rank_level <- c("Genus")

# Store table in iterable structure
tables <- list(Genus = Genus_abs_t)


# ==========================================================
# 5. Define Target Subgroups
# ==========================================================

# Each subgroup processed independently
unique_group <- list(
  Akha_CM  = c("Akha-CM"),
  Lahu_CM  = c("Lahu-CM"),
  Khuen_CM = c("Khuen-CM")
)



# ==========================================================
# 6. Initialize Output Container
# ==========================================================

results_list <- list()



# ==========================================================
# 7. CLR Processing Loop (Genus Level)
# ==========================================================

for (data_type in rank_level) {
  for (comp_name in names(unique_group)) {
    
    # Extract subgroup label
    group_pair <- unique_group[[comp_name]]
    
    # Generate CLR-transformed genus table
    res <- prepare_clr_table(
      feature_table = tables[[rank_level]],
      sample_metadata = drop_CM,
      data_type = rank_level,
      group_col = "Group",
      group = group_pair,
      factor = "Eth"
    )
    
    # Store result using structured naming
    results_list[[paste(data_type, comp_name, sep = "_")]] <- res
    
    # Progress indicator
    message(paste("Completed:", data_type, comp_name))
  }
}




###################################################################
# CM Ethnic Group
# Functional Profile Processing Loop (Metabolites)
# Generate CLR-transformed abundance tables by subgroup
##################################################################


# ==========================================================
# 1. Define Functional Resolution Levels
# ==========================================================

# Functional hierarchy levels to process
feature_level <- c("Metabolite")



# ==========================================================
# 2. Input Data Sources
# ==========================================================

# Confirm sample order matches metabolite metadata
identical(rownames(met_05_CM), drop_CM$Sample_ID)

# Transpose table: samples → columns, metabolite → rows
met_05_CM_t <- as.data.frame(t(met_05_CM))  

# Store count tables in a named list for iteration
tables <- list(Metabolite  = met_05_CM_t)




# ==========================================================
# 3. Define Analysis Subgroups
# ==========================================================

# Each subgroup is processed independently
# Values must match the "Group" column in metadata
unique_group <- list(
  Akha_CM  = c("Akha-CM"),
  Lahu_CM  = c("Lahu-CM"),
  Khuen_CM = c("Khuen-CM")
)



# ==========================================================
# 4. Initialize Output Container
# ==========================================================

# Stores processed CLR tables with structured naming
results_list <- list()



# ==========================================================
# 5. Processing Loop
# ==========================================================

for (data_type in feature_level) {
  for (comp_name in names(unique_group)) {
    group_pair <- unique_group[[comp_name]]
    
    res <- prepare_log2_scaled_metabo(
      data_table = tables[[data_type]],
      sample_metadata = drop_CM,
      data_type = data_type,
      group_col = "Group",
      group = group_pair,
      factor = "Eth"
    )
    
    results_list[[paste(data_type, comp_name, sep = "_")]] <- res
    message(paste("Completed:", data_type, comp_name))
  }
}













############################################################
# CR location
# Genus-Level Feature Processing
# Generate CLR-transformed genus count by group
############################################################

# ==========================================================
# 1. Load Genus Abundance Data
# ==========================================================

# Metadata 
drop_CR <- sample_metadata_drop %>%
  filter(Area == "CR") %>%
  mutate(Group = paste(Ethnicity, Area, sep = "-")) %>%
  relocate(Group, .before = 2) %>%
  as.data.frame()

# Data preparation
Counts_Genus_CR <- Counts_Genus %>%
  filter(Sample_ID %in% drop_CR$Sample_ID)

# Extract numeric count columns (remove metadata)
Genus_abs_CR <- Counts_Genus_CR %>%
  as.data.frame()

# Assign Sample IDs as row names
rownames(Genus_abs_CR) <- Counts_Genus_CR$Sample_ID



# ==========================================================
# 2. Prepare metabolite data
# ==========================================================
met_05_CR <- metabo_05_t[drop_CR$Sample_ID, ]



# Confirm sample order matches metabolite table structure
identical(rownames(Genus_abs_CR), drop_CR$Sample_ID)
identical(rownames(Genus_abs_CR), rownames(met_05_CR))


# ==========================================================
# 3. Format for CLR Processing
# ==========================================================

# Transpose table: samples → columns, taxa → rows
Genus_abs_CR_t <- as.data.frame(t(Genus_abs_CR[-c(1:11)]))  



# ==========================================================
# 4. Define Analysis Inputs
# ==========================================================

# Taxonomic resolution level
rank_level <- c("Genus")

# Store table in iterable structure
tables <- list(Genus = Genus_abs_t)


# ==========================================================
# 5. Define Target Subgroups
# ==========================================================

# Each subgroup processed independently
unique_group <- list(
  Akha_CR  = c("Akha-CR"),
  Lahu_CR  = c("Lahu-CR"),
  Khuen_CR = c("Lisu-CR")
)



# ==========================================================
# 6. Initialize Output Container
# ==========================================================

results_list <- list()



# ==========================================================
# 7. CLR Processing Loop (Genus Level)
# ==========================================================

for (data_type in rank_level) {
  for (comp_name in names(unique_group)) {
    
    # Extract subgroup label
    group_pair <- unique_group[[comp_name]]
    
    # Generate CLR-transformed genus table
    res <- prepare_clr_table(
      feature_table = tables[[rank_level]],
      sample_metadata = drop_CR,
      data_type = rank_level,
      group_col = "Group",
      group = group_pair,
      factor = "Eth"
    )
    
    # Store result using structured naming
    results_list[[paste(data_type, comp_name, sep = "_")]] <- res
    
    # Progress indicator
    message(paste("Completed:", data_type, comp_name))
  }
}




###################################################################
# CR Ethnic Group
# Functional Profile Processing Loop (Metabolites)
# Generate CLR-transformed abundance tables by subgroup
##################################################################


# ==========================================================
# 1. Define Functional Resolution Levels
# ==========================================================

# Functional hierarchy levels to process
feature_level <- c("Metabolite")



# ==========================================================
# 2. Input Data Sources
# ==========================================================

# Confirm sample order matches metabolite metadata
identical(rownames(met_05_CR), drop_CR$Sample_ID)

# Transpose table: samples → columns, metabolite → rows
met_05_CR_t <- as.data.frame(t(met_05_CR))  

# Store count tables in a named list for iteration
tables <- list(Metabolite  = met_05_CR_t)




# ==========================================================
# 3. Define Analysis Subgroups
# ==========================================================

# Each subgroup is processed independently
# Values must match the "Group" column in metadata
unique_group <- list(
  Akha_CR  = c("Akha-CR"),
  Lahu_CR  = c("Lahu-CR"),
  Lisu_CR = c("Lisu-CR")
)



# ==========================================================
# 4. Initialize Output Container
# ==========================================================

# Stores processed CLR tables with structured naming
results_list <- list()



# ==========================================================
# 5. Processing Loop
# ==========================================================

for (data_type in feature_level) {
  for (comp_name in names(unique_group)) {
    group_pair <- unique_group[[comp_name]]
    
    res <- prepare_log2_scaled_metabo(
      data_table = tables[[data_type]],
      sample_metadata = drop_CR,
      data_type = data_type,
      group_col = "Group",
      group = group_pair,
      factor = "Eth"
    )
    
    results_list[[paste(data_type, comp_name, sep = "_")]] <- res
    message(paste("Completed:", data_type, comp_name))
  }
}


