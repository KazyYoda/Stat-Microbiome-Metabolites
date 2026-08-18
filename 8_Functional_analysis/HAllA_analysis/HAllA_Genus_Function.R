############################################################
# Microbiome Analysis in R: HAllA Genus+Functional profiles
############################################################


#--------------------------------------------
# Data Preparation for HAllA: Genus+Funtion
#--------------------------------------------

setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function")
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





###################################################################
# Functional Profile Processing Loop (PW / EC / KO)
# Generate CLR-transformed abundance tables by subgroup
##################################################################


############################################################
# Akha Ethnic Group
# Function Feature Processing
# Generate CLR-transformed Predicted function by group
############################################################
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo")


# ==========================================================
# 1. Data preparation 
# ==========================================================

metadata_Akha <- metadata %>%
  filter(Ethnicity == "Akha") %>%
  mutate(Group = paste(Ethnicity, Area, sep = "-")) %>%
  relocate(Group, .before = 2) %>%
  as.data.frame()

# PW
PW_Akha <- PW_abs %>%
  select(all_of(metadata_Akha$Sample_ID)) %>%
  as.data.frame()

# EC
EC_Akha <- EC_abs %>%
  select(all_of(metadata_Akha$Sample_ID)) %>%
  as.data.frame()

# KO
KO_Akha <- KO_abs %>%
  select(all_of(metadata_Akha$Sample_ID)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample_ID = colnames(EC_Akha)) %>%
  relocate(Sample_ID, .before = 1) %>%
  left_join(metadata_Akha[, c(1:2)], by = "Sample_ID") %>%
  as.data.frame()

# Filter KOs with prevalence >= XXX%
n_samples <- nrow(KO_Akha)

prev_features <- KO_Akha %>%
  pivot_longer(
    cols = -c(Sample_ID, Group),
    names_to = "Feature",
    values_to = "Abundance"
  ) %>%
  mutate(Present = Abundance > 0) %>%
  group_by(Feature) %>%
  summarise(Prevalence = sum(Present, na.rm = TRUE), .groups = "drop") %>%
  filter(Prevalence >= 1 * n_samples) %>%
  pull(Feature)

KO_Akha_filtered <- KO_Akha[, prev_features]
rownames(KO_Akha_filtered) <- KO_Akha$Sample_ID
KO_Akha_filtered <- KO_Akha_filtered %>% t() %>% as.data.frame()


# ==========================================================
# 2. Checking metadata
# ==========================================================

# Confirm sample order matches metabolite table structure
identical(colnames(PW_Akha), metadata_Akha$Sample_ID)
identical(colnames(EC_Akha), metadata_Akha$Sample_ID)
identical(colnames(KO_Akha_filtered), metadata_Akha$Sample_ID)



# ==========================================================
# 3. Define Analysis Inputs
# ==========================================================

# Functional hierarchy levels to process
functional_level <- c("PW", "EC", "KO")

# Store table in iterable structure
tables <- list( 
  PW = PW_Akha,
  EC = EC_Akha,
  KO = KO_Akha_filtered
)


# ==========================================================
# 4. Define Target Subgroups
# ==========================================================

# Each subgroup processed independently
unique_group <- list(
  Akha_CM  = c("Akha-CM"),
  Akha_CR  = c("Akha-CR")
)



# ==========================================================
# 5. Initialize Output Container
# ==========================================================

results_list <- list()



# ==========================================================
# 6. CLR Processing Loop (functional_level)
# ==========================================================

# Outer loop: functional resolution level
for (data_type in functional_level) {
  
  # Inner loop: population subgroup
  for (comp_name in names(unique_group)) {
    
    # Extract subgroup label
    group_pair <- unique_group[[comp_name]]
    
    # Generate CLR-transformed table
    res <- prepare_clr_table(
      feature_table = tables[[data_type]],
      sample_metadata = metadata_Akha,
      data_type = data_type,
      group_col = "Group",
      group = group_pair,
      factor = "Geo"
    )
    
    # Store result using structured naming
    results_list[[paste(data_type, comp_name, sep = "_")]] <- res
    
    # Progress indicator
    message(paste("Completed:", data_type, comp_name))
  }
}




############################################################
# Genus-Level Feature Processing
# Generate CLR-transformed genus abundance by subgroup
############################################################

# ==========================================================
# 1. Load Genus Abundance Data
# ==========================================================

# Data preparation
Counts_Genus_Akha <- Counts_Genus %>%
  filter(Sample_ID %in% metadata_Akha$Sample_ID) %>%
  as.data.frame()

# Assign Sample IDs as row names
rownames(Counts_Genus_Akha) <- Counts_Genus_Akha$Sample_ID



# ==========================================================
# 2. Verify Sample Order Consistency
# ==========================================================

# Confirm sample order matches functional table structure
identical(Counts_Genus_Akha$Sample_ID, colnames(PW_Akha))
identical(Counts_Genus_Akha$Sample_ID, colnames(EC_Akha))
identical(Counts_Genus_Akha$Sample_ID, colnames(KO_Akha))



# ==========================================================
# 3. Format for CLR Processing
# ==========================================================

# Transpose table: samples → columns, taxa → rows
Counts_Genus_Akha_t <- as.data.frame(t(Counts_Genus_Akha[-c(1:11)]))  



# ==========================================================
# 4. Define Analysis Inputs
# ==========================================================

# Taxonomic resolution level
rank_level <- c("Genus")

# Store table in iterable structure
tables <- list(Genus = Counts_Genus_Akha_t)



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

# Outer loop: functional resolution level
for (data_type in rank_level) {
  
  # Inner loop: population subgroup
  for (comp_name in names(unique_group)) {
    
    # Extract subgroup label
    group_pair <- unique_group[[comp_name]]
    
    # Generate CLR-transformed table
    res <- prepare_clr_table(
      feature_table = tables[[data_type]],
      sample_metadata = metadata_Akha,
      data_type = data_type,
      group_col = "Group",
      group = group_pair,
      factor = "Geo"
    )
    
    # Store result using structured naming
    results_list[[paste(data_type, comp_name, sep = "_")]] <- res
    
    # Progress indicator
    message(paste("Completed:", data_type, comp_name))
  }
}












############################################################
# Lahu Ethnic Group
# Function Feature Processing
# Generate CLR-transformed Predicted function by group
############################################################
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo")


# ==========================================================
# 1. Data preparation 
# ==========================================================

metadata_Lahu <- metadata %>%
  filter(Ethnicity == "Lahu") %>%
  mutate(Group = paste(Ethnicity, Area, sep = "-")) %>%
  relocate(Group, .before = 2) %>%
  as.data.frame()

# PW
PW_Lahu <- PW_abs %>%
  select(all_of(metadata_Lahu$Sample_ID)) %>%
  as.data.frame()

# EC
EC_Lahu <- EC_abs %>%
  select(all_of(metadata_Lahu$Sample_ID)) %>%
  as.data.frame()

# KO
KO_Lahu <- KO_abs %>%
  select(all_of(metadata_Lahu$Sample_ID)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample_ID = colnames(EC_Lahu)) %>%
  relocate(Sample_ID, .before = 1) %>%
  left_join(metadata_Lahu[, c(1:2)], by = "Sample_ID") %>%
  as.data.frame()

# Filter Genus|KOs with prevalence >= 20%
n_samples <- nrow(KO_Lahu)

prev_features <- KO_Lahu %>%
  pivot_longer(
    cols = -c(Sample_ID, Group),
    names_to = "Feature",
    values_to = "Abundance"
  ) %>%
  mutate(Present = Abundance > 0) %>%
  group_by(Feature) %>%
  summarise(Prevalence = sum(Present, na.rm = TRUE), .groups = "drop") %>%
  filter(Prevalence >= 1 * n_samples) %>%
  pull(Feature)

KO_Lahu_filtered <- KO_Lahu[, prev_features]
rownames(KO_Lahu_filtered) <- KO_Lahu$Sample_ID
KO_Lahu_filtered <- KO_Lahu_filtered %>% t() %>% as.data.frame()



# ==========================================================
# 2. Checking metadata
# ==========================================================

# Confirm sample order matches metabolite table structure
identical(colnames(PW_Lahu), metadata_Lahu$Sample_ID)
identical(colnames(EC_Lahu), metadata_Lahu$Sample_ID)
identical(colnames(KO_Lahu_filtered), metadata_Lahu$Sample_ID)



# ==========================================================
# 3. Define Analysis Inputs
# ==========================================================

# Functional hierarchy levels to process
functional_level <- c("PW", "EC", "KO")

# Store table in iterable structure
tables <- list( 
  PW = PW_Lahu,
  EC = EC_Lahu,
  KO = KO_Lahu_filtered
)


# ==========================================================
# 4. Define Target Subgroups
# ==========================================================

# Each subgroup processed independently
unique_group <- list(
  Lahu_CM  = c("Lahu-CM"),
  Lahu_CR  = c("Lahu-CR")
)



# ==========================================================
# 5. Initialize Output Container
# ==========================================================

results_list <- list()



# ==========================================================
# 6. CLR Processing Loop (functional_level)
# ==========================================================

# Outer loop: functional resolution level
for (data_type in functional_level) {
  
  # Inner loop: population subgroup
  for (comp_name in names(unique_group)) {
    
    # Extract subgroup label
    group_pair <- unique_group[[comp_name]]
    
    # Generate CLR-transformed table
    res <- prepare_clr_table(
      feature_table = tables[[data_type]],
      sample_metadata = metadata_Lahu,
      data_type = data_type,
      group_col = "Group",
      group = group_pair,
      factor = "Geo"
    )
    
    # Store result using structured naming
    results_list[[paste(data_type, comp_name, sep = "_")]] <- res
    
    # Progress indicator
    message(paste("Completed:", data_type, comp_name))
  }
}




############################################################
# Genus-Level Feature Processing
# Generate CLR-transformed genus abundance by subgroup
############################################################

# ==========================================================
# 1. Load Genus Abundance Data
# ==========================================================

# Data preparation
Counts_Genus_Lahu <- Counts_Genus %>%
  filter(Sample_ID %in% metadata_Lahu$Sample_ID) %>%
  as.data.frame()

# Assign Sample IDs as row names
rownames(Counts_Genus_Lahu) <- Counts_Genus_Lahu$Sample_ID



# ==========================================================
# 2. Verify Sample Order Consistency
# ==========================================================

# Confirm sample order matches functional table structure
identical(Counts_Genus_Lahu$Sample_ID, colnames(PW_Lahu))
identical(Counts_Genus_Lahu$Sample_ID, colnames(EC_Lahu))
identical(Counts_Genus_Lahu$Sample_ID, colnames(KO_Lahu))



# ==========================================================
# 3. Format for CLR Processing
# ==========================================================

# Transpose table: samples → columns, taxa → rows
Counts_Genus_Lahu_t <- as.data.frame(t(Counts_Genus_Lahu[-c(1:11)]))  



# ==========================================================
# 4. Define Analysis Inputs
# ==========================================================

# Taxonomic resolution level
rank_level <- c("Genus")

# Store table in iterable structure
tables <- list(Genus = Counts_Genus_Lahu_t)



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

# Outer loop: functional resolution level
for (data_type in rank_level) {
  
  # Inner loop: population subgroup
  for (comp_name in names(unique_group)) {
    
    # Extract subgroup label
    group_pair <- unique_group[[comp_name]]
    
    # Generate CLR-transformed table
    res <- prepare_clr_table(
      feature_table = tables[[data_type]],
      sample_metadata = metadata_Lahu,
      data_type = data_type,
      group_col = "Group",
      group = group_pair,
      factor = "Geo"
    )
    
    # Store result using structured naming
    results_list[[paste(data_type, comp_name, sep = "_")]] <- res
    
    # Progress indicator
    message(paste("Completed:", data_type, comp_name))
  }
}











############################################################
# CM Ethnic Group
# Function Feature Processing
# Generate CLR-transformed Predicted function by group
############################################################
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth")


# ==========================================================
# 1. Data preparation 
# ==========================================================

metadata_CM <- metadata %>%
  filter(Area == "CM") %>%
  mutate(Group = paste(Ethnicity, Area, sep = "-")) %>%
  relocate(Group, .before = 2) %>%
  as.data.frame()

# PW
PW_CM <- PW_abs %>%
  select(all_of(metadata_CM$Sample_ID)) %>%
  as.data.frame()

# EC
EC_CM <- EC_abs %>%
  select(all_of(metadata_CM$Sample_ID)) %>%
  as.data.frame()

# KO
KO_CM <- KO_abs %>%
  select(all_of(metadata_CM$Sample_ID)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample_ID = colnames(EC_CM)) %>%
  relocate(Sample_ID, .before = 1) %>%
  left_join(metadata_CM[, c(1:2)], by = "Sample_ID") %>%
  as.data.frame()

# Filter Genus|KOs with prevalence >= 20%
n_samples <- nrow(KO_CM)

prev_features <- KO_CM %>%
  pivot_longer(
    cols = -c(Sample_ID, Group),
    names_to = "Feature",
    values_to = "Abundance"
  ) %>%
  mutate(Present = Abundance > 0) %>%
  group_by(Feature) %>%
  summarise(Prevalence = sum(Present, na.rm = TRUE), .groups = "drop") %>%
  filter(Prevalence >= 1 * n_samples) %>%
  pull(Feature)

KO_CM_filtered <- KO_CM[, prev_features]
rownames(KO_CM_filtered) <- KO_CM$Sample_ID
KO_CM_filtered <- KO_CM_filtered %>% t() %>% as.data.frame()

# ==========================================================
# 2. Checking metadata
# ==========================================================

# Confirm sample order matches metabolite table structure
identical(colnames(PW_CM), metadata_CM$Sample_ID)
identical(colnames(EC_CM), metadata_CM$Sample_ID)
identical(colnames(KO_CM_filtered), metadata_CM$Sample_ID)



# ==========================================================
# 3. Define Analysis Inputs
# ==========================================================

# Functional hierarchy levels to process
functional_level <- c("PW", "EC", "KO")

# Store table in iterable structure
tables <- list( 
  PW = PW_CM,
  EC = EC_CM,
  KO = KO_CM_filtered
)


# ==========================================================
# 4. Define Target Subgroups
# ==========================================================

# Each subgroup processed independently
unique_group <- list(
  Akha_CM  = c("Akha-CM"),
  Lahu_CM  = c("Lahu-CM"),
  Khuen_CM = c("Khuen-CM")
)




# ==========================================================
# 5. Initialize Output Container
# ==========================================================

results_list <- list()



# ==========================================================
# 6. CLR Processing Loop (functional_level)
# ==========================================================

# Outer loop: functional resolution level
for (data_type in functional_level) {
  
  # Inner loop: population subgroup
  for (comp_name in names(unique_group)) {
    
    # Extract subgroup label
    group_pair <- unique_group[[comp_name]]
    
    # Generate CLR-transformed table
    res <- prepare_clr_table(
      feature_table = tables[[data_type]],
      sample_metadata = metadata_CM,
      data_type = data_type,
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




############################################################
# Genus-Level Feature Processing
# Generate CLR-transformed genus abundance by subgroup
############################################################

# ==========================================================
# 1. Load Genus Abundance Data
# ==========================================================

# Data preparation
Counts_Genus_CM <- Counts_Genus %>%
  filter(Sample_ID %in% metadata_CM$Sample_ID) %>%
  as.data.frame()

# Assign Sample IDs as row names
rownames(Counts_Genus_CM) <- Counts_Genus_CM$Sample_ID



# ==========================================================
# 2. Verify Sample Order Consistency
# ==========================================================

# Confirm sample order matches functional table structure
identical(Counts_Genus_CM$Sample_ID, colnames(PW_CM))
identical(Counts_Genus_CM$Sample_ID, colnames(EC_CM))
identical(Counts_Genus_CM$Sample_ID, colnames(KO_CM_filtered))



# ==========================================================
# 3. Format for CLR Processing
# ==========================================================

# Transpose table: samples → columns, taxa → rows
Counts_Genus_CM_t <- as.data.frame(t(Counts_Genus_CM[-c(1:11)]))  



# ==========================================================
# 4. Define Analysis Inputs
# ==========================================================

# Taxonomic resolution level
rank_level <- c("Genus")

# Store table in iterable structure
tables <- list(Genus = Counts_Genus_CM_t)



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

# Outer loop: functional resolution level
for (data_type in rank_level) {
  
  # Inner loop: population subgroup
  for (comp_name in names(unique_group)) {
    
    # Extract subgroup label
    group_pair <- unique_group[[comp_name]]
    
    # Generate CLR-transformed table
    res <- prepare_clr_table(
      feature_table = tables[[data_type]],
      sample_metadata = metadata_CM,
      data_type = data_type,
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













############################################################
# CR Ethnic Group
# Function Feature Processing
# Generate CLR-transformed Predicted function by group
############################################################
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth")


# ==========================================================
# 1. Data preparation 
# ==========================================================

metadata_CR <- metadata %>%
  filter(Area == "CR") %>%
  mutate(Group = paste(Ethnicity, Area, sep = "-")) %>%
  relocate(Group, .before = 2) %>%
  as.data.frame()

# PW
PW_CR <- PW_abs %>%
  select(all_of(metadata_CR$Sample_ID)) %>%
  as.data.frame()

# EC
EC_CR <- EC_abs %>%
  select(all_of(metadata_CR$Sample_ID)) %>%
  as.data.frame()

# KO
KO_CR <- KO_abs %>%
  select(all_of(metadata_CR$Sample_ID)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample_ID = colnames(EC_CR)) %>%
  relocate(Sample_ID, .before = 1) %>%
  left_join(metadata_CR[, c(1:2)], by = "Sample_ID") %>%
  as.data.frame()

# Filter Genus|KOs with prevalence >= 20%
n_samples <- nrow(KO_CR)

prev_features <- KO_CR %>%
  pivot_longer(
    cols = -c(Sample_ID, Group),
    names_to = "Feature",
    values_to = "Abundance"
  ) %>%
  mutate(Present = Abundance > 0) %>%
  group_by(Feature) %>%
  summarise(Prevalence = sum(Present, na.rm = TRUE), .groups = "drop") %>%
  filter(Prevalence >= 1 * n_samples) %>%
  pull(Feature)

KO_CR_filtered <- KO_CR[, prev_features]
rownames(KO_CR_filtered) <- KO_CR$Sample_ID
KO_CR_filtered <- KO_CR_filtered %>% t() %>% as.data.frame()



# ==========================================================
# 2. Checking metadata
# ==========================================================

# Confirm sample order matches metabolite table structure
identical(colnames(PW_CR), metadata_CR$Sample_ID)
identical(colnames(EC_CR), metadata_CR$Sample_ID)
identical(colnames(KO_CR_filtered), metadata_CR$Sample_ID)



# ==========================================================
# 3. Define Analysis Inputs
# ==========================================================

# Functional hierarchy levels to process
functional_level <- c("PW", "EC", "KO")

# Store table in iterable structure
tables <- list( 
  PW = PW_CR,
  EC = EC_CR,
  KO = KO_CR_filtered
)


# ==========================================================
# 4. Define Target Subgroups
# ==========================================================

# Each subgroup processed independently
unique_group <- list(
  Akha_CR  = c("Akha-CR"),
  Lahu_CR  = c("Lahu-CR"),
  Lisu_CR = c("Lisu-CR")
)




# ==========================================================
# 5. Initialize Output Container
# ==========================================================

results_list <- list()



# ==========================================================
# 6. CLR Processing Loop (functional_level)
# ==========================================================

# Outer loop: functional resolution level
for (data_type in functional_level) {
  
  # Inner loop: population subgroup
  for (comp_name in names(unique_group)) {
    
    # Extract subgroup label
    group_pair <- unique_group[[comp_name]]
    
    # Generate CLR-transformed table
    res <- prepare_clr_table(
      feature_table = tables[[data_type]],
      sample_metadata = metadata_CR,
      data_type = data_type,
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




############################################################
# Genus-Level Feature Processing
# Generate CLR-transformed genus abundance by subgroup
############################################################

# ==========================================================
# 1. Load Genus Abundance Data
# ==========================================================

# Data preparation
Counts_Genus_CR <- Counts_Genus %>%
  filter(Sample_ID %in% metadata_CR$Sample_ID) %>%
  as.data.frame()

# Assign Sample IDs as row names
rownames(Counts_Genus_CR) <- Counts_Genus_CR$Sample_ID



# ==========================================================
# 2. Verify Sample Order Consistency
# ==========================================================

# Confirm sample order matches functional table structure
identical(Counts_Genus_CR$Sample_ID, colnames(PW_CR))
identical(Counts_Genus_CR$Sample_ID, colnames(EC_CR))
identical(Counts_Genus_CR$Sample_ID, colnames(KO_CR_filtered))



# ==========================================================
# 3. Format for CLR Processing
# ==========================================================

# Transpose table: samples → columns, taxa → rows
Counts_Genus_CR_t <- as.data.frame(t(Counts_Genus_CR[-c(1:11)]))  



# ==========================================================
# 4. Define Analysis Inputs
# ==========================================================

# Taxonomic resolution level
rank_level <- c("Genus")

# Store table in iterable structure
tables <- list(Genus = Counts_Genus_CR_t)



# ==========================================================
# 5. Define Target Subgroups
# ==========================================================

# Each subgroup processed independently
unique_group <- list(
  Akha_CR  = c("Akha-CR"),
  Lahu_CR  = c("Lahu-CR"),
  Lisu_CR = c("Lisu-CR")
)





# ==========================================================
# 6. Initialize Output Container
# ==========================================================

results_list <- list()



# ==========================================================
# 7. CLR Processing Loop (Genus Level)
# ==========================================================

# Outer loop: functional resolution level
for (data_type in rank_level) {
  
  # Inner loop: population subgroup
  for (comp_name in names(unique_group)) {
    
    # Extract subgroup label
    group_pair <- unique_group[[comp_name]]
    
    # Generate CLR-transformed table
    res <- prepare_clr_table(
      feature_table = tables[[data_type]],
      sample_metadata = metadata_CR,
      data_type = data_type,
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











#====================================================
# Microbiome Functional Association Analysis
#====================================================


#----------------------------------------------------
# HAllA output analysis: Genus+Functions
#----------------------------------------------------

load("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_Genus_Fn.RData")
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis")


# Import KEGG Database Mapping
library(readr)
Mapped_EC_KEGG <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/Function_Mapping/Mapped_EC_KEGG.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

Mapped_KO_KEGG <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/Function_Mapping/Mapped_KO_KEGG.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)



# Bind feature description
# EC
EC_descp <- Mapped_EC_KEGG %>% 
  select(1:2) %>%
  as.data.frame() %>%
  filter(!duplicated(EC_id)) %>%
  as.data.frame()


EC_abs_descp <- EC_abs %>%
  mutate(feature = rownames(.)) %>%
  relocate(feature, .before = 1) %>%
  left_join(EC_descp,
            by = c("feature" = "EC_id")) %>%
  as.data.frame()

rownames(EC_abs_descp) <- EC_abs_descp$feature


# KOs
KO_descp <- Mapped_KO_KEGG %>% 
  select(1:2) %>%
  as.data.frame() %>%
  filter(!duplicated(ko_id)) %>%
  as.data.frame()


KO_abs_descp <- KO_abs %>%
  mutate(feature = rownames(.)) %>%
  relocate(feature, .before = 1) %>%
  left_join(KO_descp,
            by = c("feature" = "ko_id")) %>%
  as.data.frame()

rownames(KO_abs_descp) <- KO_abs_descp$feature









#------------------------------
# Geo HAllA analysis
#------------------------------

#------------------------------
# Import HAllA output
#------------------------------
library(readr)

#-------------- Geo: KOs ----------
KO_AkhaCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/KO_Genus_Akha_CM_output_Geo/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

KO_AkhaCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/KO_Genus_Akha_CM_output_Geo/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


KO_AkhaCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/KO_Genus_Akha_CR_output_Geo/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

KO_AkhaCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/KO_Genus_Akha_CR_output_Geo/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


KO_LahuCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/KO_Genus_Lahu_CM_output_Geo/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

KO_LahuCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/KO_Genus_Lahu_CM_output_Geo/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


KO_LahuCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/KO_Genus_Lahu_CR_output_Geo/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

KO_LahuCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/KO_Genus_Lahu_CR_output_Geo/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)





#----------------------------------------------------
#  HAllA-KO manipulation
#----------------------------------------------------

#--------- Helper function: HAllA-KO manipulation -------------
process_sig_clusters <- function(sigcluster_df, asso_df, KO_descp) {
  # Unnest cluster_X
  cluster_X_split <- sigcluster_df %>%
    mutate(cluster_X_list = strsplit(as.character(cluster_X), ";")) %>%
    unnest(cols = c(cluster_X_list)) %>%
    dplyr::select(cluster_rank, cluster_X = cluster_X_list, best_adjusted_pvalue)
  
  # Unnest cluster_Y
  cluster_Y_split <- sigcluster_df %>%
    mutate(cluster_Y_list = strsplit(as.character(cluster_Y), ";")) %>%
    unnest(cols = c(cluster_Y_list)) %>%
    dplyr::select(cluster_rank, cluster_Y = cluster_Y_list)
  
  # Join both and merge with association data
  result <- full_join(cluster_X_split, cluster_Y_split, by = "cluster_rank", relationship = "many-to-many") %>%
    left_join(asso_df, by = c("cluster_X" = "X_features", "cluster_Y" = "Y_features")) %>%
    filter((association >= 0.8 | association <= -0.8) & `q-values` < 0.05) %>%
    left_join(KO_descp, by = c("cluster_Y" = "ko_id"), relationship = "many-to-many")
  
  return(result)
}




#--- Geo ---
# Named list for association and sigcluster tables
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Geo_KO")

table_lookup <- list(
  KO_AkhaCM_Geo  = list(sig = KO_AkhaCM_sigcluster,  asso = KO_AkhaCM_asso),
  KO_AkhaCR_Geo  = list(sig = KO_AkhaCR_sigcluster,  asso = KO_AkhaCR_asso),
  KO_LahuCM_Geo  = list(sig = KO_LahuCM_sigcluster,  asso = KO_LahuCM_asso),
  KO_LahuCR_Geo  = list(sig = KO_LahuCR_sigcluster,  asso = KO_LahuCR_asso)
)

functional_level <- c("KO")

unique_group <- list(
  Akha_CM   = c("AkhaCM_Geo"),
  Akha_CR   = c("AkhaCR_Geo"),
  Lahu_CM   = c("LahuCM_Geo"),
  Lahu_CR   = c("LahuCR_Geo")
)


# Run Loop with Storage
sig_results <- list()

for (func in functional_level) {
  for (group in names(unique_group)) {
    suffix <- unique_group[[group]]
    key <- paste(func, suffix, sep = "_")  # e.g., "KO_AkhaCM"
    
    if (key %in% names(table_lookup)) {
      sigcluster_df <- table_lookup[[key]]$sig
      asso_df       <- table_lookup[[key]]$asso
      
      res <- process_sig_clusters(sigcluster_df, asso_df, KO_descp)
      sig_results[[key]] <- res
      
      message("Processed: ", key)
    } else {
      warning("No data found for key: ", key)
    }
  }
}


# Save each result to CSV
output_dir <- "sig_cluster_results"
dir.create(output_dir, showWarnings = FALSE)

for (name in names(sig_results)) {
  Export(sig_results[[name]], 
         file = file.path(output_dir, paste0(name, "_sig.xlsx")), 
         row.names = FALSE)
}


# Extract elements in sig_results
for (i in names(sig_results)) {
  assign(i, sig_results[[i]])
}


#----------------------------------------------------
# Shared Association Features: Geo
#----------------------------------------------------
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Geo_KO")


#---- Helper function: Shared Genus+KO Function HAllA for paired ----
get_shared_features_geo <- function(df1, df2) {
  
  shared_features <- df1[, c("cluster_X", "cluster_Y")] %>%
    inner_join(df2[, c("cluster_X", "cluster_Y")],
               by = c("cluster_X", "cluster_Y")) 
  
  shared_g1 <- semi_join(df1, shared_features, by = c("cluster_X", "cluster_Y"))
  shared_g2 <- semi_join(df2, shared_features, by = c("cluster_X", "cluster_Y"))
  
  list(
    shared_features = shared_features, 
    shared_g1 = shared_g1, 
    shared_g2 = shared_g2
  )
}



#---- Named list for association and significant cluster tables ----
table <- list(
  KO_AkhaCM_Geo, KO_AkhaCR_Geo,
  KO_LahuCM_Geo, KO_LahuCR_Geo
)


#---- Define Functional Levels and Shared Pairs ----
functional_level <- c("KO")

# Flatten table_ into named list of dataframes
group_labels <- c("AkhaCM_Geo", "AkhaCR_Geo", 
                  "LahuCM_Geo", "LahuCR_Geo")


names(table) <- unlist(lapply(functional_level, function(lvl) paste0(lvl, "_", group_labels)))


#---- Loop to Process All Functional Levels and Shared Pairs: paired ----
# Keep pairs separate
shared_geo_pairs <- list(
  Akha_Geo = c("AkhaCM_Geo", "AkhaCR_Geo"),
  Lahu_Geo = c("LahuCM_Geo", "LahuCR_Geo")
)

shared_geo_results <- list()

for (lvl in functional_level) {
  
  lvl_dir <- file.path("shared_geo", lvl)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_geo_pairs)) {
    
    pair <- shared_geo_pairs[[pair_name]]
    
    df1 <- table[[paste0(lvl, "_", pair[1])]]
    df2 <- table[[paste0(lvl, "_", pair[2])]]
    
    shared <- get_shared_features_geo(df1, df2)
    
    key_base <- paste0(lvl, "_", pair[1], "_", pair[2])
    
    shared_geo_results[[paste0(key_base, "_shared_features")]] <- shared$shared_features
    shared_geo_results[[paste0(key_base, "_", pair[1], "_shared")]] <- shared$shared_g1
    shared_geo_results[[paste0(key_base, "_", pair[2], "_shared")]] <- shared$shared_g2
    
    Export(shared$shared_features, file = file.path(lvl_dir, paste0(key_base, "_shared_features.xlsx")), row.names = FALSE)
    Export(shared$shared_g1, file = file.path(lvl_dir, paste0(key_base, "_", pair[1], "_shared.xlsx")), row.names = FALSE)
    Export(shared$shared_g2, file = file.path(lvl_dir, paste0(key_base, "_", pair[2], "_shared.xlsx")), row.names = FALSE)
  }
}


# Extract each result into environment
for (i in names(shared_geo_results)) {
  assign(i, shared_geo_results[[i]])
}





#---------------------------------------------
# Unique Association Features: Geo
#---------------------------------------------
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Geo_KO")

#----- Helper Function: Extract Unique Features -----
get_unique_features <- function(df1, df2_shared) {
  # Create keys
  keys1 <- paste(df1$cluster_X, df1$cluster_Y, sep = "|")
  shared_keys <- paste(df2_shared$cluster_X, df2_shared$cluster_Y, sep = "|")
  
  # Filter unique
  unique_df <- df1[!(keys1 %in% shared_keys), ]
  
  # Create merged feature name
  unique_df$merged_feature <- paste(unique_df$cluster_X, unique_df$cluster_Y, sep = "|")
  
  # Filter out "unclassified"
  #unique_df <- unique_df %>%
  #  filter(cluster_X != "unclassified")
  
  return(unique_df)
}



#----- Input Lists -----
df1_table <- list(
  KO_AkhaCM_Geo, KO_AkhaCR_Geo, 
  KO_LahuCM_Geo, KO_LahuCR_Geo
)

df2_geo <- list(
  KO_AkhaCM_Geo_AkhaCR_Geo_shared_features,
  KO_LahuCM_Geo_LahuCR_Geo_shared_features
)

functional_level <- c("KO")

unique_group <- c("AkhaCM_Geo", "AkhaCR_Geo", 
                  "LahuCM_Geo", "LahuCR_Geo")


# Keep pairs separate
shared_geo_pairs <- list(
  Akha_Geo = c("AkhaCM_Geo", "AkhaCR_Geo"),
  Lahu_Geo = c("LahuCM_Geo", "LahuCR_Geo")
)



#----- Loop -----
# 1- Iterate over each functional_level (KO)
# 2- Iterate over each shared_geo group 
# 3 - Match and pull the correct df1 from df1_table and df2 from df2_shared
# 4 - Create an output directory for each level
# 5 - Store results in a list unique_results
# 6 - Export each unique result as a .xlsx file

# Note: One unique group can appear in multiple shared sets
# So, perform looping through all relevant shared sets for each group ensures completeness.

# Unique paired
unique_geo <- list()
df1_names <- names(df1_table) <- paste0(functional_level, "_", unique_group)
df2_names <- names(df2_geo) <- paste0(functional_level, "_", names(shared_geo_pairs), "_shared")

for (level in functional_level) {
  # Create output dir
  lvl_dir <- file.path("unique_geo", level)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_geo_pairs)) {
    pair <- shared_geo_pairs[[pair_name]]
    
    # Loop through each group in the pair
    for (group in pair) {
      df1_name <- paste0(level, "_", group)
      df2_name <- paste0(level, "_", pair_name, "_shared")
      
      df1 <- df1_table[[df1_name]]
      df2 <- df2_geo[[df2_name]]
      
      if (!is.null(df1) && !is.null(df2)) {
        unique_features <- get_unique_features(df1, df2)
        result_name <- paste0(level, "_", group, "_uniquepaired")
        unique_geo[[result_name]] <- unique_features
        
        # Export to Excel
        Export(unique_features, file = file.path(lvl_dir, paste0(result_name, ".xlsx")), row.names = FALSE)
      }
    }
  }
}


# Extract each result into environment
for (i in names(unique_geo)) {
  assign(i, unique_geo[[i]])
}







#===============================================================================
# SHARED & SHIFTED FEATURES ANALYSIS: FUNCTION-LEVEL (Geo GROUPS)
#===============================================================================
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Geo_KO")

### ---- Helper Function: Shared Features and shifted features based on functions ----

#----------------- Helper function: get_shared_and_diff_geo -------------------------
get_shared_and_diff_geo <- function(data_list, pair_vector, level_prefix, KO_descp) {
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  
  # ---- Validate inputs ----
  if (length(pair_vector) != 2) {
    stop("pair_vector must contain exactly 2 group identifiers.")
  }
  
  # ---- Extract the three tables ----
  table1 <- data_list[[paste0(level_prefix, "_", pair_vector[1])]]
  table2 <- data_list[[paste0(level_prefix, "_", pair_vector[2])]]
  
  # ---- Add group labels ----
  table1$Group <- paste0(level_prefix, "_", pair_vector[1])
  table2$Group <- paste0(level_prefix, "_", pair_vector[2])
  
  # ---- Combine and create pair key ----
  sig_data <- bind_rows(table1, table2) %>%
    mutate(pair_key = paste(cluster_X, cluster_Y, sep = "::"))
  
  # ---- Identify shared pairs across all three groups ----
  group_keys <- sig_data %>%
    distinct(Group, pair_key) %>%
    group_by(pair_key) %>%
    summarise(n_groups = n(), .groups = "drop")
  
  shared_keys <- group_keys %>%
    filter(n_groups == 2) %>%
    pull(pair_key)
  
  shared_geo <- sig_data %>%
    filter(pair_key %in% shared_keys)
  
  # ---- Map cluster_Y to cluster_X per group ----
  cluster_map <- sig_data %>%
    select(Group, cluster_X, cluster_Y) %>%
    distinct()
  
  cluster_y_mapping <- cluster_map %>%
    group_by(cluster_Y, Group) %>%
    summarise(X_list = list(unique(cluster_X)), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = X_list)
  
  # ---- Identify differences between groups (pairwise logic retained) ----
  colnames_expected <- paste0(level_prefix, "_", pair_vector)
  
  diff_pairs <- cluster_y_mapping %>%
    filter(all(colnames_expected %in% colnames(.))) %>%
    filter(
      pmap_lgl(
        select(., all_of(colnames_expected)),
        # “group1 differs from all other groups”
        # group1 not equal group2
        ~ !setequal(..1, ..2) 
      ),
      lengths(.[[colnames_expected[1]]]) > 0,
      lengths(.[[colnames_expected[2]]]) > 0
    ) %>%
    mutate(across(
      where(is.list),
      ~ sapply(., function(x) {
        if (is.null(x)) NA_character_ else paste(x, collapse = "; ")
      })
    )) %>%
    as.data.frame() %>%
    left_join(KO_descp, by = c("cluster_Y" = "ko_id")) %>%
    select(cluster_Y,
           all_of(colnames_expected[1:2]),
           ko_descp)
  
  return(list(
    shared = shared_geo,
    diff   = diff_pairs
  ))
}



#---- Run loop ----
# Table and metadata
table <- list(
  KO_AkhaCM_Geo, KO_AkhaCR_Geo, 
  KO_LahuCM_Geo, KO_LahuCR_Geo
)

functional_level <- c("KO")

# Flatten table_ into named list of dataframes
group_labels <- c("AkhaCM_Geo", "AkhaCR_Geo", 
                  "LahuCM_Geo", "LahuCR_Geo")


names(table) <- paste0(functional_level, "_", group_labels)


# Keep pairs separate
shared_geo_pairs <- list(
  Akha_Geo = c("AkhaCM_Geo", "AkhaCR_Geo"),
  Lahu_Geo = c("LahuCM_Geo", "LahuCR_Geo")
)


# Output storage
results <- list()

# Loop across functional levels and group pairs: paired
for (level in functional_level) {
  
  # Create output dir
  lvl_dir <- file.path("shared_shifted_results_geo", level)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_geo_pairs)) {
    pair <- shared_geo_pairs[[pair_name]]
    
    res <- get_shared_and_diff_geo(
      data_list = table,
      pair_vector = pair,
      level_prefix = level,
      KO_descp
    )
    
    results[[paste0(level, "_", pair_name)]] <- res
    
    # Naming and storing
    key_base <- paste0(level, "_", pair[1], "_", pair[2])
    results[[paste0(key_base, "_shared")]] <- res$shared
    results[[paste0(key_base, "_shifted")]] <- res$diff
    
    # Export to Excel in respective level directory
    Export(res$shared, file = file.path(lvl_dir, paste0(key_base, "_sharedpaired.xlsx")), row.names = FALSE)
    Export(res$diff,  file = file.path(lvl_dir, paste0(key_base, "_shiftedpaired.xlsx")), row.names = FALSE)
  }
}



# Extract each result into environment
for (i in names(results)) {
  assign(i, results[[i]])
}












#------------------------------
# Eth HAllA analysis
#------------------------------

#--- Eth ---
#---------- Eth: KOs ----------
KO_AkhaCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/KO_Genus_Akha_CM_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

KO_AkhaCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/KO_Genus_Akha_CM_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


KO_LahuCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/KO_Genus_Lahu_CM_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

KO_LahuCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/KO_Genus_Lahu_CM_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


KO_KhuenCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/KO_Genus_Khuen_CM_output_Eth/sig_clusters.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)

KO_KhuenCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/KO_Genus_Khuen_CM_output_Eth/all_associations.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)



KO_AkhaCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/KO_Genus_Akha_CR_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

KO_AkhaCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/KO_Genus_Akha_CR_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


KO_LahuCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/KO_Genus_Lahu_CR_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

KO_LahuCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/KO_Genus_Lahu_CR_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


KO_LisuCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/KO_Genus_Lisu_CR_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

KO_LisuCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/KO_Genus_Lisu_CR_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


# Named list for association and sigcluster tables
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Eth_KO")

table_lookup <- list(
  KO_AkhaCM_Eth  = list(sig = KO_AkhaCM_sigcluster,  asso = KO_AkhaCM_asso),
  KO_LahuCM_Eth  = list(sig = KO_LahuCM_sigcluster,  asso = KO_LahuCM_asso),
  KO_KhuenCM_Eth  = list(sig = KO_KhuenCM_sigcluster,  asso = KO_KhuenCM_asso),
  
  KO_AkhaCR_Eth  = list(sig = KO_AkhaCR_sigcluster,  asso = KO_AkhaCR_asso),
  KO_LahuCR_Eth  = list(sig = KO_LahuCR_sigcluster,  asso = KO_LahuCR_asso),
  KO_LisuCR_Eth  = list(sig = KO_LisuCR_sigcluster,  asso = KO_LisuCR_asso)
)

functional_level <- c("KO")

unique_group <- list(
  Akha_CM   = c("AkhaCM_Eth"),
  Lahu_CM   = c("LahuCM_Eth"),
  Khuen_CM   = c("KhuenCM_Eth"),
  Akha_CR   = c("AkhaCR_Eth"),
  Lahu_CR   = c("LahuCR_Eth"),
  Lisu_CR   = c("LisuCR_Eth")
)


# Run Loop with Storage
sig_results <- list()

for (func in functional_level) {
  for (group in names(unique_group)) {
    suffix <- unique_group[[group]]
    key <- paste(func, suffix, sep = "_")  # e.g., "KO_NCD0"
    
    if (key %in% names(table_lookup)) {
      sigcluster_df <- table_lookup[[key]]$sig
      asso_df       <- table_lookup[[key]]$asso
      
      res <- process_sig_clusters(sigcluster_df, asso_df, KO_descp)
      sig_results[[key]] <- res
      
      message("Processed: ", key)
    } else {
      warning("No data found for key: ", key)
    }
  }
}


# Save each result to CSV
output_dir <- "sig_cluster_results"
dir.create(output_dir, showWarnings = FALSE)

for (name in names(sig_results)) {
  Export(sig_results[[name]], 
         file = file.path(output_dir, paste0(name, "_sig.xlsx")), 
         row.names = FALSE)
}


# Extract elements in sig_results
for (i in names(sig_results)) {
  assign(i, sig_results[[i]])
}






#----------------------------------------------------
# Shared Association Features: Eth
#----------------------------------------------------
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Eth_KO")


#---- Helper function: Shared Genus+KO Function HAllA for paired ----
get_shared_features_Eth <- function(df1, df2, df3) {
  
  shared_features <- df1[, c("cluster_X", "cluster_Y")] %>%
    inner_join(df2[, c("cluster_X", "cluster_Y")],
               by = c("cluster_X", "cluster_Y")) %>%
    inner_join(df3[, c("cluster_X", "cluster_Y")],
               by = c("cluster_X", "cluster_Y"))
  
  shared_g1 <- semi_join(df1, shared_features, by = c("cluster_X", "cluster_Y"))
  shared_g2 <- semi_join(df2, shared_features, by = c("cluster_X", "cluster_Y"))
  shared_g3 <- semi_join(df3, shared_features, by = c("cluster_X", "cluster_Y"))
  
  list(
    shared_features = shared_features, 
    shared_g1 = shared_g1, 
    shared_g2 = shared_g2,
    shared_g3 = shared_g3
  )
}



#---- Named list for association and significant cluster tables ----
table <- list(
  KO_AkhaCM_Eth, KO_LahuCM_Eth, KO_KhuenCM_Eth,
  KO_AkhaCR_Eth, KO_LahuCR_Eth, KO_LisuCR_Eth
)


#---- Define Functional Levels and Shared Pairs ----
functional_level <- c("KO")

# Flatten table_ into named list of dataframes
group_labels <- c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth", 
                  "AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")


names(table) <- paste0("KO_", group_labels)


#---- Loop to Process All Functional Levels and Shared Pairs: paired ----
# Keep pairs separate
shared_Eth_pairs <- list(
  CM_Eth = c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth"),
  CR_Eth = c("AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")
)

shared_Eth_results <- list()

for (lvl in functional_level) {
  
  lvl_dir <- file.path("shared_Eth", lvl)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_Eth_pairs)) {
    
    pair <- shared_Eth_pairs[[pair_name]]
    
    df1 <- table[[paste0(lvl, "_", pair[1])]]
    df2 <- table[[paste0(lvl, "_", pair[2])]]
    df3 <- table[[paste0(lvl, "_", pair[3])]]
    
    
    # Shared-feature extraction
    shared <- get_shared_features_Eth(df1, df2, df3)
    
    # Naming and storing
    key_base <- paste0(lvl, "_", pair[1], "_", pair[2], "_", pair[3])
    shared_Eth_results[[paste0(key_base, "_shared_features")]] <- shared$shared_features
    shared_Eth_results[[paste0(key_base, "_", pair[1], "_shared")]] <- shared$shared_g1
    shared_Eth_results[[paste0(key_base, "_", pair[2], "_shared")]] <- shared$shared_g2
    shared_Eth_results[[paste0(key_base, "_", pair[3], "_shared")]] <- shared$shared_g3
    
    
    # Export to Excel in respective level directory
    Export(shared$shared_features, file = file.path(lvl_dir, paste0(key_base, "_shared_features.xlsx")), row.names = FALSE)
    Export(shared$shared_g1,       file = file.path(lvl_dir, paste0(key_base, "_", pair[1], "_shared.xlsx")), row.names = FALSE)
    Export(shared$shared_g2,       file = file.path(lvl_dir, paste0(key_base, "_", pair[2], "_shared.xlsx")), row.names = FALSE)
    Export(shared$shared_g3,       file = file.path(lvl_dir, paste0(key_base, "_", pair[3], "_shared.xlsx")), row.names = FALSE)
  }
}


# Extract each result into environment
for (i in names(shared_Eth_results)) {
  assign(i, shared_Eth_results[[i]])
}






#---------------------------------------------
# Unique Association Features: Eth
#---------------------------------------------
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Eth_KO")

#----- Helper Function: Extract Unique Features -----
get_unique_features <- function(df1, df2_shared) {
  # Create keys
  keys1 <- paste(df1$cluster_X, df1$cluster_Y, sep = "|")
  shared_keys <- paste(df2_shared$cluster_X, df2_shared$cluster_Y, sep = "|")
  
  # Filter unique
  unique_df <- df1[!(keys1 %in% shared_keys), ]
  
  # Create merged feature name
  unique_df$merged_feature <- paste(unique_df$cluster_X, unique_df$cluster_Y, sep = "|")
  
  # Filter out "unclassified"
  #unique_df <- unique_df %>%
  #  filter(cluster_X != "unclassified")
  
  return(unique_df)
}



#----- Input Lists -----
df1_table <- list(
  KO_AkhaCM_Eth, KO_LahuCM_Eth, KO_KhuenCM_Eth,
  KO_AkhaCR_Eth, KO_LahuCR_Eth, KO_LisuCR_Eth
)

df2_Eth <- list(
  KO_AkhaCM_Eth_LahuCM_Eth_KhuenCM_Eth_shared_features,
  KO_AkhaCR_Eth_LahuCR_Eth_LisuCR_Eth_shared_features
)

functional_level <- c("KO")

unique_group <- c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth", 
                  "AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")


# Keep pairs separate
shared_Eth_pairs <- list(
  CM_Eth = c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth"),
  CR_Eth = c("AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")
)



#----- Loop -----
# 1- Iterate over each functional_level (KO)
# 2- Iterate over each shared_Eth group 
# 3 - Match and pull the correct df1 from df1_table and df2 from df2_shared
# 4 - Create an output directory for each level
# 5 - Store results in a list unique_results
# 6 - Export each unique result as a .xlsx file

# Note: One unique group can appear in multiple shared sets
# So, perform looping through all relevant shared sets for each group ensures completeness.

# Unique paired
unique_Eth <- list()
df1_names <- names(df1_table) <- paste0(functional_level, "_", unique_group)
df2_names <- names(df2_Eth) <- paste0(functional_level, "_", names(shared_Eth_pairs), "_shared")

for (level in functional_level) {
  # Create output dir
  lvl_dir <- file.path("unique_Eth", level)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_Eth_pairs)) {
    pair <- shared_Eth_pairs[[pair_name]]
    
    # Loop through each group in the pair
    for (group in pair) {
      df1_name <- paste0(level, "_", group)
      df2_name <- paste0(level, "_", pair_name, "_shared")
      
      df1 <- df1_table[[df1_name]]
      df2 <- df2_Eth[[df2_name]]
      
      if (!is.null(df1) && !is.null(df2)) {
        unique_features <- get_unique_features(df1, df2)
        result_name <- paste0(level, "_", group, "_unique_Eth")
        unique_Eth[[result_name]] <- unique_features
        
        # Export to Excel
        Export(unique_features, file = file.path(lvl_dir, paste0(result_name, ".xlsx")), row.names = FALSE)
      }
    }
  }
}


# Extract each result into environment
for (i in names(unique_Eth)) {
  assign(i, unique_Eth[[i]])
}








#===============================================================================
# SHARED & SHIFTED FEATURES ANALYSIS: FUNCTION-LEVEL (ETHNIC GROUPS)
#===============================================================================
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Eth_KO")

### ---- Helper Function: Shared Features and shifted features based on functions ----

#----------------- Helper function: get_shared_and_diff_Eth -------------------------
get_shared_and_diff_Eth <- function(data_list, pair_vector, level_prefix, KO_descp) {
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  
  # ---- Validate inputs ----
  if (length(pair_vector) != 3) {
    stop("pair_vector must contain exactly 3 group identifiers.")
  }
  
  # ---- Extract the three tables ----
  table1 <- data_list[[paste0(level_prefix, "_", pair_vector[1])]] %>% as.data.frame()
  table2 <- data_list[[paste0(level_prefix, "_", pair_vector[2])]] %>% as.data.frame()
  table3 <- data_list[[paste0(level_prefix, "_", pair_vector[3])]] %>% as.data.frame()
  
  if (any(sapply(list(table1, table2, table3), is.null))) {
    stop("One or more input tables are missing. Check naming consistency.")
  }
  
  # ---- Add group labels ----
  table1$Group <- paste0(level_prefix, "_", pair_vector[1])
  table2$Group <- paste0(level_prefix, "_", pair_vector[2])
  table3$Group <- paste0(level_prefix, "_", pair_vector[3])
  
  # ---- Combine and create pair key ----
  sig_data <- bind_rows(table1, table2, table3) %>%
    mutate(pair_key = paste(cluster_X, cluster_Y, sep = "::"))
  
  # ---- Identify shared pairs across all three groups ----
  group_keys <- sig_data %>%
    distinct(Group, pair_key) %>%
    group_by(pair_key) %>%
    summarise(n_groups = n(), .groups = "drop")
  
  shared_keys <- group_keys %>%
    filter(n_groups == 3) %>%
    pull(pair_key)
  
  shared_Eth <- sig_data %>%
    filter(pair_key %in% shared_keys)
  
  # ---- Map cluster_Y to cluster_X per group ----
  cluster_map <- sig_data %>%
    select(Group, cluster_X, cluster_Y) %>%
    distinct()
  
  cluster_y_mapping <- cluster_map %>%
    group_by(cluster_Y, Group) %>%
    summarise(X_list = list(unique(cluster_X)), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = X_list)
  
  # ---- Identify differences between groups (pairwise logic retained) ----
  colnames_expected <- paste0(level_prefix, "_", pair_vector)
  
  # DEFINE HERE (before pipe)
  group_names <- colnames_expected
  
  diff_pairs <- cluster_y_mapping %>%
    filter(all(colnames_expected %in% colnames(.))) %>%
    filter(
      pmap_lgl(
        select(., all_of(colnames_expected)),
        ~ !(setequal(..1, ..2) & setequal(..1, ..3))
      ),
      lengths(.[[colnames_expected[1]]]) > 0,
      lengths(.[[colnames_expected[2]]]) > 0,
      lengths(.[[colnames_expected[3]]]) > 0
    ) %>%
    
    # PATTERN CLASSIFICATION
    mutate(
      pattern = pmap_chr(
        select(., all_of(colnames_expected)),
        ~ {
          if (setequal(..1, ..2) && !setequal(..1, ..3)) {
            paste0(group_names[3], "_shifted")
          } else if (setequal(..1, ..3) && !setequal(..1, ..2)) {
            paste0(group_names[2], "_shifted")
          } else if (setequal(..2, ..3) && !setequal(..1, ..2)) {
            paste0(group_names[1], "_shifted")
          } else {
            "Complex" # <- “Complex” = none of the three groups share identical patterns (All three different, Partial overlap but not identical, Nested but unequal)
          }
        }
      )
    ) %>%
    
    # Converts list-columns into clean text (character) columns
    mutate(across(
      where(is.list),
      ~ purrr::map_chr(.x, function(x) {
        if (is.null(x) || length(x) == 0) {
          NA_character_
        } else {
          paste(x, collapse = "; ")
        }
      })
    )) %>%
    
    as.data.frame() %>%
    left_join(KO_descp, by = c("cluster_Y" = "ko_id")) %>%
    select(cluster_Y,
           all_of(colnames_expected[1:3]),
           pattern,
           ko_descp)
  
  return(list(
    shared = shared_Eth,
    diff   = diff_pairs
  ))
}



#---- Run loop ----
# Table and metadata
table <- list(
  KO_AkhaCM_Eth, KO_LahuCM_Eth, KO_KhuenCM_Eth,
  KO_AkhaCR_Eth, KO_LahuCR_Eth, KO_LisuCR_Eth
)

functional_level <- c("KO")

# Flatten table_ into named list of dataframes
group_labels <- c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth", 
                  "AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")


names(table) <- paste0(functional_level, "_", group_labels)


# Keep pairs separate
shared_Eth_pairs <- list(
  CM_Eth = c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth"),
  CR_Eth = c("AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")
)


# Output storage
results <- list()

# Loop across functional levels and group pairs: paired
for (level in functional_level) {
  
  # Create output dir
  lvl_dir <- file.path("shared_shifted_results_Eth", level)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_Eth_pairs)) {
    pair <- shared_Eth_pairs[[pair_name]]
    
    res <- get_shared_and_diff_Eth(
      data_list = table,
      pair_vector = pair,
      level_prefix = level,
      KO_descp
    )
    
    results[[paste0(level, "_", pair_name)]] <- res
    
    # Naming and storing
    key_base <- paste0(level, "_", pair[1], "_", pair[2], "_", pair[3])
    results[[paste0(key_base, "_shared")]] <- res$shared
    results[[paste0(key_base, "_shifted")]] <- res$diff
    
    # Export to Excel in respective level directory
    Export(res$shared, file = file.path(lvl_dir, paste0(key_base, "_sharedpaired.xlsx")), row.names = FALSE)
    Export(res$diff,  file = file.path(lvl_dir, paste0(key_base, "_shiftedpaired.xlsx")), row.names = FALSE)
  }
}



# Extract each result into environment
for (i in names(results)) {
  assign(i, results[[i]])
}














#------------------------------------------------------
# Top Association Analysis: Within top 10 cluster rank
#------------------------------------------------------

#------ Akha-CM Geo --------
# Top 10 Genus-KO cluster rank 
KO_AkhaCM_Geo_uniquepaired %>% 
  group_by(cluster_X) %>% 
  count() %>% 
  arrange(desc(n))


KO_AkhaCM_Geo_uniquepaired %>% 
  filter(cluster_rank <= 10) %>%
  group_by(cluster_X, cluster_rank) %>% 
  count() %>% 
  arrange(desc(n))



ES_AkhaCM <- KO_AkhaCM_Geo_uniquepaired %>%
  filter(cluster_X == "g__Escherichia-Shigella",
         cluster_rank <= 10) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"))

ES_AkhaCM %>%
  count(Level2, Level3) %>%
  arrange(desc(n))


rank_AkhaCM_Geo <- KO_AkhaCM_Geo_uniquepaired %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank_AkhaCM_Geo %>%
  group_by(cluster_rank) %>%
  count() %>%
  arrange(cluster_rank)




#------ Akha-CR Geo --------
# Top 10 Genus-KO cluster rank 
Bactroides_AkhaCR <- KO_AkhaCR_Geo_uniquepaired %>%
  filter(cluster_X == "g__Bacteroides") %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"))

Bactroides_AkhaCR %>%
  count(Level2, Level3) %>%
  arrange(desc(n))


Alistipes_AkhaCR <- KO_AkhaCR_Geo_uniquepaired %>%
  filter(cluster_X == "g__Alistipes",
         cluster_rank <= 10) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"))

Alistipes_AkhaCR %>%
  count(Level2, Level3) %>%
  arrange(desc(n))


Prevotella_AkhaCR <- KO_AkhaCR_Geo_uniquepaired %>%
  filter(cluster_X == "g__Prevotella",
         cluster_rank <= 10) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"))

Prevotella_AkhaCR %>%
  count(Level2, Level3) %>%
  arrange(desc(n))



Parabacteroides_AkhaCR <- KO_AkhaCR_Geo_uniquepaired %>%
  filter(cluster_X == "g__Parabacteroides",
         cluster_rank <= 10) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"))

Parabacteroides_AkhaCR %>%
  count(Level2, Level3) %>%
  arrange(desc(n))



Lachno_AkhaCR <- KO_AkhaCR_Geo_uniquepaired %>%
  filter(cluster_X == "g__Lachnospiraceae_NK4A136_group",
         cluster_rank <= 10) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"))

Lachno_AkhaCR %>%
  count(Level2, Level3) %>%
  arrange(desc(n))






rank_AkhaCR_Geo <- KO_AkhaCR_Geo_uniquepaired %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank_AkhaCR_Geo %>%
  group_by(cluster_rank) %>%
  count() %>%
  arrange(cluster_rank)



#------ Lahu-CM Geo --------
rank_LahuCM <- KO_LahuCM_Geo_uniquepaired %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank_LahuCM %>%
  group_by(cluster_rank) %>%
  count() %>%
  arrange(desc(n))

KO_LahuCM_Geo_uniquepaired %>% 
  filter(cluster_rank == 3) %>%
  distinct(cluster_X)

KO_LahuCM_Geo_uniquepaired %>% 
  filter(cluster_rank == 3) %>%
  count(cluster_X) %>%
  arrange(desc(n))


  
KO_LahuCM_Geo_uniquepaired %>%
  filter(cluster_rank == 2) %>%
  group_by(cluster_X) %>%
  count() %>%
  arrange(desc(n))
 

rank2_LahuCM <- KO_LahuCM_Geo_uniquepaired %>%
  filter(cluster_rank == 3) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank2_LahuCM %>%
  count(Level2, Level3) %>%
  arrange(desc(n))



Eubacterium_LahuCM <- KO_LahuCM_Geo_uniquepaired %>%
  filter(cluster_X == "g__[Eubacterium]_coprostanoligenes_group",
         cluster_rank == 3) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"))
Eubacterium_LahuCM %>%
  count(Level2, Level3) %>%
  arrange(desc(n))



Escherichia_LahuCM <- KO_LahuCM_Geo_uniquepaired %>%
  filter(cluster_X == "g__Escherichia-Shigella",
         cluster_rank == 3) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"))

Escherichia_LahuCM %>%
  count(Level2, Level3) %>%
  arrange(desc(n))




#------ Lahu-CR Geo --------
rank_LahuCR <- KO_LahuCR_Geo_uniquepaired %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank_LahuCR


rank10_LahuCR <- KO_LahuCR_Geo_uniquepaired %>%
  filter(cluster_rank <= 10) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank10_LahuCR %>%
  count(Level2, Level3) %>%
  arrange(desc(n))






#------ Akha-CM Eth --------
rank_AkhaCM_eth <- KO_AkhaCM_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank_AkhaCM_eth %>%
  group_by(cluster_rank) %>%
  count() %>%
  arrange(desc(n))

KO_AkhaCM_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank10_AkhaCM_Eth <- KO_AkhaCM_Eth_unique_Eth %>%
  filter(cluster_rank <= 10,
         association >= 0) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank10_AkhaCM_Eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))

Escherichia_AkhaCM_Eth <- KO_AkhaCM_Eth_unique_Eth %>%
  filter(cluster_X == "g__Escherichia-Shigella",
         cluster_rank <= 10) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"))



#------ Lahu-CM Eth --------
rank_LahuCM_eth <- KO_LahuCM_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank_LahuCM_eth %>%
  group_by(cluster_rank) %>%
  count() %>%
  arrange(desc(n))


KO_LahuCM_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


rank3_LahuCM_Eth <- KO_LahuCM_Eth_unique_Eth %>%
  filter(cluster_rank == 3) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank3_LahuCM_Eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))


Eubacterium_LahuCM_Eth <- KO_LahuCM_Eth_unique_Eth %>%
  filter(cluster_X == "g__[Eubacterium]_coprostanoligenes_group",
         cluster_rank == 3) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"))

Eubacterium_LahuCM_Eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))


Escherichia_LahuCM_Eth <- KO_LahuCM_Eth_unique_Eth %>%
  filter(cluster_X == "g__Escherichia-Shigella",
         cluster_rank == 3) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"))

Escherichia_LahuCM_Eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))

#------ Khuen-CM Eth --------
rank_KhuenCM_eth <- KO_KhuenCM_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank_KhuenCM_eth %>%
  group_by(cluster_rank) %>%
  count() %>%
  arrange(desc(n))


KO_KhuenCM_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


rank10_KhuenCM_Eth <- KO_KhuenCM_Eth_unique_Eth %>%
  filter(cluster_rank <= 10) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank10_KhuenCM_Eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))


Streptococcus_KhuenCM_Eth <- KO_KhuenCM_Eth_unique_Eth %>%
  filter(cluster_X == "g__Streptococcus",
         cluster_rank <= 10) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"))




#------ Akha-CR Eth --------
rank_AkhaCR_eth <- KO_AkhaCR_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


KO_AkhaCR_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


rank1_AkhaCR_Eth <- KO_AkhaCR_Eth_unique_Eth %>%
  filter(cluster_rank == 1) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank1_AkhaCR_Eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))

Bacteroides_AkhaCR_Eth <- KO_AkhaCR_Eth_unique_Eth %>%
  filter(cluster_X == "g__Bacteroides",
         cluster_rank <= 10) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"))






#------ Lahu-CR Eth --------
rank_LahuCR_eth <- KO_LahuCR_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


KO_LahuCR_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


rank10_LahuCR_Eth <- KO_LahuCR_Eth_unique_Eth %>%
  filter(cluster_rank <= 10) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank10_LahuCR_Eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))





#------ Lisu-CR Eth --------
rank_LisuCR_eth <- KO_LisuCR_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


KO_LisuCR_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


rank10_LisuCR_eth <- KO_LisuCR_Eth_unique_Eth %>%
  filter(cluster_rank <= 10) %>%
  left_join(Mapped_KO_KEGG,
            by = c("cluster_Y" = "ko_id",
                   "ko_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank10_LisuCR_eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))


rank10_LisuCR_eth %>%
  group_by(cluster_rank) %>%
  count()


#---------------------
# aldesx2 overlapping
#---------------------
sigaldex_Akha_Geo <- ALDEx2_KO_Akha_CM_vs_Akha_CR_filtered %>%
  filter(wi.eBH < 0.05)

sigaldex_Lahu_Geo <- ALDEx2_KO_Lahu_CM_vs_Lahu_CR_filtered %>%
  filter(wi.eBH < 0.05)

# Akha Geo
sigaldesx_ES_AkhaCM_GEO <- ES_AkhaCM %>%
  filter(cluster_Y %in% sigaldex_Akha_Geo$Feature) %>% distinct(cluster_X, cluster_Y)

sigaldex_Akha_Geo %>%
  filter(Feature %in% sigaldesx_ES_AkhaCM_GEO$cluster_Y) %>%
  filter(abs(effect) >= 1)


sigaldesx_Bactroides_AkhaCR <- Bactroides_AkhaCR %>%
  filter(cluster_Y %in% sigaldex_Akha_Geo$Feature) %>% distinct(cluster_X, cluster_Y)

sigaldex_Akha_Geo %>%
  filter(Feature %in% sigaldesx_Bactroides_AkhaCR$cluster_Y) %>%
  filter(abs(effect) >= 1) %>%
  print(n = 50)


sigaldesx_Alistipes_AkhaCR <-  Alistipes_AkhaCR %>%
  filter(cluster_Y %in% sigaldex_Akha_Geo$Feature) %>% distinct(cluster_X, cluster_Y) %>%
  as.data.frame()

sigaldesx_Alistipes_AkhaCR %>% 
  left_join(sigaldex_Akha_Geo,
            by = c("cluster_Y" = "Feature"))


sigaldex_Akha_Geo %>%
  filter(Feature %in% sigaldesx_Alistipes_AkhaCR$cluster_Y) %>%
  filter(abs(effect) >= 1) %>%
  print(n = 50)




# Lahu Geo
sigaldesx_Eubacterium_LahuCM <- Eubacterium_LahuCM %>%
  filter(cluster_Y %in% sigaldex_Lahu_Geo$Feature) %>% distinct(cluster_X, cluster_Y)

sigaldesx_Escherichia_LahuCM <- Escherichia_LahuCM %>%
  filter(cluster_Y %in% sigaldex_Lahu_Geo$Feature) %>% distinct(cluster_X, cluster_Y)

sigaldex_Lahu_Geo %>%
  filter(Feature %in% sigaldesx_Eubacterium_LahuCM$cluster_Y) %>%
  filter(abs(effect) >= 1) %>%
  print(n = 50)

sigaldex_Lahu_Geo %>%
  filter(Feature %in% sigaldesx_Escherichia_LahuCM$cluster_Y) %>%
  filter(abs(effect) >= 1) %>%
  print(n = 50)




# CM Eth
sigaldex_CM1 <- ALDEx2_KO_Akha_CM_vs_Khuen_CM_filtered %>%
  filter(wi.eBH < 0.05)

sigaldex_CM2 <- ALDEx2_KO_Akha_CM_vs_Lahu_CM_filtered %>%
  filter(wi.eBH < 0.05)

sigaldex_CM3 <- ALDEx2_KO_Lahu_CM_vs_Khuen_CM_filtered %>%
  filter(wi.eBH < 0.05)


sigaldesx_Escherichia_AkhaCM_Eth1 <- Escherichia_AkhaCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM1$Feature) 

sigaldesx_Escherichia_AkhaCM_Eth2 <- Escherichia_AkhaCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM3$Feature) 




sigaldesx_Escherichia_LahuCM_Eth1 <- Escherichia_LahuCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM1$Feature) 

sigaldesx_Escherichia_LahuCM_Eth2 <- Escherichia_LahuCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM3$Feature) 

sigaldesx_Eubacterium_LahuCM_Eth1 <- Eubacterium_LahuCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM1$Feature) 

sigaldesx_Eubacterium_LahuCM_Eth2 <- Eubacterium_LahuCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM3$Feature) 


sigaldesx_Streptococcus_KhuenCM_Eth1 <- Streptococcus_KhuenCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM1$Feature) 

sigaldesx_Streptococcus_KhuenCM_Eth2 <- Streptococcus_KhuenCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM3$Feature) 




# CR Eth
sigaldex_CR1 <- ALDEx2_KO_Akha_CR_vs_Lahu_CR_filtered %>%
  filter(wi.eBH < 0.05)

sigaldex_CR2 <- ALDEx2_KO_Akha_CR_vs_Lisu_CR_filtered %>%
  filter(wi.eBH < 0.05)

sigaldex_CR3 <- ALDEx2_KO_Lahu_CR_vs_Lisu_CR_filtered %>%
  filter(wi.eBH < 0.05)



sigaldesx_Bacteroides_AkhaCR_Eth1 <- Bacteroides_AkhaCR_Eth %>%
  filter(cluster_Y %in% sigaldex_CR1$Feature) %>% distinct(cluster_X, cluster_Y)

sigaldesx_Bacteroides_AkhaCR_Eth2 <- Bacteroides_AkhaCR_Eth %>%
  filter(cluster_Y %in% sigaldex_CR2$Feature) %>% distinct(cluster_X, cluster_Y)


sigaldex_CR1 %>%
  filter(Feature %in% sigaldesx_Bacteroides_AkhaCR_Eth1$cluster_Y) %>%
  filter(abs(effect) >= 1)


sigaldex_CR2 %>%
  filter(Feature %in% sigaldesx_Bacteroides_AkhaCR_Eth2$cluster_Y) %>%
  filter(abs(effect) >= 1)


















#-------------------------
# EC HAllA Analysis
#-------------------------
library(readr)
#-------------- Geo: ECs ----------
EC_AkhaCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/EC_Genus_Akha_CM_output_Geo/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

EC_AkhaCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/EC_Genus_Akha_CM_output_Geo/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


EC_AkhaCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/EC_Genus_Akha_CR_output_Geo/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

EC_AkhaCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/EC_Genus_Akha_CR_output_Geo/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


EC_LahuCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/EC_Genus_Lahu_CM_output_Geo/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

EC_LahuCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/EC_Genus_Lahu_CM_output_Geo/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


EC_LahuCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/EC_Genus_Lahu_CR_output_Geo/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

EC_LahuCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/EC_Genus_Lahu_CR_output_Geo/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)





#----------------------------------------------------
#  HAllA-EC manipulation
#----------------------------------------------------

#--------- Helper function: HAllA-EC manipulation -------------
process_sig_clusters <- function(sigcluster_df, asso_df, EC_descp) {
  # Unnest cluster_X
  cluster_X_split <- sigcluster_df %>%
    mutate(cluster_X_list = strsplit(as.character(cluster_X), ";")) %>%
    unnest(cols = c(cluster_X_list)) %>%
    dplyr::select(cluster_rank, cluster_X = cluster_X_list, best_adjusted_pvalue)
  
  # Unnest cluster_Y
  cluster_Y_split <- sigcluster_df %>%
    mutate(cluster_Y_list = strsplit(as.character(cluster_Y), ";")) %>%
    unnest(cols = c(cluster_Y_list)) %>%
    dplyr::select(cluster_rank, cluster_Y = cluster_Y_list)
  
  # Join both and merge with association data
  result <- full_join(cluster_X_split, cluster_Y_split, by = "cluster_rank", relationship = "many-to-many") %>%
    left_join(asso_df, by = c("cluster_X" = "X_features", "cluster_Y" = "Y_features")) %>%
    filter((association >= 0.8 | association <= -0.8) & `q-values` < 0.05) %>%
    left_join(EC_descp, by = c("cluster_Y" = "EC_id"), relationship = "many-to-many")
  
  return(result)
}




#--- Geo ---
# Named list for association and sigcluster tables
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Geo_EC")

table_lookup <- list(
  EC_AkhaCM_Geo  = list(sig = EC_AkhaCM_sigcluster,  asso = EC_AkhaCM_asso),
  EC_AkhaCR_Geo  = list(sig = EC_AkhaCR_sigcluster,  asso = EC_AkhaCR_asso),
  EC_LahuCM_Geo  = list(sig = EC_LahuCM_sigcluster,  asso = EC_LahuCM_asso),
  EC_LahuCR_Geo  = list(sig = EC_LahuCR_sigcluster,  asso = EC_LahuCR_asso)
)

functional_level <- c("EC")

unique_group <- list(
  Akha_CM   = c("AkhaCM_Geo"),
  Akha_CR   = c("AkhaCR_Geo"),
  Lahu_CM   = c("LahuCM_Geo"),
  Lahu_CR   = c("LahuCR_Geo")
)


# Run Loop with Storage
sig_results <- list()

for (func in functional_level) {
  for (group in names(unique_group)) {
    suffix <- unique_group[[group]]
    key <- paste(func, suffix, sep = "_")  # e.g., "EC_AkhaCM"
    
    if (key %in% names(table_lookup)) {
      sigcluster_df <- table_lookup[[key]]$sig
      asso_df       <- table_lookup[[key]]$asso
      
      res <- process_sig_clusters(sigcluster_df, asso_df, EC_descp)
      sig_results[[key]] <- res
      
      message("Processed: ", key)
    } else {
      warning("No data found for key: ", key)
    }
  }
}


# Save each result to CSV
output_dir <- "sig_cluster_results"
dir.create(output_dir, showWarnings = FALSE)

for (name in names(sig_results)) {
  Export(sig_results[[name]], 
         file = file.path(output_dir, paste0(name, "_sig.xlsx")), 
         row.names = FALSE)
}


# Extract elements in sig_results
for (i in names(sig_results)) {
  assign(i, sig_results[[i]])
}


#----------------------------------------------------
# Shared Association Features: Geo
#----------------------------------------------------
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Geo_EC")


#---- Helper function: Shared Genus+EC Function HAllA for paired ----
get_shared_features_geo <- function(df1, df2) {
  
  shared_features <- df1[, c("cluster_X", "cluster_Y")] %>%
    inner_join(df2[, c("cluster_X", "cluster_Y")],
               by = c("cluster_X", "cluster_Y")) 
  
  shared_g1 <- semi_join(df1, shared_features, by = c("cluster_X", "cluster_Y"))
  shared_g2 <- semi_join(df2, shared_features, by = c("cluster_X", "cluster_Y"))
  
  list(
    shared_features = shared_features, 
    shared_g1 = shared_g1, 
    shared_g2 = shared_g2
  )
}



#---- Named list for association and significant cluster tables ----
table <- list(
  EC_AkhaCM_Geo, EC_AkhaCR_Geo,
  EC_LahuCM_Geo, EC_LahuCR_Geo
)


#---- Define Functional Levels and Shared Pairs ----
functional_level <- c("EC")

# Flatten table_ into named list of dataframes
group_labels <- c("AkhaCM_Geo", "AkhaCR_Geo", 
                  "LahuCM_Geo", "LahuCR_Geo")


names(table) <- unlist(lapply(functional_level, function(lvl) paste0(lvl, "_", group_labels)))


#---- Loop to Process All Functional Levels and Shared Pairs: paired ----
# Keep pairs separate
shared_geo_pairs <- list(
  Akha_Geo = c("AkhaCM_Geo", "AkhaCR_Geo"),
  Lahu_Geo = c("LahuCM_Geo", "LahuCR_Geo")
)

shared_geo_results <- list()

for (lvl in functional_level) {
  
  lvl_dir <- file.path("shared_geo", lvl)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_geo_pairs)) {
    
    pair <- shared_geo_pairs[[pair_name]]
    
    df1 <- table[[paste0(lvl, "_", pair[1])]]
    df2 <- table[[paste0(lvl, "_", pair[2])]]
    
    shared <- get_shared_features_geo(df1, df2)
    
    key_base <- paste0(lvl, "_", pair[1], "_", pair[2])
    
    shared_geo_results[[paste0(key_base, "_shared_features")]] <- shared$shared_features
    shared_geo_results[[paste0(key_base, "_", pair[1], "_shared")]] <- shared$shared_g1
    shared_geo_results[[paste0(key_base, "_", pair[2], "_shared")]] <- shared$shared_g2
    
    Export(shared$shared_features, file = file.path(lvl_dir, paste0(key_base, "_shared_features.xlsx")), row.names = FALSE)
    Export(shared$shared_g1, file = file.path(lvl_dir, paste0(key_base, "_", pair[1], "_shared.xlsx")), row.names = FALSE)
    Export(shared$shared_g2, file = file.path(lvl_dir, paste0(key_base, "_", pair[2], "_shared.xlsx")), row.names = FALSE)
  }
}


# Extract each result into environment
for (i in names(shared_geo_results)) {
  assign(i, shared_geo_results[[i]])
}





#---------------------------------------------
# Unique Association Features: Geo
#---------------------------------------------
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Geo_EC")

#----- Helper Function: Extract Unique Features -----
get_unique_features <- function(df1, df2_shared) {
  # Create keys
  keys1 <- paste(df1$cluster_X, df1$cluster_Y, sep = "|")
  shared_keys <- paste(df2_shared$cluster_X, df2_shared$cluster_Y, sep = "|")
  
  # Filter unique
  unique_df <- df1[!(keys1 %in% shared_keys), ]
  
  # Create merged feature name
  unique_df$merged_feature <- paste(unique_df$cluster_X, unique_df$cluster_Y, sep = "|")
  
  # Filter out "unclassified"
  #unique_df <- unique_df %>%
  #  filter(cluster_X != "unclassified")
  
  return(unique_df)
}



#----- Input Lists -----
df1_table <- list(
  EC_AkhaCM_Geo, EC_AkhaCR_Geo, 
  EC_LahuCM_Geo, EC_LahuCR_Geo
)

df2_geo <- list(
  EC_AkhaCM_Geo_AkhaCR_Geo_shared_features,
  EC_LahuCM_Geo_LahuCR_Geo_shared_features
)

functional_level <- c("EC")

unique_group <- c("AkhaCM_Geo", "AkhaCR_Geo", 
                  "LahuCM_Geo", "LahuCR_Geo")


# Keep pairs separate
shared_geo_pairs <- list(
  Akha_Geo = c("AkhaCM_Geo", "AkhaCR_Geo"),
  Lahu_Geo = c("LahuCM_Geo", "LahuCR_Geo")
)



#----- Loop -----
# 1- Iterate over each functional_level (EC)
# 2- Iterate over each shared_geo group 
# 3 - Match and pull the correct df1 from df1_table and df2 from df2_shared
# 4 - Create an output directory for each level
# 5 - Store results in a list unique_results
# 6 - Export each unique result as a .xlsx file

# Note: One unique group can appear in multiple shared sets
# So, perform looping through all relevant shared sets for each group ensures completeness.

# Unique paired
unique_geo <- list()
df1_names <- names(df1_table) <- paste0(functional_level, "_", unique_group)
df2_names <- names(df2_geo) <- paste0(functional_level, "_", names(shared_geo_pairs), "_shared")

for (level in functional_level) {
  # Create output dir
  lvl_dir <- file.path("unique_geo", level)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_geo_pairs)) {
    pair <- shared_geo_pairs[[pair_name]]
    
    # Loop through each group in the pair
    for (group in pair) {
      df1_name <- paste0(level, "_", group)
      df2_name <- paste0(level, "_", pair_name, "_shared")
      
      df1 <- df1_table[[df1_name]]
      df2 <- df2_geo[[df2_name]]
      
      if (!is.null(df1) && !is.null(df2)) {
        unique_features <- get_unique_features(df1, df2)
        result_name <- paste0(level, "_", group, "_uniquepaired")
        unique_geo[[result_name]] <- unique_features
        
        # Export to Excel
        Export(unique_features, file = file.path(lvl_dir, paste0(result_name, ".xlsx")), row.names = FALSE)
      }
    }
  }
}


# Extract each result into environment
for (i in names(unique_geo)) {
  assign(i, unique_geo[[i]])
}







#===============================================================================
# SHARED & SHIFTED FEATURES ANALYSIS: FUNCTION-LEVEL (Geo GROUPS)
#===============================================================================
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Geo_EC")

### ---- Helper Function: Shared Features and shifted features based on functions ----

#----------------- Helper function: get_shared_and_diff_geo -------------------------
get_shared_and_diff_geo <- function(data_list, pair_vector, level_prefix, EC_descp) {
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  
  # ---- Validate inputs ----
  if (length(pair_vector) != 2) {
    stop("pair_vector must contain exactly 2 group identifiers.")
  }
  
  # ---- Extract the three tables ----
  table1 <- data_list[[paste0(level_prefix, "_", pair_vector[1])]]
  table2 <- data_list[[paste0(level_prefix, "_", pair_vector[2])]]
  
  # ---- Add group labels ----
  table1$Group <- paste0(level_prefix, "_", pair_vector[1])
  table2$Group <- paste0(level_prefix, "_", pair_vector[2])
  
  # ---- Combine and create pair key ----
  sig_data <- bind_rows(table1, table2) %>%
    mutate(pair_key = paste(cluster_X, cluster_Y, sep = "::"))
  
  # ---- Identify shared pairs across all three groups ----
  group_keys <- sig_data %>%
    distinct(Group, pair_key) %>%
    group_by(pair_key) %>%
    summarise(n_groups = n(), .groups = "drop")
  
  shared_keys <- group_keys %>%
    filter(n_groups == 2) %>%
    pull(pair_key)
  
  shared_geo <- sig_data %>%
    filter(pair_key %in% shared_keys)
  
  # ---- Map cluster_Y to cluster_X per group ----
  cluster_map <- sig_data %>%
    select(Group, cluster_X, cluster_Y) %>%
    distinct()
  
  cluster_y_mapping <- cluster_map %>%
    group_by(cluster_Y, Group) %>%
    summarise(X_list = list(unique(cluster_X)), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = X_list)
  
  # ---- Identify differences between groups (pairwise logic retained) ----
  colnames_expected <- paste0(level_prefix, "_", pair_vector)
  
  diff_pairs <- cluster_y_mapping %>%
    filter(all(colnames_expected %in% colnames(.))) %>%
    filter(
      pmap_lgl(
        select(., all_of(colnames_expected)),
        # “group1 differs from all other groups”
        # group1 not equal group2
        ~ !setequal(..1, ..2) 
      ),
      lengths(.[[colnames_expected[1]]]) > 0,
      lengths(.[[colnames_expected[2]]]) > 0
    ) %>%
    mutate(across(
      where(is.list),
      ~ sapply(., function(x) {
        if (is.null(x)) NA_character_ else paste(x, collapse = "; ")
      })
    )) %>%
    as.data.frame() %>%
    left_join(EC_descp, by = c("cluster_Y" = "EC_id")) %>%
    select(cluster_Y,
           all_of(colnames_expected[1:2]),
           EC_descp)
  
  return(list(
    shared = shared_geo,
    diff   = diff_pairs
  ))
}



#---- Run loop ----
# Table and metadata
table <- list(
  EC_AkhaCM_Geo, EC_AkhaCR_Geo, 
  EC_LahuCM_Geo, EC_LahuCR_Geo
)

functional_level <- c("EC")

# Flatten table_ into named list of dataframes
group_labels <- c("AkhaCM_Geo", "AkhaCR_Geo", 
                  "LahuCM_Geo", "LahuCR_Geo")


names(table) <- paste0(functional_level, "_", group_labels)


# Keep pairs separate
shared_geo_pairs <- list(
  Akha_Geo = c("AkhaCM_Geo", "AkhaCR_Geo"),
  Lahu_Geo = c("LahuCM_Geo", "LahuCR_Geo")
)


# Output storage
results <- list()

# Loop across functional levels and group pairs: paired
for (level in functional_level) {
  
  # Create output dir
  lvl_dir <- file.path("shared_shifted_results_geo", level)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_geo_pairs)) {
    pair <- shared_geo_pairs[[pair_name]]
    
    res <- get_shared_and_diff_geo(
      data_list = table,
      pair_vector = pair,
      level_prefix = level,
      EC_descp
    )
    
    results[[paste0(level, "_", pair_name)]] <- res
    
    # Naming and storing
    key_base <- paste0(level, "_", pair[1], "_", pair[2])
    results[[paste0(key_base, "_shared")]] <- res$shared
    results[[paste0(key_base, "_shifted")]] <- res$diff
    
    # Export to Excel in respective level directory
    Export(res$shared, file = file.path(lvl_dir, paste0(key_base, "_sharedpaired.xlsx")), row.names = FALSE)
    Export(res$diff,  file = file.path(lvl_dir, paste0(key_base, "_shiftedpaired.xlsx")), row.names = FALSE)
  }
}



# Extract each result into environment
for (i in names(results)) {
  assign(i, results[[i]])
}












#------------------------------
# Eth HAllA analysis
#------------------------------

#---------- Eth: ECs ----------
EC_AkhaCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/EC_Genus_Akha_CM_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

EC_AkhaCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/EC_Genus_Akha_CM_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


EC_LahuCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/EC_Genus_Lahu_CM_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

EC_LahuCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/EC_Genus_Lahu_CM_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


EC_KhuenCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/EC_Genus_Khuen_CM_output_Eth/sig_clusters.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)

EC_KhuenCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/EC_Genus_Khuen_CM_output_Eth/all_associations.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)



EC_AkhaCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/EC_Genus_Akha_CR_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

EC_AkhaCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/EC_Genus_Akha_CR_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


EC_LahuCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/EC_Genus_Lahu_CR_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

EC_LahuCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/EC_Genus_Lahu_CR_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


EC_LisuCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/EC_Genus_Lisu_CR_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

EC_LisuCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/EC_Genus_Lisu_CR_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


# Named list for association and sigcluster tables
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Eth_EC")

table_lookup <- list(
  EC_AkhaCM_Eth  = list(sig = EC_AkhaCM_sigcluster,  asso = EC_AkhaCM_asso),
  EC_LahuCM_Eth  = list(sig = EC_LahuCM_sigcluster,  asso = EC_LahuCM_asso),
  EC_KhuenCM_Eth  = list(sig = EC_KhuenCM_sigcluster,  asso = EC_KhuenCM_asso),
  
  EC_AkhaCR_Eth  = list(sig = EC_AkhaCR_sigcluster,  asso = EC_AkhaCR_asso),
  EC_LahuCR_Eth  = list(sig = EC_LahuCR_sigcluster,  asso = EC_LahuCR_asso),
  EC_LisuCR_Eth  = list(sig = EC_LisuCR_sigcluster,  asso = EC_LisuCR_asso)
)

functional_level <- c("EC")

unique_group <- list(
  Akha_CM   = c("AkhaCM_Eth"),
  Lahu_CM   = c("LahuCM_Eth"),
  Khuen_CM   = c("KhuenCM_Eth"),
  Akha_CR   = c("AkhaCR_Eth"),
  Lahu_CR   = c("LahuCR_Eth"),
  Lisu_CR   = c("LisuCR_Eth")
)


# Run Loop with Storage
sig_results <- list()

for (func in functional_level) {
  for (group in names(unique_group)) {
    suffix <- unique_group[[group]]
    key <- paste(func, suffix, sep = "_")  # e.g., "EC_NCD0"
    
    if (key %in% names(table_lookup)) {
      sigcluster_df <- table_lookup[[key]]$sig
      asso_df       <- table_lookup[[key]]$asso
      
      res <- process_sig_clusters(sigcluster_df, asso_df, EC_descp)
      sig_results[[key]] <- res
      
      message("Processed: ", key)
    } else {
      warning("No data found for key: ", key)
    }
  }
}


# Save each result to CSV
output_dir <- "sig_cluster_results"
dir.create(output_dir, showWarnings = FALSE)

for (name in names(sig_results)) {
  Export(sig_results[[name]], 
         file = file.path(output_dir, paste0(name, "_sig.xlsx")), 
         row.names = FALSE)
}


# Extract elements in sig_results
for (i in names(sig_results)) {
  assign(i, sig_results[[i]])
}






#----------------------------------------------------
# Shared Association Features: Eth
#----------------------------------------------------
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Eth_EC")


#---- Helper function: Shared Genus+EC Function HAllA for paired ----
get_shared_features_Eth <- function(df1, df2, df3) {
  
  shared_features <- df1[, c("cluster_X", "cluster_Y")] %>%
    inner_join(df2[, c("cluster_X", "cluster_Y")],
               by = c("cluster_X", "cluster_Y")) %>%
    inner_join(df3[, c("cluster_X", "cluster_Y")],
               by = c("cluster_X", "cluster_Y"))
  
  shared_g1 <- semi_join(df1, shared_features, by = c("cluster_X", "cluster_Y"))
  shared_g2 <- semi_join(df2, shared_features, by = c("cluster_X", "cluster_Y"))
  shared_g3 <- semi_join(df3, shared_features, by = c("cluster_X", "cluster_Y"))
  
  list(
    shared_features = shared_features, 
    shared_g1 = shared_g1, 
    shared_g2 = shared_g2,
    shared_g3 = shared_g3
  )
}



#---- Named list for association and significant cluster tables ----
table <- list(
  EC_AkhaCM_Eth, EC_LahuCM_Eth, EC_KhuenCM_Eth,
  EC_AkhaCR_Eth, EC_LahuCR_Eth, EC_LisuCR_Eth
)


#---- Define Functional Levels and Shared Pairs ----
functional_level <- c("EC")

# Flatten table_ into named list of dataframes
group_labels <- c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth", 
                  "AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")


names(table) <- paste0("EC_", group_labels)


#---- Loop to Process All Functional Levels and Shared Pairs: paired ----
# Keep pairs separate
shared_Eth_pairs <- list(
  CM_Eth = c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth"),
  CR_Eth = c("AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")
)

shared_Eth_results <- list()

for (lvl in functional_level) {
  
  lvl_dir <- file.path("shared_Eth", lvl)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_Eth_pairs)) {
    
    pair <- shared_Eth_pairs[[pair_name]]
    
    df1 <- table[[paste0(lvl, "_", pair[1])]]
    df2 <- table[[paste0(lvl, "_", pair[2])]]
    df3 <- table[[paste0(lvl, "_", pair[3])]]
    
    
    # Shared-feature extraction
    shared <- get_shared_features_Eth(df1, df2, df3)
    
    # Naming and storing
    key_base <- paste0(lvl, "_", pair[1], "_", pair[2], "_", pair[3])
    shared_Eth_results[[paste0(key_base, "_shared_features")]] <- shared$shared_features
    shared_Eth_results[[paste0(key_base, "_", pair[1], "_shared")]] <- shared$shared_g1
    shared_Eth_results[[paste0(key_base, "_", pair[2], "_shared")]] <- shared$shared_g2
    shared_Eth_results[[paste0(key_base, "_", pair[3], "_shared")]] <- shared$shared_g3
    
    
    # Export to Excel in respective level directory
    Export(shared$shared_features, file = file.path(lvl_dir, paste0(key_base, "_shared_features.xlsx")), row.names = FALSE)
    Export(shared$shared_g1,       file = file.path(lvl_dir, paste0(key_base, "_", pair[1], "_shared.xlsx")), row.names = FALSE)
    Export(shared$shared_g2,       file = file.path(lvl_dir, paste0(key_base, "_", pair[2], "_shared.xlsx")), row.names = FALSE)
    Export(shared$shared_g3,       file = file.path(lvl_dir, paste0(key_base, "_", pair[3], "_shared.xlsx")), row.names = FALSE)
  }
}


# Extract each result into environment
for (i in names(shared_Eth_results)) {
  assign(i, shared_Eth_results[[i]])
}






#---------------------------------------------
# Unique Association Features: Eth
#---------------------------------------------
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Eth_EC")

#----- Helper Function: Extract Unique Features -----
get_unique_features <- function(df1, df2_shared) {
  # Create keys
  keys1 <- paste(df1$cluster_X, df1$cluster_Y, sep = "|")
  shared_keys <- paste(df2_shared$cluster_X, df2_shared$cluster_Y, sep = "|")
  
  # Filter unique
  unique_df <- df1[!(keys1 %in% shared_keys), ]
  
  # Create merged feature name
  unique_df$merged_feature <- paste(unique_df$cluster_X, unique_df$cluster_Y, sep = "|")
  
  # Filter out "unclassified"
  #unique_df <- unique_df %>%
  #  filter(cluster_X != "unclassified")
  
  return(unique_df)
}



#----- Input Lists -----
df1_table <- list(
  EC_AkhaCM_Eth, EC_LahuCM_Eth, EC_KhuenCM_Eth,
  EC_AkhaCR_Eth, EC_LahuCR_Eth, EC_LisuCR_Eth
)

df2_Eth <- list(
  EC_AkhaCM_Eth_LahuCM_Eth_KhuenCM_Eth_shared_features,
  EC_AkhaCR_Eth_LahuCR_Eth_LisuCR_Eth_shared_features
)

functional_level <- c("EC")

unique_group <- c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth", 
                  "AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")


# Keep pairs separate
shared_Eth_pairs <- list(
  CM_Eth = c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth"),
  CR_Eth = c("AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")
)



#----- Loop -----
# 1- Iterate over each functional_level (EC)
# 2- Iterate over each shared_Eth group 
# 3 - Match and pull the correct df1 from df1_table and df2 from df2_shared
# 4 - Create an output directory for each level
# 5 - Store results in a list unique_results
# 6 - Export each unique result as a .xlsx file

# Note: One unique group can appear in multiple shared sets
# So, perform looping through all relevant shared sets for each group ensures completeness.

# Unique paired
unique_Eth <- list()
df1_names <- names(df1_table) <- paste0(functional_level, "_", unique_group)
df2_names <- names(df2_Eth) <- paste0(functional_level, "_", names(shared_Eth_pairs), "_shared")

for (level in functional_level) {
  # Create output dir
  lvl_dir <- file.path("unique_Eth", level)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_Eth_pairs)) {
    pair <- shared_Eth_pairs[[pair_name]]
    
    # Loop through each group in the pair
    for (group in pair) {
      df1_name <- paste0(level, "_", group)
      df2_name <- paste0(level, "_", pair_name, "_shared")
      
      df1 <- df1_table[[df1_name]]
      df2 <- df2_Eth[[df2_name]]
      
      if (!is.null(df1) && !is.null(df2)) {
        unique_features <- get_unique_features(df1, df2)
        result_name <- paste0(level, "_", group, "_unique_Eth")
        unique_Eth[[result_name]] <- unique_features
        
        # Export to Excel
        Export(unique_features, file = file.path(lvl_dir, paste0(result_name, ".xlsx")), row.names = FALSE)
      }
    }
  }
}


# Extract each result into environment
for (i in names(unique_Eth)) {
  assign(i, unique_Eth[[i]])
}








#===============================================================================
# SHARED & SHIFTED FEATURES ANALYSIS: FUNCTION-LEVEL (ETHNIC GROUPS)
#===============================================================================
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Eth_EC")

### ---- Helper Function: Shared Features and shifted features based on functions ----

#----------------- Helper function: get_shared_and_diff_Eth -------------------------
get_shared_and_diff_Eth <- function(data_list, pair_vector, level_prefix, EC_descp) {
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  
  # ---- Validate inputs ----
  if (length(pair_vector) != 3) {
    stop("pair_vector must contain exactly 3 group identifiers.")
  }
  
  # ---- Extract the three tables ----
  table1 <- data_list[[paste0(level_prefix, "_", pair_vector[1])]] %>% as.data.frame()
  table2 <- data_list[[paste0(level_prefix, "_", pair_vector[2])]] %>% as.data.frame()
  table3 <- data_list[[paste0(level_prefix, "_", pair_vector[3])]] %>% as.data.frame()
  
  if (any(sapply(list(table1, table2, table3), is.null))) {
    stop("One or more input tables are missing. Check naming consistency.")
  }
  
  # ---- Add group labels ----
  table1$Group <- paste0(level_prefix, "_", pair_vector[1])
  table2$Group <- paste0(level_prefix, "_", pair_vector[2])
  table3$Group <- paste0(level_prefix, "_", pair_vector[3])
  
  # ---- Combine and create pair key ----
  sig_data <- bind_rows(table1, table2, table3) %>%
    mutate(pair_key = paste(cluster_X, cluster_Y, sep = "::"))
  
  # ---- Identify shared pairs across all three groups ----
  group_keys <- sig_data %>%
    distinct(Group, pair_key) %>%
    group_by(pair_key) %>%
    summarise(n_groups = n(), .groups = "drop")
  
  shared_keys <- group_keys %>%
    filter(n_groups == 3) %>%
    pull(pair_key)
  
  shared_Eth <- sig_data %>%
    filter(pair_key %in% shared_keys)
  
  # ---- Map cluster_Y to cluster_X per group ----
  cluster_map <- sig_data %>%
    select(Group, cluster_X, cluster_Y) %>%
    distinct()
  
  cluster_y_mapping <- cluster_map %>%
    group_by(cluster_Y, Group) %>%
    summarise(X_list = list(unique(cluster_X)), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = X_list)
  
  # ---- Identify differences between groups (pairwise logic retained) ----
  colnames_expected <- paste0(level_prefix, "_", pair_vector)
  
  # DEFINE HERE (before pipe)
  group_names <- colnames_expected
  
  diff_pairs <- cluster_y_mapping %>%
    filter(all(colnames_expected %in% colnames(.))) %>%
    filter(
      pmap_lgl(
        select(., all_of(colnames_expected)),
        ~ !(setequal(..1, ..2) & setequal(..1, ..3))
      ),
      lengths(.[[colnames_expected[1]]]) > 0,
      lengths(.[[colnames_expected[2]]]) > 0,
      lengths(.[[colnames_expected[3]]]) > 0
    ) %>%
    
    # PATTERN CLASSIFICATION
    mutate(
      pattern = pmap_chr(
        select(., all_of(colnames_expected)),
        ~ {
          if (setequal(..1, ..2) && !setequal(..1, ..3)) {
            paste0(group_names[3], "_shifted")
          } else if (setequal(..1, ..3) && !setequal(..1, ..2)) {
            paste0(group_names[2], "_shifted")
          } else if (setequal(..2, ..3) && !setequal(..1, ..2)) {
            paste0(group_names[1], "_shifted")
          } else {
            "Complex" # <- “Complex” = none of the three groups share identical patterns (All three different, Partial overlap but not identical, Nested but unequal)
          }
        }
      )
    ) %>%
    
    # Converts list-columns into clean text (character) columns
    mutate(across(
      where(is.list),
      ~ purrr::map_chr(.x, function(x) {
        if (is.null(x) || length(x) == 0) {
          NA_character_
        } else {
          paste(x, collapse = "; ")
        }
      })
    )) %>%
    
    as.data.frame() %>%
    left_join(EC_descp, by = c("cluster_Y" = "EC_id")) %>%
    select(cluster_Y,
           all_of(colnames_expected[1:3]),
           pattern,
           EC_descp)
  
  return(list(
    shared = shared_Eth,
    diff   = diff_pairs
  ))
}



#---- Run loop ----
# Table and metadata
table <- list(
  EC_AkhaCM_Eth, EC_LahuCM_Eth, EC_KhuenCM_Eth,
  EC_AkhaCR_Eth, EC_LahuCR_Eth, EC_LisuCR_Eth
)

functional_level <- c("EC")

# Flatten table_ into named list of dataframes
group_labels <- c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth", 
                  "AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")


names(table) <- paste0(functional_level, "_", group_labels)


# Keep pairs separate
shared_Eth_pairs <- list(
  CM_Eth = c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth"),
  CR_Eth = c("AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")
)


# Output storage
results <- list()

# Loop across functional levels and group pairs: paired
for (level in functional_level) {
  
  # Create output dir
  lvl_dir <- file.path("shared_shifted_results_Eth", level)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_Eth_pairs)) {
    pair <- shared_Eth_pairs[[pair_name]]
    
    res <- get_shared_and_diff_Eth(
      data_list = table,
      pair_vector = pair,
      level_prefix = level,
      EC_descp
    )
    
    results[[paste0(level, "_", pair_name)]] <- res
    
    # Naming and storing
    key_base <- paste0(level, "_", pair[1], "_", pair[2], "_", pair[3])
    results[[paste0(key_base, "_shared")]] <- res$shared
    results[[paste0(key_base, "_shifted")]] <- res$diff
    
    # Export to Excel in respective level directory
    Export(res$shared, file = file.path(lvl_dir, paste0(key_base, "_sharedpaired.xlsx")), row.names = FALSE)
    Export(res$diff,  file = file.path(lvl_dir, paste0(key_base, "_shiftedpaired.xlsx")), row.names = FALSE)
  }
}



# Extract each result into environment
for (i in names(results)) {
  assign(i, results[[i]])
}








#------------------------------------------------------
# Top Association Analysis: Within top 10 cluster rank
#------------------------------------------------------

#------ Akha-CM Geo --------
# Top 10 Genus-EC cluster rank 
rank_AkhaCM_Geo <- EC_AkhaCM_Geo_uniquepaired %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank_AkhaCM_Geo %>%
  group_by(cluster_rank) %>%
  count() %>%
  arrange(cluster_rank)


ES_AkhaCM <- EC_AkhaCM_Geo_uniquepaired %>%
  filter(cluster_X == "g__Escherichia-Shigella",
         cluster_rank <= 10) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"))

ES_AkhaCM %>%
  count(Level2, Level3) %>%
  arrange(desc(n))







#------ Akha-CR Geo --------
# Top 10 Genus-EC cluster rank 
Bactroides_AkhaCR <- EC_AkhaCR_Geo_uniquepaired %>%
  filter(cluster_X == "g__Bacteroides", 
         cluster_rank <= 10) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"))

Bactroides_AkhaCR %>%
  count(Level1, Level2, Level3) %>%
  arrange(desc(n))


Alistipes_AkhaCR <- EC_AkhaCR_Geo_uniquepaired %>%
  filter(cluster_X == "g__Alistipes",
         cluster_rank <= 10) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"))

Alistipes_AkhaCR %>%
  count(Level2, Level3) %>%
  arrange(desc(n))


rank_AkhaCR_Geo <- EC_AkhaCR_Geo_uniquepaired %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank_AkhaCR_Geo %>%
  group_by(cluster_rank) %>%
  count() %>%
  arrange(cluster_rank)



#------ Lahu-CM Geo --------
rank_LahuCM <- EC_LahuCM_Geo_uniquepaired %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank_LahuCM %>%
  group_by(cluster_rank) %>%
  count() %>%
  arrange(desc(n))

EC_LahuCM_Geo_uniquepaired %>%
  filter(cluster_rank == 2) %>%
  group_by(cluster_X) %>%
  count() %>%
  arrange(desc(n))


rank2_LahuCM <- EC_LahuCM_Geo_uniquepaired %>%
  filter(cluster_rank == 2) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank2_LahuCM %>%
  count(Level2, Level3) %>%
  arrange(desc(n))



Eubacterium_LahuCM <- EC_LahuCM_Geo_uniquepaired %>%
  filter(cluster_X == "g__[Eubacterium]_coprostanoligenes_group",
         cluster_rank == 2) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"))
Eubacterium_LahuCM %>%
  count(Level2, Level3) %>%
  arrange(desc(n))



Escherichia_LahuCM <- EC_LahuCM_Geo_uniquepaired %>%
  filter(cluster_X == "g__Escherichia-Shigella",
         cluster_rank == 2) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"))

Escherichia_LahuCM %>%
  count(Level2, Level3) %>%
  arrange(desc(n))




#------ Lahu-CR Geo --------
rank_LahuCR <- EC_LahuCR_Geo_uniquepaired %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank_LahuCR


rank10_LahuCR <- EC_LahuCR_Geo_uniquepaired %>%
  filter(cluster_rank <= 10) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank10_LahuCR %>%
  count(Level2, Level3) %>%
  arrange(desc(n))






#------ Akha-CM Eth --------
rank_AkhaCM_eth <- EC_AkhaCM_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank_AkhaCM_eth %>%
  group_by(cluster_rank) %>%
  count() %>%
  arrange(desc(n))

EC_AkhaCM_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank10_AkhaCM_Eth <- EC_AkhaCM_Eth_unique_Eth %>%
  filter(cluster_rank <= 10,
         association >= 0) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank10_AkhaCM_Eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))

Escherichia_AkhaCM_Eth <- EC_AkhaCM_Eth_unique_Eth %>%
  filter(cluster_X == "g__Escherichia-Shigella",
         cluster_rank <= 10) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"))



#------ Lahu-CM Eth --------
rank_LahuCM_eth <- EC_LahuCM_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank_LahuCM_eth %>%
  group_by(cluster_rank) %>%
  count() %>%
  arrange(desc(n))


EC_LahuCM_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


rank3_LahuCM_Eth <- EC_LahuCM_Eth_unique_Eth %>%
  filter(cluster_rank == 3) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank3_LahuCM_Eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))


Eubacterium_LahuCM_Eth <- EC_LahuCM_Eth_unique_Eth %>%
  filter(cluster_X == "g__[Eubacterium]_coprostanoligenes_group",
         cluster_rank == 3) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"))

Eubacterium_LahuCM_Eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))


Escherichia_LahuCM_Eth <- EC_LahuCM_Eth_unique_Eth %>%
  filter(cluster_X == "g__Escherichia-Shigella",
         cluster_rank == 3) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"))

Escherichia_LahuCM_Eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))

#------ Khuen-CM Eth --------
rank_KhuenCM_eth <- EC_KhuenCM_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)

rank_KhuenCM_eth %>%
  group_by(cluster_rank) %>%
  count() %>%
  arrange(desc(n))


EC_KhuenCM_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


rank10_KhuenCM_Eth <- EC_KhuenCM_Eth_unique_Eth %>%
  filter(cluster_rank <= 10) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank10_KhuenCM_Eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))


Streptococcus_KhuenCM_Eth <- EC_KhuenCM_Eth_unique_Eth %>%
  filter(cluster_X == "g__Streptococcus",
         cluster_rank <= 10) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"))




#------ Akha-CR Eth --------
rank_AkhaCR_eth <- EC_AkhaCR_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


EC_AkhaCR_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


rank1_AkhaCR_Eth <- EC_AkhaCR_Eth_unique_Eth %>%
  filter(cluster_rank == 1) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank1_AkhaCR_Eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))

Bacteroides_AkhaCR_Eth <- EC_AkhaCR_Eth_unique_Eth %>%
  filter(cluster_X == "g__Bacteroides",
         cluster_rank <= 10) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"))






#------ Lahu-CR Eth --------
rank_LahuCR_eth <- EC_LahuCR_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


EC_LahuCR_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


rank10_LahuCR_Eth <- EC_LahuCR_Eth_unique_Eth %>%
  filter(cluster_rank <= 10) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank10_LahuCR_Eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))





#------ Lisu-CR Eth --------
rank_LisuCR_eth <- EC_LisuCR_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


EC_LisuCR_Eth_unique_Eth %>% 
  group_by(cluster_rank) %>%
  count(cluster_X) %>%
  arrange(cluster_rank)


rank10_LisuCR_eth <- EC_LisuCR_Eth_unique_Eth %>%
  filter(cluster_rank <= 10) %>%
  left_join(Mapped_EC_KEGG,
            by = c("cluster_Y" = "EC_id",
                   "EC_descp"), relationship = "many-to-many") %>%
  distinct(merged_feature, .keep_all = TRUE)

rank10_LisuCR_eth %>%
  count(Level2, Level3) %>%
  arrange(desc(n))



#---------------------
# aldesx2 overlapping
#---------------------

# Akha GEO
ALDEx2_EC_Akha_CM_vs_Akha_CR_filtered <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/ALDEx2/crossALDEx2_results/ALDEx2_EC_Akha-CM_vs_Akha-CR_filtered.tsv", 
                                                    delim = "\t", escape_double = FALSE, 
                                                    trim_ws = TRUE)

sigaldex_Akha_Geo <- ALDEx2_EC_Akha_CM_vs_Akha_CR_filtered %>%
  filter(wi.eBH < 0.05)



# Akha Geo
sigaldesx_ES_AkhaCM_GEO <- ES_AkhaCM %>%
  filter(cluster_Y %in% sigaldex_Akha_Geo$Feature) %>% distinct(cluster_X, cluster_Y)

sigaldex_Akha_Geo %>%
  filter(Feature %in% sigaldesx_ES_AkhaCM_GEO$cluster_Y) %>%
  filter(abs(effect) >= 1) %>%
  arrange(desc(abs(diff.btw))) %>%
  print(n = 50)


sigaldesx_Bactroides_AkhaCR <- Bactroides_AkhaCR %>%
  filter(cluster_Y %in% sigaldex_Akha_Geo$Feature) %>% distinct(cluster_X, cluster_Y)

sigaldex_Akha_Geo %>%
  filter(Feature %in% sigaldesx_Bactroides_AkhaCR$cluster_Y) %>%
  filter(abs(effect) >= 1) %>%
  arrange(desc(diff.btw)) %>%
  print(n = 50)


sigaldesx_Alistipes_AkhaCR <-  Alistipes_AkhaCR %>%
  filter(cluster_Y %in% sigaldex_Akha_Geo$Feature) %>% distinct(cluster_X, cluster_Y) %>%
  as.data.frame()

sigaldesx_Alistipes_AkhaCR %>% 
  left_join(sigaldex_Akha_Geo,
            by = c("cluster_Y" = "Feature"))


sigaldex_Akha_Geo %>%
  filter(Feature %in% sigaldesx_Alistipes_AkhaCR$cluster_Y) %>%
  filter(abs(effect) >= 1) %>%
  print(n = 50)




# Lahu Geo
sigaldex_Lahu_Geo <- ALDEx2_EC_Lahu_CM_vs_Lahu_CR_filtered %>%
  filter(wi.eBH < 0.05)

sigaldesx_Eubacterium_LahuCM <- Eubacterium_LahuCM %>%
  filter(cluster_Y %in% sigaldex_Lahu_Geo$Feature) %>% distinct(cluster_X, cluster_Y)

sigaldesx_Escherichia_LahuCM <- Escherichia_LahuCM %>%
  filter(cluster_Y %in% sigaldex_Lahu_Geo$Feature) %>% distinct(cluster_X, cluster_Y)

sigaldex_Lahu_Geo %>%
  filter(Feature %in% sigaldesx_Eubacterium_LahuCM$cluster_Y) %>%
  filter(abs(effect) >= 1) %>%
  print(n = 50)

sigaldex_Lahu_Geo %>%
  filter(Feature %in% sigaldesx_Escherichia_LahuCM$cluster_Y) %>%
  filter(abs(effect) >= 1) %>%
  print(n = 50)




# CM Eth
sigaldex_CM1 <- ALDEx2_EC_Akha_CM_vs_Khuen_CM_filtered %>%
  filter(wi.eBH < 0.05)

sigaldex_CM2 <- ALDEx2_EC_Akha_CM_vs_Lahu_CM_filtered %>%
  filter(wi.eBH < 0.05)

sigaldex_CM3 <- ALDEx2_EC_Lahu_CM_vs_Khuen_CM_filtered %>%
  filter(wi.eBH < 0.05)


sigaldesx_Escherichia_AkhaCM_Eth1 <- Escherichia_AkhaCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM1$Feature) 

sigaldesx_Escherichia_AkhaCM_Eth2 <- Escherichia_AkhaCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM3$Feature) 




sigaldesx_Escherichia_LahuCM_Eth1 <- Escherichia_LahuCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM1$Feature) 

sigaldesx_Escherichia_LahuCM_Eth2 <- Escherichia_LahuCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM3$Feature) 

sigaldesx_Eubacterium_LahuCM_Eth1 <- Eubacterium_LahuCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM1$Feature) 

sigaldesx_Eubacterium_LahuCM_Eth2 <- Eubacterium_LahuCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM3$Feature) 


sigaldesx_Streptococcus_KhuenCM_Eth1 <- Streptococcus_KhuenCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM1$Feature) 

sigaldesx_Streptococcus_KhuenCM_Eth2 <- Streptococcus_KhuenCM_Eth %>%
  filter(cluster_Y %in% sigaldex_CM3$Feature) 




# CR Eth
sigaldex_CR1 <- ALDEx2_EC_Akha_CR_vs_Lahu_CR_filtered %>%
  filter(wi.eBH < 0.05)

sigaldex_CR2 <- ALDEx2_EC_Akha_CR_vs_Lisu_CR_filtered %>%
  filter(wi.eBH < 0.05)

sigaldex_CR3 <- ALDEx2_EC_Lahu_CR_vs_Lisu_CR_filtered %>%
  filter(wi.eBH < 0.05)



sigaldesx_Bacteroides_AkhaCR_Eth1 <- Bacteroides_AkhaCR_Eth %>%
  filter(cluster_Y %in% sigaldex_CR1$Feature) %>% distinct(cluster_X, cluster_Y)

sigaldesx_Bacteroides_AkhaCR_Eth2 <- Bacteroides_AkhaCR_Eth %>%
  filter(cluster_Y %in% sigaldex_CR2$Feature) %>% distinct(cluster_X, cluster_Y)


sigaldex_CR1 %>%
  filter(Feature %in% sigaldesx_Bacteroides_AkhaCR_Eth1$cluster_Y) %>%
  filter(abs(effect) >= 1)


sigaldex_CR2 %>%
  filter(Feature %in% sigaldesx_Bacteroides_AkhaCR_Eth2$cluster_Y) %>%
  filter(abs(effect) >= 1)












#-------------------------
# PW HAllA Analysis
#-------------------------
library(readr)
#-------------- Geo: PWs ----------
PW_AkhaCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/PW_Genus_Akha_CM_output_Geo/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

PW_AkhaCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/PW_Genus_Akha_CM_output_Geo/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


PW_AkhaCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/PW_Genus_Akha_CR_output_Geo/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

PW_AkhaCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/PW_Genus_Akha_CR_output_Geo/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


PW_LahuCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/PW_Genus_Lahu_CM_output_Geo/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

PW_LahuCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/PW_Genus_Lahu_CM_output_Geo/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


PW_LahuCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/PW_Genus_Lahu_CR_output_Geo/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

PW_LahuCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo/PW_Genus_Lahu_CR_output_Geo/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)




#----------------------------------------------------
#  HAllA-PW manipulation
#----------------------------------------------------

#--------- Helper function: HAllA-PW manipulation -------------
process_sig_clusters <- function(sigcluster_df, asso_df) {
  # Unnest cluster_X
  cluster_X_split <- sigcluster_df %>%
    mutate(cluster_X_list = strsplit(as.character(cluster_X), ";")) %>%
    unnest(cols = c(cluster_X_list)) %>%
    dplyr::select(cluster_rank, cluster_X = cluster_X_list, best_adjusted_pvalue)
  
  # Unnest cluster_Y
  cluster_Y_split <- sigcluster_df %>%
    mutate(cluster_Y_list = strsplit(as.character(cluster_Y), ";")) %>%
    unnest(cols = c(cluster_Y_list)) %>%
    dplyr::select(cluster_rank, cluster_Y = cluster_Y_list)
  
  # Join both and merge with association data
  result <- full_join(cluster_X_split, cluster_Y_split, by = "cluster_rank", relationship = "many-to-many") %>%
    left_join(asso_df, by = c("cluster_X" = "X_features", "cluster_Y" = "Y_features")) %>%
    filter((association >= 0.8 | association <= -0.8) & `q-values` < 0.05) 
  
  return(result)
}




#--- Geo ---
# Named list for association and sigcluster tables
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Geo_PW")

table_lookup <- list(
  PW_AkhaCM_Geo  = list(sig = PW_AkhaCM_sigcluster,  asso = PW_AkhaCM_asso),
  PW_AkhaCR_Geo  = list(sig = PW_AkhaCR_sigcluster,  asso = PW_AkhaCR_asso),
  PW_LahuCM_Geo  = list(sig = PW_LahuCM_sigcluster,  asso = PW_LahuCM_asso),
  PW_LahuCR_Geo  = list(sig = PW_LahuCR_sigcluster,  asso = PW_LahuCR_asso)
)

functional_level <- c("PW")

unique_group <- list(
  Akha_CM   = c("AkhaCM_Geo"),
  Akha_CR   = c("AkhaCR_Geo"),
  Lahu_CM   = c("LahuCM_Geo"),
  Lahu_CR   = c("LahuCR_Geo")
)


# Run Loop with Storage
sig_results <- list()

for (func in functional_level) {
  for (group in names(unique_group)) {
    suffix <- unique_group[[group]]
    key <- paste(func, suffix, sep = "_")  # e.g., "PW_AkhaCM"
    
    if (key %in% names(table_lookup)) {
      sigcluster_df <- table_lookup[[key]]$sig
      asso_df       <- table_lookup[[key]]$asso
      
      res <- process_sig_clusters(sigcluster_df, asso_df)
      sig_results[[key]] <- res
      
      message("Processed: ", key)
    } else {
      warning("No data found for key: ", key)
    }
  }
}


# Save each result to CSV
output_dir <- "sig_cluster_results"
dir.create(output_dir, showWarnings = FALSE)

for (name in names(sig_results)) {
  Export(sig_results[[name]], 
         file = file.path(output_dir, paste0(name, "_sig.xlsx")), 
         row.names = FALSE)
}


# Extract elements in sig_results
for (i in names(sig_results)) {
  assign(i, sig_results[[i]])
}



#----------------------------------------------------
# Shared Association Features: Geo
#----------------------------------------------------
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Geo_PW")


#---- Helper function: Shared Genus+PW Function HAllA for paired ----
get_shared_features_geo <- function(df1, df2) {
  
  shared_features <- df1[, c("cluster_X", "cluster_Y")] %>%
    inner_join(df2[, c("cluster_X", "cluster_Y")],
               by = c("cluster_X", "cluster_Y")) 
  
  shared_g1 <- semi_join(df1, shared_features, by = c("cluster_X", "cluster_Y"))
  shared_g2 <- semi_join(df2, shared_features, by = c("cluster_X", "cluster_Y"))
  
  list(
    shared_features = shared_features, 
    shared_g1 = shared_g1, 
    shared_g2 = shared_g2
  )
}



#---- Named list for association and significant cluster tables ----
table <- list(
  PW_AkhaCM_Geo, PW_AkhaCR_Geo,
  PW_LahuCM_Geo, PW_LahuCR_Geo
)


#---- Define Functional Levels and Shared Pairs ----
functional_level <- c("PW")

# Flatten table_ into named list of dataframes
group_labels <- c("AkhaCM_Geo", "AkhaCR_Geo", 
                  "LahuCM_Geo", "LahuCR_Geo")


names(table) <- unlist(lapply(functional_level, function(lvl) paste0(lvl, "_", group_labels)))


#---- Loop to Process All Functional Levels and Shared Pairs: paired ----
# Keep pairs separate
shared_geo_pairs <- list(
  Akha_Geo = c("AkhaCM_Geo", "AkhaCR_Geo"),
  Lahu_Geo = c("LahuCM_Geo", "LahuCR_Geo")
)

shared_geo_results <- list()

for (lvl in functional_level) {
  
  lvl_dir <- file.path("shared_geo", lvl)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_geo_pairs)) {
    
    pair <- shared_geo_pairs[[pair_name]]
    
    df1 <- table[[paste0(lvl, "_", pair[1])]]
    df2 <- table[[paste0(lvl, "_", pair[2])]]
    
    shared <- get_shared_features_geo(df1, df2)
    
    key_base <- paste0(lvl, "_", pair[1], "_", pair[2])
    
    shared_geo_results[[paste0(key_base, "_shared_features")]] <- shared$shared_features
    shared_geo_results[[paste0(key_base, "_", pair[1], "_shared")]] <- shared$shared_g1
    shared_geo_results[[paste0(key_base, "_", pair[2], "_shared")]] <- shared$shared_g2
    
    Export(shared$shared_features, file = file.path(lvl_dir, paste0(key_base, "_shared_features.xlsx")), row.names = FALSE)
    Export(shared$shared_g1, file = file.path(lvl_dir, paste0(key_base, "_", pair[1], "_shared.xlsx")), row.names = FALSE)
    Export(shared$shared_g2, file = file.path(lvl_dir, paste0(key_base, "_", pair[2], "_shared.xlsx")), row.names = FALSE)
  }
}


# Extract each result into environment
for (i in names(shared_geo_results)) {
  assign(i, shared_geo_results[[i]])
}




#---------------------------------------------
# Unique Association Features: Geo
#---------------------------------------------
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Geo_PW")

#----- Helper Function: Extract Unique Features -----
get_unique_features <- function(df1, df2_shared) {
  # Create keys
  keys1 <- paste(df1$cluster_X, df1$cluster_Y, sep = "|")
  shared_keys <- paste(df2_shared$cluster_X, df2_shared$cluster_Y, sep = "|")
  
  # Filter unique
  unique_df <- df1[!(keys1 %in% shared_keys), ]
  
  # Create merged feature name
  unique_df$merged_feature <- paste(unique_df$cluster_X, unique_df$cluster_Y, sep = "|")
  
  # Filter out "unclassified"
  #unique_df <- unique_df %>%
  #  filter(cluster_X != "unclassified")
  
  return(unique_df)
}



#----- Input Lists -----
df1_table <- list(
  PW_AkhaCM_Geo, PW_AkhaCR_Geo, 
  PW_LahuCM_Geo, PW_LahuCR_Geo
)

df2_geo <- list(
  PW_AkhaCM_Geo_AkhaCR_Geo_shared_features,
  PW_LahuCM_Geo_LahuCR_Geo_shared_features
)

functional_level <- c("PW")

unique_group <- c("AkhaCM_Geo", "AkhaCR_Geo", 
                  "LahuCM_Geo", "LahuCR_Geo")


# Keep pairs separate
shared_geo_pairs <- list(
  Akha_Geo = c("AkhaCM_Geo", "AkhaCR_Geo"),
  Lahu_Geo = c("LahuCM_Geo", "LahuCR_Geo")
)



#----- Loop -----
# 1- Iterate over each functional_level (PW)
# 2- Iterate over each shared_geo group 
# 3 - Match and pull the correct df1 from df1_table and df2 from df2_shared
# 4 - Create an output directory for each level
# 5 - Store results in a list unique_results
# 6 - Export each unique result as a .xlsx file

# Note: One unique group can appear in multiple shared sets
# So, perform looping through all relevant shared sets for each group ensures completeness.

# Unique paired
unique_geo <- list()
df1_names <- names(df1_table) <- paste0(functional_level, "_", unique_group)
df2_names <- names(df2_geo) <- paste0(functional_level, "_", names(shared_geo_pairs), "_shared")

for (level in functional_level) {
  # Create output dir
  lvl_dir <- file.path("unique_geo", level)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_geo_pairs)) {
    pair <- shared_geo_pairs[[pair_name]]
    
    # Loop through each group in the pair
    for (group in pair) {
      df1_name <- paste0(level, "_", group)
      df2_name <- paste0(level, "_", pair_name, "_shared")
      
      df1 <- df1_table[[df1_name]]
      df2 <- df2_geo[[df2_name]]
      
      if (!is.null(df1) && !is.null(df2)) {
        unique_features <- get_unique_features(df1, df2)
        result_name <- paste0(level, "_", group, "_uniquepaired")
        unique_geo[[result_name]] <- unique_features
        
        # Export to Excel
        Export(unique_features, file = file.path(lvl_dir, paste0(result_name, ".xlsx")), row.names = FALSE)
      }
    }
  }
}


# Extract each result into environment
for (i in names(unique_geo)) {
  assign(i, unique_geo[[i]])
}





#===============================================================================
# SHARED & SHIFTED FEATURES ANALYSIS: FUNCTION-LEVEL (Geo GROUPS)
#===============================================================================
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Geo_PW")

### ---- Helper Function: Shared Features and shifted features based on functions ----

#----------------- Helper function: get_shared_and_diff_geo -------------------------
get_shared_and_diff_geo <- function(data_list, pair_vector, level_prefix) {
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  
  # ---- Validate inputs ----
  if (length(pair_vector) != 2) {
    stop("pair_vector must contain exactly 2 group identifiers.")
  }
  
  # ---- Extract the three tables ----
  table1 <- data_list[[paste0(level_prefix, "_", pair_vector[1])]]
  table2 <- data_list[[paste0(level_prefix, "_", pair_vector[2])]]
  
  # ---- Add group labels ----
  table1$Group <- paste0(level_prefix, "_", pair_vector[1])
  table2$Group <- paste0(level_prefix, "_", pair_vector[2])
  
  # ---- Combine and create pair key ----
  sig_data <- bind_rows(table1, table2) %>%
    mutate(pair_key = paste(cluster_X, cluster_Y, sep = "::"))
  
  # ---- Identify shared pairs across all three groups ----
  group_keys <- sig_data %>%
    distinct(Group, pair_key) %>%
    group_by(pair_key) %>%
    summarise(n_groups = n(), .groups = "drop")
  
  shared_keys <- group_keys %>%
    filter(n_groups == 2) %>%
    pull(pair_key)
  
  shared_geo <- sig_data %>%
    filter(pair_key %in% shared_keys)
  
  # ---- Map cluster_Y to cluster_X per group ----
  cluster_map <- sig_data %>%
    select(Group, cluster_X, cluster_Y) %>%
    distinct()
  
  cluster_y_mapping <- cluster_map %>%
    group_by(cluster_Y, Group) %>%
    summarise(X_list = list(unique(cluster_X)), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = X_list)
  
  # ---- Identify differences between groups (pairwise logic retained) ----
  colnames_expected <- paste0(level_prefix, "_", pair_vector)
  
  diff_pairs <- cluster_y_mapping %>%
    filter(all(colnames_expected %in% colnames(.))) %>%
    filter(
      pmap_lgl(
        select(., all_of(colnames_expected)),
        # “group1 differs from all other groups”
        # group1 not equal group2
        ~ !setequal(..1, ..2) 
      ),
      lengths(.[[colnames_expected[1]]]) > 0,
      lengths(.[[colnames_expected[2]]]) > 0
    ) %>%
    mutate(across(
      where(is.list),
      ~ sapply(., function(x) {
        if (is.null(x)) NA_character_ else paste(x, collapse = "; ")
      })
    )) %>%
    as.data.frame()
  
  return(list(
    shared = shared_geo,
    diff   = diff_pairs
  ))
}



#---- Run loop ----
# Table and metadata
table <- list(
  PW_AkhaCM_Geo, PW_AkhaCR_Geo, 
  PW_LahuCM_Geo, PW_LahuCR_Geo
)

functional_level <- c("PW")

# Flatten table_ into named list of dataframes
group_labels <- c("AkhaCM_Geo", "AkhaCR_Geo", 
                  "LahuCM_Geo", "LahuCR_Geo")


names(table) <- paste0(functional_level, "_", group_labels)


# Keep pairs separate
shared_geo_pairs <- list(
  Akha_Geo = c("AkhaCM_Geo", "AkhaCR_Geo"),
  Lahu_Geo = c("LahuCM_Geo", "LahuCR_Geo")
)


# Output storage
results <- list()

# Loop across functional levels and group pairs: paired
for (level in functional_level) {
  
  # Create output dir
  lvl_dir <- file.path("shared_shifted_results_geo", level)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_geo_pairs)) {
    pair <- shared_geo_pairs[[pair_name]]
    
    res <- get_shared_and_diff_geo(
      data_list = table,
      pair_vector = pair,
      level_prefix = level
    )
    
    results[[paste0(level, "_", pair_name)]] <- res
    
    # Naming and storing
    key_base <- paste0(level, "_", pair[1], "_", pair[2])
    results[[paste0(key_base, "_shared")]] <- res$shared
    results[[paste0(key_base, "_shifted")]] <- res$diff
    
    # Export to Excel in respective level directory
    Export(res$shared, file = file.path(lvl_dir, paste0(key_base, "_sharedpaired.xlsx")), row.names = FALSE)
    Export(res$diff,  file = file.path(lvl_dir, paste0(key_base, "_shiftedpaired.xlsx")), row.names = FALSE)
  }
}



# Extract each result into environment
for (i in names(results)) {
  assign(i, results[[i]])
}













#------------------------------
# Eth HAllA analysis
#------------------------------

#---------- Eth: PWs ----------
PW_AkhaCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/PW_Genus_Akha_CM_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

PW_AkhaCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/PW_Genus_Akha_CM_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


PW_LahuCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/PW_Genus_Lahu_CM_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

PW_LahuCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/PW_Genus_Lahu_CM_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


PW_KhuenCM_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/PW_Genus_Khuen_CM_output_Eth/sig_clusters.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)

PW_KhuenCM_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/PW_Genus_Khuen_CM_output_Eth/all_associations.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)



PW_AkhaCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/PW_Genus_Akha_CR_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

PW_AkhaCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/PW_Genus_Akha_CR_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


PW_LahuCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/PW_Genus_Lahu_CR_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

PW_LahuCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/PW_Genus_Lahu_CR_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


PW_LisuCR_sigcluster <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/PW_Genus_Lisu_CR_output_Eth/sig_clusters.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

PW_LisuCR_asso <- read_delim("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Eth/PW_Genus_Lisu_CR_output_Eth/all_associations.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


# Named list for association and sigcluster tables
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Eth_PW")

table_lookup <- list(
  PW_AkhaCM_Eth  = list(sig = PW_AkhaCM_sigcluster,  asso = PW_AkhaCM_asso),
  PW_LahuCM_Eth  = list(sig = PW_LahuCM_sigcluster,  asso = PW_LahuCM_asso),
  PW_KhuenCM_Eth  = list(sig = PW_KhuenCM_sigcluster,  asso = PW_KhuenCM_asso),
  
  PW_AkhaCR_Eth  = list(sig = PW_AkhaCR_sigcluster,  asso = PW_AkhaCR_asso),
  PW_LahuCR_Eth  = list(sig = PW_LahuCR_sigcluster,  asso = PW_LahuCR_asso),
  PW_LisuCR_Eth  = list(sig = PW_LisuCR_sigcluster,  asso = PW_LisuCR_asso)
)

functional_level <- c("PW")

unique_group <- list(
  Akha_CM   = c("AkhaCM_Eth"),
  Lahu_CM   = c("LahuCM_Eth"),
  Khuen_CM   = c("KhuenCM_Eth"),
  Akha_CR   = c("AkhaCR_Eth"),
  Lahu_CR   = c("LahuCR_Eth"),
  Lisu_CR   = c("LisuCR_Eth")
)


# Run Loop with Storage
sig_results <- list()

for (func in functional_level) {
  for (group in names(unique_group)) {
    suffix <- unique_group[[group]]
    key <- paste(func, suffix, sep = "_")  
    
    if (key %in% names(table_lookup)) {
      sigcluster_df <- table_lookup[[key]]$sig
      asso_df       <- table_lookup[[key]]$asso
      
      res <- process_sig_clusters(sigcluster_df, asso_df)
      sig_results[[key]] <- res
      
      message("Processed: ", key)
    } else {
      warning("No data found for key: ", key)
    }
  }
}


# Save each result to CSV
output_dir <- "sig_cluster_results"
dir.create(output_dir, showWarnings = FALSE)

for (name in names(sig_results)) {
  Export(sig_results[[name]], 
         file = file.path(output_dir, paste0(name, "_sig.xlsx")), 
         row.names = FALSE)
}


# Extract elements in sig_results
for (i in names(sig_results)) {
  assign(i, sig_results[[i]])
}






#----------------------------------------------------
# Shared Association Features: Eth
#----------------------------------------------------
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Eth_PW")


#---- Helper function: Shared Genus+PW Function HAllA for paired ----
get_shared_features_Eth <- function(df1, df2, df3) {
  
  shared_features <- df1[, c("cluster_X", "cluster_Y")] %>%
    inner_join(df2[, c("cluster_X", "cluster_Y")],
               by = c("cluster_X", "cluster_Y")) %>%
    inner_join(df3[, c("cluster_X", "cluster_Y")],
               by = c("cluster_X", "cluster_Y"))
  
  shared_g1 <- semi_join(df1, shared_features, by = c("cluster_X", "cluster_Y"))
  shared_g2 <- semi_join(df2, shared_features, by = c("cluster_X", "cluster_Y"))
  shared_g3 <- semi_join(df3, shared_features, by = c("cluster_X", "cluster_Y"))
  
  list(
    shared_features = shared_features, 
    shared_g1 = shared_g1, 
    shared_g2 = shared_g2,
    shared_g3 = shared_g3
  )
}



#---- Named list for association and significant cluster tables ----
table <- list(
  PW_AkhaCM_Eth, PW_LahuCM_Eth, PW_KhuenCM_Eth,
  PW_AkhaCR_Eth, PW_LahuCR_Eth, PW_LisuCR_Eth
)


#---- Define Functional Levels and Shared Pairs ----
functional_level <- c("PW")

# Flatten table_ into named list of dataframes
group_labels <- c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth", 
                  "AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")


names(table) <- paste0("PW_", group_labels)


#---- Loop to Process All Functional Levels and Shared Pairs: paired ----
# Keep pairs separate
shared_Eth_pairs <- list(
  CM_Eth = c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth"),
  CR_Eth = c("AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")
)

shared_Eth_results <- list()

for (lvl in functional_level) {
  
  lvl_dir <- file.path("shared_Eth", lvl)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_Eth_pairs)) {
    
    pair <- shared_Eth_pairs[[pair_name]]
    
    df1 <- table[[paste0(lvl, "_", pair[1])]]
    df2 <- table[[paste0(lvl, "_", pair[2])]]
    df3 <- table[[paste0(lvl, "_", pair[3])]]
    
    
    # Shared-feature extraction
    shared <- get_shared_features_Eth(df1, df2, df3)
    
    # Naming and storing
    key_base <- paste0(lvl, "_", pair[1], "_", pair[2], "_", pair[3])
    shared_Eth_results[[paste0(key_base, "_shared_features")]] <- shared$shared_features
    shared_Eth_results[[paste0(key_base, "_", pair[1], "_shared")]] <- shared$shared_g1
    shared_Eth_results[[paste0(key_base, "_", pair[2], "_shared")]] <- shared$shared_g2
    shared_Eth_results[[paste0(key_base, "_", pair[3], "_shared")]] <- shared$shared_g3
    
    
    # Export to Excel in respective level directory
    Export(shared$shared_features, file = file.path(lvl_dir, paste0(key_base, "_shared_features.xlsx")), row.names = FALSE)
    Export(shared$shared_g1,       file = file.path(lvl_dir, paste0(key_base, "_", pair[1], "_shared.xlsx")), row.names = FALSE)
    Export(shared$shared_g2,       file = file.path(lvl_dir, paste0(key_base, "_", pair[2], "_shared.xlsx")), row.names = FALSE)
    Export(shared$shared_g3,       file = file.path(lvl_dir, paste0(key_base, "_", pair[3], "_shared.xlsx")), row.names = FALSE)
  }
}


# Extract each result into environment
for (i in names(shared_Eth_results)) {
  assign(i, shared_Eth_results[[i]])
}






#---------------------------------------------
# Unique Association Features: Eth
#---------------------------------------------
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Eth_PW")

#----- Helper Function: Extract Unique Features -----
get_unique_features <- function(df1, df2_shared) {
  # Create keys
  keys1 <- paste(df1$cluster_X, df1$cluster_Y, sep = "|")
  shared_keys <- paste(df2_shared$cluster_X, df2_shared$cluster_Y, sep = "|")
  
  # Filter unique
  unique_df <- df1[!(keys1 %in% shared_keys), ]
  
  # Create merged feature name
  unique_df$merged_feature <- paste(unique_df$cluster_X, unique_df$cluster_Y, sep = "|")
  
  # Filter out "unclassified"
  #unique_df <- unique_df %>%
  #  filter(cluster_X != "unclassified")
  
  return(unique_df)
}



#----- Input Lists -----
df1_table <- list(
  PW_AkhaCM_Eth, PW_LahuCM_Eth, PW_KhuenCM_Eth,
  PW_AkhaCR_Eth, PW_LahuCR_Eth, PW_LisuCR_Eth
)

df2_Eth <- list(
  PW_AkhaCM_Eth_LahuCM_Eth_KhuenCM_Eth_shared_features,
  PW_AkhaCR_Eth_LahuCR_Eth_LisuCR_Eth_shared_features
)

functional_level <- c("PW")

unique_group <- c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth", 
                  "AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")


# Keep pairs separate
shared_Eth_pairs <- list(
  CM_Eth = c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth"),
  CR_Eth = c("AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")
)



#----- Loop -----
# 1- Iterate over each functional_level (PW)
# 2- Iterate over each shared_Eth group 
# 3 - Match and pull the correct df1 from df1_table and df2 from df2_shared
# 4 - Create an output directory for each level
# 5 - Store results in a list unique_results
# 6 - Export each unique result as a .xlsx file

# Note: One unique group can appear in multiple shared sets
# So, perform looping through all relevant shared sets for each group ensures completeness.

# Unique paired
unique_Eth <- list()
df1_names <- names(df1_table) <- paste0(functional_level, "_", unique_group)
df2_names <- names(df2_Eth) <- paste0(functional_level, "_", names(shared_Eth_pairs), "_shared")

for (level in functional_level) {
  # Create output dir
  lvl_dir <- file.path("unique_Eth", level)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_Eth_pairs)) {
    pair <- shared_Eth_pairs[[pair_name]]
    
    # Loop through each group in the pair
    for (group in pair) {
      df1_name <- paste0(level, "_", group)
      df2_name <- paste0(level, "_", pair_name, "_shared")
      
      df1 <- df1_table[[df1_name]]
      df2 <- df2_Eth[[df2_name]]
      
      if (!is.null(df1) && !is.null(df2)) {
        unique_features <- get_unique_features(df1, df2)
        result_name <- paste0(level, "_", group, "_unique_Eth")
        unique_Eth[[result_name]] <- unique_features
        
        # Export to Excel
        Export(unique_features, file = file.path(lvl_dir, paste0(result_name, ".xlsx")), row.names = FALSE)
      }
    }
  }
}


# Extract each result into environment
for (i in names(unique_Eth)) {
  assign(i, unique_Eth[[i]])
}








#===============================================================================
# SHARED & SHIFTED FEATURES ANALYSIS: FUNCTION-LEVEL (ETHNIC GROUPS)
#===============================================================================
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_analysis/Eth_PW")

### ---- Helper Function: Shared Features and shifted features based on functions ----

#----------------- Helper function: get_shared_and_diff_Eth -------------------------
get_shared_and_diff_Eth <- function(data_list, pair_vector, level_prefix) {
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  
  # ---- Validate inputs ----
  if (length(pair_vector) != 3) {
    stop("pair_vector must contain exactly 3 group identifiers.")
  }
  
  # ---- Extract the three tables ----
  table1 <- data_list[[paste0(level_prefix, "_", pair_vector[1])]] %>% as.data.frame()
  table2 <- data_list[[paste0(level_prefix, "_", pair_vector[2])]] %>% as.data.frame()
  table3 <- data_list[[paste0(level_prefix, "_", pair_vector[3])]] %>% as.data.frame()
  
  if (any(sapply(list(table1, table2, table3), is.null))) {
    stop("One or more input tables are missing. Check naming consistency.")
  }
  
  # ---- Add group labels ----
  table1$Group <- paste0(level_prefix, "_", pair_vector[1])
  table2$Group <- paste0(level_prefix, "_", pair_vector[2])
  table3$Group <- paste0(level_prefix, "_", pair_vector[3])
  
  # ---- Combine and create pair key ----
  sig_data <- bind_rows(table1, table2, table3) %>%
    mutate(pair_key = paste(cluster_X, cluster_Y, sep = "::"))
  
  # ---- Identify shared pairs across all three groups ----
  group_keys <- sig_data %>%
    distinct(Group, pair_key) %>%
    group_by(pair_key) %>%
    summarise(n_groups = n(), .groups = "drop")
  
  shared_keys <- group_keys %>%
    filter(n_groups == 3) %>%
    pull(pair_key)
  
  shared_Eth <- sig_data %>%
    filter(pair_key %in% shared_keys)
  
  # ---- Map cluster_Y to cluster_X per group ----
  cluster_map <- sig_data %>%
    select(Group, cluster_X, cluster_Y) %>%
    distinct()
  
  cluster_y_mapping <- cluster_map %>%
    group_by(cluster_Y, Group) %>%
    summarise(X_list = list(unique(cluster_X)), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = X_list)
  
  # ---- Identify differences between groups (pairwise logic retained) ----
  colnames_expected <- paste0(level_prefix, "_", pair_vector)
  
  # DEFINE HERE (before pipe)
  group_names <- colnames_expected
  
  diff_pairs <- cluster_y_mapping %>%
    filter(all(colnames_expected %in% colnames(.))) %>%
    filter(
      pmap_lgl(
        select(., all_of(colnames_expected)),
        ~ !(setequal(..1, ..2) & setequal(..1, ..3))
      ),
      lengths(.[[colnames_expected[1]]]) > 0,
      lengths(.[[colnames_expected[2]]]) > 0,
      lengths(.[[colnames_expected[3]]]) > 0
    ) %>%
    
    # PATTERN CLASSIFICATION
    mutate(
      pattern = pmap_chr(
        select(., all_of(colnames_expected)),
        ~ {
          if (setequal(..1, ..2) && !setequal(..1, ..3)) {
            paste0(group_names[3], "_shifted")
          } else if (setequal(..1, ..3) && !setequal(..1, ..2)) {
            paste0(group_names[2], "_shifted")
          } else if (setequal(..2, ..3) && !setequal(..1, ..2)) {
            paste0(group_names[1], "_shifted")
          } else {
            "Complex" # <- “Complex” = none of the three groups share identical patterns (All three different, Partial overlap but not identical, Nested but unequal)
          }
        }
      )
    ) %>%
    
    # Converts list-columns into clean text (character) columns
    mutate(across(
      where(is.list),
      ~ purrr::map_chr(.x, function(x) {
        if (is.null(x) || length(x) == 0) {
          NA_character_
        } else {
          paste(x, collapse = "; ")
        }
      })
    )) %>%
    
    as.data.frame() 
  
  return(list(
    shared = shared_Eth,
    diff   = diff_pairs
  ))
}



#---- Run loop ----
# Table and metadata
table <- list(
  PW_AkhaCM_Eth, PW_LahuCM_Eth, PW_KhuenCM_Eth,
  PW_AkhaCR_Eth, PW_LahuCR_Eth, PW_LisuCR_Eth
)

functional_level <- c("PW")

# Flatten table_ into named list of dataframes
group_labels <- c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth", 
                  "AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")


names(table) <- paste0(functional_level, "_", group_labels)


# Keep pairs separate
shared_Eth_pairs <- list(
  CM_Eth = c("AkhaCM_Eth", "LahuCM_Eth", "KhuenCM_Eth"),
  CR_Eth = c("AkhaCR_Eth", "LahuCR_Eth", "LisuCR_Eth")
)


# Output storage
results <- list()

# Loop across functional levels and group pairs: paired
for (level in functional_level) {
  
  # Create output dir
  lvl_dir <- file.path("shared_shifted_results_Eth", level)
  dir.create(lvl_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pair_name in names(shared_Eth_pairs)) {
    pair <- shared_Eth_pairs[[pair_name]]
    
    res <- get_shared_and_diff_Eth(
      data_list = table,
      pair_vector = pair,
      level_prefix = level
    )
    
    results[[paste0(level, "_", pair_name)]] <- res
    
    # Naming and storing
    key_base <- paste0(level, "_", pair[1], "_", pair[2], "_", pair[3])
    results[[paste0(key_base, "_shared")]] <- res$shared
    results[[paste0(key_base, "_shifted")]] <- res$diff
    
    # Export to Excel in respective level directory
    Export(res$shared, file = file.path(lvl_dir, paste0(key_base, "_sharedpaired.xlsx")), row.names = FALSE)
    Export(res$diff,  file = file.path(lvl_dir, paste0(key_base, "_shiftedpaired.xlsx")), row.names = FALSE)
  }
}



# Extract each result into environment
for (i in names(results)) {
  assign(i, results[[i]])
}


