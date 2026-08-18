############################################################
# Microbiome Analysis in R: HAllA Metabolite+Functions
############################################################


#--------------------------------------------------
# Data Preparation for HAllA: Metabolite+Function
#-------------------------------------------------

setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Metabolite_Function")
load("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_Genus_Fn.RData")
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
# Function Feature Processing
# Generate CLR-transformed Predicted function by group
############################################################
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Metabolite_Function/Geo")

# ==========================================================
# 1. Data preparation 
# ==========================================================

# PW
PW_Akha <- PW_abs %>%
  select(all_of(drop_Akha$Sample_ID)) %>%
  as.data.frame()

# EC
EC_Akha <- EC_abs %>%
  select(all_of(drop_Akha$Sample_ID)) %>%
  as.data.frame()

# KO
KO_Akha <- KO_abs %>%
  select(all_of(drop_Akha$Sample_ID)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample_ID = colnames(EC_Akha)) %>%
  relocate(Sample_ID, .before = 1) %>%
  left_join(drop_Akha[, c(1:2)], by = "Sample_ID") %>%
  as.data.frame()

rownames(KO_Akha) <- KO_Akha$Sample_ID



# Filter Genus|KOs with prevalence >= 20%
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
identical(colnames(PW_Akha), drop_Akha$Sample_ID)
identical(colnames(PW_Akha), rownames(met_05_Akha))

identical(colnames(EC_Akha), drop_Akha$Sample_ID)
identical(colnames(EC_Akha), rownames(met_05_Akha))

identical(colnames(KO_Akha_filtered), drop_Akha$Sample_ID)
identical(colnames(KO_Akha_filtered), rownames(met_05_Akha))



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
      sample_metadata = drop_Akha,
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
# Genus-Level Feature Processing
# Generate CLR-transformed genus count by group
############################################################

setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Metabolite_Function/Geo")

# ==========================================================
# 1. Data preparation 
# ==========================================================

# PW
PW_Lahu <- PW_abs %>%
  select(all_of(drop_Lahu$Sample_ID)) %>%
  as.data.frame()

# EC
EC_Lahu <- EC_abs %>%
  select(all_of(drop_Lahu$Sample_ID)) %>%
  as.data.frame()

# KO
KO_Lahu <- KO_abs %>%
  select(all_of(drop_Lahu$Sample_ID)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample_ID = colnames(EC_Lahu)) %>%
  relocate(Sample_ID, .before = 1) %>%
  left_join(drop_Lahu[, c(1:2)], by = "Sample_ID") %>%
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
identical(colnames(PW_Lahu), drop_Lahu$Sample_ID)
identical(colnames(PW_Lahu), rownames(met_05_Lahu))

identical(colnames(EC_Lahu), drop_Lahu$Sample_ID)
identical(colnames(EC_Lahu), rownames(met_05_Lahu))

identical(colnames(KO_Lahu_filtered), drop_Lahu$Sample_ID)
identical(colnames(KO_Lahu_filtered), rownames(met_05_Lahu))



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
      sample_metadata = drop_Lahu,
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
# CM location
# Genus-Level Feature Processing
# Generate CLR-transformed genus count by group
############################################################

setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Metabolite_Function/Eth")

# ==========================================================
# 1. Data preparation 
# ==========================================================

# PW
PW_CM <- PW_abs %>%
  select(all_of(drop_CM$Sample_ID)) %>%
  as.data.frame()

# EC
EC_CM <- EC_abs %>%
  select(all_of(drop_CM$Sample_ID)) %>%
  as.data.frame()

# KO
KO_CM <- KO_abs %>%
  select(all_of(drop_CM$Sample_ID)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample_ID = colnames(EC_CM)) %>%
  relocate(Sample_ID, .before = 1) %>%
  left_join(drop_CM[, c(1:2)], by = "Sample_ID") %>%
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
identical(colnames(PW_CM), drop_CM$Sample_ID)
identical(colnames(PW_CM), rownames(met_05_CM))

identical(colnames(EC_CM), drop_CM$Sample_ID)
identical(colnames(EC_CM), rownames(met_05_CM))

identical(colnames(KO_CM_filtered), drop_CM$Sample_ID)
identical(colnames(KO_CM_filtered), rownames(met_05_CM))


# ==========================================================
# 3. Define Analysis Inputs
# ==========================================================

# Functional hierarchy levels to process
functional_level <- c("PW", "EC", "KO")

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
# 6. CLR Processing Loop (Genus Level)
# ==========================================================

for (data_type in functional_level) {
  for (comp_name in names(unique_group)) {
    
    # Extract subgroup label
    group_pair <- unique_group[[comp_name]]
    
    # Generate CLR-transformed genus table
    res <- prepare_clr_table(
      feature_table = tables[[data_type]],
      sample_metadata = drop_CM,
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
# CR location
# Genus-Level Feature Processing
# Generate CLR-transformed genus count by group
############################################################

setwd("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Metabolite_Function/Eth")

# ==========================================================
# 1. Data preparation 
# ==========================================================

# PW
PW_CR <- PW_abs %>%
  select(all_of(drop_CR$Sample_ID)) %>%
  as.data.frame()

# EC
EC_CR <- EC_abs %>%
  select(all_of(drop_CR$Sample_ID)) %>%
  as.data.frame()

# KO
KO_CR <- KO_abs %>%
  select(all_of(drop_CR$Sample_ID)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample_ID = colnames(EC_CR)) %>%
  relocate(Sample_ID, .before = 1) %>%
  left_join(drop_CR[, c(1:2)], by = "Sample_ID") %>%
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
identical(colnames(PW_CR), drop_CR$Sample_ID)
identical(colnames(PW_CR), rownames(met_05_CR))

identical(colnames(EC_CR), drop_CR$Sample_ID)
identical(colnames(EC_CR), rownames(met_05_CR))

identical(colnames(KO_CR_filtered), drop_CR$Sample_ID)
identical(colnames(KO_CR_filtered), rownames(met_05_CR))


# ==========================================================
# 3. Define Analysis Inputs
# ==========================================================

# Functional hierarchy levels to process
functional_level <- c("PW", "EC", "KO")

tables <- list( 
  PW = PW_CR,
  EC = EC_CR,
  KO = KO_CR_filtered
)



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
# 5. Initialize Output Container
# ==========================================================

results_list <- list()



# ==========================================================
# 6. CLR Processing Loop (Genus Level)
# ==========================================================

for (data_type in functional_level) {
  for (comp_name in names(unique_group)) {
    
    # Extract subgroup label
    group_pair <- unique_group[[comp_name]]
    
    # Generate CLR-transformed genus table
    res <- prepare_clr_table(
      feature_table = tables[[data_type]],
      sample_metadata = drop_CR,
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



