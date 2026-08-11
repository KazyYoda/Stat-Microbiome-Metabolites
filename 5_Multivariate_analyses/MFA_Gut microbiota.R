############################################################
# Multiple Factor Analysis (MFA)
############################################################
# Purpose:
# Perform Multiple Factor Analysis (MFA) to jointly analyze metadata,
# anthropometric variables, and significant gut microbiota taxa.
#
# Taxonomic variables are organized into separate MFA groups according to
# taxonomic rank:
#   - Phylum
#   - Class
#   - Order
#   - Family
#   - Genus
#
# The analysis is performed separately for the Akha and Lahu ethnic groups.


# ==============================
# 1. Environment Setup
# ==============================

# Set the working directory and load the processed taxonomic dataset.
setwd("~/Documents/HillTribe_NGS/6.MFA")
load("~/Documents/HillTribe_NGS/3.Compositional_profile/Heatmap_plot.RData")


# Load required packages -------------------------------------------------------
# Packages used for data manipulation, statistical analysis, MFA, and
# visualization.
library(dplyr)
library(readxl)
library(ggplot2)
library(rio)
library(car)  
library(FactoMineR)
library(factoextra)
library(ggplot2)



# ==========================================================
# 2. MFA – Akha Ethnic Group (Significant Taxa)
# ==========================================================
# Perform MFA for the Akha ethnic group using the significant taxonomic
# features identified in the preceding taxonomic analysis.
#
# The MFA integrates:
#   1. Metadata
#   2. Anthropometric variables
#   3. Phylum-level taxa
#   4. Class-level taxa
#   5. Order-level taxa
#   6. Family-level taxa
#   7. Genus-level taxa
#
# The metadata group is treated as supplementary and therefore does not
# contribute to the construction of the principal dimensions.


# Prepare MFA dataset ----------------------------------------------------------
# Remove sample identifiers and variables not included in the MFA.
# Categorical variables are converted to factors.
# Age is relocated to maintain the intended variable order within the
# anthropometric group.
Akha_MFA <- sig_taxo_Akha %>%
  select(-Sample_ID, -Ethnicity, -BMl) %>%
  mutate(across(c(Gender, Group, Area, 
                  BMI_group, Hypertension), factor)) %>% 
  relocate(Age, .before = 7)


# Identify taxonomic columns by rank prefix -----------------------------------
# Identify columns belonging to each taxonomic rank using the corresponding
# taxonomic prefixes.
p_cols <- which(startsWith(colnames(Akha_MFA), "p__"))
c_cols <- which(startsWith(colnames(Akha_MFA), "c__"))
o_cols <- which(startsWith(colnames(Akha_MFA), "o__"))
f_cols <- which(startsWith(colnames(Akha_MFA), "f__"))
g_cols <- which(startsWith(colnames(Akha_MFA), "g__"))


# Check column ranges and counts per taxonomic level --------------------------
# Inspect the positions and number of variables assigned to each taxonomic
# rank. These checks are used to verify the taxonomic group structure before
# running the MFA.
range(p_cols)
range(c_cols)
range(o_cols)
range(f_cols)
range(g_cols)

length(p_cols)
length(c_cols)
length(o_cols)
length(f_cols)
length(g_cols)


# Run MFA ----------------------------------------------------------------------
# Integrate metadata, anthropometric variables, and taxonomic variables
# as separate MFA groups.
#
#   Group 1: Metadata
#     5 categorical variables
#
#   Group 2: Anthropometric
#     5 quantitative variables
#
#   Group 3: Phylum
#     5 quantitative variables
#
#   Group 4: Class
#     6 quantitative variables
#
#   Group 5: Order
#     9 quantitative variables
#
#   Group 6: Family
#     17 quantitative variables
#
#   Group 7: Genus
#     27 quantitative variables
#
# The metadata group is specified as supplementary using num.group.sup = 1.
res.mfa <- MFA(Akha_MFA,
               
               group = c(5,
                         5,
                         5,6,9,17,27),
               
               type = c("n",
                        "s",
                        "s","s","s","s","s"),
               
               name.group = c("Metadata",
                              "Anthropometric",
                              "Phylum", "Class", "Order", "Family", "Genus"),
               
               num.group.sup = 1,
               
               graph = FALSE
)




# ==========================================================
# 3. MFA – Lahu Ethnic Group
# ==========================================================
# Perform MFA for the Lahu ethnic group using the significant taxonomic
# features identified in the preceding taxonomic analysis.
#
# The MFA structure follows the same framework as the Akha analysis, but
# the number of significant taxa available at each taxonomic rank differs.


# Prepare MFA dataset ----------------------------------------------------------
# Remove sample identifiers and variables not included in the MFA.
# Categorical variables are converted to factors.
# Age is relocated to maintain the intended variable order within the
# anthropometric group.
Lahu_MFA <- sig_taxo_Lahu %>%
  select(-Sample_ID, -Ethnicity, -BMl) %>%
  mutate(across(c(Gender, Group, Area, 
                  BMI_group, Hypertension), factor)) %>%
  relocate(Age, .before = 7)


# Identify taxonomic ranks ----------------------------------------------------
# Identify columns belonging to each taxonomic rank using the corresponding
# taxonomic prefixes.
p_cols <- which(startsWith(colnames(Lahu_MFA), "p__"))
c_cols <- which(startsWith(colnames(Lahu_MFA), "c__"))
o_cols <- which(startsWith(colnames(Lahu_MFA), "o__"))
f_cols <- which(startsWith(colnames(Lahu_MFA), "f__"))
g_cols <- which(startsWith(colnames(Lahu_MFA), "g__"))


# Check the number of variables at each taxonomic level.
length(p_cols)
length(c_cols)
length(o_cols)
length(f_cols)
length(g_cols)


# Run MFA ----------------------------------------------------------------------
# Integrate metadata, anthropometric variables, and the available taxonomic
# ranks as separate MFA groups.
#
#   Group 1: Metadata
#     5 categorical variables
#
#   Group 2: Anthropometric
#     5 quantitative variables
#
#   Group 3: Phylum
#     3 quantitative variables
#
#   Group 4: Class
#     3 quantitative variables
#
#   Group 5: Order
#     2 quantitative variables
#
#   Group 6: Family
#     4 quantitative variables
#
# Genus-level variables are not included in the MFA group specification for
# the Lahu dataset.
#
# The metadata group is specified as supplementary using num.group.sup = 1.
res.mfa <- MFA(Lahu_MFA,
               
               group = c(5,
                         5,
                         3,3,2,4),
               
               type = c("n",
                        "s",
                        "s","s","s","s"),
               
               name.group = c("Metadata",
                              "Anthropometric",
                              "Phylum", "Class", "Order", "Family"),
               
               num.group.sup = 1,
               
               graph = FALSE
)






# ==========================================================
# 4. MFA – CM Area
# ==========================================================
# Perform MFA for samples from Chiang Mai (CM) using the significant
# taxonomic features identified in the preceding taxonomic analysis.
#
# The MFA integrates:
#   1. Metadata
#   2. Anthropometric variables
#   3. Phylum-level taxa
#   4. Class-level taxa
#   5. Order-level taxa
#   6. Family-level taxa
#   7. Genus-level taxa
#
# The metadata group is treated as supplementary and therefore does not
# contribute to the construction of the principal dimensions.


# Prepare MFA dataset ----------------------------------------------------------
# Remove sample identifiers and variables not included in the MFA.
# Categorical variables are converted to factors.
# Age is relocated to maintain the intended variable order within the
# anthropometric group.
CM_MFA <- sig_taxo_CM %>%
  select(-Sample_ID, -Area, -BMl) %>%
  mutate(across(c(Gender, Group, Ethnicity, 
                  BMI_group, Hypertension), factor)) %>%
  relocate(Age, .before = 7)


# Identify taxonomic ranks ----------------------------------------------------
# Identify columns belonging to each taxonomic rank using the corresponding
# taxonomic prefixes.
p_cols <- which(startsWith(colnames(CM_MFA), "p__"))
c_cols <- which(startsWith(colnames(CM_MFA), "c__"))
o_cols <- which(startsWith(colnames(CM_MFA), "o__"))
f_cols <- which(startsWith(colnames(CM_MFA), "f__"))
g_cols <- which(startsWith(colnames(CM_MFA), "g__"))


# Check the number of variables at each taxonomic level.
length(p_cols)
length(c_cols)
length(o_cols)
length(f_cols)
length(g_cols)


# Run MFA ----------------------------------------------------------------------
# Integrate metadata, anthropometric variables, and taxonomic variables
# as separate MFA groups.
#
#   Group 1: Metadata
#     5 categorical variables
#
#   Group 2: Anthropometric
#     5 quantitative variables
#
#   Group 3: Phylum
#     3 quantitative variables
#
#   Group 4: Class
#     4 quantitative variables
#
#   Group 5: Order
#     6 quantitative variables
#
#   Group 6: Family
#     12 quantitative variables
#
#   Group 7: Genus
#     29 quantitative variables
#
# The metadata group is specified as supplementary using num.group.sup = 1.
res.mfa <- MFA(CM_MFA,
               
               group = c(5,
                         5,
                         3,4,6,12,29),
               
               type = c("n",
                        "s",
                        "s","s","s","s","s"),
               
               name.group = c("Metadata",
                              "Anthropometric",
                              "Phylum", "Class", "Order", "Family", "Genus"),
               
               num.group.sup = 1,
               
               graph = FALSE
)


# ==========================================================
# 5. MFA – CR Area
# ==========================================================
# Perform MFA for samples from Chiang Rai (CR) using the significant
# taxonomic features identified in the preceding taxonomic analysis.
#
# The same MFA framework is applied to the CR dataset, with the number of
# significant taxa at each taxonomic rank determined separately.


# Prepare MFA dataset ----------------------------------------------------------
# Remove sample identifiers and variables not included in the MFA.
# Categorical variables are converted to factors.
# Age is relocated to maintain the intended variable order within the
# anthropometric group.
CR_MFA <- sig_taxo_CR %>%
  select(-Sample_ID, -Area, -BMl) %>%
  mutate(across(c(Gender, Group, Ethnicity, 
                  BMI_group, Hypertension), factor)) %>%
  relocate(Age, .before = 7)


# Identify taxonomic ranks ----------------------------------------------------
# Identify columns belonging to each taxonomic rank using the corresponding
# taxonomic prefixes.
p_cols <- which(startsWith(colnames(CR_MFA), "p__"))
c_cols <- which(startsWith(colnames(CR_MFA), "c__"))
o_cols <- which(startsWith(colnames(CR_MFA), "o__"))
f_cols <- which(startsWith(colnames(CR_MFA), "f__"))
g_cols <- which(startsWith(colnames(CR_MFA), "g__"))


# Check the number of variables at each taxonomic level.
length(p_cols)
length(c_cols)
length(o_cols)
length(f_cols)
length(g_cols)


# Run MFA ----------------------------------------------------------------------
# Integrate metadata, anthropometric variables, and taxonomic variables
# as separate MFA groups.
#
#   Group 1: Metadata
#     5 categorical variables
#
#   Group 2: Anthropometric
#     5 quantitative variables
#
#   Group 3: Phylum
#     2 quantitative variables
#
#   Group 4: Class
#     3 quantitative variables
#
#   Group 5: Order
#     6 quantitative variables
#
#   Group 6: Family
#     12 quantitative variables
#
#   Group 7: Genus
#     39 quantitative variables
#
# The metadata group is specified as supplementary using num.group.sup = 1.
res.mfa <- MFA(CR_MFA,
               
               group = c(5,
                         5,
                         2,3,6,12,39),
               
               type = c("n",
                        "s",
                        "s","s","s","s","s"),
               
               name.group = c("Metadata",
                              "Anthropometric",
                              "Phylum", "Class", "Order", "Family", "Genus"),
               
               num.group.sup = 1,
               
               graph = FALSE
)




