############################################################
# Microbial Abundance Analysis
#
# Purpose:
# Compare relative abundances of microbial taxa between
# geographic locations within each ethnicity across multiple
# taxonomic ranks (Phylum, Class, Order, Family, and Genus).
#
# Statistical approach:
# - Wilcoxon rank sum test for comparisons between two areas
#   within each ethnicity.
# - Kruskal–Wallis Test followed by post hoc aalysis: Dunn’s Test (BH adjusted)
# - Benjamini-Hochberg (BH) adjustment for multiple
#   comparisons across taxa within each analysis.
############################################################


# ----------------------------------------
# 1. Set Working Directory & Load Results
# ----------------------------------------

# Set the working directory for compositional-profile
# statistical analyses.
setwd("~/Documents/HillTribe_NGS/3.Compositional_profile")

# Load the previously generated compositional-profile
# tables and statistical helper functions.
load("~/Documents/HillTribe_NGS/3.Compositional_profile/Compositional_profile.RData")
load("~/Documents/HillTribe_NGS/3.Compositional_profile/Stat_microbiota/Stat_microbiota.RData")


# Load packages required for data import, manipulation,
# reshaping, statistical analysis, and result export.
library(readxl)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr)
library(rio)
library(car)


# ----------------------------------------
# Import Relative Abundance Tables
# ----------------------------------------
# Import relative abundance tables generated during
# the compositional-profile analysis.
#
# Each table contains sample metadata and taxon-level
# relative abundances for the corresponding taxonomic rank.

Phylum <- read_excel("~/Documents/HillTribe_NGS/3.Compositional_profile/1_Rel_Phylum.xlsx")
Class <- read_excel("~/Documents/HillTribe_NGS/3.Compositional_profile/2_Rel_Class.xlsx")
Order <- read_excel("~/Documents/HillTribe_NGS/3.Compositional_profile/3_Rel_Order.xlsx")
Family <- read_excel("~/Documents/HillTribe_NGS/3.Compositional_profile/4_Rel_Family.xlsx")
Genus <- read_excel("~/Documents/HillTribe_NGS/3.Compositional_profile/5_Rel_Genus.xlsx")



# ---------------------------------------------------------------
# 2. Taxon Abundance Analysis by Taxonomic Rank
# ---------------------------------------------------------------
# Set the output directory for statistical results.
setwd("~/Documents/HillTribe_NGS/3.Compositional_profile/Stat_microbiota")


# ------------------------------------------------
# Wilcoxon Rank Sum Tests: Akha
# ------------------------------------------------
# Compare taxon relative abundances between CM and CR
# within the Akha group.
#
# The analysis is repeated independently for each
# taxonomic rank from Phylum to Genus.

taxa_wilcox(Phylum, 
           rank_name = "Phylum", 
           ethnicity_group = "Akha", 
           factor_col = "Area") 

taxa_wilcox(Class, 
            rank_name = "Class", 
            ethnicity_group = "Akha", 
            factor_col = "Area") 

taxa_wilcox(Order, 
            rank_name = "Order", 
            ethnicity_group = "Akha", 
            factor_col = "Area") 

taxa_wilcox(Family, 
            rank_name = "Family", 
            ethnicity_group = "Akha", 
            factor_col = "Area") 

taxa_wilcox(Genus, 
            rank_name = "Genus", 
            ethnicity_group = "Akha", 
            factor_col = "Area") 



# ------------------------------------------------
# Wilcoxon Rank Sum Tests: Lahu
# ------------------------------------------------
# Compare taxon relative abundances between CM and CR
# within the Lahu group.
#
# The analysis is repeated independently for each
# taxonomic rank from Phylum to Genus.

taxa_wilcox(Phylum, 
            rank_name = "Phylum", 
            ethnicity_group = "Lahu", 
            factor_col = "Area") 

taxa_wilcox(Class, 
            rank_name = "Class", 
            ethnicity_group = "Lahu", 
            factor_col = "Area") 

taxa_wilcox(Order, 
            rank_name = "Order", 
            ethnicity_group = "Lahu", 
            factor_col = "Area") 

taxa_wilcox(Family, 
            rank_name = "Family", 
            ethnicity_group = "Lahu", 
            factor_col = "Area") 

taxa_wilcox(Genus, 
            rank_name = "Genus", 
            ethnicity_group = "Lahu", 
            factor_col = "Area") 




# --------------------------------------------------
# Kruskal-Wallis Test: Phylum Level (Chiang Mai)
# --------------------------------------------------
#
# Objective:
# Identify phyla that differ in abundance among hill
# tribe groups within Chiang Mai (CM).
#
# Analysis workflow:
# 1. Perform a global Kruskal-Wallis test for each phylum.
# 2. Identify phyla with significant global differences.
# 3. Perform pairwise Dunn's tests for significant phyla.
# 4. Adjust pairwise P-values using the
#    Benjamini-Hochberg (BH) method.
# 5. Extract significant pairwise comparisons.
# --------------------------------------------------


# --------------------------------------------------
# Prepare Phylum-Level Data
# --------------------------------------------------
# Rename an invalid column name to ensure compatibility
# with downstream modeling and formula-based analyses.

Phylum <- Phylum %>% 
  dplyr::rename(p__WPS_2 = `p__WPS-2`)


# --------------------------------------------------
# Global Kruskal-Wallis Test
# --------------------------------------------------
# Test for differences in phylum-level abundance
# among hill tribe groups within Chiang Mai.
#
# The Kruskal-Wallis test is a non-parametric test
# used to evaluate whether the distributions of a
# feature differ among multiple independent groups.

stat_kruskal(
  Phylum, 
  rank_name   = "Phylum",
  group_col   = "Area",
  group_value = "CM",
  factor_col  = "Group",
  factor_levels = c("Akha-CM", "Lahu-CM", "Khuen-CM")
)


# --------------------------------------------------
# Post Hoc Analysis: Dunn's Test
# --------------------------------------------------
# Dunn's test is performed only for phyla showing
# a significant global Kruskal-Wallis test.
#
# Pairwise P-values are adjusted using the
# Benjamini-Hochberg (BH) method.


# Load Kruskal-Wallis summary results
kruskal_CM_Phylum <- read.delim(
  "~/Documents/HillTribe_NGS/3.Compositional_profile/Stat_microbiota/kruskal_CM_Phylum.txt"
)

# Extract phyla with a significant global
# Kruskal-Wallis test (P < 0.05).
sig_krus_CM <- kruskal_CM_Phylum %>% 
  filter(P_sig == "Yes") %>% 
  pull(Taxa)


# --------------------------------------------------
# Prepare Data for Dunn's Test
# --------------------------------------------------
# Create a combined grouping variable representing
# ethnicity and geographic area.
#
# The data are restricted to Chiang Mai and the
# factor levels are explicitly defined to ensure
# consistent ordering of groups in pairwise comparisons.

library(FSA)

Phylum_CM <- Phylum %>%
  mutate(
    Group = paste(Ethnicity, Area, sep = "-"),
    Group = factor(
      Group,
      levels = c("Akha-CM", "Lahu-CM", "Khuen-CM")
    )
  ) %>%
  relocate(Group, .before = 12) %>%
  filter(Area == "CM")


# --------------------------------------------------
# Dunn's Test for Significant Phyla
# --------------------------------------------------
# Perform pairwise Dunn's tests for each phylum
# identified as significant by the global
# Kruskal-Wallis test.
#
# The formula is constructed dynamically because
# the taxon names are stored as character strings.
#
# P-values from each Dunn's test are adjusted using
# the Benjamini-Hochberg (BH) method.

dunn_phylum_CM <- lapply(sig_krus_CM, function(taxon) {
  dunnTest(
    as.formula(paste(taxon, "~ Group")),
    data   = Phylum_CM,
    method = "bh"  # Benjamini-Hochberg correction
  )
})

# Assign taxon names to the corresponding result list elements
names(dunn_phylum_CM) <- sig_krus_CM
dunn_phylum_CM


# --------------------------------------------------
# Tidy Dunn's Test Results
# --------------------------------------------------
# Combine pairwise Dunn's test results for all
# significant phyla into a single data frame.
#
# Each row represents one pairwise comparison
# for one taxon.

dunn_phylum_tidy_CM <- do.call(
  rbind,
  lapply(names(dunn_phylum_CM), function(taxon) {
    df <- as.data.frame(dunn_phylum_CM[[taxon]]$res)
    df$Taxon <- taxon
    df
  })
)

dunn_phylum_tidy_CM

# Export the combined Dunn's test results
Export(dunn_phylum_tidy_CM, "dunn_phylum_tidy_CM.txt")


# --------------------------------------------------
# Extract Significant Pairwise Comparisons
# --------------------------------------------------
# Retain pairwise comparisons with BH-adjusted
# P-values < 0.05.

sig_dunn_phylum_CM <- dunn_phylum_tidy_CM %>% filter(P.adj < 0.05)
sig_dunn_phylum_CM


# --------------------------------------------------
# Alternative: Use Helper Function
# --------------------------------------------------
# Extract significant pairwise comparisons using
# the custom sig_dunnT() helper function.

sig_dunnT(dunn_phylum_CM)




# --------------------------------------------------
# Class-Level Analysis
# --------------------------------------------------

#------------
# [Class]
#------------

# Inspect current column names before renaming
colnames(Class)

# Rename taxonomic columns containing hyphens to ensure
# compatibility with formula-based statistical analyses
# and downstream modeling.
Class <- Class %>% 
  dplyr::rename(c__WPS_2 = `c__WPS-2`)


# --------------------------------------------------
# Global Kruskal-Wallis Test
# --------------------------------------------------
# Test for differences in class-level abundance among
# hill tribe groups within Chiang Mai (CM).
#
# The Kruskal-Wallis test is a non-parametric test used
# to evaluate whether the distribution of a feature
# differs among multiple independent groups.

stat_kruskal(
  Class, 
  rank_name   = "Class",
  group_col   = "Area",
  group_value = "CM",
  factor_col  = "Group",
  factor_levels = c("Akha-CM", "Lahu-CM", "Khuen-CM")
)


# --------------------------------------------------
# Post Hoc Analysis: Dunn's Test
# --------------------------------------------------
# Dunn's test is performed for classes showing a
# significant global Kruskal-Wallis test.
#
# Pairwise P-values are adjusted using the
# Benjamini-Hochberg (BH) method.


# Load Kruskal-Wallis summary results
kruskal_CM_Class <- read.delim(
  "~/Documents/HillTribe_NGS/3.Compositional_profile/Stat_microbiota/kruskal_CM_Class.txt"
)

# Extract classes with a significant global
# Kruskal-Wallis test (P < 0.05).
sig_krus_CM <- kruskal_CM_Class %>% 
  filter(P_sig == "Yes") %>% 
  pull(Taxa)


# --------------------------------------------------
# Prepare Data for Dunn's Test
# --------------------------------------------------
# Create a combined grouping variable representing
# ethnicity and geographic area.
#
# The data are restricted to Chiang Mai (CM), and
# factor levels are explicitly defined to control
# group ordering in pairwise comparisons.

library(FSA)

Class_CM <- Class %>%
  mutate(
    Group = paste(Ethnicity, Area, sep = "-"),
    Group = factor(
      Group,
      levels = c("Akha-CM", "Lahu-CM", "Lisu-CM", "Khuen-CM")
    )
  ) %>%
  relocate(Group, .before = 12) %>%
  filter(Area == "CM")


# --------------------------------------------------
# Run Dunn's Test for Significant Classes
# --------------------------------------------------
# Perform pairwise Dunn's tests for each class
# identified as significant by the global
# Kruskal-Wallis test.
#
# The formula is constructed dynamically because
# taxon names are stored as character strings.
#
# P-values are adjusted using the
# Benjamini-Hochberg (BH) method.

dunn_Class_CM <- lapply(sig_krus_CM, function(taxon) {
  dunnTest(
    as.formula(paste(taxon, "~ Group")),
    data   = Class_CM,
    method = "bh"  # Benjamini-Hochberg correction
  )
})

# Assign class names to the corresponding result
# list elements.
names(dunn_Class_CM) <- sig_krus_CM
dunn_Class_CM


# --------------------------------------------------
# Tidy Dunn's Test Results
# --------------------------------------------------
# Combine pairwise Dunn's test results for all
# significant classes into a single table.
#
# Each row represents one pairwise comparison
# for one taxon.

dunn_class_tidy_CM <- do.call(
  rbind,
  lapply(names(dunn_Class_CM), function(taxon) {
    df <- as.data.frame(dunn_Class_CM[[taxon]]$res)
    df$Taxon <- taxon
    df
  })
)

dunn_class_tidy_CM

# Export the combined Dunn's test results
Export(dunn_class_tidy_CM, "dunn_class_tidy_CM.txt")


# --------------------------------------------------
# Extract Significant Pairwise Comparisons
# --------------------------------------------------
# Retain pairwise comparisons with BH-adjusted
# P-values < 0.05.

sig_dunn_class_CM <- dunn_class_tidy_CM %>% filter(P.adj < 0.05)
sig_dunn_class_CM


# --------------------------------------------------
# Order-Level Analysis
# --------------------------------------------------

#------------
# [Order]
#------------

# Inspect current column names before renaming
colnames(Order)

# Rename taxonomic columns containing hyphens to ensure
# compatibility with formula-based statistical analyses
# and downstream modeling.
Order <- Order %>% 
  dplyr::rename(o__Clostridia_UCG_014 = `o__Clostridia_UCG-014`,
                o__WPS_2 = `o__WPS-2`,
                o__WCHB1_41 = `o__WCHB1-41`,
                o__Veillonellales_Selenomonadales = `o__Veillonellales-Selenomonadales`,
                o__Peptostreptococcales_Tissierellales = `o__Peptostreptococcales-Tissierellales`
                )


# --------------------------------------------------
# Global Kruskal-Wallis Test
# --------------------------------------------------
# Test for differences in order-level abundance among
# hill tribe groups within Chiang Mai (CM).

stat_kruskal(
  Order, 
  rank_name   = "Order",
  group_col   = "Area",
  group_value = "CM",
  factor_col  = "Group",
  factor_levels = c("Akha-CM", "Lahu-CM", "Lisu-CM", "Khuen-CM")
)


# --------------------------------------------------
# Post Hoc Analysis: Dunn's Test
# --------------------------------------------------
# Dunn's test is performed for orders showing a
# significant global Kruskal-Wallis test.
#
# Pairwise P-values are adjusted using the
# Benjamini-Hochberg (BH) method.


# Load Kruskal-Wallis summary results
kruskal_CM_Order <- read.delim(
  "~/Documents/HillTribe_NGS/3.Compositional_profile/Stat_microbiota/kruskal_CM_Order.txt"
)

# Extract orders with a significant global
# Kruskal-Wallis test (P < 0.05).
sig_krus_CM <- kruskal_CM_Order %>% 
  filter(P_sig == "Yes") %>% 
  pull(Taxa)


# --------------------------------------------------
# Prepare Data for Dunn's Test
# --------------------------------------------------
# Create a combined grouping variable representing
# ethnicity and geographic area.
#
# The data are restricted to Chiang Mai (CM), and
# factor levels are explicitly defined to control
# group ordering in pairwise comparisons.

library(FSA)

Order_CM <- Order %>%
  mutate(
    Group = paste(Ethnicity, Area, sep = "-"),
    Group = factor(
      Group,
      levels = c("Akha-CM", "Lahu-CM", "Khuen-CM")
    )
  ) %>%
  relocate(Group, .before = 12) %>%
  filter(Area == "CM")


# --------------------------------------------------
# Run Dunn's Test for Significant Orders
# --------------------------------------------------
# Perform pairwise Dunn's tests for each order
# identified as significant by the global
# Kruskal-Wallis test.
#
# The formula is constructed dynamically because
# taxon names are stored as character strings.
#
# P-values are adjusted using the
# Benjamini-Hochberg (BH) method.

dunn_Order_CM <- lapply(sig_krus_CM, function(taxon) {
  dunnTest(
    as.formula(paste(taxon, "~ Group")),
    data   = Order_CM,
    method = "bh"  # Benjamini-Hochberg correction
  )
})

# Assign order names to the corresponding result
# list elements.
names(dunn_Order_CM) <- sig_krus_CM
dunn_Order_CM


# --------------------------------------------------
# Tidy Dunn's Test Results
# --------------------------------------------------
# Combine pairwise Dunn's test results for all
# significant orders into a single table.
#
# Each row represents one pairwise comparison
# for one taxon.

dunn_Order_tidy_CM <- do.call(
  rbind,
  lapply(names(dunn_Order_CM), function(taxon) {
    df <- as.data.frame(dunn_Order_CM[[taxon]]$res)
    df$Taxon <- taxon
    df
  })
)

dunn_Order_tidy_CM

# Export the combined Dunn's test results
Export(dunn_Order_tidy_CM, "dunn_Order_tidy_CM.txt")


# --------------------------------------------------
# Extract Significant Pairwise Comparisons
# --------------------------------------------------
# Retain pairwise comparisons with BH-adjusted
# P-values < 0.05.

sig_dunn_Order_CM <- dunn_Order_tidy_CM %>% filter(P.adj < 0.05)
sig_dunn_Order_CM





# --------------------------------------------------
# Family-Level Analysis
# --------------------------------------------------

#------------
# [Family]
#------------

# Inspect current column names before renaming
colnames(Family)

# Rename taxonomic columns containing characters that
# may cause problems in formula-based statistical
# analyses and downstream modeling.
Family <- Family %>% 
  dplyr::rename(f__Clostridia_UCG_014 = `f__Clostridia_UCG-014`,
                f__UCG_010 = `f__UCG-010`,
                f__WPS_2 = `f__WPS-2`,
                f__UCG_011 = `f__UCG-011`,
                f__WCHB1_41 = `f__WCHB1-41`,
                f__Eubacterium_coprostanoligenes_group = `f__[Eubacterium]_coprostanoligenes_group`,
                f__Peptostreptococcales_Tissierellales = `f__Peptostreptococcales-Tissierellales`,
                f__Clostridium_methylpentosum_group = `f__[Clostridium]_methylpentosum_group`
  )


# --------------------------------------------------
# Global Kruskal-Wallis Test
# --------------------------------------------------
# Test for differences in family-level abundance among
# hill tribe groups within Chiang Mai (CM).
#
# The Kruskal-Wallis test is a non-parametric test used
# to evaluate whether the distribution of a feature
# differs among multiple independent groups.

stat_kruskal(
  Family, 
  rank_name   = "Family",
  group_col   = "Area",
  group_value = "CM",
  factor_col  = "Group",
  factor_levels = c("Akha-CM", "Lahu-CM", "Khuen-CM")
)


# --------------------------------------------------
# Post Hoc Analysis: Dunn's Test
# --------------------------------------------------
# Dunn's test is performed for families showing a
# significant global Kruskal-Wallis test.
#
# Pairwise P-values are adjusted using the
# Benjamini-Hochberg (BH) method.


# Load Kruskal-Wallis summary results
kruskal_CM_Family <- read.delim(
  "~/Documents/HillTribe_NGS/3.Compositional_profile/Stat_microbiota/kruskal_CM_Family.txt"
)

# Extract families with a significant global
# Kruskal-Wallis test (P < 0.05).
sig_krus_CM <- kruskal_CM_Family %>% 
  filter(P_sig == "Yes") %>% 
  pull(Taxa)


# --------------------------------------------------
# Prepare Data for Dunn's Test
# --------------------------------------------------
# Create a combined grouping variable representing
# ethnicity and geographic area.
#
# The data are restricted to Chiang Mai (CM), and
# factor levels are explicitly defined to control
# group ordering in pairwise comparisons.

library(FSA)

Family_CM <- Family %>%
  mutate(
    Group = paste(Ethnicity, Area, sep = "-"),
    Group = factor(
      Group,
      levels = c("Akha-CM", "Lahu-CM", "Khuen-CM")
    )
  ) %>%
  relocate(Group, .before = 12) %>%
  filter(Area == "CM")


# --------------------------------------------------
# Run Dunn's Test for Significant Families
# --------------------------------------------------
# Perform pairwise Dunn's tests for each family
# identified as significant by the global
# Kruskal-Wallis test.
#
# The formula is constructed dynamically because
# taxon names are stored as character strings.
#
# P-values are adjusted using the
# Benjamini-Hochberg (BH) method.

dunn_Family_CM <- lapply(sig_krus_CM, function(taxon) {
  dunnTest(
    as.formula(paste(taxon, "~ Group")),
    data   = Family_CM,
    method = "bh"  # Benjamini-Hochberg correction
  )
})

# Assign family names to the corresponding result
# list elements.
names(dunn_Family_CM) <- sig_krus_CM
dunn_Family_CM


# --------------------------------------------------
# Tidy Dunn's Test Results
# --------------------------------------------------
# Combine pairwise Dunn's test results for all
# significant families into a single table.
#
# Each row represents one pairwise comparison
# for one taxon.

dunn_Family_tidy_CM <- do.call(
  rbind,
  lapply(names(dunn_Family_CM), function(taxon) {
    df <- as.data.frame(dunn_Family_CM[[taxon]]$res)
    df$Taxon <- taxon
    df
  })
)

dunn_Family_tidy_CM

# Export the combined Dunn's test results
Export(dunn_Family_tidy_CM, "dunn_Family_tidy_CM.txt")


# --------------------------------------------------
# Extract Significant Pairwise Comparisons
# --------------------------------------------------
# Retain pairwise comparisons with BH-adjusted
# P-values < 0.05.

sig_dunn_Family_CM <- dunn_Family_tidy_CM %>% filter(P.adj < 0.05)
sig_dunn_Family_CM


# --------------------------------------------------
# Genus-Level Analysis
# --------------------------------------------------

#------------
# [Genus]
#------------

# Inspect current column names before cleaning
colnames(Genus)


# --------------------------------------------------
# Helper Function: Taxon Name Cleaning
# --------------------------------------------------
# Clean taxonomic column names by removing square
# brackets and replacing hyphens with underscores.
#
# This standardizes taxon names and facilitates their
# use in formula-based statistical analyses.

prep_tax_table <- function(df) {
  
  library(dplyr)
  library(stringr)
  
  # Clean column names: remove brackets and replace hyphen with underscore
  colnames(df) <- colnames(df) %>%
    str_replace_all("\\[|\\]", "") %>% 
    str_replace_all("-", "_")
  
  return(df)
}

#-------------------------------------------------------------------

# Apply taxonomic name cleaning to the genus-level table
Genus <- prep_tax_table(Genus)

# Inspect cleaned column names
colnames(Genus)


# --------------------------------------------------
# Global Kruskal-Wallis Test
# --------------------------------------------------
# Test for differences in genus-level abundance among
# hill tribe groups within Chiang Mai (CM).

stat_kruskal(
  Genus, 
  rank_name   = "Genus",
  group_col   = "Area",
  group_value = "CM",
  factor_col  = "Group",
  factor_levels = c("Akha-CM", "Lahu-CM", "Khuen-CM")
)


# --------------------------------------------------
# Post Hoc Analysis: Dunn's Test
# --------------------------------------------------
# Dunn's test is performed for genera showing a
# significant global Kruskal-Wallis test.
#
# Pairwise P-values are adjusted using the
# Benjamini-Hochberg (BH) method.


# Load Kruskal-Wallis summary results
kruskal_CM_Genus <- read.delim(
  "~/Documents/HillTribe_NGS/3.Compositional_profile/Stat_microbiota/kruskal_CM_Genus.txt"
)

# Extract genera with a significant global
# Kruskal-Wallis test (P < 0.05).
sig_krus_CM <- kruskal_CM_Genus %>% 
  filter(P_sig == "Yes") %>% 
  pull(Taxa)


# --------------------------------------------------
# Prepare Data for Dunn's Test
# --------------------------------------------------
# Create a combined grouping variable representing
# ethnicity and geographic area.
#
# The data are restricted to Chiang Mai (CM), and
# factor levels are explicitly defined to control
# group ordering in pairwise comparisons.

library(FSA)

Genus_CM <- Genus %>%
  mutate(
    Group = paste(Ethnicity, Area, sep = "-"),
    Group = factor(
      Group,
      levels = c("Akha-CM", "Lahu-CM", "Khuen-CM")
    )
  ) %>%
  relocate(Group, .before = 12) %>%
  filter(Area == "CM")


# --------------------------------------------------
# Run Dunn's Test for Significant Genera
# --------------------------------------------------
# Perform pairwise Dunn's tests for each genus
# identified as significant by the global
# Kruskal-Wallis test.
#
# The formula is constructed dynamically because
# taxon names are stored as character strings.
#
# P-values are adjusted using the
# Benjamini-Hochberg (BH) method.

dunn_Genus_CM <- lapply(sig_krus_CM, function(taxon) {
  dunnTest(
    as.formula(paste(taxon, "~ Group")),
    data   = Genus_CM,
    method = "bh"  # Benjamini-Hochberg correction
  )
})

# Assign genus names to the corresponding result
# list elements.
names(dunn_Genus_CM) <- sig_krus_CM
dunn_Genus_CM


# --------------------------------------------------
# Tidy Dunn's Test Results
# --------------------------------------------------
# Combine pairwise Dunn's test results for all
# significant genera into a single table.
#
# Each row represents one pairwise comparison
# for one taxon.

dunn_Genus_tidy_CM <- do.call(
  rbind,
  lapply(names(dunn_Genus_CM), function(taxon) {
    df <- as.data.frame(dunn_Genus_CM[[taxon]]$res)
    df$Taxon <- taxon
    df
  })
)

dunn_Genus_tidy_CM

# Export the combined Dunn's test results
Export(dunn_Genus_tidy_CM, "dunn_Genus_tidy_CM.txt")


# --------------------------------------------------
# Extract Significant Pairwise Comparisons
# --------------------------------------------------
# Retain pairwise comparisons with BH-adjusted
# P-values < 0.05.

sig_dunn_Genus_CM <- dunn_Genus_tidy_CM %>% filter(P.adj < 0.05)
sig_dunn_Genus_CM









# --------------------------------------------------
# Kruskal-Wallis Test: Phylum Level (Chiang Rai)
# --------------------------------------------------
#
# Objective:
# Identify phyla that differ in abundance among hill
# tribe groups within Chiang Rai (CR).
#
# Analysis workflow:
# 1. Perform a global Kruskal-Wallis test for each phylum.
# 2. Identify phyla with significant global differences.
# 3. Perform pairwise Dunn's tests for significant phyla.
# 4. Adjust pairwise P-values using the
#    Benjamini-Hochberg (BH) method.
# 5. Extract significant pairwise comparisons.
# --------------------------------------------------


#------------
# [Phylum]
#------------


# --------------------------------------------------
# Global Kruskal-Wallis Test
# --------------------------------------------------
# Test for differences in phylum-level abundance among
# hill tribe groups within Chiang Rai (CR).
#
# The Kruskal-Wallis test is a non-parametric test used
# to evaluate whether the distribution of a feature
# differs among multiple independent groups.

stat_kruskal(
  Phylum, 
  rank_name   = "Phylum",
  group_col   = "Area",
  group_value = "CR",
  factor_col  = "Group",
  factor_levels = c("Akha-CR", "Lahu-CR", "Lisu-CR")
)


# --------------------------------------------------
# Post Hoc Analysis: Dunn's Test
# --------------------------------------------------
# Dunn's test is performed for phyla showing a
# significant global Kruskal-Wallis test.
#
# Pairwise P-values are adjusted using the
# Benjamini-Hochberg (BH) method.


# Load Kruskal-Wallis summary results
kruskal_CR_Phylum <- read.delim(
  "~/Documents/HillTribe_NGS/3.Compositional_profile/Stat_microbiota/kruskal_CR_Phylum.txt"
)

# Extract phyla with a significant global
# Kruskal-Wallis test (P < 0.05).
sig_krus_CR <- kruskal_CR_Phylum %>% 
  filter(P_sig == "Yes") %>% 
  pull(Taxa)


# --------------------------------------------------
# Prepare Data for Dunn's Test
# --------------------------------------------------
# Create a combined grouping variable representing
# ethnicity and geographic area.
#
# The data are restricted to Chiang Rai (CR), and
# factor levels are explicitly defined to control
# group ordering in pairwise comparisons.

library(FSA)

Phylum_CR <- Phylum %>%
  mutate(
    Group = paste(Ethnicity, Area, sep = "-"),
    Group = factor(
      Group,
      levels = c("Akha-CR", "Lahu-CR", "Lisu-CR")
    )
  ) %>%
  relocate(Group, .before = 12) %>%
  filter(Area == "CR")


# --------------------------------------------------
# Run Dunn's Test for Significant Phyla
# --------------------------------------------------
# Perform pairwise Dunn's tests for each phylum
# identified as significant by the global
# Kruskal-Wallis test.
#
# The formula is constructed dynamically because
# taxon names are stored as character strings.
#
# P-values are adjusted using the
# Benjamini-Hochberg (BH) method.

dunn_phylum_CR <- lapply(sig_krus_CR, function(taxon) {
  dunnTest(
    as.formula(paste(taxon, "~ Group")),
    data   = Phylum_CR,
    method = "bh"  # Benjamini-Hochberg correction
  )
})

# Assign phylum names to the corresponding result
# list elements.
names(dunn_phylum_CR) <- sig_krus_CR
dunn_phylum_CR


# --------------------------------------------------
# Tidy Dunn's Test Results
# --------------------------------------------------
# Combine pairwise Dunn's test results for all
# significant phyla into a single table.
#
# Each row represents one pairwise comparison
# for one taxon.

dunn_phylum_tidy_CR <- do.call(
  rbind,
  lapply(names(dunn_phylum_CR), function(taxon) {
    df <- as.data.frame(dunn_phylum_CR[[taxon]]$res)
    df$Taxon <- taxon
    df
  })
)

dunn_phylum_tidy_CR

# Export the combined Dunn's test results
Export(dunn_phylum_tidy_CR, "dunn_phylum_tidy_CR.txt")


# --------------------------------------------------
# Extract Significant Pairwise Comparisons
# --------------------------------------------------
# Retain pairwise comparisons with BH-adjusted
# P-values < 0.05.

sig_dunn_phylum_CR <- dunn_phylum_tidy_CR %>% filter(P.adj < 0.05)
sig_dunn_phylum_CR


# Alternative: Use the helper function to extract
# significant pairwise comparisons.
sig_dunnT(dunn_phylum_CR)



# --------------------------------------------------
# Class-Level Analysis
# --------------------------------------------------

#------------
# [Class]
#------------


# --------------------------------------------------
# Global Kruskal-Wallis Test
# --------------------------------------------------
# Test for differences in class-level abundance among
# hill tribe groups within Chiang Rai (CR).

stat_kruskal(
  Class, 
  rank_name   = "Class",
  group_col   = "Area",
  group_value = "CR",
  factor_col  = "Group",
  factor_levels = c("Akha-CR", "Lahu-CR", "Lisu-CR")
)


# --------------------------------------------------
# Post Hoc Analysis: Dunn's Test
# --------------------------------------------------
# Dunn's test is performed for classes showing a
# significant global Kruskal-Wallis test.
#
# Pairwise P-values are adjusted using the
# Benjamini-Hochberg (BH) method.


# Load Kruskal-Wallis summary results
kruskal_CR_Class <- read.delim(
  "~/Documents/HillTribe_NGS/3.Compositional_profile/Stat_microbiota/kruskal_CR_Class.txt"
)

# Extract classes with a significant global
# Kruskal-Wallis test (P < 0.05).
sig_krus_CR <- kruskal_CR_Class %>% 
  filter(P_sig == "Yes") %>% 
  pull(Taxa)


# --------------------------------------------------
# Prepare Data for Dunn's Test
# --------------------------------------------------
# Create a combined grouping variable representing
# ethnicity and geographic area.
#
# The data are restricted to Chiang Rai (CR), and
# factor levels are explicitly defined to control
# group ordering in pairwise comparisons.

library(FSA)

Class_CR <- Class %>%
  mutate(
    Group = paste(Ethnicity, Area, sep = "-"),
    Group = factor(
      Group,
      levels = c("Akha-CR", "Lahu-CR", "Lisu-CR")
    )
  ) %>%
  relocate(Group, .before = 12) %>%
  filter(Area == "CR")


# --------------------------------------------------
# Run Dunn's Test for Significant Classes
# --------------------------------------------------
# Perform pairwise Dunn's tests for each class
# identified as significant by the global
# Kruskal-Wallis test.
#
# The formula is constructed dynamically because
# taxon names are stored as character strings.
#
# P-values are adjusted using the
# Benjamini-Hochberg (BH) method.

dunn_Class_CR <- lapply(sig_krus_CR, function(taxon) {
  dunnTest(
    as.formula(paste(taxon, "~ Group")),
    data   = Class_CR,
    method = "bh"  # Benjamini-Hochberg correction
  )
})

# Assign class names to the corresponding result
# list elements.
names(dunn_Class_CR) <- sig_krus_CR
dunn_Class_CR


# --------------------------------------------------
# Tidy Dunn's Test Results
# --------------------------------------------------
# Combine pairwise Dunn's test results for all
# significant classes into a single table.
#
# Each row represents one pairwise comparison
# for one taxon.

dunn_class_tidy_CR <- do.call(
  rbind,
  lapply(names(dunn_Class_CR), function(taxon) {
    df <- as.data.frame(dunn_Class_CR[[taxon]]$res)
    df$Taxon <- taxon
    df
  })
)

dunn_class_tidy_CR

# Export the combined Dunn's test results
Export(dunn_class_tidy_CR, "dunn_class_tidy_CR.txt")


# --------------------------------------------------
# Extract Significant Pairwise Comparisons
# --------------------------------------------------
# Retain pairwise comparisons with BH-adjusted
# P-values < 0.05.

sig_dunn_class_CR <- dunn_class_tidy_CR %>% filter(P.adj < 0.05)
sig_dunn_class_CR




# --------------------------------------------------
# Order-Level Analysis
# --------------------------------------------------

#------------
# [Order]
#------------


# --------------------------------------------------
# Global Kruskal-Wallis Test
# --------------------------------------------------
# Test for differences in order-level abundance among
# hill tribe groups within Chiang Rai (CR).

stat_kruskal(
  Order, 
  rank_name   = "Order",
  group_col   = "Area",
  group_value = "CR",
  factor_col  = "Group",
  factor_levels = c("Akha-CR", "Lahu-CR", "Lisu-CR")
)


# --------------------------------------------------
# Post Hoc Analysis: Dunn's Test
# --------------------------------------------------
# Dunn's test is performed for orders showing a
# significant global Kruskal-Wallis test.
#
# Pairwise P-values are adjusted using the
# Benjamini-Hochberg (BH) method.


# Load Kruskal-Wallis summary results
kruskal_CR_Order <- read.delim(
  "~/Documents/HillTribe_NGS/3.Compositional_profile/Stat_microbiota/kruskal_CR_Order.txt"
)

# Extract orders with a significant global
# Kruskal-Wallis test (P < 0.05).
sig_krus_CR <- kruskal_CR_Order %>% 
  filter(P_sig == "Yes") %>% 
  pull(Taxa)


# --------------------------------------------------
# Prepare Data for Dunn's Test
# --------------------------------------------------
# Create a combined grouping variable representing
# ethnicity and geographic area.
#
# The data are restricted to Chiang Rai (CR), and
# factor levels are explicitly defined to control
# group ordering in pairwise comparisons.

library(FSA)

Order_CR <- Order %>%
  mutate(
    Group = paste(Ethnicity, Area, sep = "-"),
    Group = factor(
      Group,
      levels = c("Akha-CR", "Lahu-CR", "Lisu-CR")
    )
  ) %>%
  relocate(Group, .before = 12) %>%
  filter(Area == "CR")


# --------------------------------------------------
# Run Dunn's Test for Significant Orders
# --------------------------------------------------
# Perform pairwise Dunn's tests for each order
# identified as significant by the global
# Kruskal-Wallis test.
#
# The formula is constructed dynamically because
# taxon names are stored as character strings.
#
# P-values are adjusted using the
# Benjamini-Hochberg (BH) method.

dunn_Order_CR <- lapply(sig_krus_CR, function(taxon) {
  dunnTest(
    as.formula(paste(taxon, "~ Group")),
    data   = Order_CR,
    method = "bh"  # Benjamini-Hochberg correction
  )
})

# Assign order names to the corresponding result
# list elements.
names(dunn_Order_CR) <- sig_krus_CR
dunn_Order_CR


# --------------------------------------------------
# Tidy Dunn's Test Results
# --------------------------------------------------
# Combine pairwise Dunn's test results for all
# significant orders into a single table.
#
# Each row represents one pairwise comparison
# for one taxon.

dunn_Order_tidy_CR <- do.call(
  rbind,
  lapply(names(dunn_Order_CR), function(taxon) {
    df <- as.data.frame(dunn_Order_CR[[taxon]]$res)
    df$Taxon <- taxon
    df
  })
)

dunn_Order_tidy_CR

# Export the combined Dunn's test results
Export(dunn_Order_tidy_CR, "dunn_Order_tidy_CR.txt")


# --------------------------------------------------
# Extract Significant Pairwise Comparisons
# --------------------------------------------------
# Retain pairwise comparisons with BH-adjusted
# P-values < 0.05.

sig_dunn_Order_CR <- dunn_Order_tidy_CR %>% filter(P.adj < 0.05)
sig_dunn_Order_CR



# --------------------------------------------------
# Family-Level Analysis
# --------------------------------------------------

#------------
# [Family]
#------------


# --------------------------------------------------
# Global Kruskal-Wallis Test
# --------------------------------------------------
# Test for differences in family-level abundance among
# hill tribe groups within Chiang Rai (CR).

stat_kruskal(
  Family, 
  rank_name   = "Family",
  group_col   = "Area",
  group_value = "CR",
  factor_col  = "Group",
  factor_levels = c("Akha-CR", "Lahu-CR", "Lisu-CR")
)


# --------------------------------------------------
# Post Hoc Analysis: Dunn's Test
# --------------------------------------------------
# Dunn's test is performed for families showing a
# significant global Kruskal-Wallis test.
#
# Pairwise P-values are adjusted using the
# Benjamini-Hochberg (BH) method.


# Load Kruskal-Wallis summary results
kruskal_CR_Family <- read.delim(
  "~/Documents/HillTribe_NGS/3.Compositional_profile/Stat_microbiota/kruskal_CR_Family.txt"
)

# Extract families with a significant global
# Kruskal-Wallis test (P < 0.05).
sig_krus_CR <- kruskal_CR_Family %>% 
  filter(P_sig == "Yes") %>% 
  pull(Taxa)


# --------------------------------------------------
# Prepare Data for Dunn's Test
# --------------------------------------------------
# Create a combined grouping variable representing
# ethnicity and geographic area.
#
# The data are restricted to Chiang Rai (CR), and
# factor levels are explicitly defined to control
# group ordering in pairwise comparisons.

library(FSA)

Family_CR <- Family %>%
  mutate(
    Group = paste(Ethnicity, Area, sep = "-"),
    Group = factor(
      Group,
      levels = c("Akha-CR", "Lahu-CR", "Lisu-CR")
    )
  ) %>%
  relocate(Group, .before = 12) %>%
  filter(Area == "CR")


# --------------------------------------------------
# Run Dunn's Test for Significant Families
# --------------------------------------------------
# Perform pairwise Dunn's tests for each family
# identified as significant by the global
# Kruskal-Wallis test.
#
# The formula is constructed dynamically because
# taxon names are stored as character strings.
#
# P-values are adjusted using the
# Benjamini-Hochberg (BH) method.

dunn_Family_CR <- lapply(sig_krus_CR, function(taxon) {
  dunnTest(
    as.formula(paste(taxon, "~ Group")),
    data   = Family_CR,
    method = "bh"  # Benjamini-Hochberg correction
  )
})

# Assign family names to the corresponding result
# list elements.
names(dunn_Family_CR) <- sig_krus_CR
dunn_Family_CR


# --------------------------------------------------
# Tidy Dunn's Test Results
# --------------------------------------------------
# Combine pairwise Dunn's test results for all
# significant families into a single table.
#
# Each row represents one pairwise comparison
# for one taxon.

dunn_Family_tidy_CR <- do.call(
  rbind,
  lapply(names(dunn_Family_CR), function(taxon) {
    df <- as.data.frame(dunn_Family_CR[[taxon]]$res)
    df$Taxon <- taxon
    df
  })
)

dunn_Family_tidy_CR

# Export the combined Dunn's test results
Export(dunn_Family_tidy_CR, "dunn_Family_tidy_CR.txt")


# --------------------------------------------------
# Extract Significant Pairwise Comparisons
# --------------------------------------------------
# Retain pairwise comparisons with BH-adjusted
# P-values < 0.05.

sig_dunn_Family_CR <- dunn_Family_tidy_CR %>% filter(P.adj < 0.05)
sig_dunn_Family_CR



# --------------------------------------------------
# Genus-Level Analysis
# --------------------------------------------------

#------------
# [Genus]
#------------


# --------------------------------------------------
# Global Kruskal-Wallis Test
# --------------------------------------------------
# Test for differences in genus-level abundance among
# hill tribe groups within Chiang Rai (CR).

stat_kruskal(
  Genus, 
  rank_name   = "Genus",
  group_col   = "Area",
  group_value = "CR",
  factor_col  = "Group",
  factor_levels = c("Akha-CR", "Lahu-CR", "Lisu-CR")
)


# --------------------------------------------------
# Post Hoc Analysis: Dunn's Test
# --------------------------------------------------
# Dunn's test is performed for genera showing a
# significant global Kruskal-Wallis test.
#
# Pairwise P-values are adjusted using the
# Benjamini-Hochberg (BH) method.


# Load Kruskal-Wallis summary results
kruskal_CR_Genus <- read.delim(
  "~/Documents/HillTribe_NGS/3.Compositional_profile/Stat_microbiota/kruskal_CR_Genus.txt"
)

# Extract genera with a significant global
# Kruskal-Wallis test (P < 0.05).
sig_krus_CR <- kruskal_CR_Genus %>% 
  filter(P_sig == "Yes") %>% 
  pull(Taxa)


# --------------------------------------------------
# Prepare Data for Dunn's Test
# --------------------------------------------------
# Create a combined grouping variable representing
# ethnicity and geographic area.
#
# The data are restricted to Chiang Rai (CR), and
# factor levels are explicitly defined to control
# group ordering in pairwise comparisons.

library(FSA)

Genus_CR <- Genus %>%
  mutate(
    Group = paste(Ethnicity, Area, sep = "-"),
    Group = factor(
      Group,
      levels = c("Akha-CR", "Lahu-CR", "Lisu-CR")
    )
  ) %>%
  relocate(Group, .before = 12) %>%
  filter(Area == "CR")


# --------------------------------------------------
# Run Dunn's Test for Significant Genera
# --------------------------------------------------
# Perform pairwise Dunn's tests for each genus
# identified as significant by the global
# Kruskal-Wallis test.
#
# The formula is constructed dynamically because
# taxon names are stored as character strings.
#
# P-values are adjusted using the
# Benjamini-Hochberg (BH) method.

dunn_Genus_CR <- lapply(sig_krus_CR, function(taxon) {
  dunnTest(
    as.formula(paste(taxon, "~ Group")),
    data   = Genus_CR,
    method = "bh"  # Benjamini-Hochberg correction
  )
})

# Assign genus names to the corresponding result
# list elements.
names(dunn_Genus_CR) <- sig_krus_CR
dunn_Genus_CR


# --------------------------------------------------
# Tidy Dunn's Test Results
# --------------------------------------------------
# Combine pairwise Dunn's test results for all
# significant genera into a single table.
#
# Each row represents one pairwise comparison
# for one taxon.

dunn_Genus_tidy_CR <- do.call(
  rbind,
  lapply(names(dunn_Genus_CR), function(taxon) {
    df <- as.data.frame(dunn_Genus_CR[[taxon]]$res)
    df$Taxon <- taxon
    df
  })
)

dunn_Genus_tidy_CR

# Export the combined Dunn's test results
Export(dunn_Genus_tidy_CR, "dunn_Genus_tidy_CR.txt")


# --------------------------------------------------
# Extract Significant Pairwise Comparisons
# --------------------------------------------------
# Retain pairwise comparisons with BH-adjusted
# P-values < 0.05.

sig_dunn_Genus_CR <- dunn_Genus_tidy_CR %>% filter(P.adj < 0.05)
sig_dunn_Genus_CR








# --------------------------------------------------
# Abundance Profile Summary
# --------------------------------------------------
#
# Objective:
# Generate descriptive summaries of taxonomic abundance
# for Akha participants, stratified by geographic area.
#
# Taxonomic levels:
#   - Phylum
#   - Class
#   - Order
#   - Family
#   - Genus
#
# Summary statistics:
#   - Mean abundance
#   - Standard deviation (SD)
#   - Formatted mean ± SD
#
# Missing values are excluded when calculating the
# mean and standard deviation.
# --------------------------------------------------

setwd("~/Documents/HillTribe_NGS/3.Compositional_profile/Stat_microbiota/abundance_summary")

# Load tidyr for reshaping the summary tables
library(tidyr)


# --------------------------------------------------
# Akha: Phylum-Level Abundance
# --------------------------------------------------
# Summarize phylum-level abundance among Akha
# participants by geographic area.

sum_Phylum_Akha <-  Phylum %>%
  filter(Ethnicity == "Akha") %>%
  select(-Height:-DBP) %>%
  group_by(Area) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Area,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Phylum_Akha
Export(sum_Phylum_Akha, "sum_Phylum_Akha.txt")


# --------------------------------------------------
# Akha: Class-Level Abundance
# --------------------------------------------------
# Summarize class-level abundance among Akha
# participants by geographic area.

sum_Class_Akha <-  Class %>%
  filter(Ethnicity == "Akha") %>%
  select(-Height:-DBP) %>%
  group_by(Area) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Area,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Class_Akha
Export(sum_Class_Akha, "sum_Class_Akha.txt")


# --------------------------------------------------
# Akha: Order-Level Abundance
# --------------------------------------------------
# Summarize order-level abundance among Akha
# participants by geographic area.

sum_Order_Akha <-  Order %>%
  filter(Ethnicity == "Akha") %>%
  select(-Height:-DBP) %>%
  group_by(Area) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Area,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Order_Akha
Export(sum_Order_Akha, "sum_Order_Akha.txt")


# --------------------------------------------------
# Akha: Family-Level Abundance
# --------------------------------------------------
# Summarize family-level abundance among Akha
# participants by geographic area.

sum_Family_Akha <-  Family %>%
  filter(Ethnicity == "Akha") %>%
  select(-Height:-DBP) %>%
  group_by(Area) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Area,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Family_Akha
Export(sum_Family_Akha, "sum_Family_Akha.txt")


# --------------------------------------------------
# Akha: Genus-Level Abundance
# --------------------------------------------------
# Summarize genus-level abundance among Akha
# participants by geographic area.

sum_Genus_Akha <-  Genus %>%
  filter(Ethnicity == "Akha") %>%
  select(-Height:-DBP) %>%
  group_by(Area) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Area,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Genus_Akha
Export(sum_Genus_Akha, "sum_Genus_Akha.txt")







# --------------------------------------------------
# Abundance Profile Summary: Lahu
# --------------------------------------------------
#
# Generate descriptive summaries of taxonomic abundance
# for Lahu participants, stratified by geographic area.
#
# Summary statistics:
#   - Mean abundance
#   - Standard deviation (SD)
#   - Formatted mean ± SD
#
# Missing values are excluded when calculating the
# mean and standard deviation.
# --------------------------------------------------


# --------------------------------------------------
# Lahu: Phylum-Level Abundance
# --------------------------------------------------

sum_Phylum_Lahu <-  Phylum %>%
  filter(Ethnicity == "Lahu") %>%
  select(-Height:-DBP) %>%
  group_by(Area) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Area,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Phylum_Lahu
Export(sum_Phylum_Lahu, "sum_Phylum_Lahu.txt")


# --------------------------------------------------
# Lahu: Class-Level Abundance
# --------------------------------------------------

sum_Class_Lahu <-  Class %>%
  filter(Ethnicity == "Lahu") %>%
  select(-Height:-DBP) %>%
  group_by(Area) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Area,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Class_Lahu
Export(sum_Class_Lahu, "sum_Class_Lahu.txt")


# --------------------------------------------------
# Lahu: Order-Level Abundance
# --------------------------------------------------

sum_Order_Lahu <-  Order %>%
  filter(Ethnicity == "Lahu") %>%
  select(-Height:-DBP) %>%
  group_by(Area) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Area,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Order_Lahu
Export(sum_Order_Lahu, "sum_Order_Lahu.txt")


# --------------------------------------------------
# Lahu: Family-Level Abundance
# --------------------------------------------------

sum_Family_Lahu <-  Family %>%
  filter(Ethnicity == "Lahu") %>%
  select(-Height:-DBP) %>%
  group_by(Area) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Area,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Family_Lahu
Export(sum_Family_Lahu, "sum_Family_Lahu.txt")


# --------------------------------------------------
# Lahu: Genus-Level Abundance
# --------------------------------------------------

sum_Genus_Lahu <-  Genus %>%
  filter(Ethnicity == "Lahu") %>%
  select(-Height:-DBP) %>%
  group_by(Area) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Area,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Genus_Lahu
Export(sum_Genus_Lahu, "sum_Genus_Lahu.txt")








# --------------------------------------------------
# Abundance Profile Summary: Chiang Mai (CM)
# --------------------------------------------------
#
# Objective:
# Generate descriptive summaries of taxonomic abundance
# among ethnic groups within Chiang Mai (CM).
#
# The data are filtered to Chiang Mai and summarized
# separately for each ethnicity.
#
# Taxonomic levels:
#   - Phylum
#   - Class
#   - Order
#   - Family
#   - Genus
#
# Summary statistics:
#   - Mean abundance
#   - Standard deviation (SD)
#   - Formatted mean ± SD
#
# Missing values are excluded when calculating the
# mean and standard deviation.
# --------------------------------------------------


# --------------------------------------------------
# Chiang Mai: Phylum-Level Abundance
# --------------------------------------------------
# Summarize phylum-level abundance across ethnic groups
# within Chiang Mai.

sum_Phylum_CM <-  Phylum %>%
  filter(Area == "CM") %>%
  select(-Height:-DBP) %>%
  group_by(Ethnicity) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Ethnicity,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Phylum_CM
Export(sum_Phylum_CM, "sum_Phylum_CM.txt")


# --------------------------------------------------
# Chiang Mai: Class-Level Abundance
# --------------------------------------------------
# Summarize class-level abundance across ethnic groups
# within Chiang Mai.

sum_Class_CM <-  Class %>%
  filter(Area == "CM") %>%
  select(-Height:-DBP) %>%
  group_by(Ethnicity) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Ethnicity,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Class_CM
Export(sum_Class_CM, "sum_Class_CM.txt")


# --------------------------------------------------
# Chiang Mai: Order-Level Abundance
# --------------------------------------------------
# Summarize order-level abundance across ethnic groups
# within Chiang Mai.

sum_Order_CM <-  Order %>%
  filter(Area == "CM") %>%
  select(-Height:-DBP) %>%
  group_by(Ethnicity) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Ethnicity,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Order_CM
Export(sum_Order_CM, "sum_Order_CM.txt")


# --------------------------------------------------
# Chiang Mai: Family-Level Abundance
# --------------------------------------------------
# Summarize family-level abundance across ethnic groups
# within Chiang Mai.

sum_Family_CM <-  Family %>%
  filter(Area == "CM") %>%
  select(-Height:-DBP) %>%
  group_by(Ethnicity) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Ethnicity,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Family_CM
Export(sum_Family_CM, "sum_Family_CM.txt")


# --------------------------------------------------
# Chiang Mai: Genus-Level Abundance
# --------------------------------------------------
# Summarize genus-level abundance across ethnic groups
# within Chiang Mai.

sum_Genus_CM <-  Genus %>%
  filter(Area == "CM") %>%
  select(-Height:-DBP) %>%
  group_by(Ethnicity) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Ethnicity,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Genus_CM
Export(sum_Genus_CM, "sum_Genus_CM.txt")







# --------------------------------------------------
# Abundance Profile Summary: Chiang Rai (CR)
# --------------------------------------------------
#
# Objective:
# Generate descriptive summaries of taxonomic abundance
# among ethnic groups within Chiang Rai (CR).
#
# The data are filtered to Chiang Rai and summarized
# separately for each ethnicity.
#
# Taxonomic levels:
#   - Phylum
#   - Class
#   - Order
#   - Family
#   - Genus
#
# Summary statistics:
#   - Mean abundance
#   - Standard deviation (SD)
#   - Formatted mean ± SD
#
# Missing values are excluded when calculating the
# mean and standard deviation.
# --------------------------------------------------


# --------------------------------------------------
# Chiang Rai: Phylum-Level Abundance
# --------------------------------------------------
# Summarize phylum-level abundance across ethnic groups
# within Chiang Rai.

sum_Phylum_CR <-  Phylum %>%
  filter(Area == "CR") %>%
  select(-Height:-DBP) %>%
  group_by(Ethnicity) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Ethnicity,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Phylum_CR
Export(sum_Phylum_CR, "sum_Phylum_CR.txt")


# --------------------------------------------------
# Chiang Rai: Class-Level Abundance
# --------------------------------------------------
# Summarize class-level abundance across ethnic groups
# within Chiang Rai.

sum_Class_CR <-  Class %>%
  filter(Area == "CR") %>%
  select(-Height:-DBP) %>%
  group_by(Ethnicity) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Ethnicity,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Class_CR
Export(sum_Class_CR, "sum_Class_CR.txt")


# --------------------------------------------------
# Chiang Rai: Order-Level Abundance
# --------------------------------------------------
# Summarize order-level abundance across ethnic groups
# within Chiang Rai.

sum_Order_CR <-  Order %>%
  filter(Area == "CR") %>%
  select(-Height:-DBP) %>%
  group_by(Ethnicity) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Ethnicity,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Order_CR
Export(sum_Order_CR, "sum_Order_CR.txt")


# --------------------------------------------------
# Chiang Rai: Family-Level Abundance
# --------------------------------------------------
# Summarize family-level abundance across ethnic groups
# within Chiang Rai.

sum_Family_CR <-  Family %>%
  filter(Area == "CR") %>%
  select(-Height:-DBP) %>%
  group_by(Ethnicity) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Ethnicity,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Family_CR
Export(sum_Family_CR, "sum_Family_CR.txt")


# --------------------------------------------------
# Chiang Rai: Genus-Level Abundance
# --------------------------------------------------
# Summarize genus-level abundance across ethnic groups
# within Chiang Rai.

sum_Genus_CR <-  Genus %>%
  filter(Area == "CR") %>%
  select(-Height:-DBP) %>%
  group_by(Ethnicity) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Ethnicity,
    names_to = c("Data", "Statistic"),
    names_pattern = "(.+)_(mean|sd)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    `mean ± sd` = sprintf("%.2f ± %.2f", mean, sd)
  )

sum_Genus_CR
Export(sum_Genus_CR, "sum_Genus_CR.txt")
