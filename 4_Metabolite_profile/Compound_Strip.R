############################################################################
# Metabolite Name Cleaning Pipeline
#
# Purpose:
# Standardize metabolite/compound names prior to compound identification
# and downstream mapping. This pipeline removes GC-MS derivatization
# annotations, excludes known analytical contaminants, and prepares the
# cleaned compound list for external compound mapping.
#
# Workflow:
#   1. Load required packages and input compound list
#   2. Remove derivatization annotations (TMS/TBDMS)
#   3. Remove known GC contaminants
#   4. Export the cleaned compound list
#   5. Split the compound list into smaller files for external mapping
############################################################################


setwd("~/Documents/HillTribe_NGS/5.Metabolites/Compound_Mapping")


# ==============================
# 1. Environment Setup
# ==============================
# Load packages required for data manipulation, file import/export,
# statistical functions, and string processing.

library(dplyr)
library(rio)
library(car)  
library(stringr)


# ==========================================================
# 2. Compound List 0.05 Processing
# ==========================================================
# Load the compound list generated using the 0.05 significance threshold.
# The input file contains compound names without a predefined header.

compound_list_05 <- read.delim(
  "~/Documents/HillTribe_NGS/5.Metabolites/Compound_Mapping/compound_list_05.txt",
  header = FALSE
)

# Rename the input column for easier reference.

compound_list_05 <- compound_list_05 %>% 
  rename(Compound = V1)



# ==========================================================
# 3. Remove Derivatization Annotations
# ==========================================================
#
# GC-MS metabolite identification may report derivatized compounds
# with annotations such as TMS or TBDMS. These annotations represent
# analytical derivatization rather than the underlying compound name.
#
# Method 1:
# Remove common TMS/TBDMS derivatization annotations and standardize
# whitespace.

compound_list_05_clean <- compound_list_05 %>%
  mutate(
    Compound_clean = str_remove_all(
      Compound,
      ",\\s*\\d*TMS( derivative)?|\\b\\d+TMS\\b|\\s*TMS derivative|,\\s*\\d+TBDMS\\b"
    ) %>%
      str_squish()
  ) 

Export(compound_list_05_clean, "compound_list_05_clean.txt")


# Method 2:
# Apply a broader case-insensitive cleaning pattern to remove TMS/TBDMS
# derivatization and ester annotations.
#
# Known GC contaminants are subsequently removed from the compound list
# because they are not considered biological metabolites of interest.

df_clean <- compound_list_05 %>%
  mutate(
    Compound_clean = Compound %>%
      str_remove_all("(?i)\\s*,?\\s*\\d*\\s*(tms|tbdms)\\s*(derivative)?\\b") %>%
      str_remove_all("(?i)\\s*,?\\s*(tms|tbdms)\\s*ester\\b") %>%
      str_remove_all("(?i)\\b(tms|tbdms)\\s*derivative\\b") %>%
      str_squish()
  ) %>%
  # Filter out known GC contaminants
  filter(!Compound_clean %in% c("Bis(heptamethylcyclotetrasiloxy)siloxane",
                                "Cyclotrisiloxane, hexamethyl-",
                                "2-Vinylphenol",
                                "Dimethyl phthalate",
                                "Naphthalene, 2-methyl-",
                                "Acetonitrile, trifluoro-"))

df_clean

Export(df_clean, "compound_list_05_clean.txt")


# ==========================================================
# 4. Prepare Compound Lists for External Mapping
# ==========================================================
#
# If the cleaned compound list contains more than 350 compounds,
# split the list into smaller batches before external mapping.
#
# This helps reduce mapping errors or information loss when processing
# large compound lists through external compound-mapping resources.

cleanlist1 <- df_clean[1:350, 2]
cleanlist2 <- df_clean[351:700, 2]
cleanlist3 <- df_clean[701:nrow(df_clean), 2]

Export(cleanlist1, "cleanlist_1.txt")
Export(cleanlist2, "cleanlist_2.txt")
Export(cleanlist3, "cleanlist_3.txt")
