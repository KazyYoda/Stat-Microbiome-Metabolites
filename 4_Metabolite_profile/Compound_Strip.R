############################################################################
# Metabolite Name Cleaning Pipeline
# Purpose: Remove derivatization artifacts and standardize compound names
############################################################################


setwd("~/Documents/HillTribe_NGS/5.Metabolites/Compound_Mapping")


# ==============================
# 1. Environment Setup
# ==============================

library(dplyr)
library(rio)
library(car)  
library(stringr)


# ==========================================================
# 2. Compound List 0.05 Processing
# ==========================================================

# Load compound list
compound_list_05 <- read.delim(
  "~/Documents/HillTribe_NGS/5.Metabolites/Compound_Mapping/compound_list_05.txt",
  header = FALSE
)

# Rename column
compound_list_05 <- compound_list_05 %>% 
  rename(Compound = V1)



# ==========================================================
# 3. Remove derivatization annotations
# ==========================================================

# Method 1
compound_list_05_clean <- compound_list_05 %>%
  mutate(
    Compound_clean = str_remove_all(
      Compound,
      ",\\s*\\d*TMS( derivative)?|\\b\\d+TMS\\b|\\s*TMS derivative|,\\s*\\d+TBDMS\\b"
    ) %>%
      str_squish()
  ) 

Export(compound_list_05_clean, "compound_list_05_clean.txt")


# Method 2
df_clean <- compound_list_05 %>%
  mutate(
    Compound_clean = Compound %>%
      str_remove_all("(?i)\\s*,?\\s*\\d*\\s*(tms|tbdms)\\s*(derivative)?\\b") %>%
      str_remove_all("(?i)\\s*,?\\s*(tms|tbdms)\\s*ester\\b") %>%
      str_remove_all("(?i)\\b(tms|tbdms)\\s*derivative\\b") %>%
      str_squish()
  ) %>%
  # Filter out any GC contaminants
  filter(!Compound_clean %in% c("Bis(heptamethylcyclotetrasiloxy)siloxane",
                                "Cyclotrisiloxane, hexamethyl-",
                                "2-Vinylphenol",
                                "Dimethyl phthalate",
                                "Naphthalene, 2-methyl-",
                                "Acetonitrile, trifluoro-"))

df_clean

Export(df_clean, "compound_list_05_clean.txt")


