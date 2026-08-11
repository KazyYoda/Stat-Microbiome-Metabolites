#============================================================
# Most Enriched Group: Comparison of Significant Metabolites
#============================================================
#
# This section compares significant upregulated metabolites
# identified from geographic and ethnic analytical frameworks.
#
# The objective is to determine whether the relative strength
# of metabolite enrichment is greater when groups are compared
# geographically or ethnically.
#
# For each matched metabolite, the log2FC from the geographic
# comparison is compared with the corresponding log2FC from the
# ethnic comparison. The absolute difference between the two
# log2FC values (Log2FC_diff) is used to quantify the difference
# in enrichment strength between analytical frameworks.
#
# A larger Log2FC_diff indicates a greater difference in the
# estimated enrichment strength between the two frameworks.
#
# Only metabolites identified as upregulated in the relevant
# comparisons are considered.

load("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat/Metabolites_Akha.RData")
load("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat/Metabolites_Lahu.RData")
load("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat/Metabolites_CM.RData")
load("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat/Metabolites_CR.RData")

setwd("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat")

# Load required packages
library(dplyr)
library(tidyr)
library(car)
library(readxl)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggh4x)
library(tidytext)


#------------------------------------------------------------
# Akha-CM: Geographic vs. Ethnic Enrichment
#------------------------------------------------------------
#
# Geographic framework:
#   Akha-CM vs. Akha-CR
#
# Ethnic framework:
#   Akha-CM vs. another ethnic group within CM
#
# Shared significant upregulated metabolites are identified
# between the two analytical frameworks and their log2FC values
# are compared.


#------------------------------------------------------------
# Akha-CM: Geographic vs. Akha-CM vs. Lahu-CM
#------------------------------------------------------------

# Match significant upregulated metabolites between the
# geographic comparison and the ethnic comparison using
# Lahu-CM as the ethnic reference group.
matched_AkhaGeo_EthvsLahu <- log2FC_sum_Met05_Akha_sig %>%
  filter(Group == "Akha-CM", Direction == "Upregulated") %>%
  select(Group, Metabolite, log2FC, Compound_Clean, Direction) %>%
  inner_join(
    log2FC_sum_Met05_CM2_sig_anno_majorsubclass %>%
      filter(Direction == "Upregulated") %>%
      select(Group, Metabolite, log2FC, Major_Subclass, Direction),
    by = c("Group", "Metabolite"),
    suffix = c("_AkhaGeo", "_AkhavsLahu")
  ) %>%
  # Calculate mean log2FC difference
  mutate(
    Log2FC_strength = case_when(
      log2FC_AkhaGeo > log2FC_AkhavsLahu ~
        "Akha-CM Geo",
      log2FC_AkhavsLahu > log2FC_AkhaGeo ~
        "Akha-CM Eth (vs. Lahu-CM)",
      TRUE ~ "Equal"
    )
  ) %>%
  mutate(
    Log2FC_diff = abs(log2FC_AkhaGeo - log2FC_AkhavsLahu),
    Framework = "AkhaCM Geo vs AkhaCM Eth (Lahu contrast)"
  ) %>%
  select(-Direction_AkhaGeo,
         -Direction_AkhavsLahu) %>%
  relocate(c(Compound_Clean, Major_Subclass), .before = c(3, 4)) %>%
  filter(!Log2FC_strength == "Equal") %>%
  rename(log2FC_AkhavsOther = log2FC_AkhavsLahu)


#------------------------------------------------------------
# Akha-CM: Geographic vs. Akha-CM vs. Khuen-CM
#------------------------------------------------------------

# Match significant upregulated Akha-CM geographic metabolites
# with metabolites that are downregulated in the Khuen-CM vs.
# Akha-CM comparison.
#
# The log2FC is reversed so that the ethnic comparison is expressed
# in the same direction as the Akha-CM enrichment:
#
#   Akha-CM > Khuen-CM
#
# This allows the geographic and ethnic log2FC values to be
# directly compared.
matched_AkhaGeo_EthvsKhuen <- log2FC_sum_Met05_Akha_sig %>%
  filter(Group == "Akha-CM", Direction == "Upregulated") %>%
  select(
    Group,
    Metabolite,
    log2FC_AkhaGeo = log2FC,
    Compound_Clean,
    Direction
  ) %>%
  
  inner_join(
    log2FC_sum_Met05_CM_sig_anno_majorsubclass %>%
      filter(Group == "Khuen-CM", Direction == "Downregulated") %>%
      select(
        Metabolite,
        log2FC_KhuenvsAkha = log2FC,
        Major_Subclass
      ),
    by = "Metabolite"
  ) %>%
  
  mutate(
    # Reverse the Khuen-CM vs. Akha-CM log2FC so that the
    # comparison is expressed as Akha-CM vs. Khuen-CM.
    log2FC_AkhavsKhuen = -log2FC_KhuenvsAkha,
    
    # Identify which analytical framework shows the stronger
    # enrichment of the metabolite in Akha-CM.
    Log2FC_strength = case_when(
      log2FC_AkhaGeo > log2FC_AkhavsKhuen ~
        "Akha-CM Geo",
      
      log2FC_AkhavsKhuen > log2FC_AkhaGeo ~
        "Akha-CM Eth (vs. Khuen-CM)",
      
      TRUE ~ "Equal"
    ),
    
    # Absolute difference in enrichment strength between the
    # geographic and ethnic frameworks.
    Log2FC_diff = abs(log2FC_AkhaGeo - log2FC_AkhavsKhuen),
    Framework = "AkhaCM Geo vs AkhaCM Eth (Khuen contrast)"
  )   %>%
  select(-Direction, -log2FC_KhuenvsAkha) %>%
  relocate(c(Compound_Clean, Major_Subclass), .before = c(3, 4)) %>%
  rename(log2FC_AkhavsOther = log2FC_AkhavsKhuen)


#------------------------------------------------------------
# Combine Geographic vs. Ethnic Frameworks
#------------------------------------------------------------

# Combine the two Akha-CM comparison frameworks into a single
# table for downstream filtering and visualization.
enr_AkhaCMGeo_Eth <- rbind(
  matched_AkhaGeo_EthvsLahu,
  matched_AkhaGeo_EthvsKhuen
)


#------------------------------------------------------------
# Identify Metabolites with Meaningful Differences in Enrichment
#------------------------------------------------------------

# Retain metabolites for which the absolute difference between
# geographic and ethnic log2FC values is at least 1.
#
# A difference of 1 on the log2 scale corresponds to a two-fold
# difference in the ratio of the corresponding fold changes.
Akha <- enr_AkhaCMGeo_Eth %>%
  filter(Log2FC_diff >= 1)


#------------------------------------------------------------
# Set Factor Levels for Visualization
#------------------------------------------------------------

# Order metabolite categories and comparison types for consistent
# plotting and interpretation.
Akha <- Akha %>%
  mutate(
    Major_Subclass = factor(
      Major_Subclass,
      levels = c(
        sort(unique(Major_Subclass[Major_Subclass != "Others"])),
        "Others"
      )
    ),
    
    Log2FC_strength = factor(
      Log2FC_strength,
      levels = c(
        "Akha-CM Geo",
        "Akha-CM Eth (vs. Lahu-CM)",
        "Akha-CM Eth (vs. Khuen-CM)"
      ))
  )


#------------------------------------------------------------
# Visualization Settings
#------------------------------------------------------------

# Okabe-Ito color palette for metabolite categories.
# The palette is designed to provide good visual distinction
# and is generally suitable for color-vision deficiencies.
okabe_ito <- c(
  "Hydrocarbons" = "#E69F00",
  "Lipids & FA derivatives" = "#56B4E9",
  "Oxygenated molecules" = "#009E73",
  "Phenolic & aromatic" = "#CC79A7",
  "Terpenoids" = "#0072B2",
  "Nitrogen-containing" = "#D55E00",
  "Unknown" = "grey70",
  "Others" = "grey50"
)


# Shapes indicate which analytical framework produced the
# greater log2FC for a given metabolite.
group_shapes <- c(
  "Akha-CM Geo" = 16,   
  "Akha-CM Eth (vs. Lahu-CM)" = 2,     
  "Akha-CM Eth (vs. Khuen-CM)" = 18
)


#------------------------------------------------------------
# Plot: Geographic vs. Ethnic Enrichment Strength
#------------------------------------------------------------

# Visualize the absolute difference in log2FC between geographic
# and ethnic analytical frameworks.
#
# Bar height represents Log2FC_diff.
# Fill represents the metabolite major subclass.
# Point shape identifies the framework with the greater log2FC.
Akha %>%
  mutate(
    Framework = factor(Framework, levels = c("AkhaCM Geo vs AkhaCM Eth (Lahu contrast)",
                                             "AkhaCM Geo vs AkhaCM Eth (Khuen contrast)")),
    
    Metabolite = reorder_within(Metabolite, Log2FC_diff, Framework)
  ) %>%

  rename(GreaterDiffPair = Log2FC_strength) %>%
  ggplot(aes(x = reorder(Metabolite, desc(Log2FC_diff)),
             y = Log2FC_diff,
             fill = Major_Subclass)) +
  
  geom_col(width = 0.7) +
  
  # Add a point at the end of each bar to indicate which
  # analytical framework has the greater log2FC.
  geom_point(aes(
    y = ifelse(Log2FC_diff > 0, Log2FC_diff + 0.5, Log2FC_diff - 0.3),
    shape = GreaterDiffPair,
    color = GreaterDiffPair
  ),
  size = 1.3) +
  
  scale_x_reordered() +
  
  scale_fill_manual(values = okabe_ito) +
  scale_shape_manual(values = group_shapes) +
  scale_color_manual(values = c(
    "Akha-CM Geo" = "black",
    "Akha-CM Eth (vs. Lahu-CM)" = "black",
    "Akha-CM Eth (vs. Khuen-CM)" = "black"
  )) +
  
  coord_flip() +
  
  facet_wrap2(
    ~Framework,
    scales = "free_y",
    strip = strip_themed(
      background_x = elem_list_rect(
        fill = c("lightblue2", "darkorange")
      ),
      text_x = elem_list_text(
        colour = "black",
        face = "bold",
        size = 7)
      )
  ) +
  
  theme_bw() +
  
  labs(
    x = "Metabolite",
    y = "Absolute difference in log2FC between geographic and ethnic analytical frameworks (|log2FCdiff >= 1|",
    fill = "Metabolite category",
    title = "Akha-CM metabolites",
    subtitle = "Comparison of absolute differences in mean log2FC between geographic and ethnic analytical frameworks") +
  theme(axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 8),
        title = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 8),
        strip.text = element_text(size = 6)
  )






#------------------------------------------------------------
# Khuen-CM: Comparison of Ethnic Enrichment Strength
#------------------------------------------------------------
#
# This analysis identifies metabolites that are significantly
# upregulated in Khuen-CM and compares their enrichment strength
# between two ethnic contrasts:
#
#   1. Khuen-CM vs. Akha-CM
#   2. Khuen-CM vs. Lahu-CM
#
# Only metabolites identified as upregulated in both analytical
# frameworks are retained.
#
# For each shared metabolite, the absolute difference between the
# two log2FC values (Log2FC_diff) is calculated to quantify how
# differently the metabolite is enriched in Khuen-CM relative to
# the two comparison groups.
#
# A larger Log2FC_diff indicates a greater difference in enrichment
# strength between the two ethnic contrasts.
#------------------------------------------------------------


matched_KhuenCM_CM <- log2FC_sum_Met05_CM_sig %>%
  filter(Group == "Khuen-CM", Direction == "Upregulated") %>%
  select(Group, Metabolite, log2FC, Compound_Clean, Direction) %>%
  inner_join(
    log2FC_sum_Met05_CM2_sig_anno_majorsubclass %>%
      filter(Direction == "Upregulated") %>%
      select(Group, Metabolite, log2FC, Major_Subclass, Direction),
    by = c("Group", "Metabolite"),
    suffix = c("_KhuenvsAkha", "_KhuenvsLahu")
  ) %>%
  
  mutate(
    # Identify the ethnic comparison in which the metabolite
    # shows the greater log2FC in Khuen-CM.
    Log2FC_strength = case_when(
      log2FC_KhuenvsAkha > log2FC_KhuenvsLahu ~
        "Khuen-CM (vs. Akha-CM)",
      
      log2FC_KhuenvsLahu > log2FC_KhuenvsAkha ~
        "Khuen-CM (vs. Lahu-CM)",
      
      TRUE ~ "Equal"
    ),
    
    # Quantify the absolute difference in enrichment strength
    # between the two ethnic comparisons.
    Log2FC_diff = abs(log2FC_KhuenvsAkha - log2FC_KhuenvsLahu)
  ) %>%
  
  # Remove duplicated direction variables after the comparison.
  select(-Direction_KhuenvsAkha,
         -Direction_KhuenvsLahu) %>%
  
  # Move metabolite annotation variables next to the metabolite
  # identifier for easier inspection.
  relocate(c(Compound_Clean, Major_Subclass), .before = c(3, 4)) %>%
  
  # Retain metabolites with different enrichment strengths and
  # an absolute log2FC difference of at least 1.
  filter(!Log2FC_strength == "Equal",
         Log2FC_diff >= 1)


# Export the shared Khuen-CM metabolite comparison table.
Export(matched_KhuenCM_CM, "shared_KhuenCM.txt")



#------------------------------------------------------------
# Filter Metabolites with Large Differences in Enrichment
#------------------------------------------------------------

# Retain metabolites with a substantial difference in log2FC
# between the two ethnic comparisons.
#
# Here, Log2FC_diff >= 5 corresponds to a difference of five
# units on the log2 scale between the two enrichment estimates.
KhuenCM <- matched_KhuenCM_CM %>%
  filter(Log2FC_diff >= 5)


#------------------------------------------------------------
# Visualization Settings
#------------------------------------------------------------

# Shapes indicate which ethnic comparison has the greater
# enrichment strength in Khuen-CM.
group_shapes <- c(
  "Khuen-CM (vs. Akha-CM)" = 16,   # circle
  "Khuen-CM (vs. Lahu-CM)" = 2     # triangle
)


# Order metabolite subclasses for consistent visualization.
KhuenCM <- KhuenCM %>%
  mutate(
    Major_Subclass = factor(
      Major_Subclass,
      levels = c(
        sort(unique(Major_Subclass[Major_Subclass != "Others"])),
        "Others")
    )
  )


# Okabe-Ito color palette for metabolite categories.
okabe_ito <- c(
  "Hydrocarbons" = "#E69F00",
  "Lipids & FA derivatives" = "#56B4E9",
  "Oxygenated molecules" = "#009E73",
  "Phenolic & aromatic" = "#CC79A7",
  "Terpenoids" = "#0072B2",
  "Nitrogen-containing" = "#D55E00",
  "Unknown" = "grey70",
  "Others" = "grey50"
)


#------------------------------------------------------------
# Plot: Khuen-CM Enrichment Strength Across Ethnic Contrasts
#------------------------------------------------------------

# Visualize metabolites showing large differences in Khuen-CM
# enrichment strength between the Akha-CM and Lahu-CM contrasts.
#
# Bar height:
#   Absolute difference in log2FC (Log2FC_diff)
#
# Bar fill:
#   Major metabolite subclass
#
# Point shape:
#   Ethnic comparison with the greater Khuen-CM log2FC
KhuenCM %>%
  rename(GreaterDiffPair = Log2FC_strength) %>%
  ggplot(aes(x = reorder(Metabolite, desc(Log2FC_diff)),
             y = Log2FC_diff,
             fill = Major_Subclass)) +
  
  geom_col(width = 0.7) +
  
  # Add a point at the end of each bar to indicate which
  # ethnic comparison has the greater log2FC.
  geom_point(aes(
    y = ifelse(Log2FC_diff > 0, Log2FC_diff + 0.2, Log2FC_diff - 0.2),
    shape = GreaterDiffPair,
    color = GreaterDiffPair
  ),
  size = 1.5) +
  
  scale_fill_manual(values = okabe_ito) +
  scale_shape_manual(values = group_shapes) +
  scale_color_manual(values = c(
    "Khuen-CM (vs. Akha-CM)" = "black",
    "Khuen-CM (vs. Lahu-CM)" = "black"
  )) +
  
  coord_flip() +
  
  theme_bw() +
  
  labs(
    x = "Metabolite",
    y = "Absolute difference in log2FC within the CM group (|log2FCdiff >= 5|)",
    fill = "Metabolite category",
    title = "Khuen-CM metabolites",
    subtitle = "Comparison of absolute differences in mean log2FC within the CM group") +
  theme(axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 8),
        title = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 8)
  )






#------------------------------------------------------------
# Akha-CR: Comparison of Geographic Enrichment Strength
#------------------------------------------------------------
#
# This analysis compares the enrichment strength of metabolites
# associated with Akha-CR across two analytical frameworks:
#
#   1. Akha-CR vs. Akha-CM
#      -> geographic comparison within the Akha ethnicity
#
#   2. Akha-CR vs. Lahu-CR
#      -> ethnic comparison within the CR geographic area
#
# Shared metabolites are matched between the two analyses using
# the metabolite identifier. For each matched metabolite, the
# log2FC values from the two comparisons are compared.
#
# Log2FC_diff represents the absolute difference between the two
# log2FC estimates and therefore quantifies the difference in
# enrichment strength of Akha-CR between the two analytical
# frameworks.
#
# Log2FC_strength identifies the comparison in which the metabolite
# shows the greater log2FC for Akha-CR.
#------------------------------------------------------------


matched_AkhaCR_CR2 <- log2FC_sum_Met05_Akha_sig %>%
  select(Group, Metabolite, log2FC, Compound_Clean, Direction) %>%
  
  inner_join(
    log2FC_sum_Met05_CR2_sig_anno_majorsubclass %>%
      select(Group, Metabolite, log2FC, Compound_Clean, Major_Subclass, Direction),
    by = c("Group", "Metabolite"),
    suffix = c("_Akha", "_CR2")
  ) %>%
  
  mutate(
    # Identify the analytical framework in which the metabolite
    # shows the greater enrichment strength for Akha-CR.
    Log2FC_strength = case_when(
      log2FC_Akha > log2FC_CR2 ~
        "Akha-CR (vs. Akha-CM) enriched",
      
      log2FC_CR2 > log2FC_Akha ~
        "Akha-CR (vs. Lahu-CR) enriched",
      
      TRUE ~ "Equal"
    ),
    
    # Absolute difference in log2FC between the two analytical
    # frameworks.
    Log2FC_diff = abs(log2FC_Akha - log2FC_CR2)
  )

# No matched metabolites identified in this comparison.
# Number of matched metabolites: 0





#------------------------------------------------------------
# Lahu-CM: Comparison of Geographic and Ethnic Enrichment
#------------------------------------------------------------
#
# Compare the log2FC values of Lahu-CM metabolites between:
#
#   1. Geographic framework (Lahu-CM comparison)
#   2. Ethnic framework within the CM area
#
# Shared metabolites are matched by metabolite identifier.
# Log2FC_strength identifies the framework with the greater
# enrichment magnitude, while Log2FC_diff quantifies the
# absolute difference between the two log2FC estimates.
#------------------------------------------------------------

matched_LahuCM_CM <- log2FC_sum_Met05_Lahu_sig %>%
  select(Group, Metabolite, log2FC, Compound_Clean, Direction) %>%
  inner_join(
    log2FC_sum_Met05_CM_sig_anno_majorsubclass %>%
      select(Group, Metabolite, log2FC, Compound_Clean, Major_Subclass, Direction),
    by = c("Group", "Metabolite"),
    suffix = c("_Lahu", "_CM")
  ) %>%
  mutate(
    # Identify the analytical framework with the greater
    # enrichment strength for Lahu-CM.
    Log2FC_strength = case_when(
      log2FC_Lahu > log2FC_CM ~ "Lahu-CM (Geo) enriched",
      log2FC_CM > log2FC_Lahu ~ "Lahu-CM (Eth) enriched",
      TRUE ~ "Equal"
    ),
    
    # Absolute difference in log2FC between the two frameworks.
    Log2FC_diff = abs(log2FC_Lahu - log2FC_CM)
  )

# Number of matched metabolites: 33
Export(matched_LahuCM_CM, "shared_LahuCM.txt")



#------------------------------------------------------------
# Lahu-CR: Comparison of Geographic and Ethnic Enrichment
#------------------------------------------------------------
#
# Compare the log2FC values of Lahu-CR metabolites between:
#
#   1. Geographic framework (Lahu-CR comparison)
#   2. Ethnic framework within the CR area
#
# Shared metabolites are matched by metabolite identifier.
# Log2FC_strength identifies the framework with the greater
# enrichment magnitude, while Log2FC_diff quantifies the
# absolute difference between the two log2FC estimates.
#------------------------------------------------------------

matched_LahuCR_CR <- log2FC_sum_Met05_Lahu_sig %>%
  select(Group, Metabolite, log2FC, Compound_Clean, Direction) %>%
  inner_join(
    log2FC_sum_Met05_CR_sig_anno_majorsubclass %>%
      select(Group, Metabolite, log2FC, Compound_Clean, Major_Subclass, Direction),
    by = c("Group", "Metabolite"),
    suffix = c("_Lahu", "_CR")
  ) %>%
  mutate(
    # Identify the analytical framework with the greater
    # enrichment strength for Lahu-CR.
    Log2FC_strength = case_when(
      log2FC_Lahu > log2FC_CR ~ "Lahu-CR (Geo) enriched",
      log2FC_CR > log2FC_Lahu ~ "Lahu-CR (Eth) enriched",
      TRUE ~ "Equal"
    ),
    
    # Absolute difference in log2FC between the two frameworks.
    Log2FC_diff = abs(log2FC_Lahu - log2FC_CR)
  )

# No matched metabolites identified in this comparison.
# Number of matched metabolites: 0



#------------------------------------------------------------
# Lisu-CR: Comparison of Ethnic Enrichment Strength
#------------------------------------------------------------
#
# Compare the log2FC values of Lisu-CR metabolites between two
# ethnic contrasts within the CR area:
#
#   1. Lisu-CR vs. Akha-CR
#   2. Lisu-CR vs. Lahu-CR
#
# Shared metabolites are matched by metabolite identifier.
# Log2FC_strength identifies the ethnic comparison with the
# greater enrichment strength for Lisu-CR.
#
# Log2FC_diff represents the absolute difference between the
# two log2FC estimates.
#------------------------------------------------------------

matched_LisuCR_CR <- log2FC_sum_Met05_CR_sig %>%
  select(Group, Metabolite, log2FC, Compound_Clean, Direction) %>%
  inner_join(
    log2FC_sum_Met05_CR2_sig_anno_majorsubclass %>%
      select(Group, Metabolite, log2FC, Compound_Clean, Major_Subclass, Direction),
    by = c("Group", "Metabolite"),
    suffix = c("_Lisu", "_CR2")
  ) %>%
  mutate(
    # Identify the ethnic comparison with the greater
    # enrichment strength for Lisu-CR.
    Log2FC_strength = case_when(
      log2FC_Lisu > log2FC_CR2 ~ "Lisu-CR (vs. Akha-CR) enriched",
      log2FC_CR2 > log2FC_Lisu ~ "Lisu-CR (vs. Lahu-CR) enriched",
      TRUE ~ "Equal"
    ),
    
    # Absolute difference in log2FC between the two ethnic
    # comparisons.
    Log2FC_diff = abs(log2FC_Lisu - log2FC_CR2)
  )

# Number of matched metabolites: 62
Export(matched_LisuCR_CR, "shared_LisuCR2.txt")





# ==============================================================================
# Lisu-CR enriched metabolites
# ==============================================================================
# Identify metabolites enriched in Lisu-CR that show consistent or opposing
# enrichment directions across the geographic and ethnicity comparisons.
#
# Metabolites are retained when:
#   1. The enrichment direction is either Upregulated or Opposite.
#   2. The difference in log2 fold-change between comparisons is >= 1.
#
# The resulting dataset is annotated with metabolite category information and
# used to visualize the difference in enrichment strength between comparisons.
# ==============================================================================

# Data manipulation ------------------------------------------------------------
# Classify metabolites according to the direction of enrichment in the two
# comparisons and retain metabolites meeting the predefined log2FC difference
# threshold.
LisuCR <- matched_LisuCR_CR %>%
  mutate(Direction = case_when(
    Direction_Lisu == "Upregulated" & Direction_CR2 == "Upregulated" ~ "Upregulated",
    Direction_Lisu == "Downregulated" & Direction_CR2 == "Downregulated" ~ "Downregulated",
    Direction_Lisu == "Upregulated" & Direction_CR2 == "Downregulated" ~ "Opposite",
    Direction_Lisu == "Downregulated" & Direction_CR2 == "Upregulated" ~ "Opposite"
  )) %>%
  filter(Direction %in% c("Upregulated", "Opposite"), 
         Log2FC_diff >= 1) %>%
  select(Group, Metabolite, Log2FC_strength, Log2FC_diff, Direction, Major_Subclass) %>%
  left_join(metabo_05_descp,
            by = c("Metabolite" = "Code"))

# Define point shapes used to indicate the comparison represented by
# Log2FC_strength in the visualization.
group_shapes <- c(
  "Lisu-CR (vs. Akha-CR) enriched" = 16,   # circle
  "Lisu-CR (vs. Lahu-CR) enriched" = 17     # triangle
)

# Visualization ---------------------------------------------------------------
# Plot the difference in enrichment strength between the two comparisons.
# Bar color represents the major metabolite subclass, while point shape and
# color identify the comparison associated with the enrichment strength.
LisuCR %>%
  rename(Log2FC_level = Log2FC_strength) %>%
  ggplot(aes(x = reorder(Metabolite, desc(Log2FC_diff)),
             y = Log2FC_diff,
             fill = Major_Subclass)) +
  
  geom_col(width = 0.7) +
  
  # Add a point at the end of each bar to indicate the corresponding
  # enrichment comparison.
  geom_point(aes(
    y = ifelse(Log2FC_diff > 0, Log2FC_diff + 0.2, Log2FC_diff - 0.2),
    shape = Log2FC_level,
    color = Log2FC_level
  ),
  size = 2) +
  
  scale_fill_manual(values = okabe_ito) +
  scale_shape_manual(values = group_shapes) +
  scale_color_manual(values = c(
    "Lisu-CR (vs. Akha-CR) enriched" = "black",
    "Lisu-CR (vs. Lahu-CR) enriched" = "black"
  )) +
  
  coord_flip() +
  
  theme_bw() +
  
  labs(
    x = "Enriched metabolite",
    y = "Log2FC difference",
    fill = "Metabolite category",
    title = "Lisu-CR enriched metabolites",
    subtitle = "Differences in enrichment strength between Lisu-CR vs AKha-CR and Lisu-CR vs Lahu-CR comparisons") +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        title = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 8)
  )



# ==============================================================================
# Terpenoid metabolites: CM
# ==============================================================================
# Identify terpenoid metabolites enriched in the CM geographic and ethnicity
# analytical frameworks.
#
# Three comparison groups are considered:
#   - Akha-CM Geo vs. Akha-CR
#   - Akha-CM Eth vs. Lahu-CM
#   - Khuen-CM Eth vs. Lahu-CM
#
# The resulting datasets are combined for visualization using the metabolite
# subclass as the bar color and the analytical comparison as the facet.
# ==============================================================================

# Data manipulation ------------------------------------------------------------
# Select significantly upregulated terpenoid metabolites from the geographic
# comparison involving Akha-CM and Akha-CR.
AkhaCMGeo_up_terpene <- log2FC_sum_Met05_Akha_sig_anno_majorsubclass %>%
  filter(Group == "Akha-CM",
         Direction == "Upregulated",
         Major_Subclass == "Terpenoids") %>%
  select(Group, Metabolite, Compound_Clean, log2FC, Subclass, Major_Subclass) %>%
  mutate(Enriched_group = case_when(
    Group == "Akha-CM" ~ "Akha-CM Geo (vs. Akha-CR) enriched")
  )


# Select significantly upregulated terpenoid metabolites from the CM ethnicity
# comparisons involving Akha-CM and Khuen-CM relative to Lahu-CM.
CMEth_up_terpene <- log2FC_sum_Met05_CM2_sig_anno_majorsubclass %>%
  filter(Direction == "Upregulated",
         Major_Subclass == "Terpenoids") %>%
  select(Group, Metabolite, Compound_Clean, log2FC, Subclass, Major_Subclass) %>%
  mutate(Enriched_group = case_when(
    Group == "Akha-CM" ~ "Akha-CM Eth (vs. Lahu-CM) enriched",
    Group == "Khuen-CM" ~ "Khuen-CM Eth (vs. Lahu-CM) enriched")
  )


# Combine terpenoid metabolites from the geographic and ethnicity frameworks.
terpene_CM <- rbind(AkhaCMGeo_up_terpene,
                    CMEth_up_terpene)


# Define colors for metabolite subclasses used in the plot.
okabe_ito <- c(
  "grey50",   # black
  "#0072B2",  # blue
  "#bcbd22" # Olive
  
)


# Define facet-strip colors for each analytical comparison.
group_colors <- c(
  "Akha-CM Geo (vs. Akha-CR) enriched" = "gold1",
  "Akha-CM Eth (vs. Lahu-CM) enriched" = "lightblue2",
  "Khuen-CM Eth (vs. Lahu-CM) enriched" = "darkorange"
)

# Visualization ---------------------------------------------------------------
# Plot enriched terpenoid metabolites separately for each analytical framework.
# Bar color represents the metabolite subclass, while facet-strip colors
# distinguish the geographic and ethnicity comparisons.
terpene_CM %>%
  mutate(
    Enriched_group = factor(
      Enriched_group, levels = c("Akha-CM Geo (vs. Akha-CR) enriched",
                                 "Akha-CM Eth (vs. Lahu-CM) enriched",
                                 "Khuen-CM Eth (vs. Lahu-CM) enriched")
    )
  ) %>%
  ggplot(aes(x = reorder(Metabolite, desc(log2FC)),
             y = log2FC,
             fill = Subclass)) +
  
  geom_col(width = 0.7) +
  
  facet_wrap2(
    ~Enriched_group,
    strip = strip_themed(
      background_x = elem_list_rect(fill = group_colors),
      text_x = elem_list_text(size = 8)
    )
  ) +
  
  scale_fill_manual(values = okabe_ito) +
  
  coord_flip() +
  
  theme_bw() +
  
  labs(
    x = "Enriched metabolite",
    y = "Mean Log2FC",
    fill = "Metabolite category",
    title = "Enriched metabolites",
    subtitle = "Terpenoid-enriched metabolites across geographic and ethnicity analytical frameworks") +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        title = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 8)
  )





#========================================================
# 6. Heatmap Visualization (log2-transformed Data)
#========================================================
# Generate a heatmap of selected CM terpenoid metabolites using
# log2-transformed metabolite abundance data.
#
# The workflow includes:
#   1. Selecting terpenoid metabolites identified in the CM ethnicity
#      comparisons.
#   2. Retrieving compound annotation information.
#   3. Matching metabolites and samples to the log2-transformed matrix.
#   4. Performing a variance check before clustering.
#   5. Defining the order and colors of sample groups.
#   6. Generating and displaying the annotated heatmap.


#-----------------------------------------------------------------------------
# 6.2 Select terpenoid metabolites: CM group
#-----------------------------------------------------------------------------
# Extract unique metabolite identifiers for the selected CM terpenoid
# metabolites.
CM_terpene <- CMEth_up_terpene %>%
  pull(Metabolite) %>%
  unique()

# Retrieve compound identifiers and chemical annotations for the selected
# terpenoid metabolites. Metabolites are ordered according to their log2FC,
# and duplicate metabolite entries are removed.
CM_inchi_terpene <- CMEth_up_terpene %>%
  left_join(
    Pub_Classy_05_clean,
    by = c("Compound_Clean", "Subclass")
  ) %>%
  arrange(desc(log2FC)) %>%
  select(Metabolite, Compound_Clean, CID, IUPAC, Subclass) %>%
  distinct(Metabolite, .keep_all = TRUE) 

# Export the metabolite annotation table for reference and downstream use.
Export(CM_inchi_terpene, "CM_inchi_terpene.txt")


# Match metabolite rows --------------------------------------------------------
# Identify metabolites present in both the selected CM terpenoid list and
# the log2-transformed metabolite matrix.
valid_rows <- intersect(
  rownames(log_matrix_metabo_05_t),
  CM_terpene
)

# Match CM sample columns ------------------------------------------------------
# Identify samples present in both the log2-transformed metabolite matrix and
# the CM sample metadata.
valid_cols <- intersect(
  colnames(log_matrix_metabo_05_t),
  Met05_CM$Sample_ID
)

# Subset matrix ---------------------------------------------------------------
# Extract the selected CM terpenoid metabolites and matching CM samples.
# drop = FALSE preserves the matrix structure when only one row or column
# is present.
log_matrix_metabo_05_t_CM <-
  log_matrix_metabo_05_t[
    valid_rows,
    valid_cols,
    drop = FALSE
  ]


#--------------------------------------------------------
# 6.3 Variance Check (Pre-Clustering QC)
#--------------------------------------------------------
# Check for zero-variance metabolites. Features with zero variance do not
# contribute to Euclidean distance-based clustering.
sum(apply(log_matrix_metabo_05_t_CM, 1, sd, na.rm = TRUE) == 0)

# Check for zero-variance samples. Samples with zero variance may also affect
# distance calculations and hierarchical clustering.
sum(apply(log_matrix_metabo_05_t_CM, 2, sd, na.rm = TRUE) == 0)

# Note: Zero-variance features can affect Euclidean distance calculations
# and downstream clustering.


#--------------------------------------------------------
# 6.4 Define Group Order for Annotation
#--------------------------------------------------------
# Set the desired order of CM groups for sample annotation in the heatmap.
Met05_CM$Group <- factor(
  Met05_CM$Group,
  levels = c("Akha-CM", "Lahu-CM", "Khuen-CM")
)


#========================================================
# 6.5 Generate Heatmap
#========================================================
# Define colors for the three CM sample groups.
group_colors <- c(
  "Akha-CM" = "gold1",
  "Lahu-CM" = "lightblue2",
  "Khuen-CM" = "darkorange"
)

# Generate the heatmap using the log2-transformed metabolite matrix and
# corresponding sample metadata.
CM <- plot_heatmap(
  mat = log_matrix_metabo_05_t_CM,
  metadata = Met05_CM[1:13],
  legend_title = "Log2 (CM)",
  group_color = group_colors
)

# Draw the heatmap with both the heatmap and annotation legends positioned
# on the left side of the figure.
draw(
  CM,
  heatmap_legend_side = "left",
  annotation_legend_side = "left"
)
