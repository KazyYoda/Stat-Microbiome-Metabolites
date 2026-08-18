############################################################
# Microbiome Analysis in R: Functional Analysis
############################################################

# Set a working directory:
setwd("~/Documents/HillTribe_NGS/4.Functional_profile/ALDEx2")


# Load required packages
library(dplyr)
library(tidyr)
library(readxl)
library(car)
library(ggplot2)
library(ggh4x)
library(gridExtra)


# Define group colors
group_colors <- c(
  "Akha-CM" = "gold1",
  "Akha-CR" = "dodgerblue4"
)


group_colors <- c(
  "Lahu-CM" = "lightblue2",
  "Lahu-CR" = "coral1"
)

group_colors <- c(
  "Akha-CM" = "gold1",
  "Lahu-CM" = "lightblue2",
  "Khuen-CM" = "darkorange"
)


group_colors <- c(
  "Akha-CR" = "dodgerblue4",
  "Lahu-CR" = "coral1",
  "Lisu-CR" = "grey40"
)




#------------------------------------------
# 1. Import PICRUSt2 outputs
#------------------------------------------
load("~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/HAllA_prep.RData")


# EC and KO IDs prep
EC_ID <- data.frame(EC_ID = rownames(EC_abs))

write.table(
  EC_ID,
  file = "EC_ID.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)


KO_ID <- data.frame(KO_ID = rownames(KO_abs))

write.table(
  KO_ID,
  file = "KO_ID.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)



# Update KEGG KO descp
KO_KEGG <- Mapped_KO_KEGG %>%
  select(ko_id, ko_descp) %>%
  distinct()


KO_abs_descp <- KO_absolute %>% 
  select(-description) %>%
  left_join(KO_KEGG, by = c("KO_ID" = "ko_id")) %>%
  as.data.frame()

Export(KO_abs_descp, "KO_abs_updatedescp.xlsx")


#------------------------------------------
# 2. Perform ALDEx2 for Pathway
#------------------------------------------
library(ALDEx2)

#----------- ALDEx2 Subgroup comparisons: Helper Function with Filtering Option ---------------
run_aldex2_analysis <- function(df,        # Absolute table (e.g., PW, KO, EC)
                                metadata,  # Sample metadata table
                                fun_level, # Character string: functional name, e.g. "PW"
                                group_col, # Column in metadata to filter group comparison from (e.g. "Group")
                                group1,    # First group name (e.g., "Akha-CM")
                                group2,    # Second group name (e.g., "Akha-CR")
                                paired.test = FALSE, # whether to perform paired or cross t.test (aldex.ttest)
                                output_dir = ".",
                                filter_sig = FALSE,
                                p_thresh = 0.05,
                                effect_thresh = 0.5) {
  
  message(paste("Processing:", fun_level, "|", group1, "vs", group2))
  
  # SAFETY CHECK: group names exist
  if (!all(c(group1, group2) %in% metadata[[group_col]])) {
    stop(paste("Invalid comparison:", group1, "vs", group2,
               "- group not found in metadata"))
  }
  
  
  df <- round(df)
  
  # Add pseudocount
  df_pseudo <- df + 1
  
  # Select sample IDs for the two comparison groups
  selected_samples <- metadata %>%
    filter(!!sym(group_col) %in% c(group1, group2)) %>%
    pull(Sample_ID)
  
  df_filtered <- df_pseudo %>%
    dplyr::select(all_of(selected_samples))
  
  # Condition vector, ordered to match df_filtered
  conditions <- metadata %>%
    filter(Sample_ID %in% selected_samples) %>%
    arrange(match(Sample_ID, selected_samples)) %>%
    pull(!!sym(group_col))
  
  if (!identical(names(df_filtered), selected_samples)) {
    stop("Sample ID order mismatch.")
  }
  
  if (paired.test) {
    message("Running paired aldex2 test")
  } else {
    message("Running unpaired aldex2 test")
  }
  
  # Run ALDEx2
  set.seed(123)
  aldex_clr <- aldex.clr(df_filtered, conds = conditions,
                         mc.samples = 128, denom = "all", verbose = FALSE)
  
  aldex_ttest <- aldex.ttest(aldex_clr, paired.test = paired.test)
  aldex_effect <- aldex.effect(aldex_clr)
  
  aldex_combined <- data.frame(
    Feature = rownames(aldex_ttest),
    aldex_ttest,
    aldex_effect
  )
  
  # Filter by significance
  if (filter_sig) {
    aldex_combined <- aldex_combined %>%
      dplyr::filter(we.ep < p_thresh,
                    wi.ep < p_thresh,
                    abs(effect) >= effect_thresh)
  }
  
  # Output
  file_name <- paste0("ALDEx2_", fun_level, "_", group1, "_vs_", group2,
                      if (filter_sig) "_filtered" else "", ".tsv")
  Export(aldex_combined, file.path(output_dir, file_name))
  
  return(aldex_combined)
}





# Cross: Run Loop through count tables: PW_abs, KO_abs, EC_abs
# Input files (PICRUSt2 outputs): 
# 1. PW_abs from PW_abun_unstrat.tsv
# 2. EC_abs from EC_pred_metagenome_unstrat.tsv
# 3. KO_abs from KO_pred_metagenome_unstrat.tsv                      
                    
funtional_level <- c("PW", "KO", "EC") 
tables <- list(PW = PW_abs, KO = KO_abs, EC = EC_abs)

comparisons <- list(
  Akha   = c("Akha-CM", "Akha-CR"),
  Lahu   = c("Lahu-CM", "Lahu-CR"),
  CM_1   = c("Akha-CM", "Lahu-CM"),
  CM_2   = c("Akha-CM", "Khuen-CM"),
  CM_3   = c("Lahu-CM", "Khuen-CM"),
  CR_1   = c("Akha-CR", "Lahu-CR"),
  CR_2   = c("Akha-CR", "Lisu-CR"),
  CR_3   = c("Lahu-CR", "Lisu-CR")
)


# Loop through all combinations
results_list <- list()

for (fun_level in funtional_level) {
  for (comp in names(comparisons)) {
    group_pair <- comparisons[[comp]]
    
    res <- run_aldex2_analysis(df = tables[[fun_level]],
                               metadata = sample_metadata,
                               fun_level = fun_level,
                               group_col = "Group",
                               group1 = group_pair[1],
                               group2 = group_pair[2],
                               paired.test = FALSE,
                               filter_sig = FALSE,
                               output_dir = "crossALDEx2_results")
    
    results_list[[paste(fun_level, comp, sep = "_")]] <- res
  }
}





#------------------------------------------
# Visualization of ALDEx2 results
#------------------------------------------
#-- ALDEx2 “effect”: How large the between-group difference is relative to within-group variation.
#-- interpret the magnitude: 
#-- < 0.5 = Very small / negligible
#-- 0.5-1 = Small difference
#-- 1-2 = Moderate, biologically meaningful
# >2 = Large, strong, robust difference


setwd("~/Documents/HillTribe_NGS/4.Functional_profile/ALDEx2")
load("~/Documents/HillTribe_NGS/4.Functional_profile/ALDEx2/ALDEx2.RData")


library(forcats) 
library(ggplot2)

# ------------------------------------------------------------
# PW
# ------------------------------------------------------------
PW_aldex2_AkhaCM_AkhaCR <- read_excel("PW_aldex2_AkhaCM_AkhaCR.xlsx")
PW_aldex2_LahuCM_LahuCR <- read_excel("PW_aldex2_LahuCM_LahuCR.xlsx")
PW_aldex2_CM <- read_excel("PW_aldex2_CM.xlsx")
PW_aldex2_CR <- read_excel("PW_aldex2_CR.xlsx")


# Filter wi.eBH < 0.05 (Akha and Lahu groups)
sigPW_aldex2_Akha <- PW_aldex2_AkhaCM_AkhaCR %>% 
  filter(wi.eBH < 0.05)

sigPW_aldex2_Lahu <- PW_aldex2_LahuCM_LahuCR %>% 
  filter(wi.eBH < 0.05)


# Adjust wi.ep for multiple comparisons (CM and CR groups)
sigPW_aldex2_CM <- PW_aldex2_CM %>% 
  mutate(q_value = p.adjust(wi.ep, method = "BH")) %>%
  filter(q_value < 0.05, abs(effect) >= 0.5)

sigPW_aldex2_CR <- PW_aldex2_CR %>% 
  mutate(q_value = p.adjust(wi.ep, method = "BH")) %>%
  filter(q_value < 0.05, abs(effect) >= 0.5)



# Show Level as label and reorder features by diff.btw:
sigPW_aldex2_Akha %>%
  filter(
    wi.eBH < 0.05 &
      abs(effect) >= 1.5) %>%
  mutate(
    Feature = fct_reorder(Feature, diff.btw),
    label_text = round(wi.eBH, 3),
  ) %>%
  ggplot(aes(x = diff.btw, y = Feature)) +
  geom_col(aes(fill = effect)) +
  geom_text(
    aes(
      x = ifelse(diff.btw >= 0, diff.btw + 0.1, diff.btw - 0.1),
      label = label_text,
      hjust = ifelse(diff.btw >= 0, 0, 1)
    ),
    color = "black",
    size = 1.5
  ) +
  scale_fill_gradient2(
    low = "steelblue3",
    mid = "ivory1",
    high = "#FC766A",
    midpoint = 0,
    name = "Effect size"
  ) +
  labs(
    x = "Difference between Groups (diff.btw: ALDEx2)",
    y = "Feature",
    title = "Differential abundance identified by ALDEx2 (wi.eBH < 0.05, abs(effect) >= 1.5 )"
  ) +
  facet_wrap(~ Comparison, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 7),
    title = element_text(size = 7, face = "bold"),
    panel.grid = element_blank(),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 8),
    plot.margin = margin(10, 100, 10, 10)
  ) +
  coord_cartesian(clip = "off") +
  expand_limits(x = c(min(sigPW_aldex2_Akha$diff.btw) - 1.2, max(sigPW_aldex2_Akha$diff.btw) + 1.2))





# ------------------------------------------------------------
# EC
# ------------------------------------------------------------
EC_aldex2_AkhaCM_AkhaCR <- read_excel("EC_aldex2_AkhaCM_AkhaCR.xlsx")
EC_aldex2_LahuCM_LahuCR <- read_excel("EC_aldex2_LahuCM_LahuCR.xlsx")
EC_aldex2_CM <- read_excel("EC_aldex2_CM.xlsx")
EC_aldex2_CR <- read_excel("EC_aldex2_CR.xlsx")


# Filter wi.eBH < 0.05 (Akha and Lahu groups)
sigEC_aldex2_Akha <- EC_aldex2_AkhaCM_AkhaCR %>% 
  filter(wi.eBH < 0.05)

sigEC_aldex2_Lahu <- EC_aldex2_LahuCM_LahuCR %>% 
  filter(wi.eBH < 0.05)


# Adjust wi.ep for multiple comparisons (CM and CR groups)
sigEC_aldex2_CM <- EC_aldex2_CM %>% 
  mutate(q_value = p.adjust(wi.ep, method = "BH")) %>%
  filter(q_value < 0.05, abs(effect) >= 0.5)

sigEC_aldex2_CR <- EC_aldex2_CR %>% 
  mutate(q_value = p.adjust(wi.ep, method = "BH")) %>%
  filter(q_value < 0.05, abs(effect) >= 0.5)



# Show Level as label and reorder features by diff.btw:
sigEC_aldex2_Akha %>%
  filter(
    wi.eBH < 0.05 &
      abs(effect) >= 1.5) %>%
  mutate(
    Feature = fct_reorder(Feature, diff.btw),
    label_text = round(wi.eBH, 3),
  ) %>%
  ggplot(aes(x = diff.btw, y = Feature)) +
  geom_col(aes(fill = effect)) +
  geom_text(
    aes(
      x = ifelse(diff.btw >= 0, diff.btw + 0.1, diff.btw - 0.1),
      label = label_text,
      hjust = ifelse(diff.btw >= 0, 0, 1)
    ),
    color = "black",
    size = 1.5
  ) +
  scale_fill_gradient2(
    low = "steelblue3",
    mid = "ivory1",
    high = "#FC766A",
    midpoint = 0,
    name = "Effect size"
  ) +
  labs(
    x = "Difference between Groups (diff.btw: ALDEx2)",
    y = "Feature",
    title = "Differential abundance identified by ALDEx2 (wi.eBH < 0.05, abs(effect) >= 1.5 )"
  ) +
  facet_wrap(~ Comparison, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 7),
    title = element_text(size = 7, face = "bold"),
    panel.grid = element_blank(),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 8),
    plot.margin = margin(10, 100, 10, 10)
  ) +
  coord_cartesian(clip = "off") +
  expand_limits(x = c(min(sigEC_aldex2_Akha$diff.btw) - 1.2, max(sigEC_aldex2_Akha$diff.btw) + 1.2))






# ------------------------------------------------------------
# KO
# ------------------------------------------------------------
KO_aldex2_AkhaCM_AkhaCR <- read_excel("KO_aldex2_AkhaCM_AkhaCR.xlsx")
KO_aldex2_LahuCM_LahuCR <- read_excel("KO_aldex2_LahuCM_LahuCR.xlsx")
KO_aldex2_CM <- read_excel("KO_aldex2_CM.xlsx")
KO_aldex2_CR <- read_excel("KO_aldex2_CR.xlsx")


# Filter wi.eBH < 0.05 (Akha and Lahu groups)
sigKO_aldex2_Akha <- KO_aldex2_AkhaCM_AkhaCR %>% 
  filter(wi.eBH < 0.05)

sigKO_aldex2_Lahu <- KO_aldex2_LahuCM_LahuCR %>% 
  filter(wi.eBH < 0.05)


# Adjust wi.ep for multiple comparisons (CM and CR groups)
sigKO_aldex2_CM <- KO_aldex2_CM %>% 
  mutate(q_value = p.adjust(wi.ep, method = "BH")) %>%
  filter(q_value < 0.05, abs(effect) >= 0.5)

sigKO_aldex2_CR <- KO_aldex2_CR %>% 
  mutate(q_value = p.adjust(wi.ep, method = "BH")) %>%
  filter(q_value < 0.05, abs(effect) >= 0.5)



# Show Level as label and reorder features by diff.btw:
sigKO_aldex2_Akha %>%
  filter(
    wi.eBH < 0.05 &
      abs(effect) >= 1.6) %>%
  mutate(
    Feature = fct_reorder(Feature, diff.btw),
    label_text = round(wi.eBH, 3),
  ) %>%
  ggplot(aes(x = diff.btw, y = Feature)) +
  geom_col(aes(fill = effect)) +
  geom_text(
    aes(
      x = ifelse(diff.btw >= 0, diff.btw + 0.1, diff.btw - 0.1),
      label = label_text,
      hjust = ifelse(diff.btw >= 0, 0, 1)
    ),
    color = "black",
    size = 1.5
  ) +
  scale_fill_gradient2(
    low = "steelblue3",
    mid = "ivory1",
    high = "#FC766A",
    midpoint = 0,
    name = "Effect size"
  ) +
  labs(
    x = "Difference between Groups (diff.btw: ALDEx2)",
    y = "Feature",
    title = "Differential abundance identified by ALDEx2 (wi.eBH < 0.05, abs(effect) >= 1.5 )"
  ) +
  facet_wrap(~ Comparison, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 7),
    title = element_text(size = 7, face = "bold"),
    panel.grid = element_blank(),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 8),
    plot.margin = margin(10, 100, 10, 10)
  ) +
  coord_cartesian(clip = "off") +
  expand_limits(x = c(min(sigKO_aldex2_Akha$diff.btw) - 1.2, max(sigKO_aldex2_Akha$diff.btw) + 1.2))






#---- Summarize ALDEx2 output
# ------------------------------------------------------------
# PW
# ------------------------------------------------------------
# Akha
sigPW_aldex2_Akha %>%
  filter(wi.eBH < 0.05) %>%
  summarise(max_abs_effect = max(abs(effect), na.rm = TRUE)) 

sigPW_aldex2_Akha %>% 
  group_by(Comparison) %>%
  filter(wi.eBH < 0.05,
         abs(effect) >= 1.5) %>%
  count()


# Lahu
sigPW_aldex2_Lahu %>%
  filter(wi.eBH < 0.05) %>%
  summarise(max_abs_effect = max(abs(effect), na.rm = TRUE)) 

sigPW_aldex2_Lahu %>% 
  group_by(Comparison) %>%
  filter(wi.eBH < 0.05,
         abs(effect) >= 1.5) %>%
  count()

# CM
sigPW_aldex2_CM %>%
  filter(wi.eBH < 0.05) %>%
  summarise(max_abs_effect = max(abs(effect), na.rm = TRUE)) 

sigPW_aldex2_CM %>% 
  group_by(Comparison) %>%
  filter(wi.eBH < 0.05,
         abs(effect) >= 1.5) %>%
  count()


# CR
sigPW_aldex2_CR %>%
  filter(wi.eBH < 0.05) %>%
  summarise(max_abs_effect = max(abs(effect), na.rm = TRUE)) 

sigPW_aldex2_CR %>% 
  group_by(Comparison) %>%
  filter(wi.eBH < 0.05,
         abs(effect) >= 1.5) %>%
  count()


# Bind rows
PW_2effect <- bind_rows(
  sigPW_aldex2_Akha %>%
    filter(wi.eBH < 0.05, abs(effect) >= 1) %>%
    count(Comparison, name = "n") %>%
    mutate(Group_factor = "Akha"),
  
  sigPW_aldex2_Lahu %>%
    filter(wi.eBH < 0.05, abs(effect) >= 1) %>%
    count(Comparison, name = "n") %>%
    mutate(Group_factor = "Lahu"),
  
  sigPW_aldex2_CM %>%
    filter(wi.eBH < 0.05, abs(effect) >= 1) %>%
    count(Comparison, name = "n") %>%
    mutate(Group_factor = "CM"),
  
  sigPW_aldex2_CR %>%
    filter(wi.eBH < 0.05, abs(effect) >= 1) %>%
    count(Comparison, name = "n") %>%
    mutate(Group_factor = "CR")
)


group_color <- c(
  "Akha"  = "indianred3",
  "Lahu" = "darkgoldenrod2",
  "CM"  = "deepskyblue4",
  "CR" = "grey40"
)


PW_2effectplot <- ggplot(PW_2effect, aes(x = Comparison, y = n, fill = Group_factor)) +
  geom_col(width = 0.7) +
  coord_flip()+
  scale_fill_manual(values = group_color) +
  labs(
    x = "ALDEx2 Comparison",
    y = "Number of features (KOs)",
    fill = "Group",
    title = "Number of differentially abundant PWs identified by ALDEx2\n(wi.eBH < 0.05, abs(effect) ≥ 1)") +
  theme_bw() +
  theme(
    title = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(angle = 0),
    panel.grid = element_blank()
  )

PW_2effectplot







# ------------------------------------------------------------
# EC
# ------------------------------------------------------------
# Akha
sigEC_aldex2_Akha %>%
  filter(wi.eBH < 0.05) %>%
  summarise(max_abs_effect = max(abs(effect), na.rm = TRUE)) 

sigEC_aldex2_Akha %>% 
  group_by(Comparison) %>%
  filter(wi.eBH < 0.05,
         abs(effect) >= 1.5) %>%
  count()


# Lahu
sigEC_aldex2_Lahu %>%
  filter(wi.eBH < 0.05) %>%
  summarise(max_abs_effect = max(abs(effect), na.rm = TRUE)) 

sigEC_aldex2_Lahu %>% 
  group_by(Comparison) %>%
  filter(wi.eBH < 0.05,
         abs(effect) >= 1.5) %>%
  count()

# CM
sigEC_aldex2_CM %>%
  filter(wi.eBH < 0.05) %>%
  summarise(max_abs_effect = max(abs(effect), na.rm = TRUE)) 

sigEC_aldex2_CM %>% 
  group_by(Comparison) %>%
  filter(wi.eBH < 0.05,
         abs(effect) >= 1.5) %>%
  count()


# CR
sigEC_aldex2_CR %>%
  filter(wi.eBH < 0.05) %>%
  summarise(max_abs_effect = max(abs(effect), na.rm = TRUE)) 

sigEC_aldex2_CR %>% 
  group_by(Comparison) %>%
  filter(wi.eBH < 0.05,
         abs(effect) >= 1.5) %>%
  count()



# Bind rows
EC_2effect <- bind_rows(
  sigEC_aldex2_Akha %>%
    filter(wi.eBH < 0.05, abs(effect) >= 1) %>%
    count(Comparison, name = "n") %>%
    mutate(Group_factor = "Akha"),
  
  sigEC_aldex2_Lahu %>%
    filter(wi.eBH < 0.05, abs(effect) >= 1) %>%
    count(Comparison, name = "n") %>%
    mutate(Group_factor = "Lahu"),
  
  sigEC_aldex2_CM %>%
    filter(wi.eBH < 0.05, abs(effect) >= 1) %>%
    count(Comparison, name = "n") %>%
    mutate(Group_factor = "CM"),
  
  sigEC_aldex2_CR %>%
    filter(wi.eBH < 0.05, abs(effect) >= 1) %>%
    count(Comparison, name = "n") %>%
    mutate(Group_factor = "CR")
)


group_color <- c(
  "Akha"  = "indianred3",
  "Lahu" = "darkgoldenrod2",
  "CM"  = "deepskyblue4",
  "CR" = "grey40"
)


EC_2effectplot <- ggplot(EC_2effect, aes(x = Comparison, y = n, fill = Group_factor)) +
  geom_col(width = 0.7) +
  coord_flip()+
  scale_fill_manual(values = group_color) +
  labs(
    x = "ALDEx2 Comparison",
    y = "Number of features (KOs)",
    fill = "Group",
    title = "Number of differentially abundant ECs identified by ALDEx2\n(wi.eBH < 0.05, abs(effect) ≥ 1)") +
  theme_bw() +
  theme(
    title = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(angle = 0),
    panel.grid = element_blank()
  )

EC_2effectplot









# ------------------------------------------------------------
# KO
# ------------------------------------------------------------
# Akha
sigKO_aldex2_Akha %>%
  filter(wi.eBH < 0.05) %>%
  summarise(max_abs_effect = max(abs(effect), na.rm = TRUE)) 

sigKO_aldex2_Akha %>% 
  group_by(Comparison) %>%
  filter(wi.eBH < 0.05,
         abs(effect) >= 1.5) %>%
  count()


# Lahu
sigKO_aldex2_Lahu %>%
  filter(wi.eBH < 0.05) %>%
  summarise(max_abs_effect = max(abs(effect), na.rm = TRUE)) 

sigKO_aldex2_Lahu %>% 
  group_by(Comparison) %>%
  filter(wi.eBH < 0.05,
         abs(effect) >= 1.5) %>%
  count()

# CM
sigKO_aldex2_CM %>%
  filter(wi.eBH < 0.05) %>%
  summarise(max_abs_effect = max(abs(effect), na.rm = TRUE)) 

sigKO_aldex2_CM %>% 
  group_by(Comparison) %>%
  filter(wi.eBH < 0.05,
         abs(effect) >= 1.5) %>%
  count()


# CR
sigKO_aldex2_CR %>%
  filter(wi.eBH < 0.05) %>%
  summarise(max_abs_effect = max(abs(effect), na.rm = TRUE)) 

sigKO_aldex2_CR %>% 
  group_by(Comparison) %>%
  filter(wi.eBH < 0.05,
         abs(effect) >= 1.5) %>%
  count()



# Bind rows
KO_2effect <- bind_rows(
  sigKO_aldex2_Akha %>%
    filter(wi.eBH < 0.05, abs(effect) >= 1) %>%
    count(Comparison, name = "n") %>%
    mutate(Group_factor = "Akha"),
  
  sigKO_aldex2_Lahu %>%
    filter(wi.eBH < 0.05, abs(effect) >= 1) %>%
    count(Comparison, name = "n") %>%
    mutate(Group_factor = "Lahu"),
  
  sigKO_aldex2_CM %>%
    filter(wi.eBH < 0.05, abs(effect) >= 1) %>%
    count(Comparison, name = "n") %>%
    mutate(Group_factor = "CM"),
  
  sigKO_aldex2_CR %>%
    filter(wi.eBH < 0.05, abs(effect) >= 1) %>%
    count(Comparison, name = "n") %>%
    mutate(Group_factor = "CR")
)


group_color <- c(
  "Akha"  = "indianred3",
  "Lahu" = "darkgoldenrod2",
  "CM"  = "deepskyblue4",
  "CR" = "grey40"
)


KO_2effectplot <- ggplot(KO_2effect, aes(x = Comparison, y = n, fill = Group_factor)) +
  geom_col(width = 0.7) +
  coord_flip()+
  scale_fill_manual(values = group_color) +
  labs(
    x = "ALDEx2 Comparison",
    y = "Number of features (KOs)",
    fill = "Group",
    title = "Number of differentially abundant KOs identified by ALDEx2\n(wi.eBH < 0.05, abs(effect) ≥ 1)") +
  theme_bw() +
  theme(
    title = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(angle = 0),
    panel.grid = element_blank()
  )

KO_2effectplot






# CM and CR shared significantly different features
# ------------------------------------------------------------
# CM
# ------------------------------------------------------------
# Max effect size Feature
sigKO_aldex2_CM %>%
  filter(wi.eBH < 0.05) %>%
  arrange(desc(abs(effect))) %>%
  slice(1:5) %>%
  left_join(KO_absolute, by = c("Feature" = "KO_ID")) %>%
  select(Comparison, Feature, wi.eBH, effect, description) %>%
  as.data.frame()


# Max effect size Feature
sigKO_aldex2_CM %>%
  filter(wi.eBH < 0.05) %>%
  arrange(desc(abs(effect))) %>%
  slice(1:5) %>%
  as.data.frame()


# 1 Is there any duplicated (shared) Features across comparisons
sigPW_aldex2_CM %>% 
  group_by(Comparison, Feature) %>%
  filter(
    wi.eBH < 0.05,
    abs(effect) >= 1
  ) %>%
  count(Feature, sort = TRUE) 

sigEC_aldex2_CM %>% 
  group_by(Comparison, Feature) %>%
  filter(
    wi.eBH < 0.05,
    abs(effect) >= 1
  ) %>%
  count(Feature, sort = TRUE) 

sigKO_aldex2_CM %>% 
  group_by(Comparison, Feature) %>%
  filter(
    q_value < 0.05,
    abs(effect) >= 1
  ) %>%
  count(Feature, sort = TRUE) 



# ------------------------------------------------------------
# CR
# ------------------------------------------------------------
# Max effect size Feature
sigKO_aldex2_CR %>%
  filter(wi.eBH < 0.05) %>%
  arrange(desc(abs(effect))) %>%
  slice(1:5) %>%
  left_join(KO_absolute, by = c("Feature" = "KO_ID")) %>%
  select(Comparison, Feature, wi.eBH, effect, description) %>%
  as.data.frame()


# Max effect size Feature
sigKO_aldex2_CR %>%
  filter(wi.eBH < 0.05) %>%
  arrange(desc(abs(effect))) %>%
  slice(1:5) %>%
  as.data.frame()


# 1 Is there any duplicated (shared) Features across comparisons
sigPW_aldex2_CR %>% 
  group_by(Comparison, Feature) %>%
  filter(
    wi.eBH < 0.05,
    abs(effect) >= 1
  ) %>%
  count(Feature, sort = TRUE) 

sigEC_aldex2_CR %>% 
  group_by(Comparison, Feature) %>%
  filter(
    wi.eBH < 0.05,
    abs(effect) >= 1
  ) %>%
  count(Feature, sort = TRUE) 

sigKO_aldex2_CR %>% 
  group_by(Comparison, Feature) %>%
  filter(
    wi.eBH < 0.05,
    abs(effect) >= 1
  ) %>%
  count(Feature, sort = TRUE) 
