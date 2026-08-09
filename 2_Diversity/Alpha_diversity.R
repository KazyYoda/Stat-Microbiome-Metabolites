############################################################
# Alpha Diversity Analysis (Gut Microbiota)
############################################################

# This script estimates alpha diversity of the gut microbiota,
# prepares sample metadata, performs statistical comparisons
# between geographic areas and hill tribe groups, visualizes
# alpha diversity distributions, and exports the results.
#
# Alpha diversity metrics included:
#   - Observed features: number of observed microbial features
#   - Shannon: diversity incorporating richness and evenness
#   - Chao1: estimated species/features richness


# -----------------------------------------------
# Set working directory and load phyloseq object
# -----------------------------------------------

# Set the working directory for alpha diversity analyses.
setwd("~/Documents/HillTribe_NGS/2.Diversity")

# Load the previously constructed phyloseq object and associated
# objects saved from the phyloseq-building workflow.
load("~/Documents/HillTribe_NGS/1.Raw_data/Building_phyloseq.RData")


# ------------------------------
# Load required packages
# ------------------------------

# Load packages required for microbiome data handling,
# data manipulation, statistical analysis, and file export.
library(phyloseq)
library(dplyr)
library(rio)
library(car)


# ------------------------------
# 1. Alpha Diversity Estimation
# ------------------------------

# Calculate standard alpha diversity metrics from the rarefied
# phyloseq object.
#
# Observed  = number of observed microbial features
# Shannon   = diversity index accounting for richness and evenness
# Chao1     = estimator of feature richness
alpha_div <- estimate_richness(
  ps_rare,
  measures = c("Observed", "Shannon", "Chao1")
)

# Inspect the calculated alpha diversity metrics.
head(alpha_div)


# ------------------------------
# 2. Metadata Preparation
# ------------------------------

# Extract sample metadata from the rarefied phyloseq object.
meta <- as(sample_data(ps_rare), "data.frame")

# Create a combined group variable from ethnicity and study area.
# The resulting labels are used for group-wise comparisons
# and visualization.
meta <- meta %>%
  mutate(Group = paste(Ethnicity, Area, sep = "-"))

# Define Area as a factor with a fixed order for analysis
# and visualization.
meta$Area <- factor(
  meta$Area,
  levels = c("CM", "CR")
)

# Define Ethnicity as a factor with a fixed order for analysis
# and visualization.
meta$Ethnicity <- factor(
  meta$Ethnicity,
  levels = c("Akha", "Lahu", "Lisu", "Khuen")
)


# -----------------------------------------
# 3. Merge Alpha Diversity with Metadata
# -----------------------------------------

# Combine sample metadata with the corresponding alpha diversity
# measurements.
alpha_div_meta <- cbind(meta, alpha_div)

# Inspect the combined dataset.
head(alpha_div_meta)

# Rename the Observed metric for clarity and remove the standard
# error associated with the Chao1 estimator, which is not used
# in the downstream analyses.
alpha_div_meta <- alpha_div_meta %>%
  rename(Observed_features = Observed) %>%
  select(-se.chao1)


# ------------------------------
# 4. Export Alpha Diversity Table
# ------------------------------

# Export the complete alpha diversity table together with
# sample metadata for record keeping and downstream analysis.
Export(alpha_div_meta, "Alpha_diversity.txt")



# ----------------------------------------------------
# 5. Statistical Testing: Within-Ethnicity Comparison
# ----------------------------------------------------

# Objective:
# Compare alpha diversity between the two study areas (CM vs CR)
# separately within each hill tribe using the Wilcoxon rank-sum test.
#
# The analysis is performed independently for each alpha diversity
# metric using the custom alpha_wilcox() function.


# ------------------------------
# Akha group
# ------------------------------

# Check the number of samples available in each area for the Akha group.
alpha_div_meta %>% 
  filter(Ethnicity == "Akha") %>%
  count(Area)

# Perform Wilcoxon rank-sum tests comparing CM and CR within Akha.
alpha_wilcox(alpha_div_meta, 
             ethnicity_group = "Akha",
             factor_col = "Area",
             metric = "Alpha Diversity") 


# ------------------------------
# Lahu group
# ------------------------------

# Check the number of samples available in each area for the Lahu group.
alpha_div_meta %>% 
  filter(Ethnicity == "Lahu") %>%
  count(Area)

# Perform Wilcoxon rank-sum tests comparing CM and CR within Lahu.
alpha_wilcox(alpha_div_meta, 
             ethnicity_group = "Lahu",
             factor_col = "Area",
             metric = "Alpha Diversity") 



# -----------------------------------------------
# 6. Visualization: Alpha Diversity Boxplots
# -----------------------------------------------

# Visualize alpha diversity metrics for Akha samples,
# comparing the two study areas.
alpha_boxplot(
  data = alpha_div_meta %>% filter(Ethnicity == "Akha"),
  group_col = "Group",
  metrics = c("Observed_features", "Chao1", "Shannon"),
  fill_colors = c(
    "Akha-CM" = "gold1",
    "Akha-CR" = "dodgerblue4"
  )
)



# Visualize alpha diversity metrics for Lahu samples,
# comparing the two study areas.
alpha_boxplot(
  data = alpha_div_meta %>% filter(Ethnicity == "Lahu"),
  group_col = "Group",
  metrics = c("Observed_features", "Chao1", "Shannon"),
  fill_colors = c(
    "Lahu-CM" = "lightblue2",
    "Lahu-CR" = "coral1"
  )
)





# ---------------------------------------------------
# 7. Statistical Testing: Between-Group Differences
# ---------------------------------------------------

# Objective:
# Compare alpha diversity among hill tribe groups separately
# within each study area (CM and CR).
#
# A Kruskal–Wallis test is used to evaluate overall differences
# among the groups because alpha diversity measures are analyzed
# using non-parametric methods.


# ------------------------------
# CM Area
# ------------------------------

# Check the number of samples available in each group within CM.
alpha_div_meta %>%
  filter(Area == "CM") %>%
  count(Group)

# Perform Kruskal–Wallis tests across hill tribe groups within CM.
alpha_kruskal(
  alpha_div_meta,
  group_col = "Area",
  group_value = "CM",
  factor_col = "Group",
  factor_levels = c("Akha-CM", "Lahu-CM", "Khuen-CM")
)


# Visualize alpha diversity across hill tribe groups within CM.
alpha_boxplot(
  data = alpha_div_meta %>% filter(Area == "CM"),
  group_col = "Group",
  metrics = c("Observed_features", "Chao1", "Shannon"),
  fill_colors = c(
    "Akha-CM" = "gold1",
    "Lahu-CM" = "lightblue2",
    "Khuen-CM" = "darkorange"
  )
)



# ------------------------------
# CR Area
# ------------------------------

# Check the number of samples available in each group within CR.
alpha_div_meta %>%
  filter(Area == "CR") %>%
  count(Group)

# Perform Kruskal–Wallis tests across hill tribe groups within CR.
alpha_kruskal(
  alpha_div_meta,
  group_col = "Area",
  group_value = "CR",
  factor_col = "Group",
  factor_levels = c("Akha-CR", "Lahu-CR", "Lisu-CR")
)

# Visualize alpha diversity across hill tribe groups within CR.
alpha_boxplot(
  data = alpha_div_meta %>% filter(Area == "CR"),
  group_col = "Group",
  metrics = c("Observed_features", "Chao1", "Shannon"),
  fill_colors = c(
    "Akha-CR" = "dodgerblue4",
    "Lahu-CR" = "coral1",
    "Lisu-CR" = "grey40"
    )
)



# ---------------------------------------------------
# 8. Post hoc Analysis: Dunn's Test
# ---------------------------------------------------

# Perform pairwise Dunn's tests only for alpha diversity metrics
# with a significant overall Kruskal–Wallis test.
#
# Benjamini–Hochberg (BH) adjustment is applied to the pairwise
# comparisons within each Dunn's test.


# -----------------------------------------------
# Identify Significant Kruskal–Wallis Results (CR)
# -----------------------------------------------

# Import the Kruskal–Wallis results for the CR analysis.
CR_krus <- read.delim("~/Documents/HillTribe_NGS/2.Diversity/kruskal_CR.txt")

# Retain alpha diversity metrics with a statistically significant
# Kruskal–Wallis test (P < 0.05) for subsequent post hoc testing.
sig_CR_krus <- CR_krus %>%
  filter(P_Value < 0.05)

# Inspect the significant Kruskal–Wallis results.
print(sig_CR_krus)


# ------------------------------------
# Dunn's Post hoc Test (BH-adjusted)
# ------------------------------------

# Run Dunn's pairwise test for each alpha diversity metric that
# showed a significant Kruskal–Wallis result.
#
# The Benjamini–Hochberg method is used to adjust the pairwise
# P-values for multiple comparisons within each metric.
library(FSA)

dunn_test <- setNames(
  lapply(sig_CR_krus$Metric, function(metric) {
    
    dunn_res <- dunnTest(
      as.formula(paste(metric, "~ Group")),
      data   = alpha_div_meta %>% filter(Area == "CR"),
      method = "bh"
    )
    
    as.data.frame(dunn_res$res)
  }),
  sig_CR_krus$Metric
)

# Inspect the Dunn's test results.
dunn_test


# ------------------------------------------
# Extract Significant Pairwise Comparisons
# ------------------------------------------

# Combine the Dunn's test results for all significant metrics
# into a single data frame and retain only pairwise comparisons
# with BH-adjusted P < 0.05.
dunn_sig <- do.call(
  rbind,
  lapply(names(dunn_test), function(Metric) {
    df <- dunn_test[[Metric]]. # -> Extract that metric's Dunn result table
    df$Metric <- Metric # -> Add a Metric column
    df <- df[df$P.adj < 0.05, ] # -> Keep only significant comparisons (P.adj < 0.05)
    
    if (nrow(df) > 0) {
      return(df)
    } else {
      return(NULL) # -> NULL if no significant comparisons
    }
  })
)


# Convert the combined results into a final data frame with
# standard row numbering.
dunn_sig_alpha <- data.frame(dunn_sig, row.names = NULL)

# Inspect the final table of significant pairwise comparisons.
print(dunn_sig_alpha)


# ------------------------------
# Export Results
# ------------------------------

# Export significant Dunn's test results for the CR analysis.
Export(dunn_sig_alpha, "1_dunn_sig_alpha.txt")
