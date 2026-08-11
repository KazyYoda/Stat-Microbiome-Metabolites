############################################################
# Metabolite Analysis
############################################################

# Set the working directory and load the metabolite analysis
# workspace containing the required input objects.
setwd("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat")
load("~/Documents/HillTribe_NGS/5.Metabolites/Metabolites.RData")

# Load required packages
library(dplyr)
library(tidyr)
library(car)
library(readxl)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggh4x)


#--------------------------------------------------------
# Visualization Settings
#--------------------------------------------------------

# Define colors for the geographic groups used in the
# metabolite comparison.
group_colors <- c(
  "Akha-CM" = "gold1",
  "Akha-CR" = "dodgerblue4"
)


#========================================================
# 4. Compute log2 Fold-Change (log2FC)
#    Relative to Defined Contrast Group
#========================================================

# This section calculates metabolite log2 fold-changes
# relative to a predefined contrast group and generates
# summary tables for downstream statistical interpretation
# and visualization.


#--------------------------------------------------------
# 4.1 Merge Metadata with Metabolite Matrix
#--------------------------------------------------------

# Confirm that the sample order in the metabolite matrix
# matches the sample order in the corresponding metadata.
# This check helps prevent incorrect sample-to-metadata
# assignments during data integration.
identical(rownames(metabo_05_t), Met05_match$Sample_ID)

# Combine sample metadata with the metabolite abundance matrix.
Met05 <- cbind(Met05_match, metabo_05_t)

# Apply a log2 transformation to metabolite abundance values.
# Columns 14 onward contain the metabolite measurements.
Met05_log2 <- Met05 %>%
  mutate(across(14:ncol(.), log2))



#--------------------------------------------------------
# 4.2 Subset Data by Ethnicity and Area
#--------------------------------------------------------

# Retain samples from the Akha ethnicity for the
# Akha-CM versus Akha-CR comparison.
Met05_Akha <- Met05_log2 %>% filter(Ethnicity == "Akha")


# Identify metabolite columns for downstream summarization.
met_cols <- grep("^Met", names(Met05_Akha), value = TRUE)

# Convert the metabolite data to long format and calculate
# the mean and standard deviation for each metabolite within
# each geographic group.
summary_Akha_long <- Met05_Akha %>%
  group_by(Group) %>%
  summarise(
    across(
      all_of(met_cols),
      list(
        Mean = ~mean(.x, na.rm = TRUE),
        SD   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Group,
    names_to = c("Metabolite", ".value"),
    names_sep = "_"
  ) 



# Wide format
# Format the group-specific mean ± SD values into a wide
# summary table and merge metabolite annotations and
# statistical results.
summary_Akha_wide <- summary_Akha %>%
  mutate(
    Mean_SD = sprintf("%.2f ± %.2f", Mean, SD)
  ) %>%
  select(Metabolite, Group, Mean_SD) %>%
  pivot_wider(
    names_from = Group,
    values_from = Mean_SD
  ) %>%
  left_join(
    metabo_05_descp,
    by = c("Metabolite" = "Code")
  ) %>%
  select(-Compound) %>%
  left_join(
    wilcox_metabo_Akha,
    by = "Metabolite"
  ) %>% 
  rename(
    q_value = p_adj
  ) %>%
  mutate(`-log10(q-value)` = -log10(q_value),
         Test = "wilcox") %>%
  left_join(
    log2FC_sum_Met05_Akha %>%
      select(Metabolite, log2FC, Direction),
    by = "Metabolite"
  ) %>%
  mutate(
    across(
      Comparison:Direction,
      ~ ifelse(is.na(P_Value), NA, .)
    )
  )

# Export the descriptive and statistical summary table.
Export(summary_Akha_wide, "summary_Akha_wide.txt")


# Count metabolites showing an absolute log2FC of at least 1
# and an adjusted p-value below 0.05, summarized by direction.
summary_Akha_wide %>%
  filter(abs(log2FC) >=1 & q_value < 0.05) %>%
  group_by(Direction) %>%
  count()


# Count metabolite-group combinations for which the mean
# abundance is equal to zero. These may represent
# non-detected or zero-abundance measurements.
summary_Akha_long %>%
  group_by(Group) %>%
  summarise(
    Not_Detected = sum(Mean == 0, na.rm = TRUE)
  )


#========================================================
# 4.3 Log2FC Analysis — Akha (Akha-CM vs Akha-CR)
#========================================================

#--------------------------------------------------------
# 4.3.1 Compute log2FC Relative to Contrast Group
#--------------------------------------------------------

# Calculate metabolite log2 fold-changes using Akha-CR as
# the reference (contrast) group.
#
# Positive log2FC values indicate higher metabolite abundance
# in the comparison group relative to Akha-CR, whereas
# negative values indicate lower abundance.
log2FC_Met05_Akha <- log2FC_group_summary_relative_to_contrast(
  metabo_data = Met05_Akha[-c(1:13)], 
  metadata = Met05_Akha[1:13], 
  contrast = "Akha-CR", 
  Group = "Group"
)

# Transpose the log2FC matrix and retain metabolite codes
# as an explicit column for downstream data manipulation.
log2FC_Met05_Akha_t <- log2FC_Met05_Akha %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Code")


#--------------------------------------------------------
# 4.3.2 Summarize log2FC by Group
#--------------------------------------------------------

# Summarize metabolite-level log2FC values while using
# Akha-CR as the reference group.
log2FC_sum_Met05_Akha <- summarize_log2FC(
  log2FC_Met05_Akha, 
  filter_group = "Akha-CR"
)

# Add metabolite annotation information to the log2FC summary.
log2FC_sum_Met05_Akha <- log2FC_sum_Met05_Akha %>%
  left_join(
    metabo_05_descp, 
    by = c("Metabolite" = "Code")
  )

# Export the metabolite log2FC summary table.
Export(log2FC_sum_Met05_Akha, "log2FC_sum_Met05_Akha.txt")



#========================================================
# 5. Log2FC Visualization — Barplot
#========================================================

# Generate a grouped barplot showing metabolite log2FC
# values for the Akha-CM versus Akha-CR comparison.
ggFC_Met05_Akha <- plot_log2FC_bar(
  summary_df = log2FC_sum_Met05_Akha,
  title = "Mean Log2FC of metabolites (Akha-CM vs. Akha-CR)",
  xlab = "Metabolite",
  fill_colors = c("Upregulated" = "steelblue3",
                  "Downregulated" = "#FC766A"),
  strip_text_color = "black",
  log2FC_cutoff = 3,
  group_color = group_colors
)

# Display the log2FC barplot.
print(ggFC_Met05_Akha)



#--------------------------------------------------------
# 5.1 Export Complete Akha log2FC Results
#--------------------------------------------------------

# Merge log2FC results with Wilcoxon test results to create
# a combined table containing effect sizes, statistical
# significance, direction, and metabolite annotations.
log2FC_Akha_all <- log2FC_sum_Met05_Akha %>% 
  left_join(wilcox_metabo_Akha,
            by = "Metabolite") %>%
  select(Group, Metabolite, Compound_Clean, log2FC, P_Value, p_adj, Direction) %>%
  filter(!is.na(P_Value)) %>%
  mutate(`-log10(q-value)` = -log10(p_adj),
         Test = "wilcox",
         Comparison = "Akha-CM vs AKha-CR")

# Export the complete Akha log2FC results.
Export(log2FC_Akha_all, "log2FC_Akha_all.txt")



#--------------------------------------------------------
# 5.2 Visualization of Strongly Differential Metabolites
#--------------------------------------------------------

# Reload the Wilcoxon test results used to annotate the
# metabolite log2FC results with adjusted p-values.
wilcox_metabo_Akha <- read_delim("wilcox_metabo_Akha.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)


# Select metabolites with large absolute log2FC values and
# combine their effect sizes with statistical significance
# and metabolite annotations.
log2FC_Akha_sig <- log2FC_sum_Met05_Akha_sig %>% 
  left_join(wilcox_metabo_Akha,
            by = "Metabolite") %>%
  select(Group, Metabolite, log2FC, p_adj, Direction, Compound_Clean, Subclass) %>%
  filter(abs(log2FC) >= 16) %>%
  mutate(`-log10(q-value)` = -log10(p_adj))


# Load packages required for the dot plot.
library(ggplot2)
library(dplyr)
library(stringr)

# Define colors according to the direction of metabolite
# change relative to the reference group.
fill_colors <- c(
  "Downregulated" = "steelblue3",
  "Upregulated" = "#FC766A"
)

# Generate a dot plot for metabolites with large absolute
# log2FC values. Point size represents statistical
# significance, while fill indicates the direction of change.
log2FC_Akha_sig %>%
  mutate(
    Metabolite = reorder(Metabolite, log2FC)
  ) %>%
  ggplot(aes(
    x = Metabolite,
    y = log2FC
  )) +
  
  # Circle size represents -log10(q-value), and fill color
  # represents the direction of metabolite change.
  geom_point(aes(
    size = `-log10(q-value)`,
    fill = Direction
  ),
  shape = 21,
  color = "grey30",
  stroke = 0.5
  ) +
  
  coord_flip() +
  
  scale_fill_manual(values = fill_colors) +
  
  scale_size_continuous(
    name = expression(-log[10](q-value)),
    range = c(1, 3)
  ) +
  
  labs(
    x = "Metabolite",
    y = "Mean Log2FC",
    fill = "Direction",
    title = "Differentially enriched metabolites\n(Akha-CM vs Akha-CR)"
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 4),
    axis.text.x = element_text(size = 6),
    axis.title = element_text(size = 7),
    plot.title = element_text(size = 6, face = "bold"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    panel.grid = element_blank()
  )



#--------------------------------------------------------
# Log2FC Interpretation Guide
#--------------------------------------------------------

# Log2FC values are interpreted relative to the defined
# contrast/reference group (Akha-CR).
#
# Upregulation:
# log2FC ≥ 1  → ≥ 2-fold increase
# log2FC ≥ 2  → ≥ 4-fold increase
# log2FC ≥ 3  → ≥ 8-fold increase
# log2FC ≥ 4  → ≥ 16-fold increase
#
# Fold change = 2^(log2FC)
#
# Downregulation:
# log2FC ≤ −1 → ≤ 0.5-fold
# log2FC ≤ −2 → ≤ 0.25-fold
# log2FC ≤ −3 → ≤ 0.125-fold
# log2FC ≤ −4 → ≤ 0.0625-fold





#--------------------------------------------------------
# Volcano Plot
#--------------------------------------------------------

# Load packages required for data manipulation,
# visualization, and automatic feature labeling.
library(dplyr)
library(ggplot2)
library(ggrepel)

#--------------------------------------------------------
# Prepare Data for Volcano Plot
#--------------------------------------------------------

# Combine metabolite log2FC results with Wilcoxon test
# results using the metabolite identifier as the key.
vol_Akha <- log2FC_sum_Met05_Akha %>%
  left_join(Akha_wilcox, by = "Metabolite") 

# Calculate -log10(adjusted p-value) for visualization
# and classify metabolites according to statistical
# significance and direction of change.
vol_Akha_plot <- vol_Akha %>%
  mutate(
    logP = -log10(p_adj),   # or we.eBH depending on your choice
    class = case_when(
      p_adj < 0.05 ~ "Significant",
      TRUE ~ "Not Significant"
    )
  ) %>%
  filter(
    !is.na(Comparison)
  ) %>%
  mutate(
    Significant = case_when(
      class == "Significant" & Direction == "Downregulated" ~ "Downregulated",
      class == "Significant" & Direction == "Upregulated" ~ "Upregulated",
      class == "Not Significant" ~ "Not Significant"
    )
  )


#--------------------------------------------------------
# Select Features for Labeling
#--------------------------------------------------------

# Select the 30 significant metabolites with the smallest
# adjusted p-values for labeling on the volcano plot.
top_features <- vol_Akha_plot %>%
  filter(class == "Significant") %>%
  arrange(p_adj) %>%
  slice(1:30)


#--------------------------------------------------------
# Generate Volcano Plot
#--------------------------------------------------------

# Plot metabolite log2FC against -log10(adjusted p-value).
# Positive and negative log2FC values represent higher and
# lower metabolite abundance, respectively, relative to
# the defined contrast group.
ggplot(vol_Akha_plot, aes(x = log2FC, y = logP)) +
  geom_point(aes(color = Significant), alpha = 0.7, size = 1) +
  
  # Highlight the selected top significant metabolites.
  geom_point(
    data = top_features,
    aes(color = Significant),
    size = 1.5
  ) +
  
  # Label the selected top significant metabolites.
  geom_text_repel(
    data = top_features,
    aes(label = Metabolite),
    size = 1,
    max.overlaps = 100,
    segment.size = 0.2
  ) +
  
  # Assign colors according to statistical significance
  # and direction of metabolite change.
  scale_color_manual(values = c("Upregulated" = "#FC766A",
                                "Downregulated" = "#89ABE3",
                                "Not Significant" = "grey70")) +
  
  # Vertical lines indicate the predefined absolute
  # log2FC threshold of 1 (2-fold change).
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "darkgray") + # abs(log2FC) >=1
  
  # Horizontal line indicates the adjusted p-value
  # significance threshold of 0.05.
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "darkgray") + # significance cutoff
  
  labs(
    title = "Volcano plot of metabolites (Akha-CM vs Akha-CR)",
    x = "Mean Log2FC",
    y = "-log10(q-value)"
  ) +
  
  theme_bw()+
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    axis.text = element_text(size = 6),
    title = element_text(size = 6)
  )







#========================================================
# 6. Heatmap Visualization (Log2-Transformed Data)
#========================================================

# Generate a heatmap of selected metabolites using
# log2-transformed metabolite abundance values.
# Metabolites are displayed as rows and samples as columns.


#--------------------------------------------------------
# 6.1 Prepare Log2 Matrix
#--------------------------------------------------------

# Log2-transform the complete metabolite abundance matrix.
# The transformation is applied before heatmap visualization
# to reduce the influence of large abundance differences.
log_matrix_metabo_05 <- as.matrix(log2(metabo_05_t))

# Transpose the matrix so that metabolites are represented
# as rows and samples are represented as columns.
log_matrix_metabo_05_t <- t(log_matrix_metabo_05)



#--------------------------------------------------------
# 6.2 Select Significant Metabolites
#--------------------------------------------------------

# Select metabolites based on the predefined absolute log2FC
# threshold for inclusion in the heatmap.
log2FC_X_Akha <- log2FC_sum_Met05_Akha_sig %>%
  filter(abs(log2FC) >= 10)

# Identify metabolites present in both the log2-transformed
# abundance matrix and the selected metabolite list.
valid_rows <- intersect(
  rownames(log_matrix_metabo_05_t),
  log2FC_X_Akha$Metabolite
)

# Identify Akha samples present in both the metabolite matrix
# and the corresponding metadata.
valid_cols <- intersect(
  colnames(log_matrix_metabo_05_t),
  Met05_Akha$Sample_ID
)

# Subset the abundance matrix to the selected metabolites
# and Akha samples while preserving the matrix structure.
log_matrix_metabo_05_t_Akha <-
  log_matrix_metabo_05_t[
    valid_rows,
    valid_cols,
    drop = FALSE
  ]



#--------------------------------------------------------
# 6.3 Variance Check (Pre-Clustering QC)
#--------------------------------------------------------

# Check for metabolites with zero variance across samples.
# Zero-variance features provide no information for clustering.
sum(apply(log_matrix_metabo_05_t_Akha, 1, sd, na.rm = TRUE) == 0)

# Check for samples with zero variance across the selected
# metabolites.
sum(apply(log_matrix_metabo_05_t_Akha, 2, sd, na.rm = TRUE) == 0)

# Note: zero-variance features or samples can affect
# distance-based clustering, including Euclidean clustering.



#--------------------------------------------------------
# 6.4 Define Group Order for Annotation
#--------------------------------------------------------

# Set the order of the geographic groups for consistent
# display of sample annotations in the heatmap.
Met05_Akha$Group <- factor(
  Met05_Akha$Group,
  levels = c("Akha-CM", "Akha-CR")
)



#========================================================
# 6.5 Generate Heatmap
#========================================================

# Generate the heatmap using the selected log2-transformed
# metabolite matrix and corresponding sample metadata.
Akha <- plot_heatmap(
  mat = log_matrix_metabo_05_t_Akha,
  metadata = Met05_Akha[1:12],
  legend_title = "Log2 (Akha)",
  group_color = group_colors
)

# Draw the heatmap with the heatmap and annotation legends
# positioned on the left side.
draw(
  Akha,
  heatmap_legend_side = "left",
  annotation_legend_side = "left"
)






#========================================================
# 7. PCA Analysis — Akha Metabolites
#========================================================

library(FactoMineR)
library(factoextra)

#--------------------------------------------------------
# 7.1 Select Differential Metabolites (|log2FC| ≥ 1)
#--------------------------------------------------------

log2FC_1_Akha <- log2FC_sum_Met05_Akha_sig %>%
  filter(abs(log2FC) >= 1)


#--------------------------------------------------------
# 7.2 Subset Log2 Matrix for PCA
#--------------------------------------------------------

# Match Akha sample IDs
valid_cols <- intersect(
  colnames(log_matrix_metabo_05_t),
  Met05_Akha$Sample_ID
)

# Match selected metabolite rows
valid_rows <- intersect(
  rownames(log_matrix_metabo_05_t),
  log2FC_1_Akha$Metabolite
)

# Subset matrix (metabolites × samples)
log_matrix_metabo_05_t_Akha_pca <-
  log_matrix_metabo_05_t[
    valid_rows,
    valid_cols,
    drop = FALSE
  ]


#--------------------------------------------------------
# 7.3 Scale Data Prior to PCA
#--------------------------------------------------------

# Standardize metabolites (mean = 0, SD = 1)
pca_Akha <- scale(log_matrix_metabo_05_t_Akha_pca)
head(pca_Akha)
sum(is.nan(pca_Akha))

# If there is Nan then replace NaN values (from zero SD features) with 0
pca_Akha[is.nan(pca_Akha)] <- 0


#--------------------------------------------------------
# 7.4 Perform PCA
#--------------------------------------------------------

# Transpose: samples as rows for PCA
res_pca <- PCA(
  t(pca_Akha),
  scale.unit = FALSE,  # already scaled
  graph = FALSE
)


#--------------------------------------------------------
# 7.5 Prepare Metadata for Visualization
#--------------------------------------------------------

sample_names <- colnames(log_matrix_metabo_05_t_Akha_pca)
group_labels <- Met05_Akha$Group

metadata <- data.frame(
  Sample = sample_names,
  Group = group_labels
)

rownames(metadata) <- sample_names


#--------------------------------------------------------
# 7.6 PCA Plot with Group Ellipses
#--------------------------------------------------------

# Define group color palette
group_colors <- c(
  "Akha-CM" = "gold1",
  "Akha-CR" = "dodgerblue4"
)

fviz_pca_ind(
  res_pca,
  geom.ind = "point",
  col.ind = metadata$Group,     # Color by group
  palette = group_colors,
  addEllipses = TRUE,           # Add 95% confidence ellipses
  ellipse.type = "confidence",
  ellipse.level = 0.95,
  repel = TRUE                  # Avoid overlapping labels
) +
  labs(title = "PCA of Akha Metabolites (|Log2FC| >= 1)") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 6),
    title = element_text(size = 6)
  )







#========================================================
# 8. Normality and Homogeneity of Variance Testing
#========================================================

# Assess distributional characteristics of the metabolite
# data prior to selecting or interpreting statistical tests.
# Normality and homogeneity of variance are evaluated for
# the Akha metabolite dataset.

setwd("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat")


#--------------------------------------------------------
# 8.1 Subset Akha Metabolite Matrix
#--------------------------------------------------------

# Identify metabolite columns that are present in both the
# log2-transformed metabolite matrix and the Akha log2FC results.
valid_cols <- intersect(
  colnames(log_matrix_metabo_05),
  log2FC_sum_Met05_Akha$Metabolite
)

# Identify Akha samples that are present in both the
# metabolite matrix and the corresponding sample metadata.
valid_rows <- intersect(
  rownames(log_matrix_metabo_05),
  Met05_Akha$Sample_ID
)

# Subset the log2-transformed metabolite matrix to Akha
# samples and the metabolites included in the analysis.
# The resulting matrix contains samples as rows and
# metabolites as columns.
log_matrix_metabo_05_Akha <-
  log_matrix_metabo_05[
    valid_rows,
    valid_cols,
    drop = FALSE
  ]


#--------------------------------------------------------
# 8.2 Identify and Remove Zero-Variance Metabolites
#--------------------------------------------------------

# Count metabolites with zero variance across Akha samples.
# Zero-variance features cannot be evaluated for distributional
# characteristics and may cause errors in downstream analyses.
sum(apply(Met05_Akha[-c(1:12)], 2, sd, na.rm = TRUE) == 0)

# List the metabolites with zero variance.
colnames(Met05_Akha[-c(1:12)])[
  apply(Met05_Akha[-c(1:12)], 2, sd, na.rm = TRUE) == 0
]

# Remove zero-variance metabolites to create a dataset
# containing only features with non-zero variability.
log_matrix_metabo_05_Akha_nonzero <-
  Met05_Akha[-c(1:12)][,
                       apply(Met05_Akha[-c(1:12)], 2, sd, na.rm = TRUE) != 0,
                       drop = FALSE
  ]


#--------------------------------------------------------
# 8.3 Normality and Variance Homogeneity Tests
#--------------------------------------------------------

# Evaluate normality and homogeneity of variance for each
# metabolite using the predefined norm_var_check() function.
# The resulting table contains the corresponding test
# statistics and p-values for each metabolite.
Akha_norm_var <- norm_var_check(
  Met05_Akha, 
  metabo_05_descp,
  tag = "metabo05_Akha"
)

# Inspect the normality and variance-homogeneity results.
Akha_norm_var

# Count metabolites for which the variance-homogeneity
# p-value is missing.
Akha_norm_var %>%
  count(is.na(Var_p.value))

# Export the normality and variance-homogeneity test results.
Export(Akha_norm_var, "Akha_norm_var.txt")








#========================================================
# 9. Statistical Analysis — Group Comparison
#========================================================

# Compare metabolite abundances between the two Akha
# geographic groups (Akha-CM and Akha-CR) using the
# Wilcoxon rank-sum test.
#
# Metabolites are subsequently filtered based on the
# FDR-adjusted p-value and log2FC to identify features
# showing both statistical significance and a predefined
# magnitude of difference.


#--------------------------------------------------------
# 9.1 Wilcoxon Test (Akha-CM vs Akha-CR)
#--------------------------------------------------------

# Perform the Wilcoxon rank-sum test for each metabolite
# to compare abundance between Akha-CM and Akha-CR.
#
# The norm_var results are supplied to the helper function
# to incorporate the previously evaluated distributional
# characteristics and variance homogeneity.
Akha_wilcox <- run_wilcoxon_cross(
  metabo_data = Met05_Akha, 
  metabo_norm_var = Akha_norm_var, 
  group_label = "Group",
  comparison_label = "Akha-CM vs Akha-CR",
  tag = "metabo_Akha"
)

# Inspect the first results from the Wilcoxon analysis.
head(Akha_wilcox)


#--------------------------------------------------------
# 9.2 Filter Significant Metabolites
#--------------------------------------------------------

# Identify metabolites showing statistically significant
# differences after false discovery rate (FDR) adjustment.
sig_Akha_05 <- Akha_wilcox %>% 
  filter(p_adj < 0.05) %>%
  pull(Metabolite)

# Retain metabolites that meet both criteria:
# 1. FDR-adjusted p-value < 0.05
# 2. Absolute log2FC >= 1, corresponding to at least
#    a 2-fold difference between groups.
#
# Metabolite classification information is then added
# using the metabolite annotation table.
log2FC_sum_Met05_Akha_sig <- log2FC_sum_Met05_Akha %>%
  filter(Metabolite %in% sig_Akha_05) %>%
  filter(abs(log2FC) >= 1) %>%
  left_join(Pub_Classy_05_clean,
            by = c("Compound",
                   "Compound_Clean",
                   "Metabolite" = "Code"))

# Export the significant metabolite results.
Export(log2FC_sum_Met05_Akha_sig, "Akha_sig_log2FC_1.txt")


#========================================================
# 10. Visualization — Significant Metabolites + Annotation
#========================================================

# The following optional filters can be used to exclude
# selected metabolite superclass or subclass categories
# from the annotated visualization.

#drop_superclass <- c("Hydrocarbons",
#                     "Organometallic compounds",
#                     "Organohalogen compounds")

#drop_subclass <- c("Organic cyanides")

#log2FC_sum_Met05_Akha_sig_anno <- log2FC_sum_Met05_Akha_sig %>%
#  filter(!Superclass %in% drop_superclass &
#           !Subclass %in% drop_subclass)

#log2FC_sum_Met05_Akha_sig_anno %>%
#  distinct(Subclass) %>%
#  print(n = 40)


# Generate an annotated barplot of significant metabolites.
# Metabolites are displayed according to their log2FC and
# classified by direction of change.
#
# Note: the statistical filtering above uses |log2FC| >= 1
# and q < 0.05. The log2FC cutoff supplied to the plotting
# function controls the visual display threshold.
ggFC_sig_Met05_Akha_anno <- plot_log2FC_bar_anno(
  summary_df = log2FC_sum_Met05_Akha_sig,
  title = "Log2FC of metabolites (Akha-CM vs. Akha-CR)",
  xlab = "Significant Metabolite\n(|Log2FC| >= 20, q < 0.05)",
  fill_colors = c("Upregulated" = "steelblue3",
                  "Downregulated" = "#FC766A"),
  log2FC_cutoff = 20,
  strip_text_color = "black",
  group_color =  group_colors
)

# Display the annotated metabolite barplot.
print(ggFC_sig_Met05_Akha_anno)






#=========================================================
# 11. Visualization — Major Subclass Annotation
#=========================================================

# Load the metabolite annotation and analysis workspace
# required for the major-subclass classification and
# downstream annotation checks.
load("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat/Metabolites_Akha.RData")


#---------------------------------------------------------
# 11.1 Assign Major Metabolite Subclasses
#---------------------------------------------------------

# Group individual metabolite subclasses into broader
# chemical categories to facilitate interpretation and
# visualization of significant metabolites.
#
# The classification is based on the original Subclass
# annotation. Subclasses that are not explicitly assigned
# to a major category are labeled as "Unknown".
log2FC_sum_Met05_Akha_sig_anno_majorsubclass <- 
  log2FC_sum_Met05_Akha_sig %>%
  mutate(Major_Subclass = case_when(
    
    # Lipids and fatty acid derivatives
    Subclass %in% c(
      "Fatty acid esters",
      "Fatty acids and conjugates",
      "Fatty aldehydes",
      "Fatty alcohols",
      "Dicarboxylic acids and derivatives"
    ) ~ "Lipids & FA derivatives",
    
    # Terpenoids
    Subclass %in% c(
      "Sesquiterpenoids",
      "Diterpenoids",
      "Monoterpenoids"
    ) ~ "Terpenoids",
    
    # Phenolic and aromatic compounds
    Subclass %in% c(
      "Phenylpropanes",
      "Benzoic acids and derivatives",
      "Benzoyl derivatives",
      "Benzenediols",
      "Biphenyls and derivatives",
      "Toluenes",
      "Anthracenecarboxylic acids and derivatives",
      "Phenoxy compounds"
    ) ~ "Phenolic & aromatic",
    
    # Nitrogen-containing compounds
    Subclass %in% c(
      "Amines",
      "Amino acids, peptides, and analogues",
      "Ureas",
      "Morpholines",
      "Substituted pyrroles",
      "Carbazoles",
      "Isoindolines",
      "Imidazopyridinones",
      "Phenylquinolines",
      "Phenylmethylamines",
      "Organic cyanides"
    ) ~ "Nitrogen-containing",
    
    # Oxygenated molecules
    Subclass %in% c(
      "Alcohols and polyols",
      "Alpha hydroxy acids and derivatives",
      "Ethers",
      "Carbonyl compounds",
      "Gamma butyrolactones",
      "Carbohydrates and carbohydrate conjugates"
    ) ~ "Oxygenated molecules",
    
    # Hydrocarbons
    Subclass %in% c(
      "Alkanes",
      "Cycloalkanes",
      "Unsaturated aliphatic hydrocarbons",
      "Branched unsaturated hydrocarbons",
      "Olefins"
    ) ~ "Hydrocarbons",
    
    # Sulfur-containing compounds
    Subclass %in% c(
      "Dialkylthioethers"
    ) ~ "Sulfur-containing",
    
    # Silicon-containing compounds
    Subclass %in% c(
      "Organosilicon compounds"
    ) ~ "Silicon-containing",
    
    # Retain unclassified subclasses as "Unknown".
    TRUE ~ "Unknown"
    
  ))


#---------------------------------------------------------
# 11.2 Identify Top Enriched Metabolites
#---------------------------------------------------------

# Identify strongly upregulated metabolites based on the
# predefined log2FC threshold and rank them by effect size.
topenrich_Akha <- log2FC_sum_Met05_Akha_sig_anno_majorsubclass %>% 
  filter(Direction == "Upregulated" & log2FC >= 20) %>% 
  arrange(desc(log2FC)) %>% 
  select(Group, Metabolite, Compound_Clean, 
         log2FC, Direction, Subclass, Major_Subclass) %>%
  rename("Code" = "Metabolite",
         "Compound" = "Compound_Clean")

# Export the list of highly enriched metabolites.
Export(topenrich_Akha, "topenrich_Akha.txt")

# Identify metabolites that overlap between the highly
# enriched metabolite list and the top features selected
# for volcano-plot labeling.
intersect(topenrich_Akha$Code, top_features$Metabolite)


#---------------------------------------------------------
# 11.3 Inspect Major Subclass Categories
#---------------------------------------------------------

# List the unique major subclasses represented among the
# significant metabolites.
log2FC_sum_Met05_Akha_sig_anno_majorsubclass %>%
  distinct(Major_Subclass)

# Count significant metabolites without a successfully
# matched compound identifier (CID).
log2FC_sum_Met05_Akha_sig_anno_majorsubclass %>% 
  filter(CID == "Not found") %>%
  count()


#---------------------------------------------------------
# 11.4 Visualize Metabolites by Major Subclass
#---------------------------------------------------------

# Generate a subclass-specific barplot.
# The example below displays metabolites assigned to the
# Terpenoids major subclass.
ggFC_sig_Met05_Akha_subclass <- plot_log2FC_bar_subclass(
  
  summary_df = log2FC_sum_Met05_Akha_sig_anno_majorsubclass %>% 
    filter(Major_Subclass == "Terpenoids"),
  
  title = "Log2FC of metabolites (Akha-CM vs. Akha-CR)\nxxx",
  
  xlab = "Significant Metabolite\n(|Log2FC| >= 1, q < 0.05)",
  
  fill_colors = c("Upregulated" = "#FC766A",
                  "Downregulated" = "steelblue3"),
  
  log2FC_cutoff = 1,
  
  strip_text_color = "black",
  
  group_color =  group_colors
)

# Display the subclass-specific metabolite barplot.
print(ggFC_sig_Met05_Akha_subclass)




#=========================================================
# 12. Statistical Summary and Feature-Level QC
#=========================================================

# Summarize the number of significant metabolites according
# to their direction of change.
log2FC_sum_Met05_Akha_sig_anno_majorsubclass %>%
  count(Direction)

# Summarize significant metabolites by direction of change
# and major chemical subclass.
log2FC_sum_Met05_Akha_sig_anno_majorsubclass %>% 
  group_by(Direction, Major_Subclass) %>%
  count() %>%
  arrange(desc(n)) %>%
  print(n = 50)



#---------------------------------------------------------
# 12.1 Identify Zero-Variance Metabolites
#---------------------------------------------------------

# Load packages required for data manipulation.
library(dplyr)
library(tidyr)

# Define a helper function to identify metabolites with
# zero standard deviation within a given dataset.
get_zero_var <- function(df) {

  df %>%
    summarise(across(starts_with("Met"), ~ sd(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Metabolite", values_to = "SD") %>%
    filter(SD == 0) %>%
    pull(Metabolite)
}


# Identify metabolites with zero variance across all Akha
# samples, regardless of geographic group.
Met05_Akha_zero <- Met05_Akha %>%
  filter(Ethnicity == "Akha") %>%
  get_zero_var()

# Identify metabolites with zero variance specifically
# within the Akha-CM group.
Met05_AkhaCM_zero <- Met05_Akha %>%
  filter(Group == "Akha-CM") %>%
  get_zero_var()

# Identify metabolites with zero variance specifically
# within the Akha-CR group.
Met05_AkhaCR_zero <- Met05_Akha %>%
  filter(Group == "Akha-CR") %>%
  get_zero_var()


# Identify metabolites with zero variance in Akha-CM but
# non-zero variance across the complete Akha dataset.
zero_AkhaCM <- setdiff(Met05_AkhaCM_zero, Met05_Akha_zero)

# Identify metabolites with zero variance in Akha-CR but
# non-zero variance across the complete Akha dataset.
zero_AkhaCR <- setdiff(Met05_AkhaCR_zero, Met05_Akha_zero)


#---------------------------------------------------------
# 12.2 Check Zero-Variance Metabolites Among Significant
#      Features
#---------------------------------------------------------

# Combine metabolites with group-specific zero variance
# and identify those that are also statistically significant.
diff_CM_CR <- union(zero_AkhaCM, zero_AkhaCR)
sig_zero_var <- intersect(diff_CM_CR, sig_Akha_05)

# Create a table of significant metabolites that have
# zero variance in either Akha-CM or Akha-CR.
sigzero_either_Akha <- bind_rows(
  tibble(Metabolite = zero_AkhaCM, Group = "Akha-CM"),
  tibble(Metabolite = zero_AkhaCR, Group = "Akha-CR")
) %>%
  filter(Metabolite %in% sig_zero_var) %>%
  rename(zero_var = Group) %>%
  left_join(log2FC_sum_Met05_Akha_sig_anno_majorsubclass,
            by = "Metabolite")

# Summarize these significant zero-variance metabolites
# by group and original metabolite subclass.
sigzero_either_Akha %>% 
  group_by(Group, Subclass) %>%
  count() %>%
  arrange(desc(n)) %>%
  print(n = 50)

# Summarize significant zero-variance metabolites according
# to the group in which zero variance was observed.
sigzero_either_Akha %>% 
  group_by(zero_var) %>%
  count()
  


#---------------------------------------------------------
# 12.3 Identify Low-Detection Metabolites
#---------------------------------------------------------

# Identify metabolites detected in fewer than three
# Akha-CM samples but detected in at least one sample.
# Detection is defined as an abundance greater than zero.
Met05_AkhaCM_low_detect <- Met05_Akha %>% 
  filter(Group == "Akha-CM") %>%
  summarise(across(starts_with("Met"), ~ sum(.x > 0, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Metabolite", values_to = "Detected_n") %>%
  filter(Detected_n > 0, Detected_n < 3) %>%
  pull(Metabolite)

# Identify metabolites detected in fewer than three
# Akha-CR samples but detected in at least one sample.
Met05_AkhaCR_low_detect <- Met05_Akha %>% 
  filter(Group == "Akha-CR") %>%
  summarise(across(starts_with("Met"), ~ sum(.x > 0, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Metabolite", values_to = "Detected_n") %>%
  filter(Detected_n > 0, Detected_n < 3) %>%
  pull(Metabolite)


# Identify metabolites with low detection in both groups.
intersect(Met05_AkhaCM_low_detect, Met05_AkhaCR_low_detect)

# Identify metabolites with low detection specifically
# in Akha-CM but not Akha-CR.
setdiff(Met05_AkhaCM_low_detect, Met05_AkhaCR_low_detect)


#---------------------------------------------------------
# 12.4 Identify Significant Metabolites with Low Detection
#---------------------------------------------------------

# Check whether low-detection metabolites are included
# among the statistically significant metabolites.
intersect(Met05_AkhaCM_low_detect, sig_Akha_05)
intersect(Met05_AkhaCR_low_detect, sig_Akha_05)

# Combine low-detection metabolites identified in either
# geographic group.
low_detect_Akha <- union(Met05_AkhaCM_low_detect, Met05_AkhaCR_low_detect)

# Remove low-detection metabolites from the list of
# statistically significant metabolites to generate a
# cleaned feature list for downstream analysis.
sig_Akha_05_clean <- setdiff(sig_Akha_05, low_detect_Akha)
