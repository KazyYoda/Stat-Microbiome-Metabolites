############################################################
# Metabolite Analysis — Chiang Mai (CM)
############################################################

# This section compares metabolite profiles among the three
# ethnic groups sampled in Chiang Mai (CM):
#   - Akha-CM
#   - Lahu-CM
#   - Khuen-CM
#
# Metabolite abundances are log2-transformed before summary
# statistics and fold-change calculations.
#
# Pairwise log2FC values are calculated relative to defined
# reference groups to facilitate interpretation of the
# between-group differences.


#------------------------------------------------------------
# 1. Working Directory and Data
#------------------------------------------------------------

setwd("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat")

# Load the metabolite analysis workspace.
load("~/Documents/HillTribe_NGS/5.Metabolites/Metabolites.RData")


#------------------------------------------------------------
# 2. Load Required Packages
#------------------------------------------------------------

library(dplyr)
library(tidyr)
library(car)
library(readxl)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggh4x)


#------------------------------------------------------------
# 3. Define Group Colors
#------------------------------------------------------------

# Define colors used consistently for the three CM groups
# throughout subsequent visualizations.
group_colors <- c(
  "Akha-CM" = "gold1",
  "Lahu-CM" = "lightblue2",
  "Khuen-CM" = "darkorange"
)


#========================================================
# 4. Compute log2 Fold-Change (log2FC)
#    Relative to Defined Contrast Groups
#========================================================

#--------------------------------------------------------
# 4.1 Merge Metadata with Metabolite Matrix
#--------------------------------------------------------

# Confirm that the sample order in the metabolite matrix
# matches the sample order in the metadata before combining
# the two datasets.
identical(rownames(metabo_05_t), Met05_match$Sample_ID)

# Combine sample metadata with metabolite abundance data.
Met05 <- cbind(Met05_match, metabo_05_t)

# Apply log2 transformation to all metabolite abundance
# columns. The first 13 columns contain sample metadata.
Met05_log2 <- Met05 %>%
  mutate(across(14:ncol(.), log2))



#--------------------------------------------------------
# 4.2 Subset Data by Area
#--------------------------------------------------------

# Select samples collected from Chiang Mai (CM).
Met05_CM <- Met05_log2 %>% filter(Area == "CM")


#--------------------------------------------------------
# 4.2.1 Calculate Group-Level Summary Statistics
#--------------------------------------------------------

# Identify metabolite columns.
met_cols <- grep("^Met", names(Met05_CM), value = TRUE)

# Calculate the mean and standard deviation of each metabolite
# within each CM group.
summary_CM_long <- Met05_CM %>%
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



#--------------------------------------------------------
# 4.2.2 Create Wide-Format Summary Table
#--------------------------------------------------------

# Combine the mean and SD into a single formatted value
# for each metabolite and group.
summary_CM_wide <- summary_CM_long %>%
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
  
  # Add Kruskal-Wallis results for the three CM groups.
  left_join(
    CM_kruskal_metabo,
    by = "Metabolite"
  ) %>% 
  
  # Add log2FC for Akha-CM relative to Lahu-CM.
  left_join(
    log2FC_sum_Met05_CM2 %>%
      filter(Group == "Akha-CM") %>%
      select(Metabolite, Log2FC_AkhavsLahu = log2FC),
    by = "Metabolite"
  ) %>%
  
  # Add log2FC for Khuen-CM relative to Lahu-CM.
  left_join(
    log2FC_sum_Met05_CM2 %>%
      filter(Group == "Khuen-CM") %>%
      select(Metabolite, Log2FC_KhuenvsLahu = log2FC),
    by = "Metabolite"
  ) %>%
  
  # Add log2FC for Khuen-CM relative to Akha-CM.
  left_join(
    log2FC_sum_Met05_CM %>%
      filter(Group == "Khuen-CM") %>%
      select(Metabolite, Log2FC_KhuenvsAkha = log2FC),
    by = "Metabolite"
  )
    


# Export the complete CM metabolite summary table.
Export(summary_CM_wide, "summary_CM_wide.txt")


# Count significant metabolites according to direction
# using the specified log2FC and q-value thresholds.
summary_CM_wide %>%
  filter(abs(log2FC) >=1 & q_value < 0.05) %>%
  group_by(Direction) %>%
  count()


# Count metabolite-group combinations with a group mean
# equal to zero, indicating no detected signal after
# summarization.
summary_CM_long %>%
  group_by(Group) %>%
  summarise(
    Not_Detected = sum(Mean == 0, na.rm = TRUE)
  )




#========================================================
# 4.3 Log2FC Analysis — CM
#     Akha-CM, Lahu-CM, and Khuen-CM
#========================================================

# Two sets of pairwise log2FC values are calculated using
# different reference groups:
#
#   log2FC_Met05_CM:
#       Reference = Akha-CM
#
#   log2FC_Met05_CM2:
#       Reference = Lahu-CM
#
# This allows the CM groups to be compared using the
# desired reference group for each downstream analysis.


#--------------------------------------------------------
# 4.3.1 Compute log2FC Relative to Contrast Groups
#--------------------------------------------------------

# Calculate group-level log2FC values using Akha-CM
# as the reference group.
log2FC_Met05_CM <- log2FC_group_summary_relative_to_contrast(
  metabo_data = Met05_CM[-c(1:13)], 
  metadata = Met05_CM[1:13], 
  contrast = "Akha-CM", 
  Group = "Group"
)

# Calculate group-level log2FC values using Lahu-CM
# as the reference group.
log2FC_Met05_CM2 <- log2FC_group_summary_relative_to_contrast(
  metabo_data = Met05_CM[-c(1:13)], 
  metadata = Met05_CM[1:13], 
  contrast = "Lahu-CM", 
  Group = "Group"
)


#--------------------------------------------------------
# 4.3.2 Transpose log2FC Results
#--------------------------------------------------------

# Transpose the log2FC matrices for downstream formatting
# and export.
log2FC_Met05_CM_t <- log2FC_Met05_CM %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Code")

log2FC_Met05_CM2_t <- log2FC_Met05_CM2 %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Code")



#--------------------------------------------------------
# 4.3.3 Summarize log2FC Using Akha-CM as Reference
#--------------------------------------------------------

# Summarize log2FC values for comparisons relative to
# Akha-CM.
log2FC_sum_Met05_CM <- summarize_log2FC(
  log2FC_Met05_CM, 
  filter_group = "Akha-CM"
)

# Add metabolite annotation.
log2FC_sum_Met05_CM <- log2FC_sum_Met05_CM %>%
  left_join(
    metabo_05_descp, 
    by = c("Metabolite" = "Code")
  )

# Export the log2FC summary.
Export(log2FC_sum_Met05_CM, "log2FC_sum_Met05_CM.txt")



#--------------------------------------------------------
# 4.3.4 Summarize log2FC Using Lahu-CM as Reference
#--------------------------------------------------------

# Summarize log2FC values for comparisons relative to
# Lahu-CM.
log2FC_sum_Met05_CM2 <- summarize_log2FC(
  log2FC_Met05_CM2, 
  filter_group = "Lahu-CM"
)

# Add metabolite annotation.
log2FC_sum_Met05_CM2 <- log2FC_sum_Met05_CM2 %>%
  left_join(
    metabo_05_descp, 
    by = c("Metabolite" = "Code")
  )

# Export the log2FC summary.
Export(log2FC_sum_Met05_CM2, "log2FC_sum_Met05_CM2.txt")






#========================================================
# 5. Log2FC Visualization — Barplot
#========================================================

# Generate a grouped barplot showing metabolite log2FC values
# for comparisons against Akha-CM.
#
# The log2FC summary contains comparisons of:
#   - Lahu-CM vs Akha-CM
#   - Khuen-CM vs Akha-CM
#
# Akha-CM is therefore used as the reference group for
# interpreting the direction of the log2FC.
ggFC_Met05_CM <- plot_log2FC_bar(
  summary_df = log2FC_sum_Met05_CM,
  title = "Log2FC of metabolites (Khuen-CM or Lahu-CM vs. Akha-CM)",
  xlab = "Metabolite",
  fill_colors = c("Upregulated" = "steelblue3",
                  "Downregulated" = "#FC766A"),
  strip_text_color = "black",
  log2FC_cutoff = 15,
  group_color = group_colors
)

# Display the log2FC barplot.
print(ggFC_Met05_CM)





#========================================================
# 5.1 Pairwise Dunn-Test Results
#     Akha-CM vs Lahu-CM and Khuen-CM vs Lahu-CM
#========================================================

# Extract the two pairwise comparisons in which Lahu-CM
# is the reference group.
#
# The original Dunn-test results contain unadjusted and
# adjusted p-values for each pairwise metabolite comparison.
# Comparison labels are simplified for downstream plotting.
CM2_dunn <- CM_dunn %>%
  filter(Comparison %in% c(
    "Akha-CM - Lahu-CM",
    "Khuen-CM - Lahu-CM"
  )) %>%
  mutate(
    Comparison_CM2 = case_when(
      Comparison == "Akha-CM - Lahu-CM" ~ "Akha vs Lahu",
      Comparison == "Khuen-CM - Lahu-CM" ~ "Khuen vs Lahu"
    ),
    Group = case_when(
      Comparison == "Akha-CM - Lahu-CM" ~ "Akha-CM",
      Comparison == "Khuen-CM - Lahu-CM" ~ "Khuen-CM"
    )
  ) %>%
  select(Group, Comparison_CM2, Metabolite, Compound_Clean, P.unadj, P.adj)


#--------------------------------------------------------
# 5.2 Prepare Pairwise Differential Metabolites for
#     Dot-Plot Visualization
#--------------------------------------------------------

# Combine log2FC values with the corresponding pairwise
# Dunn-test results.
#
# The current filter retains Khuen-CM metabolites with
# an absolute log2FC >= 19 and a non-missing Dunn-test
# adjusted p-value.
log2FC_CM2_sig <- log2FC_sum_Met05_CM2_sig %>% 
  filter(Group == "Khuen-CM") %>%
  left_join(CM2_dunn %>% filter(Group == "Khuen-CM"),
            by = c("Metabolite", "Group", "Compound_Clean")) %>%
  select(Group, Comparison_CM2, Metabolite, Compound_Clean, Subclass, log2FC, P.adj, Direction) %>%
  filter(abs(log2FC) >= 19) %>%
  mutate(`-log10(q-value)` = -log10(P.adj)) %>%
  filter(!is.na(P.adj))



#--------------------------------------------------------
# 5.3 Dot-Plot of Pairwise Differential Metabolites
#--------------------------------------------------------

library(ggplot2)
library(dplyr)
library(stringr)

# Define colors according to the direction of the log2FC.
fill_colors <- c(
  "Downregulated" = "steelblue3",
  "Upregulated" = "#FC766A"
)


# Generate a dot plot in which:
#   - x-axis = metabolite
#   - y-axis = mean log2FC
#   - point size = -log10(Dunn-test adjusted p-value)
#   - point fill = direction of change
log2FC_CM2_sig %>%
  mutate(
    Metabolite = reorder(Metabolite, log2FC)
  ) %>%
  ggplot(aes(
    x = Metabolite,
    y = log2FC
  )) +
  
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
    title = "Differentially enriched metabolites\n(Khuen-CM vs Lahu-CM)"
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
# 5.4 Log2FC Interpretation Guide
#--------------------------------------------------------

# Log2FC is interpreted relative to the specified reference
# group. Positive values indicate higher metabolite abundance
# in the comparison group, whereas negative values indicate
# lower abundance relative to the reference group.

# Upregulation
# log2FC >= 1  -> >= 2-fold increase
# log2FC >= 2  -> >= 4-fold increase
# log2FC >= 3  -> >= 8-fold increase
# log2FC >= 4  -> >= 16-fold increase

# Fold change = 2^(log2FC)

# Downregulation
# log2FC <= -1 -> <= 0.5-fold
# log2FC <= -2 -> <= 0.25-fold
# log2FC <= -3 -> <= 0.125-fold
# log2FC <= -4 -> <= 0.0625-fold




#========================================================
# 5.5 Volcano Plot — CM Metabolite Analysis
#========================================================

# Prepare log2FC results using Akha-CM and Lahu-CM as
# separate reference groups.
log2FC_CM1 <- log2FC_sum_Met05_CM %>% mutate(Comp = "vs Akha-CM")
log2FC_CM2 <- log2FC_sum_Met05_CM2 %>% mutate(Comp = "vs Lahu-CM")

# Combine the two reference-based log2FC tables.
log2FC_CM <- rbind(log2FC_CM1, log2FC_CM2)

# Add the overall Kruskal-Wallis results to the log2FC
# table for the three-group CM comparison.
vol_CM <- log2FC_CM %>%
  left_join(CM_kruskal_metabo, by = "Metabolite") 


#--------------------------------------------------------
# 5.6 Prepare Volcano-Plot Data
#--------------------------------------------------------

# Calculate -log10(p-value) from the Kruskal-Wallis test
# and classify metabolites according to the overall
# significance threshold.
vol_CM_plot <- vol_CM %>%
  mutate(
    logP = -log10(kruskal_p),
    class = case_when(
      kruskal_p < 0.05 ~ "Significant",
      TRUE ~ "Not Significant"
    )
  ) %>%
  filter(
    is.finite(logP),
    !Direction == "No change"
  ) %>%
  mutate(
    Significant = case_when(
      class == "Significant" & Direction == "Downregulated" ~ "Downregulated",
      class == "Significant" & Direction == "Upregulated" ~ "Upregulated",
      class == "Not Significant" ~ "Not Significant"
    )
  )


#--------------------------------------------------------
# 5.7 Select Top Features for Annotation
#--------------------------------------------------------

# Select the 30 metabolites with the smallest overall
# Kruskal-Wallis p-values for labeling on the volcano plot.
top_features <- vol_CM_plot %>%
  filter(class == "Significant") %>%
  arrange(kruskal_p) %>%
  slice(1:30)


#--------------------------------------------------------
# 5.8 Generate Volcano Plot
#--------------------------------------------------------

# The volcano plot displays:
#   - x-axis = mean log2FC
#   - y-axis = -log10(Kruskal-Wallis p-value)
#   - point color = direction and overall significance
#
# The vertical dashed lines indicate |log2FC| = 1.
# The horizontal dashed line indicates p = 0.05.
ggplot(vol_CM_plot, aes(x = log2FC, y = logP)) +
  geom_point(aes(color = Significant), alpha = 0.7, size = 1) +
  
  # Highlight the selected top features.
  geom_point(
    data = top_features,
    aes(color = Significant),
    size = 1
  ) +
  
  # Label the selected top features.
  geom_text_repel(
    data = top_features,
    aes(label = Metabolite),
    size = 1,
    max.overlaps = 100,
    segment.size = 0.2
  ) +
  
  scale_color_manual(values = c("Upregulated" = "#FC766A",
                                "Downregulated" = "#89ABE3",
                                "Not Significant" = "grey70")) +
  
  # Effect-size reference lines.
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "darkgray") +
  
  # Statistical significance reference line.
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "darkgray") +
  
  labs(
    title = "Volcano plot of metabolites (CM group)",
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
# 6. Heatmap Visualization (Log2-Transformed Metabolite Data)
#========================================================

#--------------------------------------------------------
# 6.1 Prepare Log2-Transformed Metabolite Matrix
#--------------------------------------------------------

# Apply log2 transformation to the complete metabolite abundance matrix.
# The original matrix contains samples as rows and metabolites as columns.
log_matrix_metabo_05 <- as.matrix(log2(metabo_05_t))

# Transpose the matrix so that metabolites are represented as rows
# and samples are represented as columns, which is the format required
# for heatmap visualization.
log_matrix_metabo_05_t <- t(log_matrix_metabo_05)



#-----------------------------------------------------------------------------
# 6.2 Select Metabolites for Heatmap Visualization
#     Based on the CM2 Comparison
#-----------------------------------------------------------------------------
# Identify metabolites showing a significant overall group difference
# according to the Kruskal-Wallis test.
sig_kruskal <- CM_kruskal_metabo %>%
  filter(kruskal_p < 0.05) %>%
  pull(Metabolite)

# Retain metabolites that are significant in the Kruskal-Wallis test
# and are also present in the CM2 log2FC results.
sig_across <- log2FC_sum_Met05_CM2_sig %>%
  filter(Metabolite %in% sig_kruskal)

# Apply the predefined absolute log2FC threshold to select metabolites
# with a relatively large magnitude of change for heatmap visualization.
log2FC_X_CM <- log2FC_sum_Met05_CM2_sig %>%
  filter(abs(log2FC) >= 15)

# Match the selected metabolites to the rows available in the
# log2-transformed metabolite matrix.
valid_rows <- intersect(
  rownames(log_matrix_metabo_05_t),
  log2FC_X_CM$Metabolite
)

# Match CM sample IDs to the columns available in the metabolite matrix.
valid_cols <- intersect(
  colnames(log_matrix_metabo_05_t),
  Met05_CM$Sample_ID
)

# Subset the matrix to the selected metabolites and CM samples.
# drop = FALSE preserves the matrix structure even when only one
# metabolite or sample is retained.
log_matrix_metabo_05_t_CM <-
  log_matrix_metabo_05_t[
    valid_rows,
    valid_cols,
    drop = FALSE
  ]



#--------------------------------------------------------
# 6.3 Variance Check Before Hierarchical Clustering
#--------------------------------------------------------

# Count metabolites with zero variance across the selected CM samples.
# Zero-variance features provide no information for distance-based
# clustering and should therefore be identified before interpretation.
sum(apply(log_matrix_metabo_05_t_CM, 1, sd, na.rm = TRUE) == 0)

# Count samples with zero variance across the selected metabolites.
# This provides an additional quality-control check before clustering.
sum(apply(log_matrix_metabo_05_t_CM, 2, sd, na.rm = TRUE) == 0)

# Note:
# Zero-variance features can affect Euclidean distance calculations
# and hierarchical clustering.



#--------------------------------------------------------
# 6.4 Define Group Order for Sample Annotation
#--------------------------------------------------------

# Convert Group to a factor and explicitly define the order of groups
# displayed in the heatmap annotation.
Met05_CM$Group <- factor(
  Met05_CM$Group,
  levels = c("Akha-CM", "Lahu-CM", "Khuen-CM")
)



#========================================================
# 6.5 Generate Heatmap
#========================================================

# Generate the heatmap using the selected log2-transformed metabolite
# matrix and corresponding sample metadata.
#
# Matrix structure:
#   Rows    = selected metabolites
#   Columns = CM samples
#   Values  = log2-transformed metabolite abundance
#
# Group colors are used to annotate samples according to their
# population/ethnic group.
CM <- plot_heatmap(
  mat = log_matrix_metabo_05_t_CM,
  metadata = Met05_CM[1:13],
  legend_title = "Log2 (CM)",
  group_color = group_colors
)

# Draw the heatmap and place both the heatmap and annotation legends
# on the left side of the figure.
draw(
  CM,
  heatmap_legend_side = "left",
  annotation_legend_side = "left"
)






#========================================================
# 7. Principal Component Analysis (PCA) — CM Metabolites
#========================================================

# PCA is performed to visualize the overall structure of the
# metabolite profiles across CM groups.
#
# The analysis uses metabolites meeting the predefined absolute
# log2FC threshold. Metabolites are standardized before PCA so
# that each metabolite contributes on a comparable scale.
#
# Groups included:
#   - Akha-CM
#   - Lahu-CM
#   - Khuen-CM


library(FactoMineR)
library(factoextra)


#--------------------------------------------------------
# 7.1 Select Metabolites for PCA
#--------------------------------------------------------

# Select metabolites with an absolute log2FC of at least 1
# from the CM2 comparison for inclusion in PCA.
log2FC_1_CM <- log2FC_sum_Met05_CM2_sig %>%
  filter(abs(log2FC) >= 1)


#--------------------------------------------------------
# 7.2 Subset Log2-Transformed Matrix for PCA
#--------------------------------------------------------

# Match CM sample IDs between the metabolite matrix and metadata.
valid_cols <- intersect(
  colnames(log_matrix_metabo_05_t),
  Met05_CM$Sample_ID
)

# Match selected metabolite identifiers between the differential
# metabolite results and the log2-transformed metabolite matrix.
valid_rows <- intersect(
  rownames(log_matrix_metabo_05_t),
  log2FC_1_CM$Metabolite
)

# Subset the matrix to the selected metabolites and CM samples.
#
# Matrix structure at this stage:
#   Rows    = metabolites
#   Columns = samples
#   Values  = log2-transformed metabolite abundance
log_matrix_metabo_05_t_CM_pca <-
  log_matrix_metabo_05_t[
    valid_rows,
    valid_cols,
    drop = FALSE
  ]


#--------------------------------------------------------
# 7.3 Standardize Metabolite Abundance Before PCA
#--------------------------------------------------------

# Standardize each metabolite to mean = 0 and SD = 1.
# This ensures that metabolites are placed on a comparable scale
# before PCA.
pca_CM <- scale(log_matrix_metabo_05_t_CM_pca)

# Inspect the standardized matrix and check for NaN values.
head(pca_CM)
sum(is.nan(pca_CM))

# Replace NaN values with zero.
# NaN values can occur when a metabolite has zero standard deviation
# during standardization.
pca_CM[is.nan(pca_CM)] <- 0


#--------------------------------------------------------
# 7.4 Perform Principal Component Analysis
#--------------------------------------------------------

# Transpose the matrix so that samples are rows and metabolites
# are columns, as required by FactoMineR::PCA.
#
# scale.unit = FALSE is used because the metabolite variables
# have already been standardized in the previous step.
res_pca <- PCA(
  t(pca_CM),
  scale.unit = FALSE,
  graph = FALSE
)


#--------------------------------------------------------
# 7.5 Prepare Sample Metadata for Visualization
#--------------------------------------------------------

# Extract sample identifiers corresponding to the PCA input matrix.
sample_names <- colnames(log_matrix_metabo_05_t_CM_pca)

# Extract the corresponding CM group assignments.
group_labels <- Met05_CM$Group

# Create a metadata table for assigning group information to PCA
# observations.
metadata <- data.frame(
  Sample = sample_names,
  Group = group_labels
)

# Use sample IDs as row names to facilitate correspondence between
# metadata and PCA observations.
rownames(metadata) <- sample_names


#--------------------------------------------------------
# 7.6 Visualize PCA with Group Confidence Ellipses
#--------------------------------------------------------

# Define colors for the three CM groups.
group_colors <- c(
  "Akha-CM" = "gold1",
  "Lahu-CM" = "lightblue2",
  "Khuen-CM" = "darkorange"
)

# Plot individual samples in PCA space and color them according
# to their CM group.
#
# A 95% confidence ellipse is added for each group to visualize
# the approximate distribution of samples within PCA space.
fviz_pca_ind(
  res_pca,
  geom.ind = "point",
  col.ind = metadata$Group,
  palette = group_colors,
  addEllipses = TRUE,
  ellipse.type = "confidence",
  ellipse.level = 0.95,
  repel = TRUE
) +
  labs(title = "PCA of CM Metabolites (|Log2FC| >= 1)") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 6),
    title = element_text(size = 6)
  )






#========================================================
# 8. Normality and Homogeneity of Variance Testing
#========================================================

# This section prepares the metabolite data for statistical testing
# and evaluates distributional and variance-related characteristics
# of the metabolites across the CM groups.

setwd("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat")


#--------------------------------------------------------
# 8.1 Prepare CM Metabolite Matrix
#--------------------------------------------------------

# Match metabolite columns between the log2-transformed metabolite
# matrix and the metabolites included in the CM log2FC results.
valid_cols <- intersect(
  colnames(log_matrix_metabo_05),
  log2FC_sum_Met05_CM$Metabolite
)

# Match CM sample IDs between the metabolite matrix and metadata.
valid_rows <- intersect(
  rownames(log_matrix_metabo_05),
  Met05_CM$Sample_ID
)

# Subset the matrix to the CM samples and selected metabolites.
#
# Matrix structure:
#   Rows    = samples
#   Columns = metabolites
log_matrix_metabo_05_CM <-
  log_matrix_metabo_05[
    valid_rows,
    valid_cols,
    drop = FALSE
  ]


#--------------------------------------------------------
# 8.2 Identify and Remove Zero-Variance Metabolites
#--------------------------------------------------------

# Count metabolites with zero variance across the CM samples.
# Zero-variance metabolites contain no within-dataset variability
# and may not be suitable for statistical or distance-based analyses.
sum(apply(Met05_CM[-c(1:13)], 2, sd, na.rm = TRUE) == 0)

# Identify metabolites with zero variance.
colnames(Met05_CM[-c(1:13)])[
  apply(Met05_CM[-c(1:13)], 2, sd, na.rm = TRUE) == 0
]

# Create a metabolite matrix containing only non-zero-variance features.
log_matrix_metabo_05_CM_nonzero <-
  Met05_CM[-c(1:13)][,
                       apply(Met05_CM[-c(1:13)], 2, sd, na.rm = TRUE) != 0,
                       drop = FALSE
  ]


#--------------------------------------------------------
# 8.3 Normality and Homogeneity of Variance Tests
#--------------------------------------------------------

# Evaluate distributional characteristics and homogeneity of variance
# for each metabolite using the custom quality-control function.
#
# The resulting table contains the corresponding test statistics
# and p-values used to inform downstream statistical testing.
CM_norm_var <- norm_var_check(
  Met05_CM,
  metabo_05_descp,
  tag = "metabo05_CM"
)

CM_norm_var

# Count metabolites for which the variance-test p-value is missing.
# Missing values may occur when a test cannot be performed because
# of insufficient variability or other data limitations.
CM_norm_var %>%
  count(is.na(Var_p.value))

# Export normality and variance-homogeneity test results.
Export(CM_norm_var, "CM_norm_var.txt")



#========================================================
# 9. Statistical Analysis — Group Comparison
#========================================================

# This section performs group-level statistical comparisons
# among the three CM groups:
#
#   - Akha-CM
#   - Lahu-CM
#   - Khuen-CM
#
# The primary workflow uses:
#
#   Kruskal-Wallis test
#           ↓
#       Dunn's post-hoc test
#           ↓
#   Multiple-testing-adjusted p-values
#           ↓
#   Effect-size filtering using log2FC
#
# Individual ANOVA/Tukey analyses shown below are performed
# specifically for Met26.


#--------------------------------------------------------
# 9.0 Example Parametric Analysis — Met26
#--------------------------------------------------------

# Prepare the Met26 data for a one-way ANOVA.
Met26_CM <- Met05_CM %>% 
  select(Group, Met26) %>%
  as.data.frame()

# Perform one-way ANOVA to test for differences in Met26
# across the three CM groups.
aov_result <- aov(Met26 ~ Group, data = Met26_CM)
summary(aov_result)

# Perform Tukey's HSD post-hoc test following ANOVA to identify
# pairwise group differences in Met26.
tukey_result <- TukeyHSD(aov_result)
tukey_result



#--------------------------------------------------------
# 9.1 Kruskal-Wallis Test — CM Groups
#--------------------------------------------------------

# Perform the Kruskal-Wallis test for each metabolite to evaluate
# whether metabolite abundance differs among the three CM groups.
#
# The custom function uses the normality/variance assessment and
# returns the corresponding omnibus test results.
CM_kruskal_metabo <- kruskal_metabo(
  metabo_norm_var = CM_norm_var,
  metabo_data = log_matrix_metabo_05_CM,
  metadata = Met05_CM[1:13], 
  tag = "metabo_CM"
)

head(CM_kruskal_metabo)



#----------------------------------------------------------
# 9.2 Dunn's Test — Post-hoc Pairwise Comparisons
#----------------------------------------------------------

# Perform Dunn's post-hoc test for metabolites showing evidence
# of group differences in the Kruskal-Wallis analysis.
#
# Dunn's test evaluates the pairwise differences among the CM groups
# following the omnibus Kruskal-Wallis test.
CM_dunn <- dunnTest_metabo(
  krus_pvalue_descp = CM_kruskal_metabo, 
  metabo_data = log_matrix_metabo_05_CM, 
  metadata = Met05_CM[1:13], 
  metabo_descp = metabo_05_descp, 
  tag = "metabo_CM"
)

head(CM_dunn)


#--------------------------------------------------------
# 9.3 Identify Significant Metabolites Relative to Akha-CM
#--------------------------------------------------------

# Extract metabolites showing significant pairwise differences
# involving Akha-CM based on the FDR-adjusted Dunn's test p-value.
sig_CM_05 <- CM_dunn %>% 
  filter(Comparison %in% c("Akha-CM - Lahu-CM",
                           "Akha-CM - Khuen-CM") &
         P.adj < 0.05) %>%
  pull(Metabolite)

# Retain metabolites identified as significant in the selected
# pairwise comparisons and apply an additional absolute log2FC
# threshold to retain features with a predefined magnitude of change.
#
# Metabolite class annotations are then added for downstream
# interpretation and visualization.
log2FC_sum_Met05_CM_sig <- log2FC_sum_Met05_CM %>%
  filter(Metabolite %in% sig_CM_05) %>%
  filter(abs(log2FC) >= 1) %>%
  left_join(Pub_Classy_05_clean,
            by = c("Compound",
                   "Compound_Clean",
                   "Metabolite" = "Code"))

# Export the significant metabolite results relative to Akha-CM.
Export(log2FC_sum_Met05_CM_sig, "relAkhaCM_sig_log2FC_1.txt")



#--------------------------------------------------------
# 9.4 Identify Significant Metabolites Relative to Lahu-CM
#--------------------------------------------------------

# Extract metabolites showing significant pairwise differences
# involving Lahu-CM based on the FDR-adjusted Dunn's test p-value.
sig_CM2_05 <- CM_dunn %>% 
  filter(Comparison %in% c("Akha-CM - Lahu-CM",
                           "Khuen-CM - Lahu-CM") &
           P.adj < 0.05) %>%
  pull(Metabolite)

# Retain metabolites identified as significant in the selected
# pairwise comparisons and apply the predefined absolute log2FC
# threshold.
#
# Metabolite class annotations are added for downstream analysis.
log2FC_sum_Met05_CM2_sig <- log2FC_sum_Met05_CM2 %>%
  filter(Metabolite %in% sig_CM2_05) %>%
  filter(abs(log2FC) >= 1) %>%
  left_join(Pub_Classy_05_clean,
            by = c("Compound",
                   "Compound_Clean",
                   "Metabolite" = "Code"))

# Export the significant metabolite results relative to Lahu-CM.
Export(log2FC_sum_Met05_CM2_sig, "relLahuCM_sig_log2FC_1.txt")







#========================================================
# 10. Visualization — Significant Metabolites
#========================================================

# Generate a grouped barplot of significantly different metabolites.
# Metabolites are displayed when they meet the specified statistical
# significance criterion (q < 0.05) and log2FC visualization threshold.
# The direction of change is defined relative to the designated
# contrast group used in the log2FC calculation.
ggFC_sig_Met05_CM <- plot_log2FC_bar(
  summary_df = log2FC_sum_Met05_CM_sig,
  title = "Log2FC of metabolites (Khuen-CM or CM-CM vs. Akha-CM)",
  xlab = "Significant Metabolite\n(|Log2FC| >= 2, q < 0.05)",
  fill_colors = c("Upregulated" = "steelblue3",
                  "Downregulated" = "#FC766A"),
  strip_text_color = "black",
  log2FC_cutoff = 2,
  group_color = group_colors
)

# Display the barplot.
print(ggFC_sig_Met05_CM)



#=========================================================
# 11. Visualization — Major Subclass Annotation
#=========================================================
#
# Metabolites are further grouped into broader chemical categories
# based on their HMDB/PubChem-derived Subclass annotations.
#
# The major subclass categories are intended to facilitate biological
# interpretation and visualization by grouping related chemical
# subclasses into broader compound classes.
#
# Categories include:
#   - Lipids & fatty-acid derivatives
#   - Steroids
#   - Terpenoids
#   - Phenolic & aromatic compounds
#   - Nitrogen-containing compounds
#   - Oxygenated molecules
#   - Hydrocarbons
#   - Others
#   - Unknown
#
# The same classification scheme is applied independently to the two
# CM log2FC datasets, corresponding to their respective reference
# groups.
#=========================================================

log2FC_sum_Met05_CM_sig_anno_majorsubclass <- 
  log2FC_sum_Met05_CM_sig %>%
  mutate(Major_Subclass = case_when(
    
    # Lipids and fatty-acid related compounds
    Subclass %in% c(
      "Fatty acid esters",
      "Fatty acids and conjugates",
      "Fatty aldehydes",
      "Fatty alcohols",
      "Dicarboxylic acids and derivatives"
    ) ~ "Lipids & FA derivatives",
    
    # Steroid compounds
    Subclass %in% c(
      "Cholestane steroids"
    ) ~ "Steroids",
    
    # Terpenoid compounds
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
      "Anthraquinones",
      "Phenoxy compounds",
      "Styrenes",
      "Cyclohexylphenols"
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
      "Indoles",
      "Phenylazides"
    ) ~ "Nitrogen-containing",
    
    # Oxygen-containing compounds
    Subclass %in% c(
      "Alcohols and polyols",
      "Alpha hydroxy acids and derivatives",
      "Ethers",
      "Carbonyl compounds",
      "Gamma butyrolactones",
      "Carbohydrates and carbohydrate conjugates",
      "Carbonic acid diesters"
    ) ~ "Oxygenated molecules",
    
    # Hydrocarbon compounds
    Subclass %in% c(
      "Alkanes",
      "Cycloalkanes",
      "Unsaturated aliphatic hydrocarbons",
      "Branched unsaturated hydrocarbons",
      "Olefins"
    ) ~ "Hydrocarbons",
    
    # Organosilicon compounds and other compounds not assigned
    # to the major chemical categories above
    Subclass %in% c(
      "Organosilicon compounds"
    ) ~ "Others",
    
    # Retain a separate category for subclasses not explicitly
    # assigned to one of the predefined major categories
    TRUE ~ "Unknown"
    
  ))


# Apply the same major-subclass classification to the CM2
# log2FC dataset.
log2FC_sum_Met05_CM2_sig_anno_majorsubclass <- 
  log2FC_sum_Met05_CM2_sig %>%
  mutate(Major_Subclass = case_when(
    
    # Lipids and fatty-acid related compounds
    Subclass %in% c(
      "Fatty acid esters",
      "Fatty acids and conjugates",
      "Fatty aldehydes",
      "Fatty alcohols",
      "Dicarboxylic acids and derivatives"
    ) ~ "Lipids & FA derivatives",
    
    # Steroid compounds
    Subclass %in% c(
      "Cholestane steroids"
    ) ~ "Steroids",
    
    # Terpenoid compounds
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
      "Anthraquinones",
      "Phenoxy compounds",
      "Styrenes",
      "Cyclohexylphenols"
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
      "Indoles",
      "Phenylazides"
    ) ~ "Nitrogen-containing",
    
    # Oxygen-containing compounds
    Subclass %in% c(
      "Alcohols and polyols",
      "Alpha hydroxy acids and derivatives",
      "Ethers",
      "Carbonyl compounds",
      "Gamma butyrolactones",
      "Carbohydrates and carbohydrate conjugates",
      "Carbonic acid diesters"
    ) ~ "Oxygenated molecules",
    
    # Hydrocarbon compounds
    Subclass %in% c(
      "Alkanes",
      "Cycloalkanes",
      "Unsaturated aliphatic hydrocarbons",
      "Branched unsaturated hydrocarbons",
      "Olefins"
    ) ~ "Hydrocarbons",
    
    # Organosilicon compounds and other compounds not assigned
    # to the major chemical categories above
    Subclass %in% c(
      "Organosilicon compounds"
    ) ~ "Others",
    
    # Retain a separate category for subclasses not explicitly
    # assigned to one of the predefined major categories
    TRUE ~ "Unknown"
    
  ))



#---------------------------------------------------------
# 11.1 Identify Highly Enriched Metabolites
#---------------------------------------------------------

# Identify metabolites classified as upregulated with a large
# positive log2FC in the primary CM comparison.
#
# A log2FC threshold of 10 corresponds to a >1,000-fold abundance
# ratio on the original scale (2^10 = 1,024).
topenrich_CM <- log2FC_sum_Met05_CM_sig_anno_majorsubclass %>% 
  filter(Direction == "Upregulated" & log2FC >= 10) %>% 
  arrange(desc(log2FC)) %>% 
  select(Group, Metabolite, Compound_Clean, 
         log2FC, Direction, Subclass, Major_Subclass) %>%
  rename("Code" = "Metabolite",
         "Compound" = "Compound_Clean")

# Export the highly enriched metabolite list.
Export(topenrich_CM, "topenrich_CM.txt")

# Identify metabolites shared between the highly enriched CM list
# and the previously selected top features.
intersect(topenrich_CM$Code, top_features$Metabolite)



# Identify highly enriched metabolites in the CM2 comparison.
#
# Here, a more stringent log2FC threshold of 20 is used.
# A log2FC of 20 corresponds to approximately a 1,048,576-fold
# abundance ratio on the original scale (2^20).
topenrich_CM2 <- log2FC_sum_Met05_CM2_sig_anno_majorsubclass %>% 
  filter(Direction == "Upregulated" & log2FC >= 20) %>% 
  arrange(desc(log2FC)) %>% 
  select(Group, Metabolite, Compound_Clean, 
         log2FC, Direction, Subclass, Major_Subclass) %>%
  rename("Code" = "Metabolite",
         "Compound" = "Compound_Clean")

# Export the highly enriched CM2 metabolite list.
Export(topenrich_CM2, "topenrich_CM2.txt")

# Identify metabolites shared between the highly enriched CM2 list
# and the previously selected top features.
intersect(topenrich_CM2$Code, top_features$Metabolite)



#---------------------------------------------------------
# 11.2 Review Major Subclass Distribution and Annotation
#---------------------------------------------------------

# List the major chemical categories represented among the
# significant metabolites.
log2FC_sum_Met05_CM_sig_anno_majorsubclass %>%
  distinct(Major_Subclass)

# Count significant metabolites without a resolved CID annotation.
# These records are retained in the analysis but flagged for
# annotation-quality assessment.
log2FC_sum_Met05_CM_sig_anno_majorsubclass %>% 
  filter(CID == "Not found") %>%
  count()



#---------------------------------------------------------
# 11.3 Visualize a Selected Major Subclass
#---------------------------------------------------------

# Generate a barplot for significant metabolites assigned to the
# "Oxygenated molecules" category.
#
# Only metabolites meeting the specified log2FC threshold are
# displayed by the plotting function. Statistical significance
# is represented by the q < 0.05 criterion used to define the
# input significant-metabolite dataset.
ggFC_sig_Met05_CM_subclass <- plot_log2FC_bar_subclass(
  
  summary_df = log2FC_sum_Met05_CM_sig_anno_majorsubclass %>% 
    filter(Major_Subclass == "Oxygenated molecules"),
  
  title = "Log2FC of metabolites (CM-CM vs. CM-CR)\nOxygenated molecules",
  
  xlab = "Significant Metabolite\n(|Log2FC| >= 1, q < 0.05)",
  
  fill_colors = c("Upregulated" = "steelblue3",
                  "Downregulated" = "#FC766A"),
  
  log2FC_cutoff = 1,
  
  strip_text_color = "black",
  
  group_color =  group_colors
)

# Display the major-subclass-specific barplot.
print(ggFC_sig_Met05_CM_subclass)









#=========================================================
# 12. Statistical Summary and Zero-Variance Assessment
#=========================================================
# Summarize the number of significant metabolites by
# experimental group and direction of change.
log2FC_sum_Met05_CM_sig_anno_majorsubclass %>%
  count(Group, Direction)

# Summarize significant metabolites by experimental group,
# direction of change, and metabolite major subclass.
# The results are ordered by group and descending number
# of significant metabolites within each category.
log2FC_sum_Met05_CM_sig_anno_majorsubclass %>% 
  group_by(Group, Direction, Major_Subclass) %>%
  count() %>%
  arrange(Group, desc(n)) %>%
  print(n = 50)



# Repeat the summary for the second metabolite dataset.
log2FC_sum_Met05_CM2_sig_anno_majorsubclass %>%
  count(Group, Direction)

# Summarize significant metabolites by experimental group,
# direction of change, and metabolite major subclass.
log2FC_sum_Met05_CM2_sig_anno_majorsubclass %>% 
  group_by(Group, Direction, Major_Subclass) %>%
  count() %>%
  arrange(Group, desc(n)) %>%
  print(n = 50)



#---------------------------------------------------------
# Identify Metabolites with Zero Variance
#---------------------------------------------------------
# Load packages used for data manipulation and reshaping.
library(dplyr)
library(tidyr)

# Identify metabolites with zero variance within a dataset.
# Standard deviation is calculated for all columns beginning
# with "Met". Metabolites with an SD of zero are returned.
get_zero_var <- function(df) {
  
  df %>%
    summarise(across(starts_with("Met"), ~ sd(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Metabolite", values_to = "SD") %>%
    filter(SD == 0) %>%
    pull(Metabolite)
}


# Identify metabolites with zero variance across all CM samples.
# This represents global zero variance within the CM dataset.
Met05_CM_zero <- Met05_CM %>%
  filter(Area == "CM") %>%
  get_zero_var()

# Identify metabolites with zero variance within the Akha-CM group.
Met05_Akha_zero <- Met05_CM %>%
  filter(Group == "Akha-CM") %>%
  get_zero_var()

# Identify metabolites with zero variance within the Lahu-CM group.
Met05_Lahu_zero <- Met05_CM %>%
  filter(Group == "Lahu-CM") %>%
  get_zero_var()

# Identify metabolites with zero variance within the Khuen-CM group.
Met05_Khuen_zero <- Met05_CM %>%
  filter(Group == "Khuen-CM") %>%
  get_zero_var()


#---------------------------------------------------------
# Identify Group-Specific Zero-Variance Metabolites
#---------------------------------------------------------
# Combine zero-variance metabolites from the two groups that
# are not being evaluated as group-specific zero variance.
zero_Lahu_Khuen <- c(Met05_Lahu_zero, Met05_Khuen_zero) # -> used to identify zero variance only in Akha-CM
zero_Akha_Khuen <- c(Met05_Akha_zero, Met05_Khuen_zero) # -> used to identify zero variance only in Lahu-CM
zero_Akha_Lahu <- c(Met05_Akha_zero, Met05_Lahu_zero) # -> used to identify zero variance only in Khuen-CM

# Identify metabolites with zero variance only in Akha-CM.
zero_Akha <- setdiff(Met05_Akha_zero, zero_Lahu_Khuen)

# Identify metabolites with zero variance only in Lahu-CM.
zero_Lahu <- setdiff(Met05_Lahu_zero, zero_Akha_Khuen)

# Identify metabolites with zero variance only in Khuen-CM.
zero_Khuen <- setdiff(Met05_Khuen_zero, zero_Akha_Lahu)


#---------------------------------------------------------
# Check Whether Group-Specific Zero-Variance Metabolites
# Are Present Among Significant Metabolites
#---------------------------------------------------------
# Identify metabolites with group-specific zero variance
# for each pairwise group comparison.
diff_Akha_Lahu <- union(zero_Akha, zero_Lahu)
diff_Akha_Khuen <- union(zero_Akha, zero_Khuen)
diff_Khuen_Lahu <- union(zero_Khuen, zero_Lahu)

# Combine all metabolites with group-specific zero variance
# across the CM groups.
diff_CM <- c(diff_Akha_Lahu, diff_Akha_Khuen, diff_Khuen_Lahu)

# Identify metabolites that are both group-specific zero-variance
# features and significant metabolites in the CM analysis.
sig_zero_var <- intersect(diff_CM, sig_CM_05)

# Create a summary table of significant metabolites that have
# zero variance in a specific CM group. The group in which zero
# variance occurs is recorded in the "zero_var" column, and
# metabolite annotation and statistical summary information are
# added by joining with the significant-metabolite table.
sigzero_any_CM <- bind_rows(
  tibble(Metabolite = zero_Akha, Group = "Akha-CM"),
  tibble(Metabolite = zero_Lahu, Group = "Lahu-CM"),
  tibble(Metabolite = zero_Khuen, Group = "Khuen-CM")
) %>%
  filter(Metabolite %in% sig_zero_var) %>%
  rename(zero_var = Group) %>%
  left_join(log2FC_sum_Met05_CM_sig_anno_majorsubclass,
            by = "Metabolite")

# Summarize significant zero-variance metabolites by
# experimental group and metabolite subclass.
sigzero_any_CM %>% 
  group_by(Group, Subclass) %>%
  count() %>%
  arrange(desc(n)) %>%
  print(n = 50)

# Count significant zero-variance metabolites according to
# the CM group in which zero variance was observed.
sigzero_any_CM %>% 
  group_by(zero_var) %>%
  count()
