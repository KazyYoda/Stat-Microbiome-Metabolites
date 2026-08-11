############################################################
# Metabolite Analysis - Chiang Rai (CR)
############################################################
# This section performs metabolite-level statistical analysis,
# including data preparation, log2 transformation, summary
# statistics, and log2 fold-change (log2FC) calculations for
# the CR study area.

setwd("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat")
load("~/Documents/HillTribe_NGS/5.Metabolites/Metabolites.RData")


#--------------------------------------------------------
# Load Required Packages
#--------------------------------------------------------
# Load packages used for data manipulation, statistical
# analysis, visualization, and heatmap generation.
library(dplyr)
library(tidyr)
library(car)
library(readxl)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggh4x)


#--------------------------------------------------------
# Define Group Colors
#--------------------------------------------------------
# Define consistent colors for the three ethnic groups
# included in the CR study area.
group_colors <- c(
  "Akha-CR" = "dodgerblue4",
  "Lahu-CR" = "coral1",
  "Lisu-CR" = "grey40"
)


#========================================================
# 4. Compute log2 Fold-Change (log2FC)
# Relative to Defined Contrast Group
#========================================================

#--------------------------------------------------------
# 4.1 Merge Metadata with Metabolite Matrix
#--------------------------------------------------------

# Confirm that the sample order in the metabolite matrix
# matches the sample order in the corresponding metadata.
# This check helps ensure that metadata are correctly aligned
# with metabolite measurements before combining the datasets.
identical(rownames(metabo_05_t), Met05_match$Sample_ID)

# Combine sample metadata with the metabolite abundance matrix.
Met05 <- cbind(Met05_match, metabo_05_t)

# Apply a log2 transformation to the metabolite abundance
# measurements. Columns 14 onward contain the metabolite data.
Met05_log2 <- Met05 %>%
  mutate(across(14:ncol(.), log2))



#--------------------------------------------------------
# 4.2 Subset Data by Ethnicity and Area
#--------------------------------------------------------

# Select samples from the Chiang Rai (CR) study area.
Met05_CR <- Met05_log2 %>% filter(Area == "CR")


#--------------------------------------------------------
# Summary Statistics
#--------------------------------------------------------
# Identify metabolite columns for calculation of summary
# statistics.
met_cols <- grep("^Met", names(Met05_CR), value = TRUE)

# Calculate the mean and standard deviation of each metabolite
# within each ethnic group and reshape the results to long format.
summary_CR_long <- Met05_CR %>%
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



# Reshape the summary statistics to wide format, with each
# ethnic group represented as a separate column.
#
# The table is then annotated with metabolite descriptions,
# Kruskal-Wallis results, and log2FC values from the specified
# pairwise contrasts.
summary_CR_wide <- summary_CR_long %>%
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
    CR_kruskal_metabo,
    by = "Metabolite"
  ) %>% 
  left_join(
    log2FC_sum_Met05_CR %>%
      filter(Group == "Akha-CR") %>%
      select(Metabolite, Log2FC_AkhavsLisu = log2FC),
    by = "Metabolite"
  ) %>%
  left_join(
    log2FC_sum_Met05_CR %>%
      filter(Group == "Lahu-CR") %>%
      select(Metabolite, Log2FC_LahuvsLisu = log2FC),
    by = "Metabolite"
  ) %>%
  left_join(
    log2FC_sum_Met05_CR2 %>%
      filter(Group == "Akha-CR") %>%
      select(Metabolite, Log2FC_AkhavsLahu = log2FC),
    by = "Metabolite"
  )
    

# Export the annotated CR summary table.
Export(summary_CR_wide, "summary_CR_wide.txt")


# Count significant metabolites based on an absolute log2FC
# threshold of 1 and a BH-adjusted q-value below 0.05.
# Results are summarized according to the direction of change.
summary_CR_wide %>%
  filter(abs(log2FC) >=1 & q_value < 0.05) %>%
  group_by(Direction) %>%
  count()


# Count the number of metabolite-group combinations for which
# the mean abundance on the log2 scale is zero.
# Here, a mean of zero indicates that the corresponding
# metabolite was not detected at the transformed scale.
summary_CR_long %>%
  group_by(Group) %>%
  summarise(
    Not_Detected = sum(Mean == 0, na.rm = TRUE)
  )

#========================================================
# 4.3 Log2FC Analysis — CR (Akha, Lahu, Lisu)
#========================================================

#--------------------------------------------------------
# 4.3.1 Compute log2FC Relative to Contrast Group
#--------------------------------------------------------

# Calculate group-level log2FC values using Lisu-CR as the
# contrast/reference group.
log2FC_Met05_CR <- log2FC_group_summary_relative_to_contrast(
  metabo_data = Met05_CR[-c(1:13)], 
  metadata = Met05_CR[1:13], 
  contrast = "Lisu-CR", 
  Group = "Group"
)

# Transpose the log2FC results to facilitate downstream
# formatting and annotation.
log2FC_Met05_CR_t <- log2FC_Met05_CR %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Code")



# Calculate group-level log2FC values using Lahu-CR as the
# contrast/reference group.
log2FC_Met05_CR2 <- log2FC_group_summary_relative_to_contrast(
  metabo_data = Met05_CR[-c(1:13)], 
  metadata = Met05_CR[1:13], 
  contrast = "Lahu-CR", 
  Group = "Group"
)

# Transpose the log2FC results to facilitate downstream
# formatting and annotation.
log2FC_Met05_CR2_t <- log2FC_Met05_CR2 %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Code")

#--------------------------------------------------------
# 4.3.2 Summarize log2FC by Group
#--------------------------------------------------------

# Summarize the log2FC results while retaining the
# comparison groups relative to Lisu-CR.
log2FC_sum_Met05_CR <- summarize_log2FC(
  log2FC_Met05_CR, 
  filter_group = "Lisu-CR"
)

# Add metabolite annotation to the log2FC summary table.
log2FC_sum_Met05_CR <- log2FC_sum_Met05_CR %>%
  left_join(
    metabo_05_descp, 
    by = c("Metabolite" = "Code")
  )

# Export the log2FC summary table for the Lisu-CR-referenced
# analysis.
Export(log2FC_sum_Met05_CR, "log2FC_sum_Met05_CR.txt")



# Summarize the log2FC results while retaining the
# comparison groups relative to Lahu-CR.
log2FC_sum_Met05_CR2 <- summarize_log2FC(
  log2FC_Met05_CR2, 
  filter_group = "Lahu-CR"
)

# Add metabolite annotation to the log2FC summary table.
log2FC_sum_Met05_CR2 <- log2FC_sum_Met05_CR2 %>%
  left_join(
    metabo_05_descp, 
    by = c("Metabolite" = "Code")
  )

# Export the log2FC summary table for the Lahu-CR-referenced
# analysis.
Export(log2FC_sum_Met05_CR2, "log2FC_sum_Met05_CR2.txt")








#========================================================
# 5. Log2FC Visualization — Barplot
#========================================================

# Generate a grouped barplot showing metabolite log2FC values
# for comparisons relative to the defined contrast groups.
ggFC_Met05_CR <- plot_log2FC_bar(
  summary_df = log2FC_sum_Met05_CR,
  title = "Log2FC of metabolites (Lisu-CR or Lahu-CR vs. Akha-CR)",
  xlab = "Metabolite",
  fill_colors = c("Upregulated" = "steelblue3",
                  "Downregulated" = "#FC766A"),
  strip_text_color = "black",
  log2FC_cutoff = 3,
  group_color = group_colors
)

# Display the log2FC barplot.
print(ggFC_Met05_CR)



#--------------------------------------------------------
# Log2FC Interpretation Guide
#--------------------------------------------------------
# Log2FC represents the log2-transformed fold change between
# the comparison group and the defined contrast/reference group.
#
# Upregulation:
# log2FC ≥ 1  → ≥ 2-fold increase
# log2FC ≥ 2  → ≥ 4-fold increase
# log2FC ≥ 3  → ≥ 8-fold increase
# log2FC ≥ 4  → ≥ 16-fold increase
#
# Fold change can be obtained from log2FC as:
# Fold change = 2^(log2FC)

# Downregulation:
# log2FC ≤ −1 → ≤ 0.5-fold of the contrast group
# log2FC ≤ −2 → ≤ 0.25-fold of the contrast group
# log2FC ≤ −3 → ≤ 0.125-fold of the contrast group
# log2FC ≤ −4 → ≤ 0.0625-fold of the contrast group



#--------------------------------------------------------
# Volcano Plot for CR
#--------------------------------------------------------
# Load ggrepel for non-overlapping metabolite labels.
library(ggrepel)

# Add comparison labels to the two log2FC result tables.
log2FC_CR1 <- log2FC_sum_Met05_CR %>% mutate(Comp = "vs Lisu-CR")
log2FC_CR2 <- log2FC_sum_Met05_CR2 %>% mutate(Comp = "vs Lahu-CR")

# Combine the two comparison tables.
log2FC_CR <- rbind(log2FC_CR1, log2FC_CR2)

# Add Kruskal-Wallis statistical results to the log2FC table.
vol_CR <- log2FC_CR %>%
  left_join(CR_kruskal_metabo, by = "Metabolite") 

# Prepare data for the volcano plot.
# The y-axis is calculated from the Kruskal-Wallis p-value.
# Metabolites are classified as significant when kruskal_p < 0.05.
vol_CR_plot <- vol_CR %>%
  mutate(
    logP = -log10(kruskal_p),   # Can be replaced with an adjusted p-value if required
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


# Select the 30 most statistically significant metabolites
# based on the smallest Kruskal-Wallis p-values for labeling.
top_features <- vol_CR_plot %>%
  filter(class == "Significant") %>%
  arrange(kruskal_p) %>%
  slice(1:30)

# Generate the volcano plot.
ggplot(vol_CR_plot, aes(x = log2FC, y = logP)) +
  geom_point(aes(color = Significant), alpha = 0.7, size = 1) +
  
  # Highlight the selected top significant metabolites.
  geom_point(
    data = top_features,
    aes(color = Significant),
    size = 1
  ) +
  
  # Label the selected top significant metabolites.
  geom_text_repel(
    data = top_features,
    aes(label = Metabolite),
    size = 2,
    max.overlaps = 100,
    segment.size = 0.2
  ) +
  
  # Define colors for metabolite direction and significance.
  scale_color_manual(values = c("Upregulated" = "#FC766A",
                                "Downregulated" = "#89ABE3",
                                "Not Significant" = "grey70")) +
  
  # Add log2FC thresholds corresponding to ±1.
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "darkgray") + # abs(log2FC) >=1
  
  # Add the p-value significance threshold at p = 0.05.
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "darkgray") + # significance cutoff
  
  labs(
    title = "Volcano plot of metabolites (CR group)",
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
# 6. Heatmap Visualization (log2-transformed Data)
#========================================================

#--------------------------------------------------------
# 6.1 Prepare Log2 Matrix
#--------------------------------------------------------

# Apply a log2 transformation to the complete metabolite matrix.
log_matrix_metabo_05 <- as.matrix(log2(metabo_05_t))

# Transpose the matrix so that metabolites are rows
# and samples are columns.
log_matrix_metabo_05_t <- t(log_matrix_metabo_05)



#--------------------------------------------------------
# 6.2 Select Significant Metabolites Based on CR Comparison
#--------------------------------------------------------
# Identify metabolites showing a significant overall
# difference among the CR groups based on the Kruskal-Wallis test.
sig_kruskal <- CR_kruskal_metabo %>%
  filter(kruskal_p < 0.05) %>%
  pull(Metabolite)

# Retain significant metabolites that are also present
# in the log2FC summary table.
sig_across <- log2FC_sum_Met05_CR_sig %>%
  filter(Metabolite %in% sig_kruskal)

# Select metabolites with an absolute log2FC of at least 1.
log2FC_X_CR <- log2FC_sum_Met05_CR_sig %>%
  filter(abs(log2FC) >= 1)

# Identify metabolites present in both the log2-transformed
# metabolite matrix and the selected log2FC results.
valid_rows <- intersect(
  rownames(log_matrix_metabo_05_t),
  log2FC_X_CR$Metabolite
)

# Identify CR samples present in both the metabolite matrix
# and the CR metadata.
valid_cols <- intersect(
  colnames(log_matrix_metabo_05_t),
  Met05_CR$Sample_ID
)

# Subset the metabolite matrix to the selected metabolites
# and CR samples while preserving matrix structure.
log_matrix_metabo_05_t_CR <-
  log_matrix_metabo_05_t[
    valid_rows,
    valid_cols,
    drop = FALSE
  ]



#--------------------------------------------------------
# 6.3 Variance Check (Pre-Clustering QC)
#--------------------------------------------------------

# Check the number of metabolites with zero variance.
# Zero-variance metabolites do not provide information for
# distance-based clustering.
sum(apply(log_matrix_metabo_05_t_CR, 1, sd, na.rm = TRUE) == 0)

# Check the number of samples with zero variance across
# the selected metabolites.
sum(apply(log_matrix_metabo_05_t_CR, 2, sd, na.rm = TRUE) == 0)

# Note: zero-variance features can affect Euclidean-distance
# calculations and downstream hierarchical clustering.



#--------------------------------------------------------
# 6.4 Define Group Order for Annotation
#--------------------------------------------------------

# Set the desired order of CR ethnic groups for heatmap
# annotations and visualization.
Met05_CR$Group <- factor(
  Met05_CR$Group,
  levels = c("Akha-CR", "Lahu-CR", "Lisu-CR")
)



#========================================================
# 6.5 Generate Heatmap
#========================================================

# Generate the heatmap using the selected log2-transformed
# metabolite matrix and corresponding sample metadata.
CR <- plot_heatmap(
  mat = log_matrix_metabo_05_t_CR,
  metadata = Met05_CR[1:13],
  legend_title = "Log2 (CR)",
  group_color = group_colors
)

# Draw the heatmap with both heatmap and annotation legends
# positioned on the left.
draw(
  CR,
  heatmap_legend_side = "left",
  annotation_legend_side = "left"
)







#========================================================
# 7. PCA Analysis — CR Metabolites
#========================================================

# Load packages for principal component analysis and
# PCA visualization.
library(FactoMineR)
library(factoextra)

#--------------------------------------------------------
# 7.1 Select Differential Metabolites (|log2FC| ≥ 1)
#--------------------------------------------------------

# Select metabolites showing an absolute log2FC of at least 1
# in the CR log2FC analysis for use as PCA features.
log2FC_1_CR <- log2FC_sum_Met05_CR %>%
  filter(abs(log2FC) >= 1)


#--------------------------------------------------------
# 7.2 Subset Log2 Matrix for PCA
#--------------------------------------------------------

# Identify CR samples present in both the log2-transformed
# metabolite matrix and the CR metadata.
valid_cols <- intersect(
  colnames(log_matrix_metabo_05_t),
  Met05_CR$Sample_ID
)

# Identify selected metabolites present in both the
# log2-transformed matrix and the differential-metabolite list.
valid_rows <- intersect(
  rownames(log_matrix_metabo_05_t),
  log2FC_1_CR$Metabolite
)

# Subset the matrix to the selected metabolites and CR samples.
# The resulting matrix has metabolites as rows and samples
# as columns.
log_matrix_metabo_05_t_CR_pca <-
  log_matrix_metabo_05_t[
    valid_rows,
    valid_cols,
    drop = FALSE
  ]


#--------------------------------------------------------
# 7.3 Scale Data Prior to PCA
#--------------------------------------------------------

# Standardize each metabolite to mean = 0 and SD = 1
# before PCA so that metabolites are placed on a comparable scale.
pca_CR <- scale(log_matrix_metabo_05_t_CR_pca)
head(pca_CR)

# Replace NaN values generated by scaling features with zero
# variance with 0.
pca_CR[is.nan(pca_CR)] <- 0


#--------------------------------------------------------
# 7.4 Perform PCA
#--------------------------------------------------------

# Transpose the matrix so that samples are rows and metabolites
# are columns, as required for PCA.
#
# scale.unit = FALSE is used because the metabolite data have
# already been standardized in the previous step.
res_pca <- PCA(
  t(pca_CR),
  scale.unit = FALSE,  # Data already standardized
  graph = FALSE
)


#--------------------------------------------------------
# 7.5 Prepare Metadata for Visualization
#--------------------------------------------------------

# Extract sample identifiers corresponding to the PCA observations.
sample_names <- colnames(log_matrix_metabo_05_t_CR_pca)

# Extract ethnic-group labels for each CR sample.
group_labels <- Met05_CR$Group

# Create a metadata table linking each sample to its
# corresponding ethnic group.
metadata <- data.frame(
  Sample = sample_names,
  Group = group_labels
)

# Use sample IDs as row names for the metadata table.
rownames(metadata) <- sample_names


#--------------------------------------------------------
# 7.6 PCA Plot with Group Ellipses
#--------------------------------------------------------

# Define the color palette for the three CR ethnic groups.
group_colors <- c(
  "Akha-CR" = "dodgerblue4",
  "Lahu-CR" = "coral1",
  "Lisu-CR" = "grey40"
)


# Visualize individual samples in PCA space and display
# 95% confidence ellipses for each ethnic group.
fviz_pca_ind(
  res_pca,
  geom.ind = "point",
  col.ind = metadata$Group,     # Color samples by ethnic group
  palette = group_colors,
  addEllipses = TRUE,           # Add 95% confidence ellipses
  ellipse.type = "confidence",
  ellipse.level = 0.95,
  repel = TRUE                  # Reduce overlap among labels
) +
  labs(title = "PCA of CR Metabolites (|Log2FC| >= 1)") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    title = element_text(size = 8)
  )








#========================================================
# 8. Normality and Homogeneity of Variance Testing
#========================================================

setwd("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat")

#--------------------------------------------------------
# 8.1 Subset CR Metabolite Matrix
#--------------------------------------------------------

# Identify metabolite columns present in both the log2-transformed
# metabolite matrix and the CR log2FC results.
valid_cols <- intersect(
  colnames(log_matrix_metabo_05),
  log2FC_sum_Met05_CR$Metabolite
)

# Identify CR samples present in both the log2-transformed
# metabolite matrix and the CR metadata.
valid_rows <- intersect(
  rownames(log_matrix_metabo_05),
  Met05_CR$Sample_ID
)

# Subset the matrix to CR samples and the selected metabolites.
# The resulting matrix contains samples as rows and metabolites
# as columns.
log_matrix_metabo_05_CR <-
  log_matrix_metabo_05[
    valid_rows,
    valid_cols,
    drop = FALSE
  ]


#--------------------------------------------------------
# 8.2 Identify and Remove Zero-Variance Metabolites
#--------------------------------------------------------

# Count metabolites with zero variance across the CR samples.
sum(apply(Met05_CR[-c(1:13)], 2, sd, na.rm = TRUE) == 0)

# List metabolites with zero variance across the CR samples.
colnames(Met05_CR[-c(1:13)])[
  apply(Met05_CR[-c(1:13)], 2, sd, na.rm = TRUE) == 0
]

# Remove metabolites with zero variance from the metabolite
# dataset. Zero-variance features cannot be evaluated using
# variance-based statistical procedures.
log_matrix_metabo_05_CR_nonzero <-
  Met05_CR[-c(1:13)][,
                     apply(Met05_CR[-c(1:13)], 2, sd, na.rm = TRUE) != 0,
                     drop = FALSE
  ]


#--------------------------------------------------------
# 8.3 Normality and Variance Homogeneity Tests
#--------------------------------------------------------

# Assess distributional characteristics and homogeneity of
# variance for the CR metabolite features.
#
# The resulting table is used to document the assumption checks
# before selecting the appropriate group-comparison procedure.
CR_norm_var <- norm_var_check(
  Met05_CR,
  metabo_05_descp,
  tag = "metabo05_CR"
)

CR_norm_var

# Count metabolites for which the variance-test p-value is missing.
CR_norm_var %>%
  count(is.na(Var_p.value))

# Export the normality and variance-homogeneity test results.
Export(CR_norm_var, "CR_norm_var.txt")



#========================================================
# 9. Statistical Analysis — Group Comparison
#========================================================

#--------------------------------------------------------
# 9.1 Kruskal-Wallis Test (CR)
#--------------------------------------------------------

# Perform the Kruskal-Wallis test to assess differences in
# metabolite abundance among the three CR ethnic groups.
#
# The normality and variance-assessment results are supplied
# to the analysis function together with the CR metabolite
# matrix and corresponding sample metadata.
CR_kruskal_metabo <- kruskal_metabo(metabo_norm_var = CR_norm_var,
                                    metabo_data = log_matrix_metabo_05_CR,
                                    metadata = Met05_CR[1:13], 
                                    tag = "metabo_CR")

head(CR_kruskal_metabo)



#----------------------------------------------------------
# 9.2 Dunn's Test — Post-hoc Test for Kruskal-Wallis (CR)
#----------------------------------------------------------

# Perform Dunn's post-hoc test for pairwise comparisons among
# CR ethnic groups following the Kruskal-Wallis analysis.
CR_dunn <- dunnTest_metabo(krus_pvalue_descp = CR_kruskal_metabo, 
                           metabo_data = log_matrix_metabo_05_CR, 
                           metadata = Met05_CR[1:13], 
                           metabo_descp = metabo_05_descp, 
                           tag = "metabo_CR")

head(CR_dunn)

#--------------------------------------------------------
# 9.2 Filter Significant Metabolites
#--------------------------------------------------------

# Extract metabolites showing significant pairwise differences
# for the specified comparisons based on adjusted Dunn's
# test p-values (P.adj < 0.05).
sig_CR_05 <- CR_dunn %>% 
  filter(Comparison %in% c("Akha-CR - Lisu-CR",
                           "Lahu-CR - Lisu-CR") &
           P.adj < 0.05) %>%
  pull(Metabolite)

# Retain significant metabolites with an absolute log2FC of
# at least 1 and add metabolite classification information.
log2FC_sum_Met05_CR_sig <- log2FC_sum_Met05_CR %>%
  filter(Metabolite %in% sig_CR_05) %>%
  filter(abs(log2FC) >= 1) %>%
  left_join(Pub_Classy_05_clean,
            by = c("Compound",
                   "Compound_Clean",
                   "Metabolite" = "Code"))

# Export the significant metabolite results for comparisons
# relative to Lisu-CR.
Export(log2FC_sum_Met05_CR_sig, "relLisuCR_sig_log2FC_1.txt")



# Extract metabolites showing significant pairwise differences
# for the specified comparisons based on adjusted Dunn's
# test p-values (P.adj < 0.05).
sig_CR2_05 <- CR_dunn %>% 
  filter(Comparison %in% c("Akha-CR - Lahu-CR",
                           "Lahu-CR - Lisu-CR") &
           P.adj < 0.05) %>%
  pull(Metabolite)

# Retain significant metabolites with an absolute log2FC of
# at least 1 and add metabolite classification information.
log2FC_sum_Met05_CR2_sig <- log2FC_sum_Met05_CR2 %>%
  filter(Metabolite %in% sig_CR2_05) %>%
  filter(abs(log2FC) >= 1) %>%
  left_join(Pub_Classy_05_clean,
            by = c("Compound",
                   "Compound_Clean",
                   "Metabolite" = "Code"))

# Export the significant metabolite results for comparisons
# relative to Lahu-CR.
Export(log2FC_sum_Met05_CR2_sig, "relLahuCR_sig_log2FC_1.txt")


#========================================================
# 10. Visualization — Significant Metabolites
#========================================================

# Generate a barplot of significant metabolites selected using
# adjusted Dunn's test p-value < 0.05 and absolute log2FC >= 1.
ggFC_sig_Met05_CR <- plot_log2FC_bar(
  summary_df = log2FC_sum_Met05_CR_sig,
  title = "Log2FC of metabolites (Lahu-CR or Lisu-CR vs. Akha-CR)",
  xlab = "Significant Metabolite\n(|Log2FC| >= 2, P.adj < 0.05)",
  fill_colors = c("Upregulated" = "steelblue3",
                  "Downregulated" = "#FC766A"),
  strip_text_color = "black",
  log2FC_cutoff = 2,
  group_color = group_colors
)

# Display the significant-metabolite barplot.
print(ggFC_sig_Met05_CR)





#=========================================================
# 11. Visualization — Major Subclass Annotation
#=========================================================
# Load the CR metabolite analysis workspace containing
# metabolite annotations and previously generated results.
load("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat/Metabolites_CR.RData")


#---------------------------------------------------------
# 11.1 Assign Major Metabolite Subclasses
#---------------------------------------------------------
# Group detailed metabolite subclasses into broader chemical
# categories to facilitate interpretation and visualization.
log2FC_sum_Met05_CR_sig_anno_majorsubclass <- 
  log2FC_sum_Met05_CR_sig %>%
  mutate(Major_Subclass = case_when(
    
    # Lipids and fatty acid derivatives
    Subclass %in% c(
      "Fatty acid esters",
      "Fatty acids and conjugates",
      "Carboxylic acid derivatives"
    ) ~ "Lipids & FA derivatives",
    
    # Terpenoid compounds
    Subclass %in% c(
      "Sesquiterpenoids",
      "Diterpenoids"
    ) ~ "Terpenoids",
    
    # Phenolic and aromatic compounds
    Subclass %in% c(
      "Benzoic acids and derivatives",
      "Benzoyl derivatives",
      "Benzenediols",
      "Cresols"
    ) ~ "Phenolic & aromatic",
    
    # Nitrogen-containing compounds
    Subclass %in% c(
      "Amines",
      "Ureas",
      "Morpholines",
      "Organic cyanides"
    ) ~ "Nitrogen-containing",
    
    # Oxygen-containing compounds
    Subclass %in% c(
      "Alcohols and polyols",
      "Carbonyl compounds",
      "Gamma butyrolactones",
      "Carbonic acid diesters"
    ) ~ "Oxygenated molecules",
    
    # Hydrocarbon compounds
    Subclass %in% c(
      "Alkanes",
      "Cycloalkanes",
      "Olefins"
    ) ~ "Hydrocarbons",
    
    # Silicon-containing compounds
    Subclass %in% c(
      "Organosilicon compounds"
    ) ~ "Silicon-containing",
    
    # Retain an explicit category for subclasses that are not
    # included in the predefined major-subclass classification.
    TRUE ~ "Unknown"
    
  ))



# Apply the same major-subclass classification to the
# alternative log2FC comparison using Lahu-CR as the
# contrast/reference group.
log2FC_sum_Met05_CR2_sig_anno_majorsubclass <- 
  log2FC_sum_Met05_CR2_sig %>%
  mutate(Major_Subclass = case_when(
    
    # Lipids and fatty acid derivatives
    Subclass %in% c(
      "Fatty acid esters",
      "Fatty acids and conjugates",
      "Carboxylic acid derivatives"
    ) ~ "Lipids & FA derivatives",
    
    # Terpenoid compounds
    Subclass %in% c(
      "Sesquiterpenoids",
      "Diterpenoids"
    ) ~ "Terpenoids",
    
    # Phenolic and aromatic compounds
    Subclass %in% c(
      "Benzoic acids and derivatives",
      "Benzoyl derivatives",
      "Benzenediols",
      "Cresols"
    ) ~ "Phenolic & aromatic",
    
    # Nitrogen-containing compounds
    Subclass %in% c(
      "Amines",
      "Ureas",
      "Morpholines",
      "Organic cyanides"
    ) ~ "Nitrogen-containing",
    
    # Oxygen-containing compounds
    Subclass %in% c(
      "Alcohols and polyols",
      "Carbonyl compounds",
      "Gamma butyrolactones",
      "Carbonic acid diesters"
    ) ~ "Oxygenated molecules",
    
    # Hydrocarbon compounds
    Subclass %in% c(
      "Alkanes",
      "Cycloalkanes",
      "Olefins"
    ) ~ "Hydrocarbons",
    
    # Silicon-containing compounds
    Subclass %in% c(
      "Organosilicon compounds"
    ) ~ "Silicon-containing",
    
    # Retain an explicit category for subclasses that are not
    # included in the predefined major-subclass classification.
    TRUE ~ "Unknown"
    
  ))



#---------------------------------------------------------
# 11.2 Identify Highly Enriched Metabolites
#---------------------------------------------------------
# Identify strongly upregulated metabolites with log2FC >= 10.
# Results are ordered from the largest to smallest log2FC and
# include metabolite annotation and major-subclass information.
topenrich_CR <- log2FC_sum_Met05_CR_sig_anno_majorsubclass %>% 
  filter(Direction == "Upregulated" & log2FC >= 10) %>% 
  arrange(desc(log2FC)) %>% 
  select(Group, Metabolite, Compound_Clean, 
         log2FC, Direction, Subclass, Major_Subclass) %>%
  rename("Code" = "Metabolite",
         "Compound" = "Compound_Clean")

# Export the highly enriched metabolite list.
Export(topenrich_CR, "topenrich_CR.txt")

# Check which highly enriched metabolites are also present
# among the previously selected top features.
intersect(topenrich_CR$Code, top_features$Metabolite)


# Repeat the identification of highly enriched metabolites
# for the alternative CR comparison.
topenrich_CR2 <- log2FC_sum_Met05_CR2_sig_anno_majorsubclass %>% 
  filter(Direction == "Upregulated" & log2FC >= 10) %>% 
  arrange(desc(log2FC)) %>% 
  select(Group, Metabolite, Compound_Clean, 
         log2FC, Direction, Subclass, Major_Subclass) %>%
  rename("Code" = "Metabolite",
         "Compound" = "Compound_Clean")

# Export the highly enriched metabolite list for the
# alternative comparison.
Export(topenrich_CR2, "topenrich_CR2.txt")

# Check which highly enriched metabolites are also present
# among the previously selected top features.
intersect(topenrich_CR2$Code, top_features$Metabolite)



#---------------------------------------------------------
# 11.3 Inspect Major Subclass Categories and Annotations
#---------------------------------------------------------

# List the unique major-subclass categories represented
# among the significant CR metabolites.
log2FC_sum_Met05_CR_sig_anno_majorsubclass %>%
  distinct(Major_Subclass)

# Count significant metabolites for which a CID annotation
# was not found.
log2FC_sum_Met05_CR_sig_anno_majorsubclass %>% 
  filter(CID == "Not found") %>%
  count()


#---------------------------------------------------------
# 11.4 Visualize a Selected Major Subclass
#---------------------------------------------------------

# Generate a barplot for metabolites classified as
# "Oxygenated molecules".
#
# The input table contains metabolites previously identified
# as significant and filtered using the specified log2FC
# criterion in the preceding analysis.
ggFC_sig_Met05_CR_subclass <- plot_log2FC_bar_subclass(
  
  summary_df = log2FC_sum_Met05_CR_sig_anno_majorsubclass %>% 
    filter(Major_Subclass == "Oxygenated molecules"),
  
  title = "Log2FC of metabolites (xx-CR vs. xx-CR)\nOxygenated molecules",
  
  xlab = "Significant Metabolite\n(|Log2FC| >= 1, P.adj < 0.05)",
  
  fill_colors = c("Upregulated" = "steelblue3",
                  "Downregulated" = "#FC766A"),
  
  log2FC_cutoff = 1,
  
  strip_text_color = "black",
  
  group_color =  group_colors
)

# Display the major-subclass barplot.
print(ggFC_sig_Met05_CR_subclass)





#=========================================================
# 12. Stat-Summary
#=========================================================

#---------------------------------------------------------
# 12.1 Summarize Significant Metabolites
#---------------------------------------------------------

# Count significant metabolites by group and direction of change
# for the Lisu-CR-referenced analysis.
log2FC_sum_Met05_CR_sig %>%
  count(Group, Direction)

# Summarize significant metabolites by group, direction of change,
# and metabolite subclass.
log2FC_sum_Met05_CR_sig %>% 
  group_by(Group, Direction, Subclass) %>%
  count() %>%
  arrange(Group, desc(n)) %>%
  print(n = 50)


# Count significant metabolites by group and direction of change
# for the Lahu-CR-referenced analysis.
log2FC_sum_Met05_CR2_sig %>%
  count(Group, Direction)

# Summarize significant metabolites by group, direction of change,
# and metabolite subclass for the alternative comparison.
log2FC_sum_Met05_CR2_sig %>% 
  group_by(Group, Direction, Subclass) %>%
  count() %>%
  arrange(Group, desc(n)) %>%
  print(n = 50)



#---------------------------------------------------------
# 12.2 Identify Metabolites with Zero Variance
#---------------------------------------------------------

# Load packages used for data manipulation and reshaping.
library(dplyr)
library(tidyr)

# Identify metabolites with zero variance within a dataset.
# Standard deviation is calculated for all columns beginning
# with "Met". Metabolites with SD = 0 are returned.
get_zero_var <- function(df) {
  
  df %>%
    summarise(across(starts_with("Met"), ~ sd(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Metabolite", values_to = "SD") %>%
    filter(SD == 0) %>%
    pull(Metabolite)
}


# Identify metabolites with zero variance across all CR groups.
# These represent globally invariant metabolites within the
# CR dataset.
Met05_CR_zero <- Met05_CR %>%
  filter(Area == "CR") %>%
  get_zero_var()

# Identify metabolites with zero variance within Akha-CR.
Met05_Akha_zero <- Met05_CR %>%
  filter(Group == "Akha-CR") %>%
  get_zero_var()

# Identify metabolites with zero variance within Lahu-CR.
Met05_Lahu_zero <- Met05_CR %>%
  filter(Group == "Lahu-CR") %>%
  get_zero_var()

# Identify metabolites with zero variance within Lisu-CR.
Met05_Lisu_zero <- Met05_CR %>%
  filter(Group == "Lisu-CR") %>%
  get_zero_var()


#---------------------------------------------------------
# 12.3 Identify Group-Specific Zero-Variance Metabolites
#---------------------------------------------------------

# Combine zero-variance metabolite lists from pairs of groups.
# These vectors are used to identify metabolites with zero
# variance specifically within one group.
zero_Lahu_Lisu <- c(Met05_Lahu_zero, Met05_Lisu_zero) # -> used to identify zero variance only in Akha-CR
zero_Akha_Lisu <- c(Met05_Akha_zero, Met05_Lisu_zero) # -> used to identify zero variance only in Lahu-CR
zero_Akha_Lahu <- c(Met05_Akha_zero, Met05_Lahu_zero) # -> used to identify zero variance only in Lisu-CR

# Identify metabolites with zero variance in Akha-CR but
# not zero variance across all CR groups.
zero_Akha <- setdiff(Met05_Akha_zero, Met05_CR_zero)

# Identify metabolites with zero variance in Lahu-CR but
# not zero variance across all CR groups.
zero_Lahu <- setdiff(Met05_Lahu_zero, Met05_CR_zero)

# Identify metabolites with zero variance in Lisu-CR but
# not zero variance across all CR groups.
zero_Lisu <- setdiff(Met05_Lisu_zero, Met05_CR_zero)


#---------------------------------------------------------
# 12.4 Identify Significant Metabolites with Group-Specific
# Zero Variance
#---------------------------------------------------------

# Identify the union of group-specific zero-variance metabolites
# for each pair of CR groups.
diff_Akha_Lahu <- union(zero_Akha, zero_Lahu)
diff_Akha_Lisu <- union(zero_Akha, zero_Lisu)
diff_Lisu_Lahu <- union(zero_Lisu, zero_Lahu)

# Combine group-specific zero-variance metabolites across
# all CR group pairs.
diff_CR <- c(diff_Akha_Lahu, diff_Akha_Lisu, diff_Lisu_Lahu)

# Identify metabolites that have group-specific zero variance
# and are also significant in the Lisu-CR-referenced analysis.
sig_zero_var <- intersect(diff_CR, sig_CR_05)

# Create a table identifying the group in which zero variance
# occurs and attach metabolite annotation from the significant
# metabolite results.
sigzero_any_CR <- bind_rows(
  tibble(Metabolite = zero_Akha, Group = "Akha-CR"),
  tibble(Metabolite = zero_Lahu, Group = "Lahu-CR"),
  tibble(Metabolite = zero_Lisu, Group = "Lisu-CR")
) %>%
  filter(Metabolite %in% sig_zero_var) %>%
  rename(zero_var = Group) %>%
  left_join(log2FC_sum_Met05_CR_sig_anno_majorsubclass,
            by = "Metabolite",
            relationship = "many-to-many")

# Summarize significant zero-variance metabolites by group
# and metabolite subclass.
sigzero_any_CR %>% 
  group_by(Group, Subclass) %>%
  count() %>%
  arrange(desc(n)) %>%
  print(n = 50)

# Count significant zero-variance metabolites according to
# the group in which zero variance was observed.
sigzero_any_CR %>% 
  group_by(zero_var) %>%
  count()



#---------------------------------------------------------
# 12.5 Repeat Zero-Variance Assessment for the Alternative
# Significant Metabolite Set
#---------------------------------------------------------

# Create a corresponding table for significant metabolites from
# the Lahu-CR-referenced analysis.
sigzero_any_CR2 <- bind_rows(
  tibble(Metabolite = zero_Akha, Group = "Akha-CR"),
  tibble(Metabolite = zero_Lahu, Group = "Lahu-CR"),
  tibble(Metabolite = zero_Lisu, Group = "Lisu-CR")
) %>%
  filter(Metabolite %in% sig_zero_var) %>%
  rename(zero_var = Group) %>%
  left_join(log2FC_sum_Met05_CR2_sig_anno_majorsubclass,
            by = "Metabolite",
            relationship = "many-to-many")

# Summarize significant zero-variance metabolites by group
# and metabolite subclass.
sigzero_any_CR2 %>% 
  group_by(Group, Subclass) %>%
  count() %>%
  arrange(desc(n)) %>%
  print(n = 50)

# Count significant zero-variance metabolites according to
# the group in which zero variance was observed.
sigzero_any_CR2 %>% 
  group_by(zero_var) %>%
  count()


#---------------------------------------------------------
# 12.6 Create Combined Zero-Variance Reference Table
#---------------------------------------------------------

# Create group-specific data frames identifying metabolites
# with zero variance in each CR group.
zero_Akha_df <- tibble(Metabolite = zero_Akha, zero_var = "Akha-CR")
zero_Lahu_df <- tibble(Metabolite = zero_Lahu, zero_var = "Lahu-CR")
zero_Lisu_df <- tibble(Metabolite = zero_Lisu, zero_var = "Lisu-CR")

# Combine group-specific zero-variance metabolites into a
# single reference table for downstream analysis.
Zero_CR_df <- rbind(zero_Akha_df,
                    zero_Lahu_df, 
                    zero_Lisu_df)
