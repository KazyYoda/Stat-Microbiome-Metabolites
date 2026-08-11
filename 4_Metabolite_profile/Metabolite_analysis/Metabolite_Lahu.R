############################################################
# Metabolite Analysis — Lahu
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

# Define colors for the Lahu geographic groups used in
# downstream metabolite visualizations.
group_colors <- c(
  "Lahu-CM" = "lightblue2",
  "Lahu-CR" = "coral1"
)


#========================================================
# 4. Compute log2 Fold-Change (log2FC)
#    Relative to Defined Contrast Group
#========================================================

# This section calculates metabolite log2 fold-changes
# relative to a predefined contrast group and generates
# descriptive summary tables for the Lahu-CM versus
# Lahu-CR comparison.


#--------------------------------------------------------
# 4.1 Merge Metadata with Metabolite Matrix
#--------------------------------------------------------

# Confirm that the sample order in the metabolite matrix
# matches the corresponding sample identifiers in the
# metadata before combining the datasets.
identical(rownames(metabo_05_t), Met05_match$Sample_ID)

# Combine sample metadata with the metabolite abundance matrix.
Met05 <- cbind(Met05_match, metabo_05_t)

# Apply a log2 transformation to the metabolite abundance
# columns. Columns 14 onward contain the metabolite
# measurements.
Met05_log2 <- Met05 %>%
  mutate(across(14:ncol(.), log2))



#--------------------------------------------------------
# 4.2 Subset Data by Ethnicity and Area
#--------------------------------------------------------

# Retain samples belonging to the Lahu ethnicity for the
# Lahu-CM versus Lahu-CR comparison.
Met05_Lahu <- Met05_log2 %>% filter(Ethnicity == "Lahu")


# Identify metabolite columns for group-wise summarization.
met_cols <- grep("^Met", names(Met05_Lahu), value = TRUE)

# Calculate the mean and standard deviation of each
# log2-transformed metabolite within each geographic group.
# The results are converted to long format for downstream
# table construction.
summary_Lahu_long <- Met05_Lahu %>%
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
# 4.2.1 Generate Wide Summary Table
#--------------------------------------------------------

# Format group-specific mean and standard deviation values
# as mean ± SD and reshape the table so that each geographic
# group is represented by a separate column.
#
# Metabolite annotations, Wilcoxon test results, and log2FC
# information are then added to create a comprehensive
# metabolite summary table.
summary_Lahu_wide <- summary_Lahu_long %>%
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
    Lahu_wilcox,
    by = "Metabolite"
  ) %>% 
  rename(
    q_value = p_adj
  ) %>%
  mutate(`-log10(q-value)` = -log10(q_value),
         Test = "wilcox") %>%
  left_join(
    log2FC_sum_Met05_Lahu %>%
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
Export(summary_Lahu_wide, "summary_Lahu_wide.txt")


# Count significant metabolites based on an absolute log2FC
# of at least 1 and an FDR-adjusted p-value below 0.05,
# summarized according to the direction of change.
summary_Lahu_wide %>%
  filter(abs(log2FC) >=1 & q_value < 0.05) %>%
  group_by(Direction) %>%
  count()


# Count the number of metabolite-group combinations with
# a mean log2 abundance of zero. These may represent
# metabolites that were not detected in a group.
summary_Lahu_long %>%
  group_by(Group) %>%
  summarise(
    Not_Detected = sum(Mean == 0, na.rm = TRUE)
  )


#========================================================
# 4.3 Log2FC Analysis — Lahu (Lahu-CM vs Lahu-CR)
#========================================================

#--------------------------------------------------------
# 4.3.1 Compute log2FC Relative to Contrast Group
#--------------------------------------------------------

# Calculate metabolite log2 fold-changes using Lahu-CR
# as the reference (contrast) group.
#
# Positive log2FC values indicate higher metabolite
# abundance in the comparison group relative to Lahu-CR,
# whereas negative values indicate lower abundance.
log2FC_Met05_Lahu <- log2FC_group_summary_relative_to_contrast(
  metabo_data = Met05_Lahu[-c(1:13)], 
  metadata = Met05_Lahu[1:13], 
  contrast = "Lahu-CR", 
  Group = "Group"
)

# Transpose the log2FC results and retain metabolite codes
# as an explicit column for downstream data manipulation.
log2FC_Met05_Lahu_t <- log2FC_Met05_Lahu %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Code")


#--------------------------------------------------------
# 4.3.2 Summarize log2FC by Group
#--------------------------------------------------------

# Summarize metabolite-level log2FC values using Lahu-CR
# as the reference group.
log2FC_sum_Met05_Lahu <- summarize_log2FC(
  log2FC_Met05_Lahu, 
  filter_group = "Lahu-CR"
)

# Add metabolite annotation information to the log2FC
# summary table.
log2FC_sum_Met05_Lahu <- log2FC_sum_Met05_Lahu %>%
  left_join(
    metabo_05_descp, 
    by = c("Metabolite" = "Code")
  )

# Export the metabolite log2FC summary table.
Export(log2FC_sum_Met05_Lahu, "log2FC_sum_Met05_Lahu.txt")






#========================================================
# 5. Log2FC Visualization — Barplot
#========================================================

# Generate a grouped barplot showing metabolite log2FC
# values for the Lahu-CM versus Lahu-CR comparison.
ggFC_Met05_Lahu <- plot_log2FC_bar(
  summary_df = log2FC_sum_Met05_Lahu,
  title = "Log2FC of metabolites (Lahu-CM vs. Lahu-CR)",
  xlab = "Metabolite",
  fill_colors = c("Upregulated" = "steelblue3",
                  "Downregulated" = "#FC766A"),
  group_color = "Lahu-CM", 
  strip_text_color = "black",
  log2FC_cutoff = 2
)

# Display the log2FC barplot.
print(ggFC_Met05_Lahu)



#--------------------------------------------------------
# Log2FC Interpretation Guide
#--------------------------------------------------------

# Log2FC values are interpreted relative to the defined
# contrast/reference group (Lahu-CR).
#
# Positive log2FC values indicate higher metabolite
# abundance in Lahu-CM relative to Lahu-CR.
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



#========================================================
# 6. Heatmap Visualization (Log2-Transformed Data)
#========================================================

# Generate a heatmap of selected metabolites using
# log2-transformed metabolite abundance values.
# Metabolites are displayed as rows and Lahu samples
# as columns.


#--------------------------------------------------------
# 6.1 Prepare Log2 Matrix
#--------------------------------------------------------

# Log2-transform the complete metabolite abundance matrix.
log_matrix_metabo_05 <- as.matrix(log2(metabo_05_t))

# Transpose the matrix so that metabolites are represented
# as rows and samples are represented as columns.
log_matrix_metabo_05_t <- t(log_matrix_metabo_05)



#--------------------------------------------------------
# 6.2 Select Significant Metabolites (|log2FC| ≥ 1)
#--------------------------------------------------------

# Select significant metabolites meeting the predefined
# absolute log2FC threshold for inclusion in the heatmap.
log2FC_X_Lahu <- log2FC_sum_Met05_Lahu_sig %>%
  filter(abs(log2FC) >= 1)

# Identify metabolites present in both the log2-transformed
# abundance matrix and the selected metabolite list.
valid_rows <- intersect(
  rownames(log_matrix_metabo_05_t),
  log2FC_X_Lahu$Metabolite
)

# Identify Lahu samples present in both the metabolite
# matrix and the corresponding metadata.
valid_cols <- intersect(
  colnames(log_matrix_metabo_05_t),
  Met05_Lahu$Sample_ID
)

# Subset the matrix to the selected metabolites and Lahu
# samples while preserving the matrix structure.
log_matrix_metabo_05_t_Lahu <-
  log_matrix_metabo_05_t[
    valid_rows,
    valid_cols,
    drop = FALSE
  ]



#--------------------------------------------------------
# 6.3 Variance Check (Pre-Clustering QC)
#--------------------------------------------------------

# Check for metabolites with zero variance across the
# selected Lahu samples.
sum(apply(log_matrix_metabo_05_t_Lahu, 1, sd, na.rm = TRUE) == 0)

# Check for samples with zero variance across the
# selected metabolites.
sum(apply(log_matrix_metabo_05_t_Lahu, 2, sd, na.rm = TRUE) == 0)

# Note: zero-variance features or samples can affect
# distance-based clustering, including Euclidean clustering.



#--------------------------------------------------------
# 6.4 Define Group Order for Annotation
#--------------------------------------------------------

# Set the order of the geographic groups to ensure
# consistent ordering of sample annotations in the heatmap.
Met05_Lahu$Group <- factor(
  Met05_Lahu$Group,
  levels = c("Lahu-CM", "Lahu-CR")
)



#========================================================
# 6.5 Generate Heatmap
#========================================================

# Generate the heatmap using the selected log2-transformed
# metabolite matrix and corresponding Lahu sample metadata.
Lahu <- plot_heatmap(
  mat = log_matrix_metabo_05_t_Lahu,
  metadata = Met05_Lahu[1:13],
  legend_title = "Log2 (Lahu)",
  group_color = group_colors
)

# Draw the heatmap with the heatmap and annotation legends
# positioned on the left.
draw(
  Lahu,
  heatmap_legend_side = "left",
  annotation_legend_side = "left"
)







#========================================================
# 7. PCA Analysis — Lahu Metabolites
#========================================================

# Perform principal component analysis (PCA) on selected
# Lahu metabolite features to visualize overall sample
# variation and separation between Lahu-CM and Lahu-CR.


# Load packages required for PCA and PCA visualization.
library(FactoMineR)
library(factoextra)


#--------------------------------------------------------
# 7.1 Select Differential Metabolites (|log2FC| ≥ 1)
#--------------------------------------------------------

# Select metabolites with an absolute log2 fold-change of
# at least 1 for inclusion in the PCA.
#
# This corresponds to a minimum 2-fold difference between
# Lahu-CM and Lahu-CR.
log2FC_1_Lahu <- log2FC_sum_Met05_Lahu %>%
  filter(abs(log2FC) >= 1)


#--------------------------------------------------------
# 7.2 Subset Log2 Matrix for PCA
#--------------------------------------------------------

# Identify Lahu samples that are present in both the
# log2-transformed metabolite matrix and the metadata.
valid_cols <- intersect(
  colnames(log_matrix_metabo_05_t),
  Met05_Lahu$Sample_ID
)

# Identify selected metabolites that are present in both
# the log2-transformed matrix and the differential
# metabolite list.
valid_rows <- intersect(
  rownames(log_matrix_metabo_05_t),
  log2FC_1_Lahu$Metabolite
)

# Subset the matrix to the selected metabolites and Lahu
# samples. The resulting matrix contains metabolites as
# rows and samples as columns.
log_matrix_metabo_05_t_Lahu_pca <-
  log_matrix_metabo_05_t[
    valid_rows,
    valid_cols,
    drop = FALSE
  ]


#--------------------------------------------------------
# 7.3 Scale Data Prior to PCA
#--------------------------------------------------------

# Standardize each metabolite to mean = 0 and SD = 1
# before PCA so that metabolites with larger variance do
# not disproportionately influence the analysis.
pca_Lahu <- scale(log_matrix_metabo_05_t_Lahu_pca)
head(pca_Lahu)

# Replace NaN values generated by features with zero
# standard deviation with 0 after scaling.
pca_Lahu[is.nan(pca_Lahu)] <- 0


#--------------------------------------------------------
# 7.4 Perform PCA
#--------------------------------------------------------

# Transpose the matrix so that samples are represented
# as rows and metabolites as columns, as required by
# FactoMineR::PCA().
#
# scale.unit = FALSE is used because the metabolite
# variables were standardized before PCA.
res_pca <- PCA(
  t(pca_Lahu),
  scale.unit = FALSE,  # already scaled
  graph = FALSE
)


#--------------------------------------------------------
# 7.5 Prepare Metadata for Visualization
#--------------------------------------------------------

# Extract sample identifiers corresponding to the PCA
# input matrix.
sample_names <- colnames(log_matrix_metabo_05_t_Lahu_pca)

# Retrieve the corresponding geographic group labels.
group_labels <- Met05_Lahu$Group

# Create a metadata table for PCA visualization.
metadata <- data.frame(
  Sample = sample_names,
  Group = group_labels
)

# Use sample identifiers as row names for the metadata.
rownames(metadata) <- sample_names


#--------------------------------------------------------
# 7.6 PCA Plot with Group Ellipses
#--------------------------------------------------------

# Define the color palette for the two Lahu geographic
# groups.
group_colors <- c(
  "Lahu-CM" = "lightblue2",
  "Lahu-CR" = "coral1"
)

# Visualize individual samples in PCA space and add
# 95% confidence ellipses for each geographic group.
fviz_pca_ind(
  res_pca,
  geom.ind = "point",
  col.ind = metadata$Group,     # Color samples by geographic group
  palette = group_colors,
  addEllipses = TRUE,           # Add 95% confidence ellipses
  ellipse.type = "confidence",
  ellipse.level = 0.95,
  repel = TRUE                  # Reduce overlap among labels
) +
  labs(title = "PCA of Lahu Metabolites (|Log2FC| >= 1)") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    title = element_text(size = 8)
  )






#========================================================
# 8. Normality and Homogeneity of Variance Testing
#========================================================

# Assess distributional normality and homogeneity of
# variance for Lahu metabolites prior to group comparison.
#
# These tests are used to characterize the distributional
# properties of the metabolite data and to support the
# selection of an appropriate statistical test.
setwd("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat")


#--------------------------------------------------------
# 8.1 Subset Lahu Metabolite Matrix
#--------------------------------------------------------

# Identify metabolite columns that are present in both the
# log2-transformed metabolite matrix and the Lahu log2FC
# summary table.
valid_cols <- intersect(
  colnames(log_matrix_metabo_05),
  log2FC_sum_Met05_Lahu$Metabolite
)

# Identify Lahu samples that are present in both the
# metabolite matrix and the Lahu metadata.
valid_rows <- intersect(
  rownames(log_matrix_metabo_05),
  Met05_Lahu$Sample_ID
)

# Subset the matrix to Lahu samples and the corresponding
# metabolite features.
# The resulting matrix contains samples as rows and
# metabolites as columns.
log_matrix_metabo_05_Lahu <-
  log_matrix_metabo_05[
    valid_rows,
    valid_cols,
    drop = FALSE
  ]


#--------------------------------------------------------
# 8.2 Identify and Remove Zero-Variance Metabolites
#--------------------------------------------------------

# Count metabolites with zero standard deviation across
# the Lahu samples.
sum(apply(Met05_Lahu[-c(1:12)], 2, sd, na.rm = TRUE) == 0)

# List metabolites with zero variance.
colnames(Met05_Lahu[-c(1:12)])[
  apply(Met05_Lahu[-c(1:12)], 2, sd, na.rm = TRUE) == 0
]

# Create a dataset containing only metabolites with
# non-zero variance for downstream analyses where
# variance is required.
log_matrix_metabo_05_Lahu_nonzero <-
  Met05_Lahu[-c(1:12)][,
                       apply(Met05_Lahu[-c(1:12)], 2, sd, na.rm = TRUE) != 0,
                       drop = FALSE
  ]


#--------------------------------------------------------
# 8.3 Normality and Variance Homogeneity Tests
#--------------------------------------------------------

# Assess normality and homogeneity of variance for each
# metabolite using the predefined norm_var_check() function.
#
# The resulting table contains the corresponding test
# statistics and p-values used to characterize the
# distributional properties of each metabolite.
Lahu_norm_var <- norm_var_check(
  Met05_Lahu,
  metabo_05_descp,
  tag = "metabo05_Lahu"
)

# Inspect the normality and variance-homogeneity results.
Lahu_norm_var

# Count metabolites for which the variance-homogeneity
# test did not return a p-value.
Lahu_norm_var %>%
  count(is.na(Var_p.value))

# Export the normality and variance-homogeneity test results.
Export(Lahu_norm_var, "Lahu_norm_var.txt")







#========================================================
# 9. Statistical Analysis — Group Comparison
#========================================================

# Perform statistical comparisons between the two Lahu
# geographic groups (Lahu-CM and Lahu-CR).
#
# The primary group-comparison workflow uses the Wilcoxon
# test implemented through run_wilcoxon_cross(). A separate
# Student's t-test is shown below as a test for a specific
# metabolite.


#--------------------------------------------------------
# 9.0 Example Student's t-Test — Met125
#--------------------------------------------------------

# Perform an independent two-sample Student's t-test for
# Met125 between Lahu-CM and Lahu-CR.
#
# var.equal = TRUE specifies the equal-variance version
# of the Student's t-test.
t.test(Met05_Lahu$Met125 ~ Met05_Lahu$Group,
       var.equal = TRUE)



#--------------------------------------------------------
# 9.1 Wilcoxon Test (Lahu-CM vs Lahu-CR)
#--------------------------------------------------------

# Perform metabolite-wise group comparisons between
# Lahu-CM and Lahu-CR using the predefined
# run_wilcoxon_cross() function.
#
# The function uses the normality and variance-assessment
# results generated in Section 8 and returns the
# corresponding statistical results and adjusted p-values.
Lahu_wilcox <- run_wilcoxon_cross(
  metabo_data = Met05_Lahu, 
  metabo_norm_var = Lahu_norm_var, 
  group_label = "Group",
  comparison_label = "Lahu-CM vs Lahu-CR",
  tag = "metabo_Lahu"
)

# Inspect the first results.
head(Lahu_wilcox)


#--------------------------------------------------------
# 9.2 Filter Significant Metabolites
#--------------------------------------------------------

# Identify metabolites showing statistically significant
# differences after false discovery rate (FDR) adjustment.
#
# Metabolites with p_adj < 0.05 are retained for the
# subsequent effect-size filtering step.
sig_Lahu_05 <- Lahu_wilcox %>% 
  filter(p_adj < 0.05) %>%
  pull(Metabolite)

# Retain statistically significant metabolites with an
# absolute log2FC of at least 1, corresponding to at least
# a 2-fold difference between groups.
#
# Metabolite classification information is then added
# using the curated metabolite annotation table.
log2FC_sum_Met05_Lahu_sig <- log2FC_sum_Met05_Lahu %>%
  filter(Metabolite %in% sig_Lahu_05) %>%
  filter(abs(log2FC) >= 1) %>%
  left_join(Pub_Classy_05_clean,
            by = c("Compound",
                   "Compound_Clean",
                   "Metabolite" = "Code"))

# Export the final list of significant metabolites meeting
# both the FDR and effect-size criteria.
Export(log2FC_sum_Met05_Lahu_sig, "Lahu_sig_log2FC_1.txt")


#========================================================
# 10. Visualization — Significant Metabolites
#========================================================

# Generate a barplot of metabolites retained after
# statistical significance and effect-size filtering.
#
# The input table contains metabolites with:
#   - FDR-adjusted p-value < 0.05
#   - absolute log2FC >= X
ggFC_sig_Met05_Lahu <- plot_log2FC_bar(
  summary_df = log2FC_sum_Met05_Lahu_sig,
  title = "Log2FC of metabolites (Lahu-CM vs. Lahu-CR)",
  xlab = "Significant Metabolite\n(|Log2FC| >= 2, q < 0.05)",
  fill_colors = c("Upregulated" = "steelblue3",
                  "Downregulated" = "#FC766A"),
  group_color = "Lahu-CM", 
  strip_text_color = "black"
)

# Display the significant-metabolite barplot.
print(ggFC_sig_Met05_Lahu)







#=========================================================
# 11. Visualization — Major Subclass Annotation
#=========================================================

# Load the Lahu-specific metabolite analysis workspace
# containing objects required for annotation and downstream
# visualization.
load("~/Documents/HillTribe_NGS/5.Metabolites/Metabolites_Lahu.RData")


#---------------------------------------------------------
# 11.1 Assign Major Metabolite Subclasses
#---------------------------------------------------------

# Group individual metabolite subclasses into broader
# chemical categories to facilitate interpretation of the
# significant metabolites.
#
# The classification is based on the original metabolite
# Subclass annotation. Subclasses that are not explicitly
# assigned to a major category are labeled as "Unknown".
log2FC_sum_Met05_Lahu_sig_anno_majorsubclass <- 
  log2FC_sum_Met05_Lahu_sig %>%
  mutate(Major_Subclass = case_when(
    
    # Lipids and fatty-acid related compounds
    Subclass %in% c(
      "Fatty acid esters",
      "Carboxylic acid derivatives"
    ) ~ "Lipids & FA derivatives",
    
    # Terpenoids
    Subclass %in% c(
      "Diterpenoids",
      "Sesquiterpenoids"
    ) ~ "Terpenoids",
    
    # Phenolic and aromatic compounds
    Subclass %in% c(
      "Benzoic acids and derivatives"
    ) ~ "Phenolic & aromatic",
    
    # Nitrogen-containing compounds
    Subclass %in% c(
      "Amines"
    ) ~ "Nitrogen-containing",
    
    # Oxygenated molecules
    Subclass %in% c(
      "Gamma butyrolactones",
      "Carbonyl compounds"
    ) ~ "Oxygenated molecules",
    
    # Hydrocarbons
    Subclass %in% c(
      "Alkanes",
      "Olefins"
    ) ~ "Hydrocarbons",
    
    # Silicon-containing compounds
    # These compounds may represent analytical artifacts or
    # contamination commonly associated with GC-MS analysis.
    Subclass %in% c(
      "Organosilicon compounds"
    ) ~ "Silicon-containing",
    
    # Retain subclasses not included above as "Unknown".
    TRUE ~ "Unknown"
    
  ))


#---------------------------------------------------------
# 11.2 Identify Top Enriched Metabolites
#---------------------------------------------------------

# Identify strongly upregulated metabolites based on the
# predefined log2FC threshold and rank them according to
# their effect size.
topenrich_Lahu <- log2FC_sum_Met05_Lahu_sig_anno_majorsubclass %>% 
  filter(Direction == "Upregulated" & log2FC >= 10) %>% 
  arrange(desc(log2FC)) %>% 
  select(Group, Metabolite, Compound_Clean, 
         log2FC, Direction, Subclass, Major_Subclass) %>%
  rename("Code" = "Metabolite",
         "Compound" = "Compound_Clean")

# Export the list of highly enriched metabolites.
Export(topenrich_Lahu, "topenrich_Lahu.txt")


#---------------------------------------------------------
# 11.3 Inspect Major Subclass and Annotation Coverage
#---------------------------------------------------------

# List the major chemical categories represented among
# significant Lahu metabolites.
log2FC_sum_Met05_Lahu_sig_anno_majorsubclass %>%
  distinct(Major_Subclass)

# Count significant metabolites for which a CID could not
# be identified in the annotation database.
log2FC_sum_Met05_Lahu_sig_anno_majorsubclass %>% 
  filter(CID == "Not found") %>%
  count()


#---------------------------------------------------------
# 11.4 Visualize Significant Metabolites by Major Subclass
#---------------------------------------------------------

# Generate a barplot for significant metabolites belonging
# to the Oxygenated molecules category.
ggFC_sig_Met05_Lahu_subclass <- plot_log2FC_bar_subclass(
  
  summary_df = log2FC_sum_Met05_Lahu_sig_anno_majorsubclass %>% 
    filter(Major_Subclass == "Oxygenated molecules"),
  
  title = "Log2FC of metabolites (Lahu-CM vs. Lahu-CR)\nOxygenated molecules",
  
  xlab = "Significant Metabolite\n(|Log2FC| >= 1, q < 0.05)",
  
  fill_colors = c("Upregulated" = "steelblue3",
                  "Downregulated" = "#FC766A"),
  
  log2FC_cutoff = 1,
  
  strip_text_color = "black",
  
  group_color =  group_colors
)

# Display the subclass-specific metabolite barplot.
print(ggFC_sig_Met05_Lahu_subclass)








#=========================================================
# 12. Statistical Summary and Feature-Level QC
#=========================================================

# Summarize the number of significant metabolites according
# to their direction of change.
log2FC_sum_Met05_Lahu_sig_anno_majorsubclass %>%
  count(Direction)

# Summarize significant metabolites by direction of change
# and major chemical subclass, and rank subclasses according
# to the number of significant metabolites.
log2FC_sum_Met05_Lahu_sig_anno_majorsubclass %>% 
  group_by(Direction, Major_Subclass) %>%
  count() %>%
  arrange(desc(n)) %>%
  print(n = 50)



#---------------------------------------------------------
# 12.1 Identify Zero-Variance Metabolites
#---------------------------------------------------------

# Load packages required for the zero-variance QC procedure.
library(dplyr)
library(tidyr)

# Define a helper function to identify metabolites with
# zero standard deviation within a specified dataset.
get_zero_var <- function(df) {
  
  df %>%
    summarise(across(starts_with("Met"), ~ sd(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Metabolite", values_to = "SD") %>%
    filter(SD == 0) %>%
    pull(Metabolite)
}


# Identify metabolites with zero variance across all
# Lahu samples, irrespective of geographic group.
Met05_Lahu_zero <- Met05_Lahu %>%
  filter(Ethnicity == "Lahu") %>%
  get_zero_var()

# Identify metabolites with zero variance specifically
# within the Lahu-CM group.
Met05_LahuCM_zero <- Met05_Lahu %>%
  filter(Group == "Lahu-CM") %>%
  get_zero_var()

# Identify metabolites with zero variance specifically
# within the Lahu-CR group.
Met05_LahuCR_zero <- Met05_Lahu %>%
  filter(Group == "Lahu-CR") %>%
  get_zero_var()


#---------------------------------------------------------
# 12.2 Identify Group-Specific Zero-Variance Metabolites
#---------------------------------------------------------

# Identify metabolites with zero variance in Lahu-CM but
# non-zero variance across the complete Lahu dataset.
zero_LahuCM <- setdiff(Met05_LahuCM_zero, Met05_Lahu_zero)

# Identify metabolites with zero variance in Lahu-CR but
# non-zero variance across the complete Lahu dataset.
zero_LahuCR <- setdiff(Met05_LahuCR_zero, Met05_Lahu_zero)


#---------------------------------------------------------
# 12.3 Assess Zero-Variance Metabolites Among Significant
#       Features
#---------------------------------------------------------

# Combine metabolites with group-specific zero variance
# and identify those that are also statistically significant.
diff_CM_CR <- union(zero_LahuCM, zero_LahuCR)
sig_zero_var <- intersect(diff_CM_CR, sig_Lahu_05)

# Create a table of significant metabolites that have
# zero variance in either Lahu-CM or Lahu-CR.
sigzero_either_Lahu <- bind_rows(
  tibble(Metabolite = zero_LahuCM, Group = "Lahu-CM"),
  tibble(Metabolite = zero_LahuCR, Group = "Lahu-CR")
) %>%
  filter(Metabolite %in% sig_zero_var) %>%
  rename(zero_var = Group) %>%
  left_join(log2FC_sum_Met05_Lahu_sig_anno_majorsubclass,
            by = "Metabolite")

# Summarize significant zero-variance metabolites by
# geographic group and original metabolite subclass.
sigzero_either_Lahu %>% 
  group_by(Group, Subclass) %>%
  count() %>%
  arrange(desc(n)) %>%
  print(n = 50)

# Summarize the number of significant zero-variance
# metabolites according to the group in which zero
# variance was observed.
sigzero_either_Lahu %>% 
  group_by(zero_var) %>%
  count()
