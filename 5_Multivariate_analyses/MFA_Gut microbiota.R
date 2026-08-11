############################################################
# Multiple Factor Analysis (MFA)
############################################################
# Purpose:
# Perform Multiple Factor Analysis (MFA) to jointly analyze metadata,
# anthropometric variables, and significant gut microbiota taxa.
#
# Taxonomic variables are organized into separate MFA groups according to
# taxonomic rank:
#   - Phylum
#   - Class
#   - Order
#   - Family
#   - Genus
#
# The analysis is performed separately for the Akha and Lahu ethnic groups.


# ==============================
# 1. Environment Setup
# ==============================

# Set the working directory and load the processed taxonomic dataset.
setwd("~/Documents/HillTribe_NGS/6.MFA")
load("~/Documents/HillTribe_NGS/3.Compositional_profile/Heatmap_plot.RData")


# Load required packages -------------------------------------------------------
# Packages used for data manipulation, statistical analysis, MFA, and
# visualization.
library(dplyr)
library(readxl)
library(ggplot2)
library(rio)
library(car)  
library(FactoMineR)
library(factoextra)
library(ggplot2)



# ==========================================================
# 2. MFA – Akha Ethnic Group (Significant Taxa)
# ==========================================================
# Perform MFA for the Akha ethnic group using the significant taxonomic
# features identified in the preceding taxonomic analysis.
#
# The MFA integrates:
#   1. Metadata
#   2. Anthropometric variables
#   3. Phylum-level taxa
#   4. Class-level taxa
#   5. Order-level taxa
#   6. Family-level taxa
#   7. Genus-level taxa
#
# The metadata group is treated as supplementary and therefore does not
# contribute to the construction of the principal dimensions.


# Prepare MFA dataset ----------------------------------------------------------
# Remove sample identifiers and variables not included in the MFA.
# Categorical variables are converted to factors.
# Age is relocated to maintain the intended variable order within the
# anthropometric group.
Akha_MFA <- sig_taxo_Akha %>%
  select(-Sample_ID, -Ethnicity, -BMl) %>%
  mutate(across(c(Gender, Group, Area, 
                  BMI_group, Hypertension), factor)) %>% 
  relocate(Age, .before = 7)


# Identify taxonomic columns by rank prefix -----------------------------------
# Identify columns belonging to each taxonomic rank using the corresponding
# taxonomic prefixes.
p_cols <- which(startsWith(colnames(Akha_MFA), "p__"))
c_cols <- which(startsWith(colnames(Akha_MFA), "c__"))
o_cols <- which(startsWith(colnames(Akha_MFA), "o__"))
f_cols <- which(startsWith(colnames(Akha_MFA), "f__"))
g_cols <- which(startsWith(colnames(Akha_MFA), "g__"))


# Check column ranges and counts per taxonomic level --------------------------
# Inspect the positions and number of variables assigned to each taxonomic
# rank. These checks are used to verify the taxonomic group structure before
# running the MFA.
range(p_cols)
range(c_cols)
range(o_cols)
range(f_cols)
range(g_cols)

length(p_cols)
length(c_cols)
length(o_cols)
length(f_cols)
length(g_cols)


# Run MFA ----------------------------------------------------------------------
# Integrate metadata, anthropometric variables, and taxonomic variables
# as separate MFA groups.
#
#   Group 1: Metadata
#     5 categorical variables
#
#   Group 2: Anthropometric
#     5 quantitative variables
#
#   Group 3: Phylum
#     5 quantitative variables
#
#   Group 4: Class
#     6 quantitative variables
#
#   Group 5: Order
#     9 quantitative variables
#
#   Group 6: Family
#     17 quantitative variables
#
#   Group 7: Genus
#     27 quantitative variables
#
# The metadata group is specified as supplementary using num.group.sup = 1.
res.mfa <- MFA(Akha_MFA,
               
               group = c(5,
                         5,
                         5,6,9,17,27),
               
               type = c("n",
                        "s",
                        "s","s","s","s","s"),
               
               name.group = c("Metadata",
                              "Anthropometric",
                              "Phylum", "Class", "Order", "Family", "Genus"),
               
               num.group.sup = 1,
               
               graph = FALSE
)




# ==========================================================
# 3. MFA – Lahu Ethnic Group
# ==========================================================
# Perform MFA for the Lahu ethnic group using the significant taxonomic
# features identified in the preceding taxonomic analysis.
#
# The MFA structure follows the same framework as the Akha analysis, but
# the number of significant taxa available at each taxonomic rank differs.


# Prepare MFA dataset ----------------------------------------------------------
# Remove sample identifiers and variables not included in the MFA.
# Categorical variables are converted to factors.
# Age is relocated to maintain the intended variable order within the
# anthropometric group.
Lahu_MFA <- sig_taxo_Lahu %>%
  select(-Sample_ID, -Ethnicity, -BMl) %>%
  mutate(across(c(Gender, Group, Area, 
                  BMI_group, Hypertension), factor)) %>%
  relocate(Age, .before = 7)


# Identify taxonomic ranks ----------------------------------------------------
# Identify columns belonging to each taxonomic rank using the corresponding
# taxonomic prefixes.
p_cols <- which(startsWith(colnames(Lahu_MFA), "p__"))
c_cols <- which(startsWith(colnames(Lahu_MFA), "c__"))
o_cols <- which(startsWith(colnames(Lahu_MFA), "o__"))
f_cols <- which(startsWith(colnames(Lahu_MFA), "f__"))
g_cols <- which(startsWith(colnames(Lahu_MFA), "g__"))


# Check the number of variables at each taxonomic level.
length(p_cols)
length(c_cols)
length(o_cols)
length(f_cols)
length(g_cols)


# Run MFA ----------------------------------------------------------------------
# Integrate metadata, anthropometric variables, and the available taxonomic
# ranks as separate MFA groups.
#
#   Group 1: Metadata
#     5 categorical variables
#
#   Group 2: Anthropometric
#     5 quantitative variables
#
#   Group 3: Phylum
#     3 quantitative variables
#
#   Group 4: Class
#     3 quantitative variables
#
#   Group 5: Order
#     2 quantitative variables
#
#   Group 6: Family
#     4 quantitative variables
#
# Genus-level variables are not included in the MFA group specification for
# the Lahu dataset.
#
# The metadata group is specified as supplementary using num.group.sup = 1.
res.mfa <- MFA(Lahu_MFA,
               
               group = c(5,
                         5,
                         3,3,2,4),
               
               type = c("n",
                        "s",
                        "s","s","s","s"),
               
               name.group = c("Metadata",
                              "Anthropometric",
                              "Phylum", "Class", "Order", "Family"),
               
               num.group.sup = 1,
               
               graph = FALSE
)






# ==========================================================
# 4. MFA – CM Area
# ==========================================================
# Perform MFA for samples from Chiang Mai (CM) using the significant
# taxonomic features identified in the preceding taxonomic analysis.
#
# The MFA integrates:
#   1. Metadata
#   2. Anthropometric variables
#   3. Phylum-level taxa
#   4. Class-level taxa
#   5. Order-level taxa
#   6. Family-level taxa
#   7. Genus-level taxa
#
# The metadata group is treated as supplementary and therefore does not
# contribute to the construction of the principal dimensions.


# Prepare MFA dataset ----------------------------------------------------------
# Remove sample identifiers and variables not included in the MFA.
# Categorical variables are converted to factors.
# Age is relocated to maintain the intended variable order within the
# anthropometric group.
CM_MFA <- sig_taxo_CM %>%
  select(-Sample_ID, -Area, -BMl) %>%
  mutate(across(c(Gender, Group, Ethnicity, 
                  BMI_group, Hypertension), factor)) %>%
  relocate(Age, .before = 7)


# Identify taxonomic ranks ----------------------------------------------------
# Identify columns belonging to each taxonomic rank using the corresponding
# taxonomic prefixes.
p_cols <- which(startsWith(colnames(CM_MFA), "p__"))
c_cols <- which(startsWith(colnames(CM_MFA), "c__"))
o_cols <- which(startsWith(colnames(CM_MFA), "o__"))
f_cols <- which(startsWith(colnames(CM_MFA), "f__"))
g_cols <- which(startsWith(colnames(CM_MFA), "g__"))


# Check the number of variables at each taxonomic level.
length(p_cols)
length(c_cols)
length(o_cols)
length(f_cols)
length(g_cols)


# Run MFA ----------------------------------------------------------------------
# Integrate metadata, anthropometric variables, and taxonomic variables
# as separate MFA groups.
#
#   Group 1: Metadata
#     5 categorical variables
#
#   Group 2: Anthropometric
#     5 quantitative variables
#
#   Group 3: Phylum
#     3 quantitative variables
#
#   Group 4: Class
#     4 quantitative variables
#
#   Group 5: Order
#     6 quantitative variables
#
#   Group 6: Family
#     12 quantitative variables
#
#   Group 7: Genus
#     29 quantitative variables
#
# The metadata group is specified as supplementary using num.group.sup = 1.
res.mfa <- MFA(CM_MFA,
               
               group = c(5,
                         5,
                         3,4,6,12,29),
               
               type = c("n",
                        "s",
                        "s","s","s","s","s"),
               
               name.group = c("Metadata",
                              "Anthropometric",
                              "Phylum", "Class", "Order", "Family", "Genus"),
               
               num.group.sup = 1,
               
               graph = FALSE
)


# ==========================================================
# 5. MFA – CR Area
# ==========================================================
# Perform MFA for samples from Chiang Rai (CR) using the significant
# taxonomic features identified in the preceding taxonomic analysis.
#
# The same MFA framework is applied to the CR dataset, with the number of
# significant taxa at each taxonomic rank determined separately.


# Prepare MFA dataset ----------------------------------------------------------
# Remove sample identifiers and variables not included in the MFA.
# Categorical variables are converted to factors.
# Age is relocated to maintain the intended variable order within the
# anthropometric group.
CR_MFA <- sig_taxo_CR %>%
  select(-Sample_ID, -Area, -BMl) %>%
  mutate(across(c(Gender, Group, Ethnicity, 
                  BMI_group, Hypertension), factor)) %>%
  relocate(Age, .before = 7)


# Identify taxonomic ranks ----------------------------------------------------
# Identify columns belonging to each taxonomic rank using the corresponding
# taxonomic prefixes.
p_cols <- which(startsWith(colnames(CR_MFA), "p__"))
c_cols <- which(startsWith(colnames(CR_MFA), "c__"))
o_cols <- which(startsWith(colnames(CR_MFA), "o__"))
f_cols <- which(startsWith(colnames(CR_MFA), "f__"))
g_cols <- which(startsWith(colnames(CR_MFA), "g__"))


# Check the number of variables at each taxonomic level.
length(p_cols)
length(c_cols)
length(o_cols)
length(f_cols)
length(g_cols)


# Run MFA ----------------------------------------------------------------------
# Integrate metadata, anthropometric variables, and taxonomic variables
# as separate MFA groups.
#
#   Group 1: Metadata
#     5 categorical variables
#
#   Group 2: Anthropometric
#     5 quantitative variables
#
#   Group 3: Phylum
#     2 quantitative variables
#
#   Group 4: Class
#     3 quantitative variables
#
#   Group 5: Order
#     6 quantitative variables
#
#   Group 6: Family
#     12 quantitative variables
#
#   Group 7: Genus
#     39 quantitative variables
#
# The metadata group is specified as supplementary using num.group.sup = 1.
res.mfa <- MFA(CR_MFA,
               
               group = c(5,
                         5,
                         2,3,6,12,39),
               
               type = c("n",
                        "s",
                        "s","s","s","s","s"),
               
               name.group = c("Metadata",
                              "Anthropometric",
                              "Phylum", "Class", "Order", "Family", "Genus"),
               
               num.group.sup = 1,
               
               graph = FALSE
)





# ==========================================================
# 6. Eigenvalues and Scree Plot
# ==========================================================
# Examine the MFA results and identify the amount of variance explained
# by each dimension.
#
# The eigenvalues are used to assess the importance of the MFA dimensions,
# while the scree plot provides a visual summary of the variance explained
# across dimensions.

# Print the complete MFA object and its main results.
print(res.mfa) # The output of the MFA function

# Extract eigenvalues and the percentage of variance explained by each
# MFA dimension.
eig.val <- get_eigenvalue(res.mfa)
head(eig.val)

# Visualize the eigenvalues using a scree plot.
fviz_screeplot(res.mfa)


# ==========================================================
# 7. Group-Level Results
# ==========================================================
# Extract MFA results at the group level.
#
# Group coordinates indicate the position of each variable group on the
# MFA dimensions.
#
# cos2 indicates the quality of representation of each group on a given
# dimension. Higher values indicate that the dimension represents the
# group more strongly.
#
# Contribution indicates the relative contribution of each group to the
# construction of a given MFA dimension.

group <- get_mfa_var(res.mfa, "group")

head(group$coord)     # Group coordinates
head(group$cos2)      # Quality of representation
head(group$contrib)   # Contribution to dimensions


# Plot group factor maps -------------------------------------------------------
# Visualize the positions of the variable groups on Dimensions 1 and 2.
# The coordinates indicate the association of each group with the two
# displayed dimensions.
fviz_mfa_var(res.mfa, "group", repel = TRUE, 
             axes = c(1,2)) + 
  theme_bw()+
  theme(legend.position = "right", legend.text = element_text(size = 12, colour = "black"),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14))

# Visualize the positions of the variable groups on Dimensions 3 and 4.
fviz_mfa_var(res.mfa, "group", repel = TRUE, 
             axes = c(3,4)) + 
  theme_bw()+
  theme(legend.position = "right", legend.text = element_text(size = 12, colour = "black"),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14))


# Group contribution plots -----------------------------------------------------
# Visualize the contribution of each variable group to Dimension 1.
# Groups with contributions above the average contribution have a greater
# influence on the construction of the dimension.
fviz_contrib(res.mfa, "group", axes = 1) +
  theme_bw() +
  theme(legend.position = "right", legend.text = element_text(size = 12, colour = "black"),
        axis.text = element_text(color = "black", size = 12),
        axis.title.x = element_blank())

# Visualize the contribution of each variable group to Dimension 2.
fviz_contrib(res.mfa, "group", axes = 2) +
  theme_bw() +
  theme(legend.position = "right", legend.text = element_text(size = 12, colour = "black"),
        axis.text = element_text(color = "black", size = 12),
        axis.title.x = element_blank())




# ==========================================================
# 8. Individual Factor Map
# ==========================================================
# Extract individual-level MFA results for examining the position of
# individual samples in the MFA factor space.
ind <- get_mfa_ind(res.mfa)
ind


# Define color palettes according to grouping variable ------------------------
# The following palettes correspond to different MFA analyses and grouping
# structures used throughout the project. The active palette is determined
# by the most recent assignment to group_colors.

group_colors <- c(
  "Akha-CM" = "gold1",
  "Akha-CR" = "dodgerblue4",
  "NHT" = "grey30",
  "HT" = "firebrick3"
)


group_colors <- c(
  "Lahu-CM" = "lightblue2",
  "Lahu-CR" = "coral1",
  "NHT" = "grey30",
  "HT" = "firebrick3"
)

group_colors <- c(
  "NHT" = "grey30",
  "HT" = "firebrick3"
)

group_colors <- c(
  "Akha-CM" = "gold1",
  "Lahu-CM" = "lightblue2",
  "Khuen-CM" = "darkorange",
  "NHT" = "grey30",
  "HT" = "firebrick3"
)


group_colors <- c(
  "Akha-CR" = "dodgerblue4",
  "Lahu-CR" = "coral1",
  "Lisu-CR" = "grey40",
  "NHT" = "grey30",
  "HT" = "firebrick3"
)

group_colors <- c(
  "UW" = "grey70",
  "N" = "grey20",
  "OW" = "gold1",
  "OB" = "dodgerblue4"
)


# Individual factor map (colored by group factor) -----------------------------
# Visualize individual samples on Dimensions 1 and 2.
#
# Samples are colored according to the selected grouping variables, with
# confidence ellipses used to visualize the distribution of observations
# within groups.
#
# The supplementary qualitative variable is excluded from the plot, and
# individual labels are omitted to improve readability.
fviz_mfa_ind(res.mfa, axes = c(1,2), 
             labelsize = 2, 
             invisible = "quali.var", #omit  the supplementary qualitative variable
             habillage = c("Group"), # colored by groups or c("Group", "Hypertension").
             palette = group_colors, 
             label = "none",
             addEllipses = TRUE, 
             ellipse.type = "confidence", 
             repel = TRUE) + # Avoid text overlapping)
  theme(aspect.ratio = 1,
        axis.title.y = element_text(face = "plain", size = rel(1),
                                    margin = margin(t = 0, r = 15, b = 0, l = 0)
                                    , hjust = 0.5), #label of y
        axis.title.x = element_text(face = "plain", size = rel(1),
                                    margin = margin(t = 15, r = 0, b = 0, l = 0)
                                    , vjust = 0.5),
        axis.text.y = element_text(face = "plain", size = rel(1), colour = "black"), #scale y label
        axis.text.x = element_text(face = "plain", size = rel(1), colour = "black",
                                   angle = 0,
                                   margin = margin(t = 0, r = 0, b = 0, l = 0),
                                   vjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.ticks.length.x = unit(0.2, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.spacing.y = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "cm"),
        legend.position = "right",
        legend.text = element_text(face = "plain", size = rel(1))
  )


# ==========================================================
# 8.1 Quantitative Variable Representation (cos²)
# ==========================================================
# Calculate the combined cos² of quantitative variables across Dimensions 1
# and 2. The resulting value represents the proportion of each variable's
# representation captured by the two-dimensional factor map.
#
# Higher combined cos² values indicate that a variable is better represented
# in the Dim1-Dim2 factor space.
cos2_vals <- res.mfa$quanti.var$cos2[, 1] + res.mfa$quanti.var$cos2[, 2]

# Select the 30 quantitative variables with the highest combined cos² values.
# These variables are the best represented in the Dim1-Dim2 factor map.
top_vars <- names(sort(cos2_vals, decreasing = TRUE))[1:30]

# Plot the selected quantitative variables on Dimensions 1 and 2.
# Color intensity represents the cos² value, allowing variables with higher
# quality of representation to be readily identified.
fviz_mfa_var(res.mfa, "quanti.var",
             axes = c(1,2),
             select.var = list(name = top_vars),  # filter by variable names
             col.var = "cos2",          
             gradient.cols = c("gray54", "#E7B800", "#FC4E07"), 
             repel = TRUE,
             labelsize = 3,
             geom = c("point", "text"),
             legend = "right") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10, colour = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 10)
  )


# ==========================================================
# 8.2 Overall cos² Visualization
# ==========================================================
# Visualize the overall quality of representation of quantitative variables
# across the MFA factor space using cos² values.
#
# Variables with higher cos² values are better represented by the displayed
# factor map. For a two-dimensional map, the sum of cos² values for
# Dimensions 1 and 2 approaches 1 when the variable is very well represented
# by these two dimensions.
fviz_mfa_var(res.mfa, col.var = "cos2",
             gradient.cols = c("gray54", "#E7B800", "#FC4E07"), 
             col.var.sup = "violet", repel = TRUE,
             geom = c("point", "text"),
             labelsize = 4,
             top = 10) +
  theme_bw()+
  theme(legend.position = "right", legend.text = element_text(size = 12, colour = "black"),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14))




# ==========================================================
# 10. Correlation with Dimensions
# ==========================================================
# Examine the relationships between quantitative variables and the MFA
# dimensions.
#
# The MFA variable coordinates and correlations can be used to assess how
# quantitative variables are positioned and associated with each dimension.
# Dimensions 1-2 and Dimensions 3-4 are visualized separately.


# Correlation between quantitative variables and Dimensions 1 and 2 -----------

# Visualize quantitative variables on the Dim1-Dim2 factor map.
# The position of each variable reflects its coordinates on the two dimensions.
fviz_mfa_var(res.mfa, "quanti.var", palette = "jco", axes = c(1,2), 
             col.var.sup = "violet", repel = TRUE,
             geom = c("point"), legend = "bottom") + 
  theme_bw()+
  theme(legend.position = "right", legend.text = element_text(size = 12, colour = "black"),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14))

# Extract variable coordinates on Dimensions 1 and 2.
res.mfa$quanti.var$coord[, 1:2]

# Extract correlations between quantitative variables and Dimensions 1 and 2.
# Correlations indicate the strength and direction of the association between
# each quantitative variable and the corresponding MFA dimension.
res.mfa$quanti.var$cor[, 1:2]


# Correlation between quantitative variables and Dimensions 3 and 4 -----------

# Visualize quantitative variables on the Dim3-Dim4 factor map.
fviz_mfa_var(res.mfa, "quanti.var", palette = "jco", axes = c(3,4), 
             col.var.sup = "violet", repel = TRUE,
             geom = c("point"), legend = "bottom") + 
  theme_bw()+
  theme(legend.position = "right", legend.text = element_text(size = 12, colour = "black"),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14))


# ==========================================================
# 11. Variable Contributions
# ==========================================================
# Identify the quantitative variables that contribute most strongly to the
# construction of each MFA dimension.
#
# Contribution measures the relative influence of each variable on a
# dimension. Variables with higher contributions have a greater influence
# on the formation of that dimension.
#
# The plots display the highest-contributing variables for each dimension.


# Contributions to Dimension 1 ------------------------------------------------
fviz_contrib(res.mfa, choice = "quanti.var", axes = 1, top = 20,
             palette = "jco") +
  theme_bw()+
  theme(legend.position = "right", legend.text = element_text(size = 9, colour = "black"),
        axis.text = element_text(color = "black", size = 8, angle = 45, hjust = 1),
        axis.title.x = element_blank())


# Contributions to Dimension 2 ------------------------------------------------
fviz_contrib(res.mfa, choice = "quanti.var", axes = 2, top = 20,
             palette = "jco") +
  theme_bw()+
  theme(legend.position = "right", legend.text = element_text(size = 9, colour = "black"),
        axis.text = element_text(color = "black", size = 8, angle = 45, hjust = 1),
        axis.title.x = element_blank())


# Contributions to Dimension 3 ------------------------------------------------
fviz_contrib(res.mfa, choice = "quanti.var", axes = 3, top = 50,
             palette = "jco") +
  theme_bw()+
  theme(legend.position = "right", legend.text = element_text(size = 9, colour = "black"),
        axis.text = element_text(color = "black", size = 7, angle = 45, hjust = 1),
        axis.title.x = element_blank())


# Contributions to Dimension 4 ------------------------------------------------
fviz_contrib(res.mfa, choice = "quanti.var", axes = 4, top = 50,
             palette = "jco")+  theme_bw()+
  theme(legend.position = "right", legend.text = element_text(size = 9, colour = "black"),
        axis.text = element_text(color = "black", size = 7, angle = 45, hjust = 1),
        axis.title.x = element_blank())


# ==========================================================
# 12. Dimension Description
# ==========================================================
# Identify variables significantly associated with each MFA dimension.
#
# dimdesc() provides statistical descriptions of the dimensions, helping
# identify quantitative variables and qualitative variables that are
# associated with the corresponding factor.
#
# Dimensions are examined individually to facilitate interpretation of the
# biological and clinical characteristics represented by each dimension.

dimdesc(res.mfa, axes = 1)
dimdesc(res.mfa, axes = 2)
dimdesc(res.mfa, axes = 3)
dimdesc(res.mfa, axes = 4)
