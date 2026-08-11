############################################################
# Redundancy Analysis (RDA)
############################################################
# Purpose:
# Evaluate how much variation in metabolite and gut microbiota profiles
# can be explained by predefined sample groups using Redundancy Analysis
# (RDA).
#
# The workflow includes:
#   1. Preparation of phylum- and genus-level microbial count data
#   2. CLR transformation of microbial count data
#   3. Preparation and quality control of scaled metabolite data
#   4. RDA using metabolites alone or combined metabolite + microbiota data
#   5. Permutation-based significance testing
#   6. RDA visualization using scaling 1 and scaling 2
#   7. Identification of features with large RDA vector lengths
#   8. Pairwise comparison of feature directions using cosine similarity
#
# The analysis below demonstrates the workflow for the Akha ethnic group.


# ==========================================================
# 1. Environment Setup
# ==========================================================

# Set the working directory for RDA analysis.
setwd("~/Documents/HillTribe_NGS/7.RDA")


# Load previously processed datasets required for the RDA workflow.
#
# MFA_composition.RData:
#   Contains processed microbiota/compositional data and associated metadata.
#
# Metabolites.RData:
#   Contains processed metabolite data.
#
# PERMANOVA_metabolite.RData:
#   Contains metabolite-level statistical results and objects generated
#   during the PERMANOVA analysis.
load("~/Documents/HillTribe_NGS/6.MFA/MFA_composition.RData")
load("~/Documents/HillTribe_NGS/5.Metabolites/Metabolites.RData")
load("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat/PERMANOVA_metabolite/PERMANOVA_metabolite.RData")


# Import microbial count tables at the phylum and genus levels.
library(readxl)
Counts_Phylum <- read_excel("~/Documents/HillTribe_NGS/3.Compositional_profile/1_Counts_Phylum.xlsx")
Counts_Genus <- read_excel("~/Documents/HillTribe_NGS/3.Compositional_profile/5_Counts_Genus.xlsx")


# Join the phylum- and genus-level count tables.
# Confirm that both tables contain samples in the same order before combining.
identical(Counts_Phylum$Sample_ID, Counts_Genus$Sample_ID)
# [TRUE]


# Combine the phylum and genus count tables while retaining only one copy
# of the shared metadata columns.
count_asv <- bind_cols(
  Counts_Phylum,
  Counts_Genus %>% select(-(1:13))
)


# Confirm that sample order was retained after combining the tables.
all(count_asv$Sample_ID == Counts_Phylum$Sample_ID)


# Create a combined geographic/ethnic grouping variable.
# The resulting Group variable follows the format:
#   Ethnicity-Area
# e.g., Akha-CM and Akha-CR.
#
# Samples listed in drop_sampleID are excluded before RDA.
count_asv <- count_asv %>%
  mutate(Group = paste(Ethnicity, Area, sep = "-")) %>%
  relocate(Group, .before = 2) %>%
  filter(!Sample_ID %in% drop_sampleID) %>%
  as.data.frame()


# Use Sample_ID as the row names so that microbial count matrices can be
# matched directly to the corresponding sample metadata.
rownames(count_asv) <- count_asv$Sample_ID



# ==========================================================
# 2. Load Packages for RDA and Data Processing
# ==========================================================
# dplyr and tibble:
#   Data manipulation and table handling.
#
# compositions:
#   CLR transformation of compositional microbiome data.
#
# vegan:
#   Redundancy Analysis and permutation-based significance testing.
library(dplyr)
library(tibble)
library(compositions)
library(vegan)



# ==========================================================
# 3. RDA – Akha Ethnic Group
# ==========================================================
# The Akha-specific RDA is performed using the significant taxonomic
# features identified in the previous taxonomic analysis.
#
# Microbial features are retained at the phylum and genus levels.
# Metabolite features are taken from the previously processed metabolite
# dataset.


# Extract microbial count data for the Akha ethnic group.
# The selected columns correspond to the significant taxonomic features
# identified in sig_taxo_Akha.
count_Akha <- count_asv %>%
  filter(Ethnicity == "Akha") %>%
  select(any_of(names(sig_taxo_Akha)[13:ncol(sig_taxo_Akha)])) %>%
  as.matrix()


# Confirm that microbial sample order matches the Akha metadata.
# Matching row order is essential because the response matrix and explanatory
# variables must refer to the same samples.
identical(rownames(count_Akha), rownames(metadata_Akha))


# ----------------------------------------------------------
# 3.1 CLR Transformation of Microbial Counts
# ----------------------------------------------------------
# Microbiome count data are compositional and therefore are transformed
# using the centered log-ratio (CLR) transformation.
#
# A pseudocount of 1 is added before CLR transformation to avoid zero values.
Akha_bac_clr <- clr(count_Akha + 1)

head(Akha_bac_clr)
rownames(Akha_bac_clr)


# ----------------------------------------------------------
# 3.2 Metabolite Data Quality Control
# ----------------------------------------------------------
# Check for NaN values in the scaled metabolite matrix.
colSums(is.nan((met_scale_Akha)))


# Remove metabolite features containing missing values.
# Features with at least one NA are excluded from the RDA response matrix.
met_scale_Akha_clean <- met_scale_Akha[
  , colSums(is.na(met_scale_Akha)) == 0
]


# Confirm that no missing values remain after filtering.
anyNA(met_scale_Akha_clean)


# Check the dimensions of the metabolite and microbiota matrices.
# Both matrices should contain the same number and order of samples.
dim(met_scale_Akha_clean)
dim(Akha_bac_clr)


# Data transformation summary:
#
#   Metabolites:
#     Previously log-transformed and scaled.
#
#   Gut microbiota:
#     Raw counts + pseudocount (+1), followed by CLR transformation.
#
# These transformed datasets are subsequently combined for the
# multi-omics RDA.


# ==========================================================
# 4. RDA Using Metabolite Data
# ==========================================================
# This model evaluates the relationship between the metabolite profile
# and the predefined sample grouping variable.
#
# Response:
#   Scaled metabolite matrix.
#
# Explanatory variable:
#   Group.
#
# The constrained proportion represents the fraction of variation in the
# response matrix explained by the explanatory variable, whereas the
# unconstrained proportion represents residual variation.
#
# Note: This model is retained as a separate analysis before the integrated
# metabolite + microbiota RDA.
rda_model <- rda(met_scale_Akha_clean ~ Group, data = metadata_Akha)

summary(rda_model)
# Constrained Proportion: proportion of variation in the response matrix
# explained by the explanatory variable.
#
# Unconstrained Proportion: residual variation not explained by the model.



# ==========================================================
# 5. RDA Using Combined Metabolite + Gut Microbiota Data
# ==========================================================
# This model integrates metabolite and gut microbiota features into a
# single response matrix.
#
# Response matrix:
#   - Scaled metabolite features
#   - CLR-transformed microbial features
#
# Explanatory variable:
#   Group
#
# The analysis evaluates how strongly the predefined groups explain
# variation across the integrated multi-omics profile.


# Combine the transformed metabolite and microbiota matrices.
response_matrix <- cbind(met_scale_Akha_clean, Akha_bac_clr)
class(response_matrix)

dim(met_scale_Akha_clean)
dim(Akha_bac_clr)


# Run RDA with Group as the explanatory variable.
#
# RDA model:
#   Integrated metabolite + gut microbiota profile ~ Group
#
# The constrained variation represents the component of multi-omics
# variation associated with Group.
rda_model <- rda(response_matrix ~ Group, data = metadata_Akha)

summary(rda_model)
# Constrained Proportion: proportion of variation in the integrated response
# matrix explained by Group.
#
# Unconstrained Proportion: residual variation not explained by Group.


# Extract canonical coefficients and adjusted R-squared.
#
# Canonical coefficients describe the fitted relationship between the
# explanatory variables and constrained axes.
#
# Adjusted R-squared provides an adjustment for model complexity and sample
# size when estimating the proportion of variation explained by the model.
coef(rda_model)
RsquareAdj(rda_model)


# ==========================================================
# 6. Permutation-Based Significance Testing
# ==========================================================
# Statistical significance of the RDA model is evaluated using permutation
# tests with 999 permutations.
#
# set.seed() ensures reproducibility of the permutation procedure.


# Overall RDA model significance.
set.seed(123)
anova.cca(rda_model, permutations = 999)


# Significance of individual constrained axes.
set.seed(123)
anova.cca(rda_model, by = "axis", permutations = 999)


# Significance of individual model terms.
# In this analysis, Group is the primary explanatory variable.
set.seed(123)
anova.cca(rda_model, by = "terms", permutations = 999)



# ==========================================================
# 7. Base R RDA Plot
# ==========================================================
# Generate the default RDA ordination using scaling = 1.
#
# Scaling 1 emphasizes relationships among sample/site scores and
# approximates distances among sites.
plot(rda_model, scaling = 1)


# Interpretation:
#
# Scaling = 1:
#   Primarily emphasizes distances among samples/sites.
#
# Scaling = 2:
#   Primarily facilitates interpretation of relationships among response
#   variables through the angles and directions of their vectors.
#
# The two scalings serve different interpretive purposes and are therefore
# visualized separately below.



# ==========================================================
# 8. RDA Biplot – Scaling 2
# ==========================================================
# Create a publication-style RDA biplot using ggplot2.
#
# Scaling = 2 is used to facilitate interpretation of relationships among
# response variables such as microbial taxa and metabolites.


library(ggplot2)
library(ggrepel)


# Extract site and species/feature scores using scaling = 2.
rda_scores <- scores(rda_model, display = c("sites", "species"), scaling = 2)


# Convert site scores to a data frame and add the corresponding sample group.
sites_df <- as.data.frame(rda_scores$sites)
sites_df$Group <- metadata_Akha$Group


# Extract response-variable scores.
# These scores represent the positions/directions of the microbial and
# metabolite features in the ordination space.
species_scores <- as.data.frame(rda_scores$species)


# Calculate the Euclidean length of each feature vector.
# Longer vectors indicate features positioned farther from the origin in
# the displayed ordination.
arrow_lengths <- sqrt(rowSums(species_scores^2))


# Select the 50 features with the largest vector lengths.
#
# This selection is based on the magnitude of the displayed RDA vectors,
# rather than on formal statistical contribution values.
top_species <- names(sort(arrow_lengths, decreasing = TRUE))[1:50]


species_df <- as.data.frame(rda_scores$species)


# Retain only the selected features.
topX_species <- species_df[rownames(species_df) %in% top_species, ]


# Multiply feature coordinates to improve visual separation of arrows
# and labels in the plotted biplot.
scale_factor <- 2 # adjust the lenght of arrows
top_species_scaled <- topX_species * scale_factor



# ----------------------------------------------------------
# 8.1 Classify Features by Data Type
# ----------------------------------------------------------
# Assign each selected feature to its corresponding data type based on
# its variable name:
#
#   Met*  = metabolite
#   p__   = phylum
#   g__   = genus
#
# Features that do not match these naming patterns are assigned to the
# generic "feature" category.


topX_species_scaled <- top_species_scaled %>%
  mutate(
    # Assign feature category based on the naming convention.
    Feature = case_when(
      grepl("^Met", rownames(top_species_scaled)) ~ "Metabolite",
      grepl("^p__", rownames(top_species_scaled)) ~ "Phylum",
      grepl("^g__", rownames(top_species_scaled)) ~ "Genus",
      TRUE ~ "feature"
    )
  )


# Inspect the classified feature table.
head(topX_species_scaled)


# Create a final data frame containing feature names and RDA coordinates.
top50_scaling2_bac_Metabo <- data.frame(Variable = rownames(topX_species_scaled),
                                         topX_species_scaled) 

top100_scaling2_bac_Metabo <- data.frame(Variable = rownames(topX_species_scaled),
                                           topX_species_scaled) 


library(car)


# Export the selected feature coordinates for downstream analysis.
Export(top50_scaling2_bac_Metabo, "top50_scaling2_bac_Metabo_Akha.txt")
Export(top100_scaling2_bac_Metabo, "top100_scaling2_bac_Metabo_Akha.txt")


# ==========================================================
# 8.2 Variance Explained by RDA Axes
# ==========================================================
# Extract the proportion of constrained variation represented by the
# first two constrained axes.
summary(rda_model)$concont$importance
rda_var <- summary(rda_model)$cont$importance["Proportion Explained", 1:2]


# ----------------------------------------------------------
# 8.3 RDA Scaling 2 Plot
# ----------------------------------------------------------
# Define colors for the two Akha geographic groups.
col_key <-  c("Akha-CM" = "gold1", "Akha-CR" = "dodgerblue4")


# Create the RDA scaling 2 biplot.
#
# Sample points are colored by Group.
# Feature arrows are colored according to feature type.
sc2 <- ggplot() +
  
  # Sample scores
  geom_point(data = sites_df, aes(RDA1, PC1, color = Group), size = 2, alpha = 0.5) +
  
  # Feature vectors
  geom_segment(data = top50_scaling2_bac_Metabo, aes(x = 0, y = 0, xend = RDA1, 
                                               yend = PC1, color = Feature), 
               arrow = arrow(angle = 20, length = unit(0.1, "cm"), type = "closed"),
               linetype = 1, size = 0.5) +
  
  # Feature labels
  geom_text_repel(data = top50_scaling2_bac_Metabo, aes(RDA1, PC1, 
                                                  label = row.names(top50_scaling2_bac_Metabo), 
                                                  color = Feature),
                  size = 2.5, force = 10, xlim = c(-9, 6), ylim = c(-8, 5),
                  point.padding = 0.5, segment.colour = "grey", segment.alpha = 0.4,
                  max.overlaps = Inf) +
  
  # Color scale for both sample groups and feature categories
  scale_color_manual(values = c("Phylum" = "#0072B2",
                                "Genus" = "#009E73",
                                "Metabolite" = "#E69F00",
                                "Akha-CM" = "gold1", 
                                "Akha-CR" = "dodgerblue4")) +
  
  # Axis labels display the percentage of constrained variation explained
  # by the corresponding dimensions.
  labs(
    x = paste0("RDA1 (", round(100*rda_var[1],2), "%)"),
    y = paste0("PC1 (", round(100*rda_var[2],2), "%)")
  ) +
  
  # Add reference lines through the origin.
  geom_hline(yintercept = 0, linetype = 3, size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 3, size = 0.5) +
  
  theme_bw() +
  
  labs(title = "RDA Biplot (Scaling 2): Top 50 Contributing Features")+
  
  # Expand plotting limits slightly beyond the observed scores and feature
  # vectors to prevent labels and arrows from being clipped.
  coord_cartesian(
    xlim = c(min(sites_df$RDA1, top50_scaling2_bac_Metabo$RDA1) - 2.5, 
             max(sites_df$RDA1, top50_scaling2_bac_Metabo$RDA1) + 2),
    ylim = c(min(sites_df$PC1, top50_scaling2_bac_Metabo$PC1) - 2, 
             max(sites_df$PC1, top50_scaling2_bac_Metabo$PC1) + 2),
    clip = "off"
  )

sc2



# ==========================================================
# 9. RDA Biplot – Scaling 1
# ==========================================================
# Repeat the ordination using scaling = 1.
#
# Scaling 1 is used here to emphasize sample/site relationships and
# visualize the distribution of samples and group centroids.


# Extract site and feature scores using scaling = 1.
rda_scores <- scores(rda_model, display = c("sites", "species"), scaling = 1)

sites_df <- as.data.frame(rda_scores$sites)
sites_df$Group <- metadata_Akha$Group


# Extract feature scores.
species_scores <- as.data.frame(rda_scores$species)


# Calculate feature vector lengths.
arrow_lengths <- sqrt(rowSums(species_scores^2))


# Select the 50 features with the largest vector lengths.
top_species <- names(sort(arrow_lengths, decreasing = TRUE))[1:50]

species_df <- as.data.frame(rda_scores$species)


# Retain the selected features.
topX_species <- species_df[rownames(species_df) %in% top_species, ]



# ----------------------------------------------------------
# 9.1 Classify Selected Features
# ----------------------------------------------------------
# Reuse the previously scaled feature coordinates and assign each feature
# to a data type based on its naming convention.
topX_species_scaled <- top_species_scaled %>%
  mutate(
    Feature = case_when(
      grepl("^Met", rownames(top_species_scaled)) ~ "Metabolite",
      grepl("^p__", rownames(top_species_scaled)) ~ "Phylum",
      grepl("^g__", rownames(top_species_scaled)) ~ "Genus",
      TRUE ~ "feature"
    )
  )


# Inspect the feature classification.
head(topX_species_scaled)



# ----------------------------------------------------------
# 9.2 Calculate Group Centroids
# ----------------------------------------------------------
# Calculate the mean RDA coordinates for each sample group.
# These centroids summarize the average position of each group in the
# ordination space.
centroids <- aggregate(cbind(RDA1, PC1) ~ Group, 
                       data = data.frame(sites_df), FUN = mean)


# Inspect group centroid coordinates.
print(centroids)


# Extract the proportion of variation represented by the first two
# constrained axes.
summary(rda_model)$concont$importance
rda_var <- summary(rda_model)$cont$importance["Proportion Explained", 1:2]



# ----------------------------------------------------------
# 9.3 RDA Scaling 1 Plot
# ----------------------------------------------------------
# Plot:
#   - Individual samples as semi-transparent points
#   - Group centroids as triangles
#   - Selected feature vectors and labels
#
# This visualization is intended primarily to show sample distribution
# and the relative position of the predefined groups.
sc1 <- ggplot() +
  
  # Individual sample scores
  geom_point(data = sites_df, aes(RDA1, PC1, color = Group), size = 2, alpha = 0.5) +
  
  # Group centroids
  geom_point(data = centroids, aes(RDA1, PC1, color = Group), size = 3, shape = 17) +
  
  # Group colors
  scale_color_manual(values = col_key) +
  
  
  # Feature vectors
  geom_segment(data = topX_species_scaled, 
               aes(x = 0, y = 0, xend = RDA1, yend = PC1),
               color = "grey60",
               arrow = arrow(angle = 20, length = unit(0.2, "cm"), type = "closed"),
               linetype = 1, size = 0.5, alpha = 0.1) +
  
  # Feature labels
  geom_text_repel(data = topX_species_scaled, 
                  aes(RDA1, PC1, label = row.names(topX_species_scaled)),
                  color = "grey60",  alpha = 0.5,  
                  size = 2.5, force = 10, xlim = c(-3, 3), ylim = c(-2.5, 2.5),
                  point.padding = 0.5, segment.colour = "grey80", segment.alpha = 0.1,
                  max.overlaps = 100) +
  
  # Axis labels
  labs(
    x = paste0("RDA1 (", round(100*rda_var[1],2), "%)"),
    y = paste0("PC1 (", round(100*rda_var[2],2), "%)")
  ) +
  
  # Reference lines at the origin
  geom_hline(yintercept = 0, linetype = 3, size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 3, size = 0.5) +
  
  theme_bw() +
  
  labs(title = "RDA Plot (Scaling 1): Sample Distribution and Group Centroids")+
  
  # Expand the plotting area to accommodate feature labels and vectors.
  coord_cartesian(
    xlim = c(min(sites_df$RDA1, topX_species_scaled$RDA1) - 1.5, 
             max(sites_df$RDA1, topX_species_scaled$RDA1) + 1.5),
    ylim = c(min(sites_df$PC1, topX_species_scaled$PC1) - 1, 
             max(sites_df$PC1, topX_species_scaled$PC1) + 1),
    clip = "off"
  )


sc1


# Display the scaling 1 and scaling 2 RDA plots together.
library(gridExtra)
grid.arrange(
  sc1, 
  sc2,
  nrow = 2
)



# ==========================================================
# 10. Directional Similarity Between RDA Feature Vectors
# ==========================================================
# Calculate cosine similarity between pairs of feature vectors based on
# their RDA1 and PC1 coordinates.
#
# Cosine similarity measures the similarity in direction between two vectors
# independently of their magnitude.
#
# Interpretation:
#   +1  = vectors point in the same direction
#    0  = vectors are orthogonal
#   -1  = vectors point in opposite directions
#
# This analysis is used to identify microbial–metabolite features that
# exhibit similar directional relationships in the RDA ordination.


# Function to calculate cosine similarity between two vectors.
cosine_similarity <- function(x1, y1, x2, y2) {
  dot_product <- x1 * x2 + y1 * y2
  magnitude1 <- sqrt(x1^2 + y1^2)
  magnitude2 <- sqrt(x2^2 + y2^2)
  return(dot_product / (magnitude1 * magnitude2))
}


# Cosine similarity is the normalized dot product of two vectors.
#
#   cos(theta) = (A · B) / (|A| |B|)
#
# Values close to +1 indicate similar directions, whereas values close
# to -1 indicate opposing directions.


# Calculate the angle between vectors in degrees.
# Values are constrained to [-1, 1] before applying acos() to avoid
# numerical rounding errors that could otherwise produce NaN values.
angle_deg_safe <- function(cos_sim) {
  acos(pmin(pmax(cos_sim, -1), 1)) * (180 / pi)
}



# ----------------------------------------------------------
# 10.1 Assign Features to RDA Quadrants
# ----------------------------------------------------------
# Classify each feature according to the signs of its RDA1 and PC1
# coordinates.
#
#   PP = RDA1 > 0, PC1 > 0
#   PN = RDA1 > 0, PC1 < 0
#   NP = RDA1 < 0, PC1 > 0
#   NN = RDA1 < 0, PC1 < 0
#
# Features located in the same quadrant generally point in broadly similar
# directions in the displayed two-dimensional ordination.


library(dplyr)
df <- top50_scaling2_bac_Metabo %>%
  mutate(quadrant = case_when(
    RDA1 > 0 & PC1 > 0 ~ "PP",
    RDA1 > 0 & PC1 < 0 ~ "PN",
    RDA1 < 0 & PC1 > 0 ~ "NP",
    RDA1 < 0 & PC1 < 0 ~ "NN"
  ))



# ----------------------------------------------------------
# 10.2 Calculate Pairwise Feature Similarity
# ----------------------------------------------------------
# Generate all unique pairs of selected features and calculate:
#
#   - Feature names
#   - Quadrant membership
#   - Cosine similarity
#   - Angular separation
#
# Only unique pairs are retained by restricting the analysis to i < j.


library(tidyr)
Akha_pairs <- expand.grid(i = 1:nrow(df), j = 1:nrow(df)) %>%
  filter(i < j) %>%
  mutate(
    Feature1 = df$Variable[i],
    Feature2 = df$Variable[j],
    
    quad1 = df$quadrant[i],
    quad2 = df$quadrant[j],
    
    cos_sim = cosine_similarity(
      df$RDA1[i], df$PC1[i],
      df$RDA1[j], df$PC1[j]
    ),
    
    angle_deg = angle_deg_safe(cos_sim)
  ) %>%
  arrange(desc(cos_sim)) %>%
  left_join(Pub_Classy_05_clean,
            by = c("Feature1" = "Code"))


head(Akha_pairs)



# ----------------------------------------------------------
# 10.3 Identify Strong Microbe–Metabolite Associations
# ----------------------------------------------------------
# Identify microbial–metabolite pairs with:
#
#   1. Both features located in the same quadrant
#   2. One feature being a phylum/genus and the other a metabolite
#   3. Cosine similarity > 0.90
#
# A cosine similarity > 0.90 indicates that the two vectors have a very
# similar direction in the displayed RDA ordination.
#
# IMPORTANT:
# These are ordination-based directional similarities and should not be
# interpreted as direct statistical correlations between the original
# microbial and metabolite measurements.


library(stringr)

Akha_strongcosine <- Akha_pairs %>%
  filter(
    quad1 == quad2,
    (grepl("^(p__|g__)", Feature1) & grepl("^Met", Feature2)) |
      (grepl("^Met", Feature1) & grepl("^(p__|g__)", Feature2)),
    cos_sim > 0.9
  ) %>%
  arrange(desc(cos_sim)) 

head(Akha_strongcosine)


# Export the complete pairwise feature comparison and the subset of
# strong microbe–metabolite directional relationships.
Export(Akha_pairs, "Akha_pairs.txt")
Export(Akha_strongcosine, "Akha_strongcosine.txt")


# Interpretation:
#
#   cos_sim:
#     Measures directional similarity between two RDA vectors.
#
#   angle_deg:
#     Measures the angular separation between two RDA vectors.
#
# Smaller angular separation corresponds to greater directional similarity.



# ============================================================
# 11. Identify Key Features Associated with Geographic Groups
# ============================================================
# Select features with relatively large absolute RDA1 coordinates.
#
# Features on the negative side of RDA1 are assigned to Akha-CM,
# whereas features on the positive side are assigned to Akha-CR.
#
# This provides a descriptive grouping of features according to their
# position along the first constrained axis.


key_Akha <- top50_scaling2_bac_Metabo %>% 
  filter(abs(RDA1) >= 0.5) %>% 
  select(!PC1) %>%
  arrange(RDA1) %>%
  mutate(Group_association = 
    case_when(
    RDA1 < 0 ~ "Akha-CM",
    RDA1 > 0 ~ "Akha-CR"
    ))


# Export features associated with the two sides of the first constrained
# axis.
Export(key_Akha, "key_Akha.txt")



# ==========================================================
# 12. Microbe–Metabolite Relationships Among Key Features
# ==========================================================
# Restrict the strong cosine-similarity pairs to features that were also
# identified as key features based on their RDA1 position.
#
# The final table contains:
#   - Feature names
#   - Quadrant information
#   - Cosine similarity
#   - Angular separation
#   - Metabolite annotation
#   - Metabolite subclass


key_Akha_cossim <- Akha_strongcosine %>%
filter(Feature1 %in% key_Akha$Variable,
       Feature2 %in% key_Akha$Variable) %>%
  select(Feature1, Feature2, quad1, quad2, cos_sim, angle_deg, 
         Compound_Clean, Subclass) %>%
  filter(cos_sim >= 0.9)


# Export the final set of high-directional-similarity relationships among
# key RDA features.
Export(key_Akha_cossim, "key_Akha_cossim.txt")







####### The same analytical workflows were applied to the Lahu, CM, and CR groups. ########


# ==========================================================
# 2. RDA – Lahu Ethnic Group 
# ==========================================================

# Count ASVs for Lahu, select significant taxa at phylum and genus levels
count_Lahu <- count_asv %>%
  filter(Ethnicity == "Lahu") %>%
  select(any_of(names(sig_taxo_Lahu)[13:ncol(sig_taxo_Lahu)])) %>%
  as.matrix()

# confirm sample order
identical(rownames(count_Lahu), rownames(metadata_Lahu))


# Apply CLR transformation to Count Phylum and Genus
Lahu_bac_clr <- clr(count_Lahu + 1)
head(Lahu_bac_clr)
rownames(Lahu_bac_clr)

# Check is there any na, NaN 
colSums(is.nan((met_scale_Lahu)))

met_scale_Lahu_clean <- met_scale_Lahu



dim(met_scale_Lahu_clean)
dim(Lahu_bac_clr)


# ✅ metabo_scaled = log-transformed + scaled
# ✅ bac count = add pseudocount (+1) + apply CLR transformation


#-------------------------------------
# RDA: Metabolites or Gut microbiota
#-------------------------------------
# Run RDA with  group as explanatory variable
# RDA model: multi-omics response explained by  group
rda_model <- rda(met_scale_Lahu_clean ~ Group, data = metadata_Lahu)

summary(rda_model)
# Constrained Proportion: variance of Y explained by X 
# Unconstrained Proportion: unexplained variance in Y 



#---------------------------------
# RDA: Metabolites+Gut microbiota
#---------------------------------
#-- RDA with multi-omics as response
# Associations between genus and metabolites across groups
# •	✅ Response matrix: both bac_clr and metabo_scaled
# •	✅ Explanatory variable: Metadata$Group 

# Combine CLR-genus and scaled metabolites
response_matrix <- cbind(met_scale_Lahu_clean, Lahu_bac_clr)
class(response_matrix)

dim(met_scale_Lahu_clean)
dim(Lahu_bac_clr)

# Run RDA with  group as explanatory variable
# RDA model: multi-omics response explained by  group
rda_model <- rda(response_matrix ~ Group, data = metadata_Lahu)

summary(rda_model)
# Constrained Proportion: variance of Y explained by X 
# Unconstrained Proportion: unexplained variance in Y 


# canonical coefficients
coef(rda_model)
RsquareAdj(rda_model)


# Test significance (PERMUTATION-BASED ANOVA)
# Test overall model
set.seed(123)
anova.cca(rda_model, permutations = 999)

# Test each constrained axis
set.seed(123)
anova.cca(rda_model, by = "axis", permutations = 999)

# Test terms (Group, metabolites, etc.)
set.seed(123)
anova.cca(rda_model, by = "terms", permutations = 999)


# Plot the RDA
plot(rda_model, scaling = 1)  


#🔬 Interpretation

# • Scaling = 1 = distances between samples matter
# •	Scaling = 2 = angles between variables (species) matter (better when you want to explore 
# relationships between features like genus or metabolites).

#This RDA will help you answer:
#•	✅ How much variance in combined genus + metabolite profiles is explained by  group?
#•	✅ Are there visible clusters or shifts in profiles across  categories?
#•	✅ Which genera and metabolites drive the separation?



## Redundancy Analysis (RDA) biplot with type 2 scaling 
library(ggplot2)
library(ggrepel)


# Creating clean RDA biplots using ggplot2
# Extract scores
rda_scores <- scores(rda_model, display = c("sites", "species"), scaling = 2)

# Site score 
sites_df <- as.data.frame(rda_scores$sites)
sites_df$Group <- metadata_Lahu$Group  # add grouping info for coloring

# Label only the top contributing species (genus, metabolites)
# Get species scores (with scaling 2)
# This declutters the plot by showing only the most influential features.
species_scores <- as.data.frame(rda_scores$species)

# Calculate vector lengths
arrow_lengths <- sqrt(rowSums(species_scores^2))

# Get top features
top_species <- names(sort(arrow_lengths, decreasing = TRUE))[1:50]

species_df <- as.data.frame(rda_scores$species)

# Keep only top species Scaling = 2
# Get the top rownames based on arrow_lengths
topX_species <- species_df[rownames(species_df) %in% top_species, ]


# Scale the species scores (e.g., multiplying by x to extend the arrows)
scale_factor <- 2
top_species_scaled <- topX_species * scale_factor




# RDA plot with ggplot2

# Create the new column with categories
topX_species_scaled <- top_species_scaled %>%
  mutate(
    # Ion type assignment based on rownames (posxx or negxx)
    Feature = case_when(
      grepl("^Met", rownames(top_species_scaled)) ~ "Metabolite",
      grepl("^p__", rownames(top_species_scaled)) ~ "Phylum",
      grepl("^g__", rownames(top_species_scaled)) ~ "Genus",
      TRUE ~ "feature"  # Default to "feature" 
    )
  )


# View the updated data frame
head(topX_species_scaled)

# create dataframe
top50_scaling2_bac_Metabo <- data.frame(Variable = rownames(topX_species_scaled),
                                         topX_species_scaled) 

top100_scaling2_bac_Metabo <- data.frame(Variable = rownames(topX_species_scaled),
                                         topX_species_scaled) 


library(car)
Export(top50_scaling2_bac_Metabo, "top50_scaling2_bac_Metabo_Lahu.txt")


# proportion of variance explained by the RDA axes
summary(rda_model)$concont$importance
rda_var <- summary(rda_model)$cont$importance["Proportion Explained", 1:2]


# RDA ggplot scaling = 2
# Define the col_key with the colors for different groups
col_key <-  c("Lahu-CM" = "lightblue2", "Lahu-CR" = "coral1")


# Plotting
sc2 <- ggplot() +
  # Points for sites (sites_df) colored based on 'expl.var' or your relevant grouping variable
  geom_point(data = sites_df, aes(RDA1, PC1, color = Group), size = 2, alpha = 0.5) +
  
  # Use category to color arrows differently (Genus and Metabolite)
  geom_segment(data = top50_scaling2_bac_Metabo, aes(x = 0, y = 0, xend = RDA1, 
                                                      yend = PC1, color = Feature), 
               arrow = arrow(angle = 20, length = unit(0.1, "cm"), type = "closed"),
               linetype = 1, size = 0.5) +
  
  # Text labels with the same color as arrows
  geom_text_repel(data = top50_scaling2_bac_Metabo, aes(RDA1, PC1, 
                                                         label = row.names(top50_scaling2_bac_Metabo), 
                                                         color = Feature),
                  size = 2.5, force = 10, xlim = c(-9, 6), ylim = c(-8, 5),
                  point.padding = 0.5, segment.colour = "grey", segment.alpha = 0.4,
                  max.overlaps = 100) +
  
  # Color scale for both points and arrows 
  scale_color_manual(values = c("Phylum" = "#0072B2",
                                "Genus" = "#009E73",
                                "Metabolite" = "#E69F00",
                                "Lahu-CM" = "lightblue2", "Lahu-CR" = "coral1")) +
  
  # Labels for the axes
  labs(
    x = paste0("RDA1 (", round(100*rda_var[1],2), "%)"),
    y = paste0("PC1 (", round(100*rda_var[2],2), "%)")
  ) +
  
  # Horizontal and vertical lines at 0
  geom_hline(yintercept = 0, linetype = 3, size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 3, size = 0.5) +
  
  theme_bw() +
  
  # Title of the plot
  labs(title = "RDA Biplot (Scaling 2): Top 50 Contributing Features")+
  
  coord_cartesian(
    xlim = c(min(sites_df$RDA1, top50_scaling2_bac_Metabo$RDA1) - 1.5, 
             max(sites_df$RDA1, top50_scaling2_bac_Metabo$RDA1) + 1.5),
    ylim = c(min(sites_df$PC1, top50_scaling2_bac_Metabo$PC1) - 1, 
             max(sites_df$PC1, top50_scaling2_bac_Metabo$PC1) + 1),
    clip = "off"
  )

sc2



## Redundancy Analysis (RDA) biplot with type 1 scaling 

# Creating clean RDA biplots using ggplot2
# Extract scores
rda_scores <- scores(rda_model, display = c("sites", "species"), scaling = 1)

sites_df <- as.data.frame(rda_scores$sites)
sites_df$Group <- metadata_Lahu$Group  # add grouping info for coloring

# Label only the top contributing species (genus, metabolites)
# Get species scores (with scaling 2)
# This declutters the plot by showing only the most influential features.
species_scores <- as.data.frame(rda_scores$species)

# Calculate vector lengths
arrow_lengths <- sqrt(rowSums(species_scores^2))

# Get top features
top_species <- names(sort(arrow_lengths, decreasing = TRUE))[1:50]

species_df <- as.data.frame(rda_scores$species)

# Keep only top species Scaling = 2
# Get the top rownames based on arrow_lengths
topX_species <- species_df[rownames(species_df) %in% top_species, ]




# RDA plot with ggplot2
# Create the new column with categories
topX_species_scaled <- top_species_scaled %>%
  mutate(
    # Ion type assignment based on rownames (posxx or negxx)
    Feature = case_when(
      grepl("^Met", rownames(top_species_scaled)) ~ "Metabolite",
      grepl("^p__", rownames(top_species_scaled)) ~ "Phylum",
      grepl("^g__", rownames(top_species_scaled)) ~ "Genus",
      TRUE ~ "feature"  # Default to "feature" 
    )
  )

# View the updated data frame
head(topX_species_scaled)



# Extract the centroids 
centroids <- aggregate(cbind(RDA1, PC1) ~ Group, 
                       data = data.frame(sites_df), FUN = mean)

# View the centroids
print(centroids)

# proportion of variance explained by the RDA axes
summary(rda_model)$concont$importance
rda_var <- summary(rda_model)$cont$importance["Proportion Explained", 1:2]


sc1 <- ggplot() +
  # Points for sites (sites_df) colored based on 'Group'
  geom_point(data = sites_df, aes(RDA1, PC1, color = Group), size = 2, alpha = 0.5) +
  
  # Points for centroids (cent) colored based on 'Group'
  geom_point(data = centroids, aes(RDA1, PC1, color = Group), size = 3, shape = 17) +
  
  # Color scale for points
  scale_color_manual(values = col_key) +
  
  
  geom_segment(data = topX_species_scaled, 
               aes(x = 0, y = 0, xend = RDA1, yend = PC1),
               color = "grey60",
               arrow = arrow(angle = 20, length = unit(0.2, "cm"), type = "closed"),
               linetype = 1, size = 0.5, alpha = 0.1) +
  
  # Text labels with the same color as arrows
  geom_text_repel(data = topX_species_scaled, 
                  aes(RDA1, PC1, label = row.names(topX_species_scaled)),
                  color = "grey60",  alpha = 0.5,  
                  size = 2.5, force = 10, xlim = c(-3, 3), ylim = c(-2.5, 2.5),
                  point.padding = 0.5, segment.colour = "grey80", segment.alpha = 0.1,
                  max.overlaps = 100) +
  
  # Labels for the axes
  labs(
    x = paste0("RDA1 (", round(100*rda_var[1],2), "%)"),
    y = paste0("PC1 (", round(100*rda_var[2],2), "%)")
  ) +
  
  # Horizontal and vertical lines at 0
  geom_hline(yintercept = 0, linetype = 3, size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 3, size = 0.5) +
  
  theme_bw() +
  
  # Title of the plot
  labs(title = "RDA Plot (Scaling 1): Sample Distribution and Group Centroids")+
  coord_cartesian(
    xlim = c(min(sites_df$RDA1, topX_species_scaled$RDA1) - 1.5, 
             max(sites_df$RDA1, topX_species_scaled$RDA1) + 1.5),
    ylim = c(min(sites_df$PC1, topX_species_scaled$PC1) - 1, 
             max(sites_df$PC1, topX_species_scaled$PC1) + 1),
    clip = "off"
  )


sc1


# Multiple plots
library(gridExtra)
grid.arrange(
  sc1, 
  sc2,
  nrow = 2
)


# Prepare the dataframe

# • cosineNP = [-x, y], RDA1 < 0 & PC1 > 0
# • cosinePP = [x, y], RDA1 > 0 & PC1 > 0
# • cosineNN = [-x, -y], RDA1 < 0 & PC1 < 0
# • cosinePN = [x, -y], RDA1 > 0 & PC1 < 0


library(dplyr)
df <- top50_scaling2_bac_Metabo %>%
  mutate(quadrant = case_when(
    RDA1 > 0 & PC1 > 0 ~ "PP",
    RDA1 > 0 & PC1 < 0 ~ "PN",
    RDA1 < 0 & PC1 > 0 ~ "NP",
    RDA1 < 0 & PC1 < 0 ~ "NN"
  ))


# Compute pairwise cosine similarity
library(tidyr)
Lahu_pairs <- expand.grid(i = 1:nrow(df), j = 1:nrow(df)) %>%
  filter(i < j) %>%
  mutate(
    Feature1 = df$Variable[i],
    Feature2 = df$Variable[j],
    
    quad1 = df$quadrant[i],
    quad2 = df$quadrant[j],
    
    cos_sim = cosine_similarity(
      df$RDA1[i], df$PC1[i],
      df$RDA1[j], df$PC1[j]
    ),
    
    angle_deg = angle_deg_safe(cos_sim)
  ) %>%
  arrange(desc(cos_sim)) %>%
  left_join(Pub_Classy_05_clean,
            by = c("Feature1" = "Code"))


head(Lahu_pairs)


# Microbe–metabolite relationships in both directions
library(stringr)

Lahu_strongcosine <- Lahu_pairs %>%
  filter(
    quad1 == quad2,
    (grepl("^(p__|g__)", Feature1) & grepl("^Met", Feature2)) |
      (grepl("^Met", Feature1) & grepl("^(p__|g__)", Feature2)),
    cos_sim > 0.9
  ) %>%
  arrange(desc(cos_sim)) 

head(Lahu_strongcosine)

Export(Lahu_pairs, "Lahu_pairs.txt")
Export(Lahu_strongcosine, "Lahu_strongcosine.txt")

# Interpretation:
# •	cos_sim → directional correlation
# •	angle_deg → angular separation

key_Lahu <- top50_scaling2_bac_Metabo %>% 
  filter(abs(RDA1) >= 0.5) %>% 
  select(!PC1) %>%
  arrange(RDA1) %>%
  mutate(Group_association = 
           case_when(
             RDA1 < 0 ~ "Lahu-CM",
             RDA1 > 0 ~ "Lahu-CR"
           ))

Export(key_Lahu, "key_Lahu.txt")


key_Lahu_cossim <- Lahu_strongcosine %>%
  filter(Feature1 %in% key_Lahu$Variable,
         Feature2 %in% key_Lahu$Variable) %>%
  select(Feature1, Feature2, quad1, quad2, cos_sim, angle_deg, 
         Compound_Clean, Subclass) %>%
  filter(cos_sim >= 0.9)

Export(key_Lahu_cossim, "key_Lahu_cossim.txt")









# ==========================================================
# 3. RDA – CM Group 
# ==========================================================

# Count ASVs for CM, select significant taxa at phylum and genus levels
count_CM <- count_asv %>%
  filter(Area == "CM") %>%
  select(any_of(names(sig_taxo_CM)[13:ncol(sig_taxo_CM)])) %>%
  as.matrix()

# confirm sample order
identical(rownames(count_CM), rownames(metadata_CM))


# Apply CLR transformation to Count Phylum and Genus
CM_bac_clr <- clr(count_CM + 1)
head(CM_bac_clr)
rownames(CM_bac_clr)

# Check is there any na, NaN 
colSums(is.nan((met_scale_CM)))

met_scale_CM_clean <- met_scale_CM



dim(met_scale_CM_clean)
dim(CM_bac_clr)


# ✅ metabo_scaled = log-transformed + scaled
# ✅ bac count = add pseudocount (+1) + apply CLR transformation


#-------------------------------------
# RDA: Metabolites or Gut microbiota
#-------------------------------------
# Run RDA with  group as explanatory variable
# RDA model: multi-omics response explained by  group
rda_model <- rda(CM_bac_clr ~ Group, data = metadata_CM)

summary(rda_model)
# Constrained Proportion: variance of Y explained by X 
# Unconstrained Proportion: unexplained variance in Y 



#---------------------------------
# RDA: Metabolites+Gut microbiota
#---------------------------------
#-- RDA with multi-omics as response
# Associations between genus and metabolites across groups
# •	✅ Response matrix: both bac_clr and metabo_scaled
# •	✅ Explanatory variable: Metadata$Group 

# Combine CLR-genus and scaled metabolites
response_matrix <- cbind(met_scale_CM_clean, CM_bac_clr)
class(response_matrix)

dim(met_scale_CM_clean)
dim(CM_bac_clr)

# Run RDA with  group as explanatory variable
# RDA model: multi-omics response explained by  group
rda_model <- rda(response_matrix ~ Group, data = metadata_CM)

summary(rda_model)
# Constrained Proportion: variance of Y explained by X 
# Unconstrained Proportion: unexplained variance in Y 


# canonical coefficients
coef(rda_model)
RsquareAdj(rda_model)


# Test significance (PERMUTATION-BASED ANOVA)
# Test overall model
set.seed(123)
anova.cca(rda_model, permutations = 999)

# Test each constrained axis
set.seed(123)
anova.cca(rda_model, by = "axis", permutations = 999)

# Test terms (Group, metabolites, etc.)
set.seed(123)
anova.cca(rda_model, by = "terms", permutations = 999)


# Plot the RDA
plot(rda_model, scaling = 1)  


#🔬 Interpretation

# • Scaling = 1 = distances between samples matter
# •	Scaling = 2 = angles between variables (species) matter (better when you want to explore 
# relationships between features like genus or metabolites).

#This RDA will help you answer:
#•	✅ How much variance in combined genus + metabolite profiles is explained by  group?
#•	✅ Are there visible clusters or shifts in profiles across  categories?
#•	✅ Which genera and metabolites drive the separation?



## Redundancy Analysis (RDA) biplot with type 2 scaling 
library(ggplot2)
library(ggrepel)


# Creating clean RDA biplots using ggplot2
# Extract scores
rda_scores <- scores(rda_model, display = c("sites", "species"), scaling = 2)

# Site score 
sites_df <- as.data.frame(rda_scores$sites)
sites_df$Group <- metadata_CM$Group  # add grouping info for coloring

# Label only the top contributing species (genus, metabolites)
# Get species scores (with scaling 2)
# This declutters the plot by showing only the most influential features.
species_scores <- as.data.frame(rda_scores$species)

# Calculate vector lengths
arrow_lengths <- sqrt(rowSums(species_scores^2))

# Get top features
top_species <- names(sort(arrow_lengths, decreasing = TRUE))[1:50]

species_df <- as.data.frame(rda_scores$species)

# Keep only top species Scaling = 2
# Get the top rownames based on arrow_lengths
topX_species <- species_df[rownames(species_df) %in% top_species, ]


# Scale the species scores (e.g., multiplying by x to extend the arrows)
scale_factor <- 2
top_species_scaled <- topX_species * scale_factor




# RDA plot with ggplot2

# Create the new column with categories
topX_species_scaled <- top_species_scaled %>%
  mutate(
    # Ion type assignment based on rownames (posxx or negxx)
    Feature = case_when(
      grepl("^Met", rownames(top_species_scaled)) ~ "Metabolite",
      grepl("^p__", rownames(top_species_scaled)) ~ "Phylum",
      grepl("^g__", rownames(top_species_scaled)) ~ "Genus",
      TRUE ~ "feature"  # Default to "feature" 
    )
  )


# View the updated data frame
head(topX_species_scaled)

# create dataframe
top50_scaling2_bac_Metabo <- data.frame(Variable = rownames(topX_species_scaled),
                                         topX_species_scaled) 

top100_scaling2_bac_Metabo <- data.frame(Variable = rownames(topX_species_scaled),
                                         topX_species_scaled) 


library(car)
Export(top50_scaling2_bac_Metabo, "top50_scaling2_bac_Metabo_CM.txt")


# proportion of variance explained by the RDA axes
summary(rda_model)$concont$importance
rda_var <- summary(rda_model)$cont$importance["Proportion Explained", 1:2]


# RDA ggplot scaling = 2
# Define the col_key with the colors for different groups
col_key <-  c("Akha-CM" = "gold1",
              "Lahu-CM" = "lightblue2",
              "Khuen-CM" = "darkorange")


# Plotting
sc2 <- ggplot() +
  # Points for sites (sites_df) colored based on 'expl.var' or your relevant grouping variable
  geom_point(data = sites_df, aes(RDA1, RDA2, color = Group), size = 2, alpha = 0.5) +
  
  # Use category to color arrows differently (Genus and Metabolite)
  geom_segment(data = top50_scaling2_bac_Metabo, aes(x = 0, y = 0, xend = RDA1, 
                                                      yend = RDA2, color = Feature), 
               arrow = arrow(angle = 20, length = unit(0.1, "cm"), type = "closed"),
               linetype = 1, size = 0.5) +
  
  # Text labels with the same color as arrows
  geom_text_repel(data = top50_scaling2_bac_Metabo, aes(RDA1, RDA2, 
                                                         label = row.names(top50_scaling2_bac_Metabo), 
                                                         color = Feature),
                  size = 2.5, force = 10, xlim = c(-9, 6), ylim = c(-8, 5),
                  point.padding = 0.5, segment.colour = "grey", segment.alpha = 0.4,
                  max.overlaps = Inf) +
  
  # Color scale for both points and arrows 
  scale_color_manual(values = c("Phylum" = "#0072B2",
                                "Genus" = "#009E73",
                                "Metabolite" = "#E69F00",
                                "Akha-CM" = "gold1",
                                "Lahu-CM" = "lightblue2",
                                "Khuen-CM" = "darkorange")) +
  
  # Labels for the axes
  labs(
    x = paste0("RDA1 (", round(100*rda_var[1],2), "%)"),
    y = paste0("RDA2 (", round(100*rda_var[2],2), "%)")
  ) +
  
  # Horizontal and vertical lines at 0
  geom_hline(yintercept = 0, linetype = 3, size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 3, size = 0.5) +
  
  theme_bw() +
  
  # Title of the plot
  labs(title = "RDA Biplot (Scaling 2): Top 50 Contributing Features")+
  
  coord_cartesian(
    xlim = c(min(sites_df$RDA1, top50_scaling2_bac_Metabo$RDA1) - 1.5, 
             max(sites_df$RDA1, top50_scaling2_bac_Metabo$RDA1) + 1.5),
    ylim = c(min(sites_df$RDA2, top50_scaling2_bac_Metabo$RDA2) - 1, 
             max(sites_df$RDA2, top50_scaling2_bac_Metabo$RDA2) + 1),
    clip = "off"
  )

sc2



## Redundancy Analysis (RDA) biplot with type 1 scaling 

# Creating clean RDA biplots using ggplot2
# Extract scores
rda_scores <- scores(rda_model, display = c("sites", "species"), scaling = 1)

sites_df <- as.data.frame(rda_scores$sites)
sites_df$Group <- metadata_CM$Group  # add grouping info for coloring

# Label only the top contributing species (genus, metabolites)
# Get species scores (with scaling 2)
# This declutters the plot by showing only the most influential features.
species_scores <- as.data.frame(rda_scores$species)

# Calculate vector lengths
arrow_lengths <- sqrt(rowSums(species_scores^2))

# Get top features
top_species <- names(sort(arrow_lengths, decreasing = TRUE))[1:50]

species_df <- as.data.frame(rda_scores$species)

# Keep only top species Scaling = 2
# Get the top rownames based on arrow_lengths
topX_species <- species_df[rownames(species_df) %in% top_species, ]




# RDA plot with ggplot2
# Create the new column with categories
topX_species_scaled <- top_species_scaled %>%
  mutate(
    # Ion type assignment based on rownames (posxx or negxx)
    Feature = case_when(
      grepl("^Met", rownames(top_species_scaled)) ~ "Metabolite",
      grepl("^p__", rownames(top_species_scaled)) ~ "Phylum",
      grepl("^g__", rownames(top_species_scaled)) ~ "Genus",
      TRUE ~ "feature"  # Default to "feature" 
    )
  )

# View the updated data frame
head(topX_species_scaled)



# Extract the centroids 
centroids <- aggregate(cbind(RDA1, RDA2) ~ Group, 
                       data = data.frame(sites_df), FUN = mean)

# View the centroids
print(centroids)

# proportion of variance explained by the RDA axes
summary(rda_model)$concont$importance
rda_var <- summary(rda_model)$cont$importance["Proportion Explained", 1:2]


sc1 <- ggplot() +
  # Points for sites (sites_df) colored based on 'Group'
  geom_point(data = sites_df, aes(RDA1, RDA2, color = Group), size = 2, alpha = 0.5) +
  
  # Points for centroids (cent) colored based on 'Group'
  geom_point(data = centroids, aes(RDA1, RDA2, color = Group), size = 3, shape = 17) +
  
  # Color scale for points
  scale_color_manual(values = col_key) +
  
  
  geom_segment(data = topX_species_scaled, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               color = "grey60",
               arrow = arrow(angle = 20, length = unit(0.2, "cm"), type = "closed"),
               linetype = 1, size = 0.5, alpha = 0.1) +
  
  # Text labels with the same color as arrows
  geom_text_repel(data = topX_species_scaled, 
                  aes(RDA1, RDA2, label = row.names(topX_species_scaled)),
                  color = "grey60",  alpha = 0.5,  
                  size = 2.5, force = 10, xlim = c(-3, 3), ylim = c(-2.5, 2.5),
                  point.padding = 0.5, segment.colour = "grey80", segment.alpha = 0.1,
                  max.overlaps = Inf) +
  
  # Labels for the axes
  labs(
    x = paste0("RDA1 (", round(100*rda_var[1],2), "%)"),
    y = paste0("RDA2 (", round(100*rda_var[2],2), "%)")
  ) +
  
  # Horizontal and vertical lines at 0
  geom_hline(yintercept = 0, linetype = 3, size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 3, size = 0.5) +
  
  theme_bw() +
  
  # Title of the plot
  labs(title = "RDA Plot (Scaling 1): Sample Distribution and Group Centroids")+
  coord_cartesian(
    xlim = c(min(sites_df$RDA1, topX_species_scaled$RDA1) - 1.5, 
             max(sites_df$RDA1, topX_species_scaled$RDA1) + 1.5),
    ylim = c(min(sites_df$RDA2, topX_species_scaled$RDA2) - 1, 
             max(sites_df$RDA2, topX_species_scaled$RDA2) + 1),
    clip = "off"
  )


sc1


# Multiple plots
library(gridExtra)
grid.arrange(
  sc1, 
  sc2,
  nrow = 2
)


# Prepare the dataframe

# • cosineNP = [-x, y], RDA1 < 0 & RDA2 > 0
# • cosinePP = [x, y], RDA1 > 0 & RDA2 > 0
# • cosineNN = [-x, -y], RDA1 < 0 & RDA2 < 0
# • cosinePN = [x, -y], RDA1 > 0 & RDA2 < 0



library(dplyr)
df <- top50_scaling2_bac_Metabo %>%
  mutate(quadrant = case_when(
    RDA1 > 0 & RDA2 > 0 ~ "PP",
    RDA1 > 0 & RDA2 < 0 ~ "PN",
    RDA1 < 0 & RDA2 > 0 ~ "NP",
    RDA1 < 0 & RDA2 < 0 ~ "NN"
  ))


# Compute pairwise cosine similarity
library(tidyr)
CM_pairs <- expand.grid(i = 1:nrow(df), j = 1:nrow(df)) %>%
  filter(i < j) %>%
  mutate(
    Feature1 = df$Variable[i],
    Feature2 = df$Variable[j],
    
    quad1 = df$quadrant[i],
    quad2 = df$quadrant[j],
    
    cos_sim = cosine_similarity(
      df$RDA1[i], df$RDA2[i],
      df$RDA1[j], df$RDA2[j]
    ),
    
    angle_deg = angle_deg_safe(cos_sim)
  ) %>%
  arrange(desc(cos_sim)) %>%
  left_join(Pub_Classy_05_clean,
            by = c("Feature1" = "Code"))


head(CM_pairs)


# Microbe–metabolite relationships in both directions
library(stringr)

CM_strongcosine <- CM_pairs %>%
  filter(
    quad1 == quad2,
    (grepl("^(p__|g__)", Feature1) & grepl("^Met", Feature2)) |
      (grepl("^Met", Feature1) & grepl("^(p__|g__)", Feature2)),
    cos_sim > 0.9
  ) %>%
  arrange(desc(cos_sim)) 

head(CM_strongcosine)

Export(CM_pairs, "CM_pairs.txt")
Export(CM_strongcosine, "CM_strongcosine.txt")

# Interpretation:
# •	cos_sim → directional correlation
# •	angle_deg → angular separation

key_CM <- top50_scaling2_bac_Metabo %>% 
  filter(abs(RDA1) >= 0.5) %>% 
  arrange(RDA1) %>%
  mutate(Group_association = case_when(
    RDA1 < 0 & RDA2 < 0 ~ "Akha-CM",
    RDA1 > 0 & RDA2 > 0 ~ "Lahu-CM",
    RDA1 > 0 & RDA2 < 0 ~ "Lahu-CM",
    RDA1 < 0 & RDA2 > 0 ~ "Khuen-CM"
  ))

Export(key_CM, "key_CM.txt")


key_CM_cossim <- CM_strongcosine %>%
  filter(Feature1 %in% key_CM$Variable,
         Feature2 %in% key_CM$Variable) %>%
  select(Feature1, Feature2, quad1, quad2, cos_sim, angle_deg, 
         Compound_Clean, Subclass) %>%
  filter(cos_sim >= 0.9)

Export(key_CM_cossim, "key_CM_cossim.txt")








# ==========================================================
# 4. RDA – CR Group 
# ==========================================================

# Count ASVs for CR, select significant taxa at phylum and genus levels
count_CR <- count_asv %>%
  filter(Area == "CR") %>%
  select(any_of(names(sig_taxo_CR)[13:ncol(sig_taxo_CR)])) %>%
  as.matrix()

# confirm sample order
identical(rownames(count_CR), rownames(metadata_CR))


# Apply CLR transformation to Count Phylum and Genus
CR_bac_clr <- clr(count_CR + 1)
head(CR_bac_clr)
rownames(CR_bac_clr)

# Check is there any na, NaN 
colSums(is.nan((met_scale_CR)))

met_scale_CR_clean <- met_scale_CR



dim(met_scale_CR_clean)
dim(CR_bac_clr)


# ✅ metabo_scaled = log-transformed + scaled
# ✅ bac count = add pseudocount (+1) + apply CLR transformation


#-------------------------------------
# RDA: Metabolites or Gut microbiota
#-------------------------------------
# Run RDA with  group as explanatory variable
# RDA model: multi-omics response explained by  group
rda_model <- rda(CR_bac_clr ~ Group, data = metadata_CR)

summary(rda_model)
# Constrained Proportion: variance of Y explained by X 
# Unconstrained Proportion: unexplained variance in Y 



#---------------------------------
# RDA: Metabolites+Gut microbiota
#---------------------------------
#-- RDA with multi-omics as response
# Associations between genus and metabolites across groups
# •	✅ Response matrix: both bac_clr and metabo_scaled
# •	✅ Explanatory variable: Metadata$Group 

# Combine CLR-genus and scaled metabolites
response_matrix <- cbind(met_scale_CR_clean, CR_bac_clr)
class(response_matrix)

dim(met_scale_CR_clean)
dim(CR_bac_clr)

# Run RDA with  group as explanatory variable
# RDA model: multi-omics response explained by  group
rda_model <- rda(response_matrix ~ Group, data = metadata_CR)

summary(rda_model)
# Constrained Proportion: variance of Y explained by X 
# Unconstrained Proportion: unexplained variance in Y 


# canonical coefficients
coef(rda_model)
RsquareAdj(rda_model)


# Test significance (PERMUTATION-BASED ANOVA)
# Test overall model
set.seed(123)
anova.cca(rda_model, permutations = 999)

# Test each constrained axis
set.seed(123)
anova.cca(rda_model, by = "axis", permutations = 999)

# Test terms (Group, metabolites, etc.)
set.seed(123)
anova.cca(rda_model, by = "terms", permutations = 999)


# Plot the RDA
plot(rda_model, scaling = 1)  


#🔬 Interpretation

# • Scaling = 1 = distances between samples matter
# •	Scaling = 2 = angles between variables (species) matter (better when you want to explore 
# relationships between features like genus or metabolites).

#This RDA will help you answer:
#•	✅ How much variance in combined genus + metabolite profiles is explained by  group?
#•	✅ Are there visible clusters or shifts in profiles across  categories?
#•	✅ Which genera and metabolites drive the separation?



## Redundancy Analysis (RDA) biplot with type 2 scaling 
library(ggplot2)
library(ggrepel)


# Creating clean RDA biplots using ggplot2
# Extract scores
rda_scores <- scores(rda_model, display = c("sites", "species"), scaling = 2)

# Site score 
sites_df <- as.data.frame(rda_scores$sites)
sites_df$Group <- metadata_CR$Group  # add grouping info for coloring

# Label only the top contributing species (genus, metabolites)
# Get species scores (with scaling 2)
# This declutters the plot by showing only the most influential features.
species_scores <- as.data.frame(rda_scores$species)

# Calculate vector lengths
arrow_lengths <- sqrt(rowSums(species_scores^2))

# Get top features
top_species <- names(sort(arrow_lengths, decreasing = TRUE))[1:50]

species_df <- as.data.frame(rda_scores$species)

# Keep only top species Scaling = 2
# Get the top rownames based on arrow_lengths
topX_species <- species_df[rownames(species_df) %in% top_species, ]


# Scale the species scores (e.g., multiplying by x to extend the arrows)
scale_factor <- 2
top_species_scaled <- topX_species * scale_factor




# RDA plot with ggplot2

# Create the new column with categories
topX_species_scaled <- top_species_scaled %>%
  mutate(
    # Ion type assignment based on rownames (posxx or negxx)
    Feature = case_when(
      grepl("^Met", rownames(top_species_scaled)) ~ "Metabolite",
      grepl("^p__", rownames(top_species_scaled)) ~ "Phylum",
      grepl("^g__", rownames(top_species_scaled)) ~ "Genus",
      TRUE ~ "feature"  # Default to "feature" 
    )
  )


# View the updated data frame
head(topX_species_scaled)

# create dataframe
top50_scaling2_bac_Metabo <- data.frame(Variable = rownames(topX_species_scaled),
                                         topX_species_scaled) 

top100_scaling2_bac_Metabo <- data.frame(Variable = rownames(topX_species_scaled),
                                         topX_species_scaled) 


library(car)
Export(top50_scaling2_bac_Metabo, "top50_scaling2_bac_Metabo_CR.txt")


# proportion of variance explained by the RDA axes
summary(rda_model)$concont$importance
rda_var <- summary(rda_model)$cont$importance["Proportion Explained", 1:2]


# RDA ggplot scaling = 2
# Define the col_key with the colors for different groups
col_key <-  c("Akha-CR" = "dodgerblue4",
              "Lahu-CR" = "coral1",
              "Lisu-CR" = "grey40")


# Plotting
sc2 <- ggplot() +
  # Points for sites (sites_df) colored based on 'expl.var' or your relevant grouping variable
  geom_point(data = sites_df, aes(RDA1, RDA2, color = Group), size = 2, alpha = 0.5) +
  
  # Use category to color arrows differently (Genus and Metabolite)
  geom_segment(data = top50_scaling2_bac_Metabo, aes(x = 0, y = 0, xend = RDA1, 
                                                      yend = RDA2, color = Feature), 
               arrow = arrow(angle = 20, length = unit(0.1, "cm"), type = "closed"),
               linetype = 1, size = 0.5) +
  
  # Text labels with the same color as arrows
  geom_text_repel(data = top50_scaling2_bac_Metabo, aes(RDA1, RDA2, 
                                                         label = row.names(top50_scaling2_bac_Metabo), 
                                                         color = Feature),
                  size = 2.5, force = 10, xlim = c(-9, 6), ylim = c(-8, 5),
                  point.padding = 0.5, segment.colour = "grey", segment.alpha = 0.4,
                  max.overlaps = Inf) +
  
  # Color scale for both points and arrows 
  scale_color_manual(values = c("Phylum" = "#0072B2",
                                "Genus" = "#009E73",
                                "Metabolite" = "#E69F00",
                                "Akha-CR" = "dodgerblue4",
                                "Lahu-CR" = "coral1",
                                "Lisu-CR" = "grey40")) +
  
  # Labels for the axes
  labs(
    x = paste0("RDA1 (", round(100*rda_var[1],2), "%)"),
    y = paste0("RDA2 (", round(100*rda_var[2],2), "%)")
  ) +
  
  # Horizontal and vertical lines at 0
  geom_hline(yintercept = 0, linetype = 3, size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 3, size = 0.5) +
  
  theme_bw() +
  
  # Title of the plot
  labs(title = "RDA Biplot (Scaling 2): Top 50 Contributing Features")+
  
  coord_cartesian(
    xlim = c(min(sites_df$RDA1, top50_scaling2_bac_Metabo$RDA1) - 2, 
             max(sites_df$RDA1, top50_scaling2_bac_Metabo$RDA1) + 2),
    ylim = c(min(sites_df$RDA2, top50_scaling2_bac_Metabo$RDA2) - 1, 
             max(sites_df$RDA2, top50_scaling2_bac_Metabo$RDA2) + 1),
    clip = "off"
  )

sc2



## Redundancy Analysis (RDA) biplot with type 1 scaling 

# Creating clean RDA biplots using ggplot2
# Extract scores
rda_scores <- scores(rda_model, display = c("sites", "species"), scaling = 1)

sites_df <- as.data.frame(rda_scores$sites)
sites_df$Group <- metadata_CR$Group  # add grouping info for coloring

# Label only the top contributing species (genus, metabolites)
# Get species scores (with scaling 2)
# This declutters the plot by showing only the most influential features.
species_scores <- as.data.frame(rda_scores$species)

# Calculate vector lengths
arrow_lengths <- sqrt(rowSums(species_scores^2))

# Get top features
top_species <- names(sort(arrow_lengths, decreasing = TRUE))[1:50]

species_df <- as.data.frame(rda_scores$species)

# Keep only top species Scaling = 2
# Get the top rownames based on arrow_lengths
topX_species <- species_df[rownames(species_df) %in% top_species, ]




# RDA plot with ggplot2
# Create the new column with categories
topX_species_scaled <- top_species_scaled %>%
  mutate(
    # Ion type assignment based on rownames (posxx or negxx)
    Feature = case_when(
      grepl("^Met", rownames(top_species_scaled)) ~ "Metabolite",
      grepl("^p__", rownames(top_species_scaled)) ~ "Phylum",
      grepl("^g__", rownames(top_species_scaled)) ~ "Genus",
      TRUE ~ "feature"  # Default to "feature" 
    )
  )

# View the updated data frame
head(topX_species_scaled)



# Extract the centroids 
centroids <- aggregate(cbind(RDA1, RDA2) ~ Group, 
                       data = data.frame(sites_df), FUN = mean)

# View the centroids
print(centroids)

# proportion of variance explained by the RDA axes
summary(rda_model)$concont$importance
rda_var <- summary(rda_model)$cont$importance["Proportion Explained", 1:2]


sc1 <- ggplot() +
  # Points for sites (sites_df) colored based on 'Group'
  geom_point(data = sites_df, aes(RDA1, RDA2, color = Group), size = 2, alpha = 0.5) +
  
  # Points for centroids (cent) colored based on 'Group'
  geom_point(data = centroids, aes(RDA1, RDA2, color = Group), size = 3, shape = 17) +
  
  # Color scale for points
  scale_color_manual(values = col_key) +
  
  
  geom_segment(data = topX_species_scaled, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               color = "grey60",
               arrow = arrow(angle = 20, length = unit(0.2, "cm"), type = "closed"),
               linetype = 1, size = 0.5, alpha = 0.1) +
  
  # Text labels with the same color as arrows
  geom_text_repel(data = topX_species_scaled, 
                  aes(RDA1, RDA2, label = row.names(topX_species_scaled)),
                  color = "grey60",  alpha = 0.5,  
                  size = 2.5, force = 10, xlim = c(-3, 3), ylim = c(-2.5, 2.5),
                  point.padding = 0.5, segment.colour = "grey80", segment.alpha = 0.1,
                  max.overlaps = Inf) +
  
  # Labels for the axes
  labs(
    x = paste0("RDA1 (", round(100*rda_var[1],2), "%)"),
    y = paste0("RDA2 (", round(100*rda_var[2],2), "%)")
  ) +
  
  # Horizontal and vertical lines at 0
  geom_hline(yintercept = 0, linetype = 3, size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 3, size = 0.5) +
  
  theme_bw() +
  
  # Title of the plot
  labs(title = "RDA Plot (Scaling 1): Sample Distribution and Group Centroids")+
  coord_cartesian(
    xlim = c(min(sites_df$RDA1, topX_species_scaled$RDA1) - 2, 
             max(sites_df$RDA1, topX_species_scaled$RDA1) + 2),
    ylim = c(min(sites_df$RDA2, topX_species_scaled$RDA2) - 1, 
             max(sites_df$RDA2, topX_species_scaled$RDA2) + 1),
    clip = "off"
  )


sc1


# Multiple plots
library(gridExtra)
grid.arrange(
  sc1, 
  sc2,
  nrow = 2
)


# Prepare the dataframe

# • cosineNP = [-x, y], RDA1 < 0 & RDA2 > 0
# • cosinePP = [x, y], RDA1 > 0 & RDA2 > 0
# • cosineNN = [-x, -y], RDA1 < 0 & RDA2 < 0
# • cosinePN = [x, -y], RDA1 > 0 & RDA2 < 0


library(dplyr)
df <- top50_scaling2_bac_Metabo %>%
  mutate(quadrant = case_when(
    RDA1 > 0 & RDA2 > 0 ~ "PP",
    RDA1 > 0 & RDA2 < 0 ~ "PN",
    RDA1 < 0 & RDA2 > 0 ~ "NP",
    RDA1 < 0 & RDA2 < 0 ~ "NN"
  ))


# Compute pairwise cosine similarity
library(tidyr)
CR_pairs <- expand.grid(i = 1:nrow(df), j = 1:nrow(df)) %>%
  filter(i < j) %>%
  mutate(
    Feature1 = df$Variable[i],
    Feature2 = df$Variable[j],
    
    quad1 = df$quadrant[i],
    quad2 = df$quadrant[j],
    
    cos_sim = cosine_similarity(
      df$RDA1[i], df$RDA2[i],
      df$RDA1[j], df$RDA2[j]
    ),
    
    angle_deg = angle_deg_safe(cos_sim)
  ) %>%
  arrange(desc(cos_sim)) %>%
  left_join(Pub_Classy_05_clean,
            by = c("Feature1" = "Code"))


head(CR_pairs)


# Microbe–metabolite relationships in both directions
library(stringr)

CR_strongcosine <- CR_pairs %>%
  filter(
    quad1 == quad2,
    (grepl("^(p__|g__)", Feature1) & grepl("^Met", Feature2)) |
      (grepl("^Met", Feature1) & grepl("^(p__|g__)", Feature2)),
    cos_sim > 0.9
  ) %>%
  arrange(desc(cos_sim)) 

head(CR_strongcosine)

Export(CR_pairs, "CR_pairs.txt")
Export(CR_strongcosine, "CR_strongcosine.txt")

# Interpretation:
# •	cos_sim → directional correlation
# •	angle_deg → angular separation

key_CR <- top50_scaling2_bac_Metabo %>% 
  filter(abs(RDA1) >= 0.5) %>% 
  arrange(RDA1) %>%
  mutate(Group_association = case_when(
    RDA1 < 0 & RDA2 < 0 ~ "Lahu-CR",
    RDA1 > 0 & RDA2 > 0 ~ "Akha-CR",
    RDA1 > 0 & RDA2 < 0 ~ "Lisu-CR",
    RDA1 < 0 & RDA2 > 0 ~ "Lahu-CR & AKha-CR"
  ))

Export(key_CR, "key_CR.txt")


key_CR_cossim <- CR_strongcosine %>%
  filter(Feature1 %in% key_CR$Variable,
         Feature2 %in% key_CR$Variable) %>%
  select(Feature1, Feature2, quad1, quad2, cos_sim, angle_deg, 
         Compound_Clean, Subclass) %>%
  filter(cos_sim >= 0.9)

Export(key_CR_cossim, "key_CR_cossim.txt")



CR_pairs %>% filter(angle_deg >= 170) %>% 
  select(Feature1, Feature2, quad1, quad2, cos_sim, angle_deg) %>% 
  arrange(Feature2)


