############################################################
# Metabolite Profiles: PERMANOVA Analysis
############################################################
# This section evaluates differences in multivariate metabolite
# profiles using PERMANOVA (adonis2). The analysis focuses on
# Akha, Lahu, CM and CR samples and evaluates the effect of geographic area
# and ethnicity while considering potential covariates, including
# gender, hypertension status, BMI group, and age.
#
# Euclidean distance is calculated from standardized log2-
# transformed metabolite abundances, and statistical significance
# is assessed using 999 permutations.

setwd("~/Documents/HillTribe_NGS/5.Metabolites/Metabo_Stat/PERMANOVA_metabolite")
load("~/Documents/HillTribe_NGS/6.MFA/Metabolite_MFA/MFA_metabolite.RData")

# Load required packages
library(dplyr)
library(tidyr)
library(car)
library(readxl)
library(ggplot2)
library(vegan)


#=========================================================
# 1. PERMANOVA — Akha
#=========================================================

#---------------------------------------------------------
# 1.1 Prepare Akha Metadata and Metabolite Matrix
#---------------------------------------------------------

# Select the metadata variables required for the PERMANOVA
# analysis from the Akha MFA metadata.
metadata_Akha <- Akha_MFA[1:6]

# Select Akha samples and metabolite features from the
# log2-transformed metabolite dataset.
#
# Each metabolite is standardized to mean = 0 and SD = 1
# before calculating Euclidean distances. The resulting
# object is converted to a numeric matrix for use with vegan.
met_scale_Akha <- Met05_log2 %>%
  filter(Ethnicity == "Akha") %>%
  select(starts_with("Met")) %>%
  mutate(across(everything(), ~ scale(.x))) %>%
  as.matrix()
                             

# Confirm that sample ordering is identical between the
# metadata and metabolite matrix before performing PERMANOVA.
identical(rownames(metadata_Akha), rownames(met_scale_Akha))

# Remove metabolite features containing NaN values, which can
# arise during standardization when a feature has zero variance.
met_scale_Akha <- met_scale_Akha[, colSums(is.nan(met_scale_Akha)) == 0]


#---------------------------------------------------------
# 1.2 PERMANOVA: Area Effect with Covariate Adjustment
#---------------------------------------------------------

# Test whether multivariate metabolite profiles differ by
# geographic area (CM vs. CR) after accounting for gender,
# hypertension status, BMI group, and age.
#
# The analysis uses Euclidean distance and 999 permutations.
set.seed(123)
adonis2(met_scale_Akha ~ Area + Gender + Hypertension + BMI_group + Age,
        data = metadata_Akha,
        method = "euclidean",
        permutations = 999)

# Check the assumption of homogeneous multivariate dispersion
# among the groups using PERMDISP.
#
# betadisper evaluates whether differences in multivariate
# dispersion could contribute to an observed PERMANOVA result.

disp_Akha <- betadisper(dist(met_scale_Akha), metadata_Akha$Group)
anova(disp_Akha)


#---------------------------------------------------------
# 1.3 Sensitivity Analysis: Restricted Permutations
#---------------------------------------------------------
# The following models evaluate the Area effect while restricting
# permutations within levels of selected covariates.
#
# Restricting permutations within a variable tests the Area effect
# while preserving the structure of that variable during permutation.
#
# These models are used to assess whether the Area effect remains
# after accounting for potential confounding variables.


#---------------------------------------------------------
# Model 1: Full Covariate-Adjusted Model
#---------------------------------------------------------
# Test whether metabolite profiles differ by Area after accounting
# for gender, hypertension status, BMI group, and age.
set.seed(123)
adonis2(met_scale_Akha ~ Area + Gender + Hypertension + BMI_group + Age,
        data = metadata_Akha,
        method = "euclidean",
        permutations = 999)

# Estimate the marginal contribution of each factor after accounting
# for the other variables in the model.
set.seed(123)
adonis2(met_scale_Akha ~ Area + Gender + Hypertension + BMI_group + Age,
        data = metadata_Akha,
        method = "euclidean",
        permutations = 999,
        by = "margin")


#---------------------------------------------------------
# Model 2: Permutations Restricted Within Gender
#---------------------------------------------------------
# Test whether metabolite profiles differ by Area after accounting
# for hypertension status, BMI group, and age, while restricting
# permutations within gender strata.
set.seed(123)
adonis2(met_scale_Akha ~ Area + Hypertension + BMI_group + Age,
        data = metadata_Akha,
        method = "euclidean",
        permutations = 999,
        strata = metadata_Akha$Gender)

# Estimate the marginal contribution of each factor while retaining
# the gender-restricted permutation scheme.
set.seed(123)
adonis2(met_scale_Akha ~ Area + Hypertension + BMI_group + Age,
        data = metadata_Akha,
        method = "euclidean",
        permutations = 999,
        strata = metadata_Akha$Gender,
        by = "margin")


#---------------------------------------------------------
# Model 3: Permutations Restricted Within Hypertension Status
#---------------------------------------------------------
# Test whether metabolite profiles differ by Area after accounting
# for gender, BMI group, and age, while restricting permutations
# within hypertension-status strata.
set.seed(123)
adonis2(met_scale_Akha ~ Area + Gender + BMI_group + Age,
        data = metadata_Akha,
        method = "euclidean",
        permutations = 999,
        strata = metadata_Akha$Hypertension)

# Estimate the marginal contribution of each factor while retaining
# the hypertension-restricted permutation scheme.
set.seed(123)
adonis2(met_scale_Akha ~ Area + Gender + BMI_group + Age,
        data = metadata_Akha,
        method = "euclidean",
        permutations = 999,
        strata = metadata_Akha$Hypertension,
        by = "margin")


#---------------------------------------------------------
# Model 4: Permutations Restricted Within BMI Group
#---------------------------------------------------------
# Test whether metabolite profiles differ by Area after accounting
# for gender, age, and hypertension status, while restricting
# permutations within BMI-group strata.
set.seed(123)
adonis2(met_scale_Akha ~ Area + Gender + Age + Hypertension,
        data = metadata_Akha,
        method = "euclidean",
        permutations = 999,
        strata = metadata_Akha$BMI_group)

# Estimate the marginal contribution of each factor while retaining
# the BMI-group-restricted permutation scheme.
set.seed(123)
adonis2(met_scale_Akha ~ Area + Gender + Age + Hypertension,
        data = metadata_Akha,
        method = "euclidean",
        permutations = 999,
        strata = metadata_Akha$BMI_group,
        by = "margin")




#=========================================================
# 2. PERMANOVA — Lahu
#=========================================================

#---------------------------------------------------------
# 2.1 Prepare Lahu Metadata and Metabolite Matrix
#---------------------------------------------------------

# Select the metadata variables required for the PERMANOVA
# analysis from the Lahu MFA metadata.
metadata_Lahu <- Lahu_MFA[1:6]

# Select Lahu metabolite features and standardize each metabolite
# to mean = 0 and SD = 1 before calculating Euclidean distances.
met_scale_Lahu <- Lahu_MFA %>%
  select(starts_with("Met")) %>%
  mutate(across(everything(), ~ scale(.x))) %>%
  as.matrix()


# Confirm that sample ordering is identical between the metadata
# and metabolite matrix before performing PERMANOVA.
identical(rownames(metadata_Lahu), rownames(met_scale_Lahu))


#---------------------------------------------------------
# 2.2 Check Homogeneity of Multivariate Dispersion
#---------------------------------------------------------

# Assess whether the multivariate dispersion differs among groups.
# This provides a complementary check for the PERMANOVA analysis,
# because significant differences in dispersion can contribute
# to significant PERMANOVA results.
disp_Lahu <- betadisper(dist(met_scale_Lahu), metadata_Lahu$Group)
anova(disp_Lahu)



#===============================================================
# 2.3 PERMANOVA and Restricted-Permutation Sensitivity Analyses
#===============================================================

# The following models evaluate whether metabolite profiles differ
# between geographic areas (CM vs. CR) among Lahu samples.
#
# The models account for potential covariates including gender,
# hypertension status, BMI group, and age.
#
# In Models 2–4, permutations are restricted within levels of a
# selected covariate using the "strata" argument. This preserves
# the structure of the specified covariate during permutation and
# provides a sensitivity analysis for potential confounding.


#---------------------------------------------------------
# Model 1: Full Covariate-Adjusted Model
#---------------------------------------------------------

# Test whether metabolite profiles differ by Area after accounting
# for gender, hypertension status, BMI group, and age.
set.seed(123)
adonis2(met_scale_Lahu ~ Area + Gender + Hypertension + BMI_group + Age,
        data = metadata_Lahu,
        method = "euclidean",
        permutations = 999)

# Estimate the marginal contribution of each factor after accounting
# for all other variables included in the model.
set.seed(123)
adonis2(met_scale_Lahu ~ Area + Gender + Hypertension + BMI_group + Age,
        data = metadata_Lahu,
        method = "euclidean",
        permutations = 999,
        by = "margin")


#---------------------------------------------------------
# Model 2: Permutations Restricted Within Gender
#---------------------------------------------------------

# Test whether metabolite profiles differ by Area after accounting
# for hypertension status, BMI group, and age, while restricting
# permutations within gender strata.
set.seed(123)
adonis2(met_scale_Lahu ~ Area + Hypertension + BMI_group + Age,
        data = metadata_Lahu,
        method = "euclidean",
        permutations = 999,
        strata = metadata_Lahu$Gender)

# Estimate the marginal contribution of each factor while retaining
# the gender-restricted permutation scheme.
set.seed(123)
adonis2(met_scale_Lahu ~ Area + Hypertension + BMI_group + Age,
        data = metadata_Lahu,
        method = "euclidean",
        permutations = 999,
        strata = metadata_Lahu$Gender,
        by = "margin")


#---------------------------------------------------------
# Model 3: Permutations Restricted Within Hypertension Status
#---------------------------------------------------------

# Test whether metabolite profiles differ by Area after accounting
# for gender, BMI group, and age, while restricting permutations
# within hypertension-status strata.
set.seed(123)
adonis2(met_scale_Lahu ~ Area + Gender + BMI_group + Age,
        data = metadata_Lahu,
        method = "euclidean",
        permutations = 999,
        strata = metadata_Lahu$Hypertension)

# Estimate the marginal contribution of each factor while retaining
# the hypertension-restricted permutation scheme.
set.seed(123)
adonis2(met_scale_Lahu ~ Area + Gender + BMI_group + Age,
        data = metadata_Lahu,
        method = "euclidean",
        permutations = 999,
        strata = metadata_Lahu$Hypertension,
        by = "margin")


#---------------------------------------------------------
# Model 4: Permutations Restricted Within BMI Group
#---------------------------------------------------------

# Test whether metabolite profiles differ by Area after accounting
# for gender, age, and hypertension status, while restricting
# permutations within BMI-group strata.
set.seed(123)
adonis2(met_scale_Lahu ~ Area + Gender + Age + Hypertension,
        data = metadata_Lahu,
        method = "euclidean",
        permutations = 999,
        strata = metadata_Lahu$BMI_group)

# Estimate the marginal contribution of each factor while retaining
# the BMI-group-restricted permutation scheme.
set.seed(123)
adonis2(met_scale_Lahu ~ Area + Gender + Age + Hypertension,
        data = metadata_Lahu,
        method = "euclidean",
        permutations = 999,
        strata = metadata_Lahu$BMI_group,
        by = "margin")






#=========================================================
# 3. PERMANOVA — CM
#=========================================================

#---------------------------------------------------------
# 3.1 Prepare CM Metadata and Metabolite Matrix
#---------------------------------------------------------

# Select the metadata variables required for the PERMANOVA
# analysis from the CM MFA metadata.
metadata_CM <- CM_MFA[1:6]

# Select CM metabolite features and standardize each metabolite
# to mean = 0 and SD = 1 before calculating Euclidean distances.
met_scale_CM <- CM_MFA %>%
  select(starts_with("Met")) %>%
  mutate(across(everything(), ~ scale(.x))) %>%
  as.matrix()


# Confirm that sample ordering is identical between the metadata
# and metabolite matrix before performing PERMANOVA.
identical(rownames(metadata_CM), rownames(met_scale_CM))


#---------------------------------------------------------
# 3.2 Check Homogeneity of Multivariate Dispersion
#---------------------------------------------------------

# Assess whether multivariate dispersion differs among CM groups.
# This provides a complementary check for the PERMANOVA analysis,
# because differences in within-group dispersion can contribute
# to significant PERMANOVA results.
disp_CM <- betadisper(dist(met_scale_CM), metadata_CM$Group)
anova(disp_CM)



#===============================================================
# 3.3 PERMANOVA and Restricted-Permutation Sensitivity Analyses
#==============================================================

# The following models evaluate whether metabolite profiles differ
# among ethnic groups within the CM study area.
#
# Ethnicity is the primary factor of interest because the analysis
# is restricted to CM samples.
#
# Potential covariates include gender, hypertension status,
# BMI group, and age.
#
# Models 2–4 additionally restrict permutations within levels of
# selected covariates using the "strata" argument. These models
# provide sensitivity analyses by preserving the structure of the
# specified covariate during permutation.


#---------------------------------------------------------
# Model 1: Full Covariate-Adjusted Model
#---------------------------------------------------------

# Test whether metabolite profiles differ by ethnicity after
# accounting for gender, hypertension status, BMI group, and age.
set.seed(123)
adonis2(met_scale_CM ~ Ethnicity + Gender + Hypertension + BMI_group + Age,
        data = metadata_CM,
        method = "euclidean",
        permutations = 999)

# Estimate the marginal contribution of each factor after accounting
# for all other variables included in the model.
set.seed(123)
adonis2(met_scale_CM ~ Ethnicity + Gender + Hypertension + BMI_group + Age,
        data = metadata_CM,
        method = "euclidean",
        permutations = 999,
        by = "margin")


#---------------------------------------------------------
# Model 2: Permutations Restricted Within Gender
#---------------------------------------------------------

# Test group differences while accounting for hypertension status,
# BMI group, and age, with permutations restricted within gender
# strata.
set.seed(123)
adonis2(met_scale_CM ~ Group + Hypertension + BMI_group + Age,
        data = metadata_CM,
        method = "euclidean",
        permutations = 999,
        strata = metadata_CM$Gender)

# Estimate the marginal contribution of each factor while retaining
# the gender-restricted permutation scheme.
set.seed(123)
adonis2(met_scale_CM ~ Ethnicity + Hypertension + BMI_group + Age,
        data = metadata_CM,
        method = "euclidean",
        permutations = 999,
        strata = metadata_CM$Gender,
        by = "margin")


#---------------------------------------------------------
# Model 3: Permutations Restricted Within Hypertension Status
#---------------------------------------------------------

# Test whether metabolite profiles differ by ethnicity after
# accounting for gender, BMI group, and age, with permutations
# restricted within hypertension-status strata.
set.seed(123)
adonis2(met_scale_CM ~ Ethnicity + Gender + BMI_group + Age,
        data = metadata_CM,
        method = "euclidean",
        permutations = 999,
        strata = metadata_CM$Hypertension)

# Estimate the marginal contribution of each factor while retaining
# the hypertension-restricted permutation scheme.
set.seed(123)
adonis2(met_scale_CM ~ Ethnicity + Gender + BMI_group + Age,
        data = metadata_CM,
        method = "euclidean",
        permutations = 999,
        strata = metadata_CM$Hypertension,
        by = "margin")


#---------------------------------------------------------
# Model 4: Permutations Restricted Within BMI Group
#---------------------------------------------------------

# Test whether metabolite profiles differ by ethnicity after
# accounting for gender, age, and hypertension status, with
# permutations restricted within BMI-group strata.
set.seed(123)
adonis2(met_scale_CM ~ Ethnicity + Gender + Age + Hypertension,
        data = metadata_CM,
        method = "euclidean",
        permutations = 999,
        strata = metadata_CM$BMI_group)

# Estimate the marginal contribution of each factor while retaining
# the BMI-group-restricted permutation scheme.
set.seed(123)
adonis2(met_scale_CM ~ Ethnicity + Gender + Age + Hypertension,
        data = metadata_CM,
        method = "euclidean",
        permutations = 999,
        strata = metadata_CM$BMI_group,
        by = "margin")





#=========================================================
# 4. PERMANOVA — CR
#=========================================================

#---------------------------------------------------------
# 4.1 Prepare CR Metadata and Metabolite Matrix
#---------------------------------------------------------

# Select the metadata variables required for the PERMANOVA
# analysis from the CR MFA metadata.
metadata_CR <- CR_MFA[1:6]

# Select CR metabolite features and standardize each metabolite
# to mean = 0 and SD = 1 before calculating Euclidean distances.
met_scale_CR <- CR_MFA %>%
  select(starts_with("Met")) %>%
  mutate(across(everything(), ~ scale(.x))) %>%
  as.matrix()


# Confirm that sample ordering is identical between the metadata
# and metabolite matrix before performing PERMANOVA.
identical(rownames(metadata_CR), rownames(met_scale_CR))


#---------------------------------------------------------
# 4.2 Check Homogeneity of Multivariate Dispersion
#---------------------------------------------------------

# Assess whether multivariate dispersion differs among CR groups.
# This provides a complementary check for the PERMANOVA analysis,
# because differences in within-group dispersion can contribute
# to significant PERMANOVA results.
disp_CR <- betadisper(dist(met_scale_CR), metadata_CR$Group)
anova(disp_CR)



#=========================================================
# 4.3 PERMANOVA and Restricted-Permutation Sensitivity Analyses
#=========================================================

# The following models evaluate whether metabolite profiles differ
# among ethnic groups within the CR study area.
#
# Ethnicity is the primary factor of interest because the analysis
# is restricted to CR samples.
#
# Potential covariates include gender, hypertension status,
# BMI group, and age.
#
# Models 2–4 additionally restrict permutations within levels of
# selected covariates using the "strata" argument. These models
# provide sensitivity analyses by preserving the structure of the
# specified covariate during permutation.


#---------------------------------------------------------
# Model 1: Full Covariate-Adjusted Model
#---------------------------------------------------------

# Test whether metabolite profiles differ by ethnicity after
# accounting for gender, hypertension status, BMI group, and age.
set.seed(123)
adonis2(met_scale_CR ~ Ethnicity + Gender + Hypertension + BMI_group + Age,
        data = metadata_CR,
        method = "euclidean",
        permutations = 999)

# Estimate the marginal contribution of each factor after accounting
# for all other variables included in the model.
set.seed(123)
adonis2(met_scale_CR ~ Ethnicity + Gender + Hypertension + BMI_group + Age,
        data = metadata_CR,
        method = "euclidean",
        permutations = 999,
        by = "margin")


#---------------------------------------------------------
# Model 2: Permutations Restricted Within Gender
#---------------------------------------------------------

# Test whether metabolite profiles differ by ethnicity after
# accounting for hypertension status, BMI group, and age, with
# permutations restricted within gender strata.
set.seed(123)
adonis2(met_scale_CR ~ Ethnicity + Hypertension + BMI_group + Age,
        data = metadata_CR,
        method = "euclidean",
        permutations = 999,
        strata = metadata_CR$Gender)

# Estimate the marginal contribution of each factor while retaining
# the gender-restricted permutation scheme.
set.seed(123)
adonis2(met_scale_CR ~ Ethnicity + Hypertension + BMI_group + Age,
        data = metadata_CR,
        method = "euclidean",
        permutations = 999,
        strata = metadata_CR$Gender,
        by = "margin")


#---------------------------------------------------------
# Model 3: Permutations Restricted Within Hypertension Status
#---------------------------------------------------------

# Test whether metabolite profiles differ by ethnicity after
# accounting for gender, BMI group, and age, with permutations
# restricted within hypertension-status strata.
set.seed(123)
adonis2(met_scale_CR ~ Ethnicity + Gender + BMI_group + Age,
        data = metadata_CR,
        method = "euclidean",
        permutations = 999,
        strata = metadata_CR$Hypertension)

# Estimate the marginal contribution of each factor while retaining
# the hypertension-restricted permutation scheme.
set.seed(123)
adonis2(met_scale_CR ~ Ethnicity + Gender + BMI_group + Age,
        data = metadata_CR,
        method = "euclidean",
        permutations = 999,
        strata = metadata_CR$Hypertension,
        by = "margin")


#---------------------------------------------------------
# Model 4: Permutations Restricted Within BMI Group
#---------------------------------------------------------

# Test whether metabolite profiles differ by ethnicity after
# accounting for gender, age, and hypertension status, with
# permutations restricted within BMI-group strata.
set.seed(123)
adonis2(met_scale_CR ~ Ethnicity + Gender + Age + Hypertension,
        data = metadata_CR,
        method = "euclidean",
        permutations = 999,
        strata = metadata_CR$BMI_group)

# Estimate the marginal contribution of each factor while retaining
# the BMI-group-restricted permutation scheme.
set.seed(123)
adonis2(met_scale_CR ~ Ethnicity + Gender + Age + Hypertension,
        data = metadata_CR,
        method = "euclidean",
        permutations = 999,
        strata = metadata_CR$BMI_group,
        by = "margin")







#=========================================================
# 5. PERMANOVA — Geographic and Ethnic Effects
#=========================================================

#---------------------------------------------------------
# 5.1 Prepare Metadata and Metabolite Matrix
#---------------------------------------------------------

# Prepare metadata containing the primary grouping variables
# (geographic area and ethnicity) and potential covariates.
#
# Categorical variables are explicitly converted to factors
# to ensure that they are treated as categorical predictors
# in the PERMANOVA models.
metadata <- Met05_log2 %>%
  select(Group, Gender, Hypertension, Ethnicity, Area, BMI_group, Age) %>%
  mutate(across(c(Gender, Group, Area, Ethnicity, 
                         BMI_group, Hypertension), factor))

# Select all metabolite features and standardize each metabolite
# to mean = 0 and SD = 1 before calculating Euclidean distances.
met_scale <- Met05_log2 %>%
  select(starts_with("Met")) %>%
  mutate(across(everything(), ~ scale(.x))) %>%
  as.matrix()


# Confirm that sample ordering is identical between the metadata
# and metabolite matrix before performing PERMANOVA.
identical(rownames(metadata), rownames(met_scale))


#---------------------------------------------------------
# 5.2 Check Homogeneity of Multivariate Dispersion
#---------------------------------------------------------

# Assess whether multivariate dispersion differs among the
# Group categories.
#
# This provides a complementary check for the PERMANOVA analysis,
# because differences in within-group dispersion can contribute
# to significant PERMANOVA results.
disp <- betadisper(dist(met_scale), metadata$Group)
anova(disp)



#=========================================================
# 5.3 PERMANOVA: Combined Geographic and Ethnic Effects
#=========================================================

#---------------------------------------------------------
# Model 1: Full Covariate-Adjusted Model
#---------------------------------------------------------

# Test whether metabolite profiles differ by geographic Area
# and Ethnicity after accounting for gender, hypertension status,
# BMI group, and age.
#
# Area and Ethnicity represent the primary geographic and
# demographic factors of interest, respectively.
set.seed(123)
adonis2(met_scale ~ Area + Ethnicity + Gender + Hypertension + BMI_group + Age,
        data = metadata,
        method = "euclidean",
        permutations = 999)

# Estimate the marginal contribution of each factor after accounting
# for all other variables included in the model.
#
# "by = margin" evaluates each term while conditioning on the
# remaining terms in the model.
set.seed(123)
adonis2(met_scale ~ Area + Ethnicity + Gender + Hypertension + BMI_group + Age,
        data = metadata,
        method = "euclidean",
        permutations = 999,
        by = "margin")


# Repeat the marginal analysis with Ethnicity listed before Area.
# This provides the marginal effect of each factor while accounting
# for the other variables in the model.
set.seed(123)
adonis2(met_scale ~ Ethnicity + Area + Gender + Hypertension + BMI_group + Age,
        data = metadata,
        method = "euclidean",
        permutations = 999,
        by = "margin")
