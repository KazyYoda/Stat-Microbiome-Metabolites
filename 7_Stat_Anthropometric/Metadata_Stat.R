############################################################
# Metadata Statistical Analysis
#
# Purpose:
#   Perform statistical analysis of demographic, anthropometric,
#   blood pressure, and sex-related metadata across ethnic,
#   geographic, and combined Group categories.
#
# Main analyses:
#   1. Data preparation and factor definition
#   2. Shapiro-Wilk tests for normality
#   3. Levene's tests for homogeneity of variance
#   4. Two-group comparisons using t-test or Wilcoxon test
#   5. Multi-group comparisons using Kruskal-Wallis or ANOVA
#   6. Post hoc pairwise comparisons using Dunn's test or Tukey HSD
#   7. Gender distribution analysis using Fisher's exact or chi-square test
#   8. Descriptive summary of continuous variables
############################################################


# ==========================================================
# 1. Set Working Directory
# ==========================================================

# Set the directory containing the metadata analysis files.
setwd("~/Documents/HillTribe_NGS/0_Stat_Metadata")



# ==========================================================
# 2. Load Required Packages
# ==========================================================

# Load packages required for data manipulation, statistical testing,
# and importing tab-delimited and Excel data.
library(dplyr)
library(car)
library(readxl)
library(readr)


# ==========================================================
# 3. Import and Prepare Metadata
# ==========================================================

# Import the metadata table.
# The input file is tab-delimited and contains participant-level
# demographic, anthropometric, geographic, and clinical information.
metadata <- read_delim("metadata.txt", delim = "\t", 
                       escape_double = FALSE, trim_ws = TRUE)


# Create a combined Group variable from ethnicity and geographic area.
# For example:
#   Akha + CM -> Akha-CM
#   Lahu + CR  -> Lahu-CR
#
# BMI is renamed from "BMl" to "BMI" for consistency.
metadata <- metadata %>%
  mutate(Group = paste(Ethnicity, Area, sep = "-")) %>%
  relocate(Group, .before = 2) %>%
  rename(BMI = BMl)


# Convert categorical metadata variables to factors.
# Explicit factor levels are defined to maintain a consistent
# ordering of groups throughout statistical analyses and plots.
metadata$Ethnicity <- factor(metadata$Ethnicity) 
metadata$Area <- factor(metadata$Area,
                        levels = c("CM", "CR")) 
metadata$Group <- factor(metadata$Group,
                        levels = c("Akha-CM", "Lahu-CM", 
                                   "Akha-CR", "Lahu-CR",
                                   "Khuen-CM", "Lisu-CR")) 
metadata$Gender <- factor(metadata$Gender,
                        levels = c("Male", "Female")) 




# ==========================================================
# 4. Test Normality of Continuous Variables
# ==========================================================

# Continuous variables evaluated for normality:
#   Age, Height, Weight, BMI, SBP, and DBP.
#
# Shapiro-Wilk tests are performed separately within each Group.
# Missing values are removed before testing.
#
# A Shapiro-Wilk test requires at least three observations, so
# groups with fewer than three non-missing observations are assigned
# NA for the corresponding p-value.
vars_to_test <- c("Age", "Height","Weight", "BMI", "SBP","DBP")

results_list <- list()

for (var in vars_to_test) {
  
  # Split the selected variable by Group and apply the Shapiro-Wilk
  # test independently to each group.
  test_result <- by(metadata[[var]], metadata$Group, function(x) {
    x_clean <- x[!is.na(x)]
    n <- length(x_clean)
    
    if (n >= 3) {
      res <- shapiro.test(x_clean)
      data.frame(n = n, p.value = res$p.value)
    } else {
      data.frame(n = n, p.value = NA)
    }
  })
  
  # Store the normality result for each Group.
  for (g in names(test_result)) { 
    results_list[[length(results_list) + 1]] <- data.frame(
      Variable = var,
      Group = g,
      N = test_result[[g]]$n,
      P_value = test_result[[g]]$p.value
    )
  }
}


# Combine all Shapiro-Wilk results into a single data frame.
normality_results <- do.call(rbind, results_list)
print(normality_results)




# ==========================================================
# 5. Test Homogeneity of Variance
# ==========================================================

# Levene's test is used to assess whether the variance of each
# continuous variable is comparable between selected two-group
# comparisons.
#
# Comparisons evaluated:
#   Akha-CM vs Akha-CR
#   Lahu-CM vs Lahu-CR
library(car)

pairs <- list(
  c("Akha-CM", "Akha-CR"),
  c("Lahu-CM", "Lahu-CR")
)

levene_results <- list()

for (var in vars_to_test) {
  for (pair in pairs) {
    
    data_sub <- metadata %>%
      dplyr::filter(Group %in% pair) %>%
      dplyr::filter(!is.na(.data[[var]]))
    
    if (length(unique(data_sub$Group)) == 2 &&
        nrow(data_sub) >= 4) {
      
      formula <- as.formula(paste(var, "~ Group"))
      test <- leveneTest(formula, data = data_sub)
      pval <- test[1, "Pr(>F)"]
      
    } else {
      pval <- NA
    }
    
    levene_results[[length(levene_results) + 1]] <- data.frame(
      Variable = var,
      Comparison = paste(pair, collapse = " vs "),
      P_value = pval
    )
  }
}


# Combine Levene's test results.
levene_df <- do.call(rbind, levene_results)
print(levene_df)



# ==========================================================
# 6. Test Homogeneity of Variance Within Each Geographic Area
# ==========================================================

# Assess variance homogeneity among groups within the CM area.
levene_CM_results <- setNames(
  lapply(vars_to_test, function(var) {
    car::leveneTest(as.formula(paste(var, "~ Group")),
                    data = metadata %>% dplyr::filter(Area == "CM"))
  }),
  vars_to_test
)

levene_CM <- do.call(rbind, levene_CM_results)
print(levene_CM)


# Assess variance homogeneity among groups within the CR area.
levene_CR_results <- setNames(
  lapply(vars_to_test, function(var) {
    car::leveneTest(as.formula(paste(var, "~ Group")),
                    data = metadata %>% dplyr::filter(Area == "CR"))
  }),
  vars_to_test
)

levene_CR <- do.call(rbind, levene_CR_results)
print(levene_CR)





# ==========================================================
# 7. Identify Variables Violating Statistical Assumptions
# ==========================================================

# Identify variable–Group combinations with evidence against
# normality (P < 0.05).
normality_results %>% 
  filter(P_value < 0.05) %>%
  pull(Group, Variable)


# Identify variable–comparison combinations with evidence of
# unequal variances (P < 0.05).
levene_df %>%
  filter(P_value < 0.05) %>%
  pull(Variable, Comparison)



# ==========================================================
# 8. Two-Group Comparison: Welch's t-test
# ==========================================================

# Example baseline comparison:
#   Lahu-CM vs Lahu-CR
#
# Welch's t-test is used because it does not require equal variances
# between the two groups.
metadata %>% 
  filter(Group %in% c("Lahu-CM", "Lahu-CR")) %>%
  t.test(Weight ~ Group, data = ., var.equal = FALSE)





# ==========================================================
# 9. Two-Group Comparisons: Student's t-test
# ==========================================================

# Perform equal-variance t-tests for selected variables within
# the Akha ethnic group:
#   Akha-CM vs Akha-CR
#
# These analyses assume approximately normally distributed data
# and homogeneous variances between the two groups.
metadata %>% 
  filter(Ethnicity == "Akha") %>%
  t.test(Age ~ Group, data = ., var.equal = TRUE)

metadata %>% 
  filter(Ethnicity == "Akha") %>%
  t.test(BMI ~ Group, data = ., var.equal = TRUE)

metadata %>% 
  filter(Ethnicity == "Akha") %>%
  t.test(SBP ~ Group, data = ., var.equal = TRUE)

metadata %>% 
  filter(Ethnicity == "Akha") %>%
  t.test(DBP ~ Group, data = ., var.equal = TRUE)


# Apply the equal-variance t-test systematically to all selected
# continuous variables for the predefined Akha and Lahu comparisons.
vars_to_test <- c("Age", "Height","Weight", "BMI", "SBP","DBP")

comparisons <- list(
  Akha   = c("Akha-CM", "Akha-CR"),
  Lahu = c("Lahu-CM", "Lahu-CR")
)

t_test_results <- list()

for (var in vars_to_test) {
  for (comp_name in names(comparisons)) {
    
    comp_groups <- comparisons[[comp_name]]
    
    df <- metadata %>%
      filter(Group %in% comp_groups) %>%
      select(Group, all_of(var)) %>%
      na.omit()
    
    # Skip the comparison if fewer than two groups remain after
    # removal of missing observations.
    if (length(unique(df$Group)) < 2) next
    
    # Set factor order so that Group1 and Group2 are consistently
    # defined across comparisons.
    df$Group <- factor(df$Group, levels = comp_groups)
    
    test <- t.test(as.formula(paste(var, "~ Group")),
                   data = df,
                   var.equal = TRUE)
    
    t_test_results[[length(t_test_results) + 1]] <- data.frame(
      Variable = var,
      Comparison = comp_name,
      Group1 = comp_groups[1],
      Group2 = comp_groups[2],
      P_value = test$p.value
    )
  }
}

t_test_df <- do.call(rbind, t_test_results)

print(t_test_df)

# Display Lahu comparisons excluding Weight.
t_test_df %>% filter(Comparison == "Lahu", 
                    !Variable == "Weight")


# ==========================================================
# 10. Two-Group Comparison: Wilcoxon Rank-Sum Test
# ==========================================================

# Use the Wilcoxon rank-sum test for selected variables that do not
# satisfy the assumptions required for the parametric t-test.
vars_to_test <- c("Height","Weight")

comparisons <- list(
  Akha   = c("Akha-CM", "Akha-CR")
)

wilcox_results <- list()

for (var in vars_to_test) {
  for (comp_name in names(comparisons)) {
    
    comp_groups <- comparisons[[comp_name]]
    
    df <- metadata %>%
      filter(Group %in% comp_groups) %>%
      select(Group, all_of(var)) %>%
      na.omit()
    
    # Skip the comparison if fewer than two groups remain after
    # removal of missing observations.
    if (length(unique(df$Group)) < 2) next
    
    # Set factor order for consistent group labeling.
    df$Group <- factor(df$Group, levels = comp_groups)
    
    test <- wilcox.test(as.formula(paste(var, "~ Group")),
                        data = df,
                        exact = FALSE)
    
    wilcox_results[[length(wilcox_results) + 1]] <- data.frame(
      Variable = var,
      Comparison = comp_name,
      Group1 = comp_groups[1],
      Group2 = comp_groups[2],
      P_value = test$p.value
    )
  }
}

wilcox_df <- do.call(rbind, wilcox_results)

print(wilcox_df)



# ==========================================================
# 11. Multi-Group Comparison: Kruskal-Wallis Test
# ==========================================================

# Perform Kruskal-Wallis tests for continuous variables across
# the three CM groups.
#
# The Kruskal-Wallis test is a non-parametric omnibus test used to
# assess whether the distributions differ among more than two groups.
kruskal.test(Age ~ Group, data = metadata %>% filter(Area == "CM"))
kruskal.test(Height ~ Group, data = metadata %>% filter(Area == "CM"))
kruskal.test(Weight ~ Group, data = metadata %>% filter(Area == "CM"))

library(FSA)

# Specify variables for post hoc Dunn's test.
var <- c("Age")
var <- c("Height")


# Perform Dunn's test with Benjamini-Hochberg adjustment for
# pairwise comparisons following the Kruskal-Wallis analysis.
dunn_results <- setNames(
  lapply(var, function(var) {
    FSA::dunnTest(
      as.formula(paste(var, "~ Group")),
      data = metadata %>% dplyr::filter(Area == "CM"),
      method = "bh"
    )$res
  }),
  var
)


# ==========================================================
# 12. Generate Compact-Letter Displays from Dunn's Test
# ==========================================================

# Convert significant/non-significant pairwise Dunn's test results
# into compact letter groupings.
#
# Group names contain internal hyphens (e.g., "Akha-CM"). The function
# temporarily removes these internal hyphens so that the comparison
# separator can be identified safely.
library(multcompView)

get_letters <- function(dunn_df) {
  
  # Extract adjusted p-values.
  pvals <- dunn_df$P.adj
  
  # Split pairwise comparison labels.
  comps <- strsplit(as.character(dunn_df$Comparison), " - ")
  
  # Clean group names and reconstruct comparison labels.
  comps_clean <- lapply(comps, function(x) {
    g1 <- gsub("-", "", x[1])
    g2 <- gsub("-", "", x[2])
    paste0(g1, "-", g2)
  })
  
  # Assign comparison names to adjusted p-values.
  names(pvals) <- unlist(comps_clean)
  
  # Generate compact-letter groupings.
  letters <- multcompView::multcompLetters(pvals)$Letters
  
  return(letters)
}


get_letters(dunn_results$Height)





# ==========================================================
# 13. Multi-Group Comparison: CR Area
# ==========================================================

# Perform the Kruskal-Wallis test for Height across CR groups,
# followed by Dunn's test with BH-adjusted pairwise P-values.
kruskal.test(Height ~ Group, data = metadata %>% filter(Area == "CR"))
dunnTest(Height ~ Group, data = metadata %>% filter(Area == "CR"), method = "bh")





# ==========================================================
# 14. Multi-Group Comparison: One-Way ANOVA
# ==========================================================

# Perform one-way ANOVA for selected continuous variables across
# the CM groups.
#
# ANOVA evaluates whether the mean values differ among groups under
# the assumptions of approximately normal residuals and homogeneous
# variances.
vars_to_test = c("BMI", "SBP", "DBP")

anova_results <- setNames(
  lapply(vars_to_test, function(var) {
    aov(as.formula(paste(var, "~ Group")),
        data = metadata %>% dplyr::filter(Area == "CM"))
  }),
  vars_to_test
)


# Extract the overall ANOVA P-value for each variable.
anova_pvalues <- sapply(anova_results, function(model) {
  summary(model)[[1]][["Pr(>F)"]][1]
})




# ==========================================================
# 15. One-Way ANOVA: CR Area
# ==========================================================

# Perform one-way ANOVA for all selected continuous metadata variables
# across the CR groups.
vars_to_test = c("Age", "Height", "Weight", "BMI", "SBP", "DBP")

anova_results <- setNames(
  lapply(vars_to_test, function(var) {
    aov(as.formula(paste(var, "~ Group")),
        data = metadata %>% dplyr::filter(Area == "CR"))
  }),
  vars_to_test
)


# Extract overall ANOVA P-values.
anova_pvalues <- sapply(anova_results, function(model) {
  summary(model)[[1]][["Pr(>F)"]][1]
})

anova_pvalues


# ==========================================================
# 16. Tukey HSD Post Hoc Analysis
# ==========================================================

# Perform Tukey's HSD test for pairwise comparisons following
# one-way ANOVA.
posthoc_results <- lapply(anova_results, function(model) {
  TukeyHSD(model)
})

# Inspect selected ANOVA and post hoc results.
anova_results$Height
anova_pvalues["Height"]
posthoc_results$Height
posthoc_results$Weight





# ==========================================================
# 17. Gender Distribution
# ==========================================================

# Summarize the number of participants of each gender within each Group.
metadata %>%
  group_by(Group) %>%
  count(Gender)


# Create subsets for different ethnicity and geographic comparisons.
gender_Akha <- metadata %>%
  filter(Ethnicity  == "Akha") 

gender_Lahu <- metadata %>%
  filter(Ethnicity  == "Lahu") 

gender_CM <- metadata %>%
  filter(Area  == "CM") 

gender_CR <- metadata %>%
  filter(Area  == "CR")


# Create a contingency table of Group by Gender for the CR area.
table_gender <- table(gender_CR$Group, gender_CR$Gender)
table_gender


# Apply Fisher's exact test for the selected two-group comparison.
gender <- fisher.test(table_gender[c(2,4), ])
gender$p.value


# Apply chi-square test for the selected comparison involving
# more than two groups.
gender <- chisq.test(table_gender[c(3,4,6), ])
gender$p.value




# ==========================================================
# 18. Descriptive Summary of Continuous Variables
# ==========================================================

# Generate descriptive statistics for demographic, anthropometric,
# and blood pressure variables.
#
# For each Group, calculate:
#   - Mean
#   - Standard deviation (SD)
#
# The results are formatted as "Mean ± SD" for reporting.
library(dplyr)
library(tidyr)
library(stringr)

vars_to_summarize <- c("Age", "BMI", "Height", "Weight", 
                       "SBP", "DBP")


# Calculate mean and standard deviation within each Group,
# then reshape the results into a publication-friendly table.
summary_table <- metadata %>%
  group_by(Group) %>%
  summarise(across(all_of(vars_to_summarize),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd = ~sd(.x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(-Group,
               names_to = c("Variable", ".value"),
               names_sep = "_") %>%
  mutate(Mean_SD = sprintf("%.2f ± %.2f", mean, sd)) %>%
  dplyr::select(Group, Variable, Mean_SD) %>%
  pivot_wider(names_from = Group, values_from = Mean_SD)


# Display the final descriptive summary table.
print(summary_table)
