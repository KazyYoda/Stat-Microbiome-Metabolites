# ==========================================================
# Helper Functions for Alpha Diversity Analysis
# ==========================================================
#
# This section defines reusable functions for:
#   1. Wilcoxon rank-sum tests for within-ethnicity comparisons
#   2. Alpha diversity visualization
#   3. Kruskal–Wallis tests for between-group comparisons
#
# The functions are used throughout the alpha diversity analysis
# to ensure a consistent statistical and visualization workflow.


# ----------------------------------------------------------
# Helper Function 1: Wilcoxon Rank-Sum Test
# ----------------------------------------------------------
#
# Purpose:
#   Compare alpha diversity metrics between two levels of a
#   grouping variable within a specified ethnicity.
#
# Statistical method:
#   Wilcoxon rank-sum test (two-sided), with exact = FALSE.
#
# Multiple-testing correction:
#   P-values across the tested alpha diversity metrics are
#   adjusted using the Benjamini–Hochberg (BH) method.
#
# Input:
#   data            : Data frame containing metadata and
#                     alpha diversity measurements.
#   ethnicity_group : Ethnicity to be analyzed.
#   metric          : Descriptive label for the analysis.
#   factor_col      : Column defining the two comparison groups.
#
# Output:
#   A data frame containing the comparison, alpha diversity
#   metric, raw P-value, and BH-adjusted P-value.
#
# The results are also exported as a text file named according
# to the specified ethnicity.


alpha_wilcox <- function(data, 
                         ethnicity_group = "Akha", 
                         metric = "Alpha Diversity",
                         factor_col = "Area") {
  
  message("\nRunning Wilcoxon tests for: ", ethnicity_group, metric)
  
  # Restrict the analysis to the specified ethnicity.
  data <- data %>%
    dplyr::filter(Ethnicity == ethnicity_group)
  
  # Helper function to perform Wilcoxon tests across
  # all alpha diversity metrics.
  wilcox_summary <- function(data, group_label) {
    
    # Extract the grouping variable and remove unused factor levels.
    grp <- droplevels(data[[factor_col]])
    group_levels <- unique(grp)
    
    # Confirm that exactly two groups are available for
    # the Wilcoxon rank-sum comparison.
    if (length(group_levels) != 2)
      stop("Wilcoxon test requires exactly 2 levels in ", factor_col)
    
    # Perform a separate Wilcoxon rank-sum test for each
    # alpha diversity variable in columns 14 through the
    # final column of the input data frame.
    tests <- lapply(data[, 14:ncol(data)], function(alpha) {
      
      group1 <- alpha[grp == group_levels[1]]
      group2 <- alpha[grp == group_levels[2]]
      
      wilcox.test(group1, group2, exact = FALSE)
    })
    
    # Summarize the test results.
    data.frame(
      Comparison = rep(group_label, length(tests)),
      alpha = colnames(data[, 14:ncol(data)]),
      P_Value = sapply(tests, function(x) x$p.value)
    )
  }
  
  # Compare the two study areas within the specified ethnicity.
  result_1 <- wilcox_summary(data, "CM vs CR")
  
  # Adjust P-values across the tested alpha diversity metrics
  # using the Benjamini–Hochberg procedure.
  result_1$p_adj <- p.adjust(result_1$P_Value, method = "BH")
  
  # Export the results for the specified ethnicity.
  Export(result_1, paste0("wilcox_alpha_", ethnicity_group, ".txt"))
  
  # Return the results to the R environment.
  return(result_1)
}








# ----------------------------------------------------------
# Helper Function 2: Alpha Diversity Boxplot Visualization
# ----------------------------------------------------------
#
# Purpose:
#   Generate faceted boxplots for multiple alpha diversity
#   metrics, with individual sample observations overlaid.
#
# Features:
#   - Accepts multiple alpha diversity metrics
#   - Reshapes the data to long format
#   - Allows a user-defined grouping variable
#   - Supports custom group colors
#   - Displays individual sample points using jitter
#
# Input:
#   data        : Data frame containing grouping and alpha
#                 diversity variables.
#   group_col   : Name of the column used to define groups.
#   metrics     : Character vector specifying alpha diversity
#                 variables to visualize.
#   fill_colors : Optional named vector of group colors.
#
# Output:
#   A ggplot object containing the alpha diversity boxplots.


alpha_boxplot <- function(data,
                          group_col,
                          metrics = c("Observed_features", "Chao1", "Shannon"),
                          fill_colors = NULL) {
  
  # Load packages required for data reshaping, visualization,
  # and tidy evaluation.
  library(tidyr)
  library(ggplot2)
  library(dplyr)
  library(rlang)
  
  # ----------------------------------------
  # Reshape data to long format
  # ----------------------------------------
  
  # Convert the selected alpha diversity metrics from wide
  # format to long format for faceted visualization.
  alpha_long <- data %>%
    pivot_longer(
      cols = all_of(metrics),
      names_to = "Index",
      values_to = "Value"
    )
  
  # Set the display order of the alpha diversity metrics.
  alpha_long$Index <- factor(alpha_long$Index, levels = metrics)
  
  # Convert the selected grouping variable to a factor.
  alpha_long[[group_col]] <- factor(alpha_long[[group_col]])
  
  # ----------------------------------------
  # Define Color Palette
  # ----------------------------------------
  
  # Use the default ggplot hue palette when custom colors
  # are not provided.
  if (is.null(fill_colors)) {
    fill_colors <- scales::hue_pal()(length(levels(alpha_long[[group_col]])))
    names(fill_colors) <- levels(alpha_long[[group_col]])
  }
  
  # ----------------------------------------
  # Generate Alpha Diversity Plot
  # ----------------------------------------
  
  # Create boxplots with individual sample observations.
  # Each alpha diversity metric is displayed in a separate facet.
  p <- ggplot(alpha_long, 
              aes(x = !!sym(group_col), y = Value, fill = !!sym(group_col))) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
    facet_wrap(~ Index, scales = "free_y", nrow = 1) +
    labs(x = group_col,
         y = "Alpha Diversity Index",
         fill = group_col) +
    scale_fill_manual(values = fill_colors) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # Return the ggplot object so that it can be displayed,
  # further modified, or saved externally.
  return(p)
}


# !!sym() is used to convert a character string containing
# a column name into a symbol for use in ggplot2 aesthetics.
# This allows group_col to be supplied dynamically.




# ----------------------------------------------------------
# Helper Function 3: Kruskal–Wallis Test
# ----------------------------------------------------------
#
# Purpose:
#   Test for differences in alpha diversity across multiple
#   groups within a specified study area.
#
# Statistical method:
#   Kruskal–Wallis rank-sum test.
#
# Input:
#   data          : Data frame containing metadata and alpha
#                   diversity measurements.
#   group_col     : Column used to subset the data.
#   group_value   : Value used to define the subset.
#   factor_col    : Column defining the comparison groups.
#   factor_levels : Desired order of the comparison groups.
#
# Output:
#   A data frame containing the Kruskal–Wallis P-value for
#   each alpha diversity metric.
#
# The results are also exported as a text file according
# to the specified group value.


alpha_kruskal <- function(data, 
                          group_col = "Area",
                          group_value = "CR",
                          factor_col = "Group",
                          factor_levels = c("Akha-CR", "Lahu-CR", "Lisu-CR", "Khuen-CR")
                          ) {
  
  # ------------------------------
  # Subset data by group
  # ------------------------------
  
  # Restrict the analysis to the specified study area.
  alpha_1 <- data %>%
    filter(!!sym(group_col) == group_value)
  
  # Set the comparison factor to the specified group order.
  alpha_1[[factor_col]] <- factor(
    alpha_1[[factor_col]],
    levels = factor_levels
  )
  
  # ------------------------------
  # Alpha Diversity Metrics
  # ------------------------------
  
  # Define the alpha diversity metrics included in the analysis.
  alpha_metrics <- c("Observed_features", "Chao1", "Shannon")
  
  # ------------------------------
  # Statistical Testing
  # ------------------------------
  
  # Apply a separate Kruskal–Wallis test to each alpha
  # diversity metric.
  kruskal_results <- lapply(alpha_metrics, function(metric) {
    kruskal.test(
      as.formula(paste(metric, "~", factor_col)),
      data = alpha_1
    )
  })
  
  # ------------------------------
  # Result Summary
  # ------------------------------
  
  # Summarize the Kruskal–Wallis P-values for all metrics.
  kruskal_summary <- data.frame(
    Metric = alpha_metrics,
    P_Value = sapply(kruskal_results, function(x) x$p.value)
  )
  
  # Display the results in the R console.
  print(kruskal_summary)
  
  # Export the results for reproducibility and record keeping.
  Export(kruskal_summary, paste0("kruskal_", group_value, ".txt"))
}
