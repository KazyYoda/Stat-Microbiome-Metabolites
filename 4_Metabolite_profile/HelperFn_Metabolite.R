##################################################################
# Metabolomics Analysis: Helper Functions and Statistical Tests
#################################################################
#
# Purpose:
# Define reusable functions for metabolomics data preparation,
# visualization, assumption checking, and statistical analysis.
#
# Main functions included:
#   1. Reorder sample metadata to match metabolite matrices
#   2. Calculate group-level log2 fold changes
#   3. Summarize log2 fold changes by group
#   4. Visualize metabolite log2 fold changes
#   5. Visualize metabolite annotations using chemical classes
#   6. Generate metabolite heatmaps
#   7. Assess normality and homogeneity of variance
#   8. Perform Wilcoxon rank-sum tests
#   9. Perform Kruskal-Wallis tests
#  10. Perform Dunn's post-hoc tests
#
# These functions are designed to be reused across different
# metabolite datasets, groups, and analytical comparisons.
############################################################



# ============================================================
# Helper Function: Reorder Metadata to Match Sample Order
# ============================================================
#
# Purpose:
# Reorder the metadata table so that its sample order exactly
# matches the order of samples in a metabolite abundance matrix.
#
# Inputs:
#   metadata     : sample metadata table
#   sample_vector: vector containing sample IDs in the desired order
#   sample_col   : metadata column containing sample IDs
#
# Returns:
#   A metadata data frame reordered according to sample_vector,
#   with Sample_ID values assigned as row names.
#
# This function also verifies that all sample IDs in the
# metabolite dataset are present in the metadata.

match_metadata_order <- function(metadata, sample_vector, sample_col = "Sample_ID") {
  # Ensure the sample_vector and metadata contain the same sample IDs
  if (!all(sample_vector %in% metadata[[sample_col]])) {
    stop("Some sample IDs in sample_vector are not found in metadata.")
  }
  
  # Reorder metadata according to the sample order in sample_vector
  metadata_ordered <- metadata %>%
    slice(match(sample_vector, .data[[sample_col]])) %>% # Reorders rows of metadata using positions in metadata.
    as.data.frame()
  
  # Set row names to sample IDs for convenient sample-level indexing
  rownames(metadata_ordered) <- metadata_ordered[[sample_col]]
  
  # Verify that the reordered metadata follows the requested sample order
  if (!identical(sample_vector, metadata_ordered[[sample_col]])) {
    warning("Sample order does not match even after reordering.")
  }
  
  return(metadata_ordered)
}





# =======================================================================
# Calculate Group-Level log2FC Relative to a Reference Group
# =======================================================================
#
# Purpose:
# Calculate the mean log2 fold change of each group relative to
# a specified reference/contrast group.
#
# IMPORTANT:
# The input metabolite table must already be log2-transformed.
#
# For log2-transformed data:
#
#   log2FC = mean(log2 abundance in comparison group) - mean(log2 abundance in reference group)
#            
#
# Inputs:
#   metabo_data: log2-transformed metabolite abundance table
#   metadata   : sample metadata containing the group variable
#   contrast   : reference group used for calculating log2FC
#   Group      : metadata column containing group labels
#
# Returns:
#   A matrix/data frame containing group-level log2FC values
#   for every metabolite.

log2FC_group_summary_relative_to_contrast <- function(metabo_data, # <- log2-transformed table
                                                      metadata, 
                                                      contrast = "Akha-CM", 
                                                      Group = "Group") {
  
  # Add group information to the metabolite abundance table
  metabo_data_with_group <- cbind(metadata[Group], metabo_data)
  colnames(metabo_data_with_group)[1] <- "Group"
  
  # Calculate mean log2 abundance for each metabolite within each group
  group_means <- aggregate(. ~ Group, data = metabo_data_with_group, 
                           FUN = function(x) mean(x, na.rm = TRUE))
  
  # Convert group labels to row names
  rownames(group_means) <- group_means$Group
  group_means$Group <- NULL
  
  # Extract the reference-group mean
  contrast_mean <- group_means[contrast, , drop = FALSE]
  
  # Calculate log2FC relative to the reference group.
  # Because the input is already log2-transformed, subtraction
  # directly gives the log2 fold change.
  ref <- as.numeric(group_means[contrast, colnames(group_means)]) # ensure same metabolite orders
  log2FC <- sweep(group_means, 2, ref, "-")
  
  return(log2FC)
}






# ============================================================
# Helper Function: Summarize log2FC by Group
# ============================================================
#
# Purpose:
# Convert a group-level log2FC matrix from wide format to long
# format and classify metabolites according to the direction
# of change.
#
# Classification:
#   log2FC > 0 : Upregulated
#   log2FC < 0 : Downregulated
#   log2FC = 0 : No change
#
# The reference/filter group is excluded from the returned
# summary because its log2FC values are zero by definition.
#
# Inputs:
#   df          : group-level log2FC table
#   filter_group: reference group to exclude
#
# Returns:
#   A long-format data frame containing Group, Metabolite,
#   log2FC, and Direction.

summarize_log2FC <- function(df,
                             filter_group = "Akha-CR") {
  
  library(dplyr)
  library(tidyr)
  library(tibble)
  
  # Convert row names containing group labels into a column
  df_long <- df %>%
    as.data.frame() %>%
    rownames_to_column("Group") %>%
    pivot_longer(
      cols = -Group,
      names_to = "Metabolite",
      values_to = "log2FC"
    )
  
  # Exclude the reference group and classify the direction
  # of metabolite change.
  summary <- df_long %>%
    filter(Group != filter_group) %>%
    mutate(Direction = case_when(
      log2FC > 0 ~ "Upregulated",
      log2FC < 0 ~ "Downregulated",
      TRUE ~ "No change"
    ))
  
  return(summary)
}





# ============================================================
# Helper Function: Plot Group-Level log2FC
# ============================================================
#
# Purpose:
# Generate faceted horizontal bar plots showing mean metabolite
# log2 fold changes for each comparison group.
#
# Only metabolites meeting the specified absolute log2FC cutoff
# are displayed.
#
# Inputs:
#   summary_df       : output from summarize_log2FC()
#   title            : plot title
#   xlab             : x-axis label
#   fill_colors      : colors for up/downregulated categories
#   Group            : group variable
#   strip_text_color : facet-strip text color
#   log2FC_cutoff    : minimum absolute log2FC for display
#   group_color      : colors used for group facet strips
#
# Returns:
#   A ggplot object.

plot_log2FC_bar <- function(summary_df, 
                            title, 
                            xlab, 
                            fill_colors, 
                            Group = "Group", 
                            strip_text_color = "white",
                            log2FC_cutoff = 2,
                            group_color = group_colors) {
  
  library(forcats)
  
  # Filter metabolites using the absolute log2FC threshold
  # and order metabolites within each group.
  summary_df <- summary_df %>%
    group_by(Group) %>%
    mutate(Metabolite = fct_reorder(Metabolite, log2FC)) %>%
    ungroup() %>%
    filter(abs(log2FC) >= log2FC_cutoff)
  
  
  # Define group colors
  group_color <- group_colors
  
  # Assign group-specific colors to facet strips
  strip_fill <- group_colors[unique(summary_df$Group)]
  
  library(ggplot2)
  library(ggh4x)
  
  # Create a horizontal bar plot of mean log2FC values.
  # Bars are colored according to the direction of change.
  p <- ggplot(summary_df, aes(x = Metabolite, y = log2FC, fill = Direction)) +
    geom_col(position = position_dodge(0.8), width = 0.7, aes(group = Group)) +
    scale_fill_manual(values = fill_colors) +
    facet_wrap2(~ Group, strip = strip_themed(
      background_x = elem_list_rect(fill = strip_fill, color = "black"),
      text_x = elem_list_text(face = "bold", color = strip_text_color)
    )) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 7, face = "bold"),
          strip.text = element_text(face = "bold"),
          title = element_text(size = 7)) +
    labs(x = xlab, y = "Mean log2 Fold Change", fill = "Direction", title = title) +
    coord_flip()
  
  return(p)
}



# ============================================================
# Helper Function: log2FC Bar Plot with Superclass Annotation
# ============================================================
#
# Purpose:
# Generate a log2FC bar plot in which each metabolite is
# additionally annotated according to its chemical superclass.
#
# A colored point is placed beside each bar to indicate the
# corresponding metabolite superclass.
#
# Inputs:
#   summary_df       : metabolite log2FC summary table
#   title            : plot title
#   xlab             : x-axis label
#   fill_colors      : colors for up/downregulated categories
#   log2FC_cutoff    : minimum absolute log2FC for display
#   strip_text_color : facet-strip text color
#   group_color      : colors used for group facet strips
#
# Returns:
#   A ggplot object.

plot_log2FC_bar_anno <- function(summary_df, 
                                 title, 
                                 xlab, 
                                 fill_colors, 
                                 log2FC_cutoff = 2,
                                 strip_text_color = "white",
                                 group_color) {
  
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(ggh4x)
  
  # Filter metabolites according to the absolute log2FC cutoff
  # and order metabolites within each comparison group.
  summary_df <- summary_df %>%
    filter(abs(log2FC) >= log2FC_cutoff) %>%
    group_by(Group) %>%
    mutate(Metabolite = fct_reorder(Metabolite, log2FC)) %>%
    ungroup()
  
  # Colorblind-friendly palette based on the Okabe-Ito palette
  # for chemical superclass annotation.
  superclass_colors <- c(
    "Lipids and lipid-like molecules" = "#E69F00",
    "Benzenoids" = "#56B4E9",
    "Organic oxygen compounds" = "#009E73",
    "Organoheterocyclic compounds" = "#F0E442",
    "Organic acids and derivatives" = "#0072B2",
    "Phenylpropanoids and polyketides" = "#D55E00",
    "Organic nitrogen compounds" = "#CC79A7",
    "Organosulfur compounds" = "#999999",
    "Not found" = "grey10"
  )
  
  # Define group colors
  group_color <- group_colors
  
  # Assign group-specific colors to facet strips
  strip_fill <- group_color[unique(summary_df$Group)]
  
  p <- ggplot(summary_df, aes(x = Metabolite, y = log2FC, fill = Direction)) +
    
    geom_col(width = 0.7) +
    
    # Add a point slightly outside each bar to indicate
    # the chemical superclass of the metabolite.
    geom_point(aes(
      y = ifelse(log2FC > 0,
                 log2FC + 0.3,
                 log2FC - 0.3),
      color = Superclass
    ),
    size = 1,
    shape = 16) +
    
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = superclass_colors, na.value = "black") +
    
    facet_wrap2(~ Group,
                strip = strip_themed(
                  background_x = elem_list_rect(fill = strip_fill, color = "black"),
                  text_x = elem_list_text(face = "bold", color = strip_text_color)
                )) +
    
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    
    theme_bw() +
    theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 8, face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "right") +
    
    labs(x = xlab,
         y = "Mean log2 Fold Change",
         fill = "Direction",
         color = "Superclass",
         title = title) +
    
    coord_flip()
  
  return(p)
}





# ============================================================
# Helper Function: log2FC Bar Plot with Subclass Annotation
# ============================================================
#
# Purpose:
# Generate a log2FC bar plot with metabolite chemical subclass
# annotation.
#
# A colored point beside each bar represents the corresponding
# metabolite subclass.
#
# The default absolute log2FC cutoff is 1.
#
# Inputs:
#   summary_df       : metabolite log2FC summary table
#   title            : plot title
#   xlab             : x-axis label
#   fill_colors      : colors for up/downregulated categories
#   log2FC_cutoff    : minimum absolute log2FC for display
#   strip_text_color : facet-strip text color
#   group_color      : colors used for group facet strips
#
# Returns:
#   A ggplot object.

plot_log2FC_bar_subclass <- function(summary_df, 
                                 title, 
                                 xlab, 
                                 fill_colors, 
                                 log2FC_cutoff = 1,
                                 strip_text_color = "white",
                                 group_color) {
  
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(ggh4x)
  
  # Filter metabolites according to the absolute log2FC cutoff
  # and order metabolites within each comparison group.
  summary_df <- summary_df %>%
    filter(abs(log2FC) >= log2FC_cutoff) %>%
    group_by(Group) %>%
    mutate(Metabolite = fct_reorder(Metabolite, log2FC)) %>%
    ungroup()
  
  # Colorblind-friendly palette for chemical subclass annotation.
  subclass_colors <- c(
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#999999",  # grey
    
    "#332288",  # dark indigo
    "#88CCEE",  # light cyan
    "#44AA99",  # teal
    "#117733",  # dark green (safe tone)
    "#AA4499",  # plum
    "#DDCC77",  # sand
    "#882255"   # wine
  )
  
  # Define group colors
  group_color <- group_colors
  
  # Assign group-specific colors to facet strips
  strip_fill <- group_color[unique(summary_df$Group)]
  
  p <- ggplot(summary_df, aes(x = Metabolite, y = log2FC, fill = Direction)) +
    
    geom_col(width = 0.7) +
    
    # Add a point slightly outside each bar to indicate
    # the chemical subclass of the metabolite.
    geom_point(aes(
      y = ifelse(log2FC > 0,
                 log2FC + 0.3,
                 log2FC - 0.3),
      color = Subclass
    ),
    size = 2,
    shape = 16) +
    
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = subclass_colors, na.value = "black") +
    
    facet_wrap2(~ Group,
                strip = strip_themed(
                  background_x = elem_list_rect(fill = strip_fill, color = "black"),
                  text_x = elem_list_text(face = "bold", color = strip_text_color)
                )) +
    
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    
    theme_bw() +
    theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          axis.title = element_text(size = 8, face = "bold"),
          strip.text = element_text(face = "bold"),
          title = element_text(size = 8),
          legend.position = "right") +
    
    labs(x = xlab,
         y = "Mean log2 Fold Change",
         fill = "Direction",
         color = "Subclass",
         title = title) +
    
    coord_flip()
  
  return(p)
}





# ============================================================
# Helper Function: Generate Metabolite Heatmap
# ============================================================
#
# Purpose:
# Generate a heatmap of metabolite abundance values with
# hierarchical clustering of metabolites and samples.
#
# Samples are annotated according to their experimental group,
# and samples are split into group-specific sections in the
# heatmap.
#
# The input matrix should contain metabolites as rows and
# samples as columns.
#
# Inputs:
#   mat          : metabolite abundance matrix with metabolites
#                  as rows and samples as columns
#   metadata     : sample metadata containing Sample_ID and
#                  Group columns
#   legend_title : title of the heatmap color legend
#   Group        : metadata column containing group labels
#   group_color  : named vector of colors corresponding to
#                  experimental groups
#
# Returns:
#   A ComplexHeatmap object.

plot_heatmap <- function(mat, 
                         metadata, 
                         legend_title = "Log2", 
                         Group = "Group",
                         group_color = group_colors) {
  
  
  # ----------------------------------------------------------
  # Ensure sample names in the matrix match the metadata
  # ----------------------------------------------------------
  
  # Assign metadata sample IDs to the columns of the matrix
  colnames(mat) <- metadata$Sample_ID
  
  # Print sample dimensions for verification
  cat("Matrix columns (samples):", ncol(mat), "\n")
  cat("Metadata samples       :", length(metadata$Group), "\n")
  
  
  # Check whether the number of samples matches
  if (ncol(mat) == length(metadata$Group)) {
    cat("✅ Sample counts matched.\n\n")
  } else {
    cat("❌ Sample counts do NOT match! Please check your data.\n\n")
  }
  
  
  # ----------------------------------------------------------
  # Create sample group annotation
  # ----------------------------------------------------------
  
  # Construct a sample annotation data frame
  sample_group <- data.frame(
    Group = metadata$Group
  )
  
  # Assign sample IDs as annotation row names
  rownames(sample_group) <- metadata$Sample_ID
  
  # Create the heatmap annotation
  ha <- HeatmapAnnotation(
    df = sample_group,
    col = list(Group = group_color)
  )
  
  
  # ----------------------------------------------------------
  # Define heatmap color scale
  # ----------------------------------------------------------
  
  # Use the minimum, mean, and maximum values of the matrix
  # to define a continuous three-color gradient.
  col_fun <- colorRamp2(
    c(min(mat, na.rm = TRUE),
      mean(mat, na.rm = TRUE),
      max(mat, na.rm = TRUE)),
    c("steelblue3", "grey96", "#FC766A")
  )
  
  
  # ----------------------------------------------------------
  # Generate heatmap
  # ----------------------------------------------------------
  
  Heatmap(
    mat,
    name = legend_title,
    
    # Hierarchical clustering of metabolites and samples
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    
    # Heatmap borders
    border = TRUE,
    border_gp = gpar(col = "grey"),
    rect_gp = gpar(col = "grey30", lwd = 0.5),
    
    # Continuous color scale
    col = col_fun,
    
    # Sample group annotation
    top_annotation = ha,
    
    # Separate samples according to experimental group
    column_split = metadata$Group,
    
    # Text formatting
    column_names_gp = gpar(fontsize = 6),
    row_names_gp = gpar(fontsize = 5),
    
    # Euclidean distance for hierarchical clustering
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    
    # Heatmap legend
    heatmap_legend_param = list(
      direction = "vertical"
    ),
    
    # Place metabolite names on the right
    row_names_side = "right"
  )
}





                           
# ===================================================
# Normality and Homogeneity of Variance Checking
# ===================================================
#
# Purpose:
# Assess distributional normality and homogeneity of variance
# for each metabolite before selecting the statistical test.
#
# Tests performed:
#   - Shapiro-Wilk test for normality
#   - Levene's test for homogeneity of variance
#
# The function identifies metabolite columns using the "Met"
# prefix and combines the test results with metabolite
# annotation information.

norm_var_check <- function(df,
                           metabo_descp, # metabo_descp: contains Compound codes and compound names
                           tag = "metabo05_Akha") {
  
  library(dplyr)
  library(car)
  
  # Identify metabolite abundance columns using the Met prefix
  metabolite_cols <- grep("^Met", names(df), value = TRUE)
  
  # -------------------------
  # Normality (Shapiro-Wilk test)
  # -------------------------
  
  shapiro_list <- lapply(metabolite_cols, function(met) {
    
    x <- df[[met]]
    
    # Skip the Shapiro-Wilk test when the metabolite has
    # two or fewer unique values.
    if (length(unique(x)) <= 2) {
      return(c(W = NA, p.value = NA))
    }
    
    test <- shapiro.test(x)
    c(W = test$statistic, p.value = test$p.value)
  })
  
  shapiro_df <- do.call(rbind, shapiro_list) %>%
    as.data.frame()
  
  shapiro_df$Metabolite <- metabolite_cols
  
  # -----------------------------------------
  # Homogeneity of variance (Levene's test)
  # -----------------------------------------
  
  levene_pvals <- sapply(metabolite_cols, function(met) {
    
    formula <- as.formula(paste0("`", met, "` ~ Group"))

    # Try to run Levene’s test; if it produces an error, 
    # don’t stop the whole loop—just return NA for this metabolite.
    test <- tryCatch(
      leveneTest(formula, data = df),
      error = function(e) NA
    )
    
    if (is.data.frame(test)) {
      test$`Pr(>F)`[1]
    } else {
      NA 
    }
  })

  # If NA produced it means a group has fewer than 2 usable observations, 
  # missing/non-finite values, or all values are identical (including all 0s).
  
  levene_df <- data.frame(
    Metabolite = metabolite_cols,
    Var_p.value = round(levene_pvals, 6)
  )

  
  # ----------------------------------------------------------
  # Combine statistical results with metabolite descriptions
  # ----------------------------------------------------------
  
  NormVar <- shapiro_df %>%
    left_join(levene_df, by = "Metabolite") %>%
    left_join(metabo_descp, by = c("Metabolite" = "Code")) %>%
    mutate(across(where(is.numeric), ~ round(.x, 6)))
  
  return(NormVar)
}

# ---------------------------------------------------------------
# Example Usage: Normality and Homogeneity of Variance Checking
# ---------------------------------------------------------------
Akha_norm_var <- norm_var_check(df,
                                metabo_descp,
                                tag = "metabo05_Akha")





# =============================================================
# Statistical Analysis: Direct Cross-Sectional Comparisons
# =============================================================
#
# Purpose:
# Perform direct statistical comparisons of metabolite
# abundances between independent groups in cross-sectional
# analyses.
#
# The workflow uses:
#
#   Normality/variance assessment
#            ↓
#   Wilcoxon rank-sum test for non-normal metabolites
#
# For comparisons involving more than two groups:
#
#   Normality/variance assessment
#            ↓
#   Kruskal-Wallis test
#            ↓
#   Dunn's post-hoc test for significant metabolites
                           
##---------------------------------------------------------
## Wilcoxon Rank-Sum Test for Cross-Sectional Comparisons
##---------------------------------------------------------

#----------------------- Wilcoxon Cross-Sectional: Helper Function --------------------------

# Purpose:
# Perform Wilcoxon rank-sum tests for two-group comparisons
# and apply Benjamini-Hochberg correction across metabolites.

run_wilcoxon_cross <- function(metabo_data, 
                               metabo_norm_var, 
                               group_label = "Group",
                               comparison_label = "Akha-CM vs Akha-CR",
                               tag = "metabo_Akha") {
  
  message("\nRunning Wilcoxon tests (cross-sectional) for: ", tag)
  
  wilcox_summary <- function(data, group_label) {

    # Apply the test after metadata columns
    tests <- lapply(data[, 14:ncol(data)], function(metabo) {
      
      grp <- data[[group_label]]
      group_levels <- unique(grp)
      
      # Continue only when exactly two groups are present
      if (length(group_levels) != 2) return(NA)
      
      group1 <- metabo[grp == group_levels[1]]
      group2 <- metabo[grp == group_levels[2]]
      
      # Wilcoxon rank-sum test for independent groups
      wilcox.test(group1, group2, exact = FALSE, paired = FALSE)
    })
    
    data.frame(
      Comparison = rep(comparison_label, length(tests)),
      Metabolite = colnames(data[, 14:ncol(data)]),
      P_Value = sapply(tests, function(x) if (is.list(x)) x$p.value else NA)
    ) %>%
      # Correct for multiple testing across metabolites
      mutate(p_adj = p.adjust(P_Value, method = "BH"))
  }
  
  # Select metabolites identified as non-normal
  nonNorm <- metabo_norm_var %>%
    filter(p.value < 0.05)
  
  # Retain metadata columns and selected non-normal metabolites
  metabo_dat <- metabo_data %>%
    select(1:13, all_of(nonNorm$Metabolite))
  
  result_1 <- wilcox_summary(metabo_dat, group_label)
  
  # Export Wilcoxon test results
  Export(result_1, paste0("wilcox_", tag, ".txt"))
  
  return(result_1)
}




# ============================================================
# Kruskal-Wallis Test: Multi-Group Metabolite Comparisons
# ============================================================
#
# Purpose:
# Perform a Kruskal-Wallis test for each metabolite when
# comparing more than two independent groups.
#
# Metabolites identified as non-normal are selected based on
# the normality assessment.
#
# The function returns the Kruskal-Wallis p-value for each
# metabolite.

kruskal_metabo <- function(metabo_norm_var,
                           metabo_data,
                           metadata, 
                           tag = "metabo_CM") {
  
  # Ensure input matrix is treated as a data frame
  metabo_data <- as.data.frame(metabo_data)
  
  # Step 1: Select metabolites that fail the normality criterion
  nonNorm <- metabo_norm_var %>% 
    filter(p.value < 0.05)
  
  # Step 2: Extract selected metabolites from the full dataset
  metabo_kruskal <- metabo_data %>% 
    select(all_of(nonNorm$Metabolite))
  
  # Step 3: Perform Kruskal-Wallis test for each metabolite
  krus_results <- lapply(metabo_kruskal, function(x) kruskal.test(x ~ factor(Group), data = metadata))
  
  # Step 4: Extract p-values and format results
  krus_pvalue <- data.frame(
    Metabolite = names(metabo_kruskal),
    kruskal_p = sapply(krus_results, function(res) res$p.value)
  ) %>%
    mutate(kruskal_p = round(kruskal_p, 6))
  
  # Step 5: Merge test results with metabolite descriptions
  krus_pvalue_descp <- krus_pvalue %>%
    left_join(nonNorm, by = "Metabolite") %>%
    select(Metabolite, kruskal_p)
  
  # Print the result to the console
  print(krus_pvalue_descp)
}

# ---------------------------------------------------------------
# Example Usage: Kruskal-Wallis Test
# ---------------------------------------------------------------
CM_kruskal_metabo <- kruskal_metabo(metabo_norm_var = CM_norm_var,
                                    metabo_data = log_matrix_metabo_05_CM,
                                    metadata = Met05_CM[1:13], 
                                    tag = "metabo_CM")

CR_kruskal_metabo <- kruskal_metabo(metabo_norm_var = CR_norm_var,
                                    metabo_data = log_matrix_metabo_05_CR,
                                    metadata = Met05_CR[1:13], 
                                    tag = "metabo_CR")



# ============================================================
# Dunn's Test: Post-Hoc Analysis Following Kruskal-Wallis
# ============================================================
#
# Purpose:
# Perform pairwise Dunn's tests for metabolites showing a
# significant global Kruskal-Wallis test.
#
# Benjamini-Hochberg correction is applied within each
# metabolite's set of pairwise comparisons.
#
# Only pairwise comparisons with adjusted P < 0.05 are retained.
#
# Metabolite annotation information is subsequently joined
# to the significant pairwise results.

dunnTest_metabo <- function(krus_pvalue_descp, 
                            metabo_data, 
                            metadata, 
                            metabo_descp, 
                            tag = "metabo_CM") {
  
  # Ensure input matrix is treated as a data frame
  metabo_data <- as.data.frame(metabo_data)
  
  
  # Load required packages
  library(FSA)
  library(dplyr)
  
  # 1. Select metabolites with significant global
  # Kruskal-Wallis p-values
  krus_sig <- krus_pvalue_descp %>% 
    filter(kruskal_p < 0.05)
  
  # 2. Run Dunn's test for each metabolite
  dunn_test <- setNames(
    lapply(krus_sig$Metabolite, function(x) {
      dunn_res <- dunnTest(
        metabo_data[[x]] ~ factor(Group),
        data = metadata,
        method = "bh"  # Benjamini-Hochberg adjustment
      )
      dunn_df <- as.data.frame(dunn_res$res) 
      return(dunn_df)
    }),
    krus_sig$Metabolite 
  )
  
  # 3. Extract pairwise comparisons with adjusted P < 0.05
  dunn_sig <- do.call(rbind, lapply(names(dunn_test), function(metabolite) {
    df <- dunn_test[[metabolite]]
    df$Metabolite <- metabolite
    df <- df[df$P.adj < 0.05, ]
    if (nrow(df) > 0) return(df)
    else return(NULL)
  }))
  
  # 4. Merge significant pairwise results with metabolite
  # annotation information, when available.
  if (!is.null(dunn_sig)) {
    dunn_sig <- data.frame(dunn_sig, row.names = NULL)
    
    # Ensure metabo_descp contains the required metabolite code
    if (!("Code" %in% colnames(metabo_descp))) {
      stop("metabo_descp must contain a 'Code' column.")
    }
    
    dunn_sig_descp <- dunn_sig %>%
      left_join(metabo_descp, by = c("Metabolite" = "Code"))
    
    return(dunn_sig_descp)
    
  } else {
    message("No significant pairwise differences found after adjustment.")
    return(NULL)
  }
}

# ---------------------------------------------------------------
# Example Usage: Dunn's Post-Hoc Test
# ---------------------------------------------------------------
CM_dunn <- dunnTest_metabo(krus_pvalue_descp, 
                           metabo_data, 
                           metadata, 
                           metabo_descp, 
                           tag = "metabo_CM")

CR_dunn <- dunnTest_metabo(krus_pvalue_descp, 
                           metabo_data, 
                           metadata, 
                           metabo_descp, 
                           tag = "metabo_CR")


# DunnT Work Flow.
Dunn results
     │
     ▼
Any P.adj < 0.05?
     │
   ┌─┴─┐
  YES  NO
   │    │
   ▼    ▼
dunn_sig  NULL
   │    │
   ▼    ▼
join    message()
annotation
   │    │
   ▼    ▼
return  return(NULL)
results
