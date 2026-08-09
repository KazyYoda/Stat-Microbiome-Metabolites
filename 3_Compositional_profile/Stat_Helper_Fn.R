# ------------------------------------------------
# Helper function: Wilcoxon Rank Sum Test
# ------------------------------------------------
# Purpose:
# Perform Wilcoxon rank sum tests to compare
# taxa diversity taxas between two locations
# within a specified group.

taxa_wilcox <- function(data, 
                        rank_name = "Phylum", 
                        ethnicity_group = "Akha", 
                        factor_col = "Area") {
  
  message("\nRunning Wilcoxon tests for: ", ethnicity_group, rank_name)
  
  # Correct subgroup filtering
  data <- data %>%
    dplyr::filter(Ethnicity == ethnicity_group)
  
  # Helper function
  wilcox_summary <- function(data, group_label) {
    
    grp <- droplevels(factor(data[[factor_col]]))
    group_levels <- unique(grp)
    
    if (length(group_levels) != 2)
      stop("Wilcoxon test requires exactly 2 levels in ", factor_col)
    
    tests <- lapply(data[, 13:ncol(data)], function(taxa) {
      
      group1 <- taxa[grp == group_levels[1]]
      group2 <- taxa[grp == group_levels[2]]
      
      wilcox.test(group1, group2, exact = FALSE)
    })
    
    data.frame(
      Comparison = rep(group_label, length(tests)),
      taxa = colnames(data[, 13:ncol(data)]),
      P_Value = sapply(tests, function(x) x$p.value)
    )
  }
  
  # Run comparison
  result_1 <- wilcox_summary(data, "CM vs CR")
  
  # Adjust p-values
  result_1$p_adj <- p.adjust(result_1$P_Value, method = "BH")
  
  # Export
  Export(result_1, paste0("wilcox_", ethnicity_group,"_", rank_name, ".txt"))
  
  return(result_1)
}






# ------------------------------------------------
# Helper function: Kruskal–Wallis test
# ------------------------------------------------
# Purpose:
# Test differences in alpha diversity taxas across
# multiple groups using a non-parataxa Kruskal–Wallis test.

stat_kruskal <- function(data, 
                         rank_name = "Phylum",
                         group_col = "Area",
                         group_value = "CM",
                         factor_col = "Group",
                         factor_levels = c("Akha-CM", "Lahu-CM", "Khuen-CM")
                         ) {
  
  message("\nRunning Kruskal-Wallis tests for: ", group_value, "_", rank_name)
  
  
  # ------------------------------
  # Subset data by group
  # ------------------------------
  # Correct subgroup filtering
  data <- data %>%
    mutate(Group = paste(Ethnicity, Area, sep = "-")) %>%
    relocate(Group, .before = 2)
    

  stat_1 <- data %>%
    filter(!!sym(group_col) == group_value)
  
  # Ensure factor column has defined level order
  stat_1[[factor_col]] <- factor(
    stat_1[[factor_col]],
    levels = factor_levels
  )
  
  # Select taxas (starting after metadata columns)
  names_taxa <- names(stat_1)[14:ncol(stat_1)]
  
  # Apply Kruskal-Wallis test
  kruskal_results <- lapply(names_taxa, function(taxa) {
    kruskal.test(as.formula(paste(taxa, "~", factor_col)), data = stat_1)
  })
  
  # Summarize p-values
  kruskal_summary <- data.frame(
    Taxa = names_taxa,
    P_Value = sapply(kruskal_results, function(x) x$p.value),
    Absence = sapply(names_taxa, function(taxa) {
      if (sum(stat_1[[taxa]]) == 0) "Yes" else "No"
    })
  )
  
  # Add significance column (vectorized)
  kruskal_summary$P_sig <- ifelse(
    kruskal_summary$P_Value < 0.05,
    "Yes",
    "No"
  )
  
  print(kruskal_summary)
  print(kruskal_summary %>% filter(P_sig == "Yes"))
  
  # Get significant taxa names
  sig_taxa <- kruskal_summary$Taxa[kruskal_summary$P_sig == "Yes"]
  
  # Get column numbers in stat_1
  sig_col_numbers <- match(sig_taxa, colnames(stat_1))
  
  # Show results
  sig_col =  data.frame(
    Taxa = sig_taxa,
    Column_Number = sig_col_numbers
  ) %>%
    filter(!is.na(Column_Number))
  
  return(sig_col)
  
  
  # Export result
  Export(kruskal_summary, paste0("kruskal_", group_value, "_", rank_name, ".txt"))
  
}


# Extract adjusted values function:
sig_dunnT <- function(dunnTest_result){
  
  library(purrr) # map function
  library(dplyr)
  ## Extracting Data in Lists
  dunn_x <- map(dunnTest_result, "res")
  
  ## combine data from the lists
  dunn_x <- do.call(rbind.data.frame, dunn_x)
  dunn_x <- data.frame(taxa = rownames(dunn_x), dunn_x)
  
  ## Filter P.adj < 0.05 
  sig_dunn_x <- dunn_x %>% 
    filter(P.adj < 0.05) %>% 
    mutate_if(is.numeric, round, 5)
  
  return(sig_dunn_x)
  
}

# Extract p.adj < 0.05
sig_dunnT(data)



#-------------Helper Function: Taxon name cleaning ----------------------
prep_tax_table <- function(df) {
  
  library(dplyr)
  library(stringr)
  
  # Clean column names: remove brackets and replace hyphen with underscore
  colnames(df) <- colnames(df) %>%
    str_replace_all("\\[|\\]", "") %>% 
    str_replace_all("-", "_")
  
  return(df)
}

