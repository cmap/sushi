# synergy related functions

# restructure l2fc data
# for all the combination treatments, add single agent l2fcs for each drug
restructure_l2fc = function(cps_l2fc, sig_cols, cell_line_cols,
                            singles_type = "trt_cp",
                            combos_type = "trt_combo",
                            l2fc_col = "median_l2fc",
                            names_prefix = c("pert1", "pert2", "combo"),
                            names_sep = "_") {
  # New column names
  new_names = paste(names_prefix, l2fc_col, sep = names_sep)
  
  # Slice data frame into single agents and combos
  singles_df = cps_l2fc %>% dplyr::filter(pert_type == singles_type)
  combos_df = cps_l2fc %>% dplyr::filter(pert_type == combos_type) %>%
    dplyr::rename(setNames(l2fc_col, new_names[3]))
  # Change the old l2fc_col name into combo + old l2fc_col
  
  # Add single agents to pert1 of the combos dataframe
  join_cols = sig_cols[sig_cols != "pert_type" & !grepl(names_prefix[2], sig_cols)]
  
  restructured_l2fc = combos_df %>%
    dplyr::left_join(singles_df %>% dplyr::rename(setNames(l2fc_col, new_names[1])),
                     by = c(cell_line_cols, join_cols), 
                     suffix = c("", ".y")) %>%
    dplyr::relocate(any_of(c(cell_line_cols, sig_cols, new_names))) %>%
    dplyr::select(!contains(".y"))
  
  # Add single agents to pert2 of the combos dataframe
  join_cols = sig_cols[!grepl("pert", sig_cols)]
  cross_join_cols = setNames(sig_cols[grepl(names_prefix[1], sig_cols)],
                             sig_cols[grepl(names_prefix[2], sig_cols)])
  
  restructured_l2fc = restructured_l2fc %>% 
    dplyr::left_join(singles_df %>% dplyr::rename(setNames(l2fc_col, new_names[2])),
                     by = c(cell_line_cols, join_cols, cross_join_cols), suffix = c("", ".y")) %>%
    dplyr::relocate(any_of(c(cell_line_cols, sig_cols, new_names))) %>% 
    dplyr::select(!contains(".y"))

  return(restructured_l2fc)
}

# Restructure 2 - grouping single agents with combinations
restructure2_l2fc = function(cps_l2fc, sig_cols, cell_line_cols,
                             singles_type = "trt_cp",
                             combos_type = "trt_combo",
                             l2fc_col = 'median_l2fc') {
  # function testing - to be removed
  cps_l2fc = test_l2fc
  sig_cols = c('x_project_id', 'pert_type', 'pert1_iname', 'pert1_dose', 'pert2_iname', 'pert2_dose')
  cell_line_cols = c('culture', 'pool_id', 'ccle_name')
  singles_type = "trt_cp"; combos_type = "trt_combo"
  l2fc_col = 'LFC'
  # -
  
  # Slice data frame into single agents and combos
  singles_df = cps_l2fc %>% dplyr::filter(pert_type == singles_type)
  combos_df = cps_l2fc %>% dplyr::filter(pert_type == combos_type)
  
  # Validate singles df
  print(unique(singles_df$pert2_iname))
  
  # Expand single agents for pert1
  join_cols = sig_cols[sig_cols != "pert_type" & !grepl('pert2_', sig_cols)]
  
  singles_pert1_df = combos_df %>% dplyr::mutate(pert2_dose = 0) %>%
    dplyr::distinct(pick(all_of(c(cell_line_cols, sig_cols)))) %>%
    dplyr::left_join(singles_df, by = c(cell_line_cols, join_cols), suffix = c('', '.y')) %>%
    dplyr::relocate(all_of(c(cell_line_cols, sig_cols, l2fc_col))) %>% dplyr::select(!contains('.y'))
  
  # Expand single agents for pert2
  join_cols = sig_cols[!grepl('pert', sig_cols)]
  cross_join_cols = setNames(sig_cols[grepl('pert1_', sig_cols)], 
                             sig_cols[grepl('pert2_', sig_cols)])
  
  singles_pert2_df = combos_df %>% dplyr::mutate(pert1_dose = 0) %>%
    dplyr::distinct(pick(all_of(c(cell_line_cols, sig_cols)))) %>%
    dplyr::left_join(singles_df, by = c(cell_line_cols, join_cols, cross_join_cols), suffix = c('', '.y')) %>%
    dplyr::relocate(all_of(c(cell_line_cols, sig_cols, l2fc_col))) %>% dplyr::select(!contains('.y'))
  
  # Full data table
  full_table = rbind(singles_pert1_df,
                     singles_pert2_df, 
                     combos_df)
  return(full_table)
}

# Calculate synergy scores
# uses the first restructuring function
# l2fc_root is set to l2fc_col
calculate_synergy = function(restructured_l2fc,
                             l2fc_root = "median_l2fc",
                             viability_cap = 1,
                             names_prefix = c("pert1", "pert2", "combo"),
                             names_sep = "_") {
  
  # Convert data frame to data.table
  if(!data.table::is.data.table(restructured_l2fc)) {
    restructured_l2fc = as.data.table(restructured_l2fc)
  }
  
  # Create vector of names used to create new columns
  l2fc_names = paste(names_prefix, l2fc_root, sep = names_sep)
  viab_names = paste(names_prefix, "viab", sep = names_sep)
  
  # Convert l2fcs to viabilities
  restructured_l2fc[, (l2fc_names) := .(pmin(2^get(l2fc_names[1]), viability_cap),
                                        pmin(2^get(l2fc_names[2]), viability_cap),
                                        pmin(2^get(l2fc_names[3]), viability_cap))]
  data.table::setnames(restructured_l2fc, l2fc_names, viab_names)
  
  # Calculate synergies
  restructured_l2fc[, hsa := pmin(get(viab_names[1]), get(viab_names[2]))]
  restructured_l2fc[, bliss := get(viab_names[1]) * get(viab_names[2])]
  restructured_l2fc[, synergy := data.table::fifelse(get(viab_names[3]) > hsa, hsa - get(viab_names[3]),
                                 data.table::fifelse(get(viab_names[3]) < bliss, bliss - get(viab_names[3]), 0))]
  
  return(restructured_l2fc)
}

# Sample DMSO l2fc values
median_sample = function(x, n_samples,
                         size = 3, seed = 2,
                         replace = TRUE, prob = NULL) {
  
  # Check that there are enough entries to sample to n
  num_pick_combinations = base::choose(length(x), size)
  print(paste0(length(x), " choose ", size, " = ", num_pick_combinations))
  
  # Stop if the number to sample to is too low
  if(num_pick_combinations < n_samples) {
    stop("ERROR: Cannot sample up to n.")
  }
  
  # Set seed
  base::set.seed(seed)
  
  # Test non loop
  resampled_values = base::replicate(n_samples, 
                                     median(base::sample(x, size, replace = replace, prob = prob)))
  
  return(resampled_values)
}

# Get cdf to change into pvalue
get_cdf_value = function(group_name, synergy, h5_file) {
  ecdf_obj = stats::ecdf(h5read(h5_file, group_name))
  out_value = ecdf_obj(synergy)
  return(out_value)
}

# smooth pvalue
smooth_pvalue <- function(pval_naive, n_samples ){
  ## use Laplace smoothing on pvalue to not output a 0 pvalue.
  return ( (pval_naive * n_samples+1)/(n_samples+2) )
}
