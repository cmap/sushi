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
  # Create vector of names used to create new columns
  l2fc_names = paste(names_prefix, l2fc_root, sep = names_sep)
  viab_names = paste(names_prefix, "viab", sep = names_sep)

  # From l2fc calculate FCs
  restructured_l2fc[[viab_names[1]]] = pmin(2^restructured_l2fc[[l2fc_names[1]]], viability_cap)
  restructured_l2fc[[viab_names[2]]] = pmin(2^restructured_l2fc[[l2fc_names[2]]], viability_cap)
  restructured_l2fc[[viab_names[3]]] = pmin(2^restructured_l2fc[[l2fc_names[3]]], viability_cap)
  
  # Use viab to calculate bliss, hsa, and synergy
  synergy_scores = restructured_l2fc %>%
    dplyr::mutate(bliss = .data[[viab_names[1]]] * .data[[viab_names[2]]]) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(hsa = min(.data[[viab_names[1]]], .data[[viab_names[2]]],
                            na.rm = FALSE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(synergy =
                    dplyr::case_when(.data[[viab_names[3]]] > hsa ~ .data[[viab_names[3]]] - hsa,
                                     bliss <= .data[[viab_names[1]]] & .data[[viab_names[1]]] <= hsa ~ 0,
                                     .data[[viab_names[3]]] < bliss ~ bliss - .data[[viab_names[3]]],
                                     .default = NA))
  
  return(synergy_scores)
}

# Sample DMSO l2fc values
sample_dmso_l2fc = function(dmso_l2fc, l2fc_col = 'l2fc',
                            grouping_cols,
                            n_samples,
                            names_prefix = c("pert1", "pert2", "combo"),
                            names_sep = "_",
                            seed = 2) {
  # Check that there are enough entries to sample to n
  min_num_reps = dmso_l2fc %>%
    dplyr::group_by(pick(all_of(grouping_cols))) %>%
    dplyr::summarize(count= dplyr::n(), .groups = "drop") %>%
    dplyr::pull(count) %>% base::min()
  # Calculate combination with the minimum
  num_pick_combinations = base::choose(min_num_reps, 3)
  print(paste0(min_num_reps, " choose ", 3, " = ", num_pick_combinations))
  
  # Stop if the number to sample to is too low
  if(num_pick_combinations < n_samples) {
    stop("ERROR: Cannot sample up to n.")
  }
    
  # Set seed
  base::set.seed(seed)
  # Create vector of names used to create new columns
  l2fc_names = paste(names_prefix, l2fc_col, sep = names_sep)
  # Create an empty list to store results
  sample_output = list()
  
  # Loop to sample each cell line on a plate
  for(i in 1:n_samples) {
    sample_output[[i]] = dmso_l2fc %>%
      dplyr::group_by(pick(all_of(grouping_cols))) %>%
      dplyr::slice_sample(n = 3) %>%
      dplyr::summarise(!!l2fc_names[1] := median(.data[[l2fc_col]]), .groups = 'drop')
    
    if(i %% 100 == 0) { print(i) }
  }
  
  # Combine sampled outputs and sample other two l2fc columns
  sample_l2fc = dplyr::bind_rows(sample_output) %>%
    dplyr::group_by(pick(all_of(grouping_cols))) %>%
    dplyr::mutate(!!l2fc_names[2] := sample(.data[[l2fc_names[1]]], n_samples, replace = TRUE),
                  !!l2fc_names[3] := sample(.data[[l2fc_names[1]]], n_samples, replace = TRUE)) %>%
    dplyr::ungroup()
  
  return(sample_l2fc)
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
