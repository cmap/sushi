# synergy related functions

# restructure l2fc data
# for all the combination treatments, add single agent l2fcs for each drug
restructure_l2fc = function(cps_l2fc, sig_cols, cell_line_cols,
                            singles_type = "trt_cp",
                            combos_type = "trt_combo",
                            l2fc_col = 'median_l2fc') {
  
  # Slice data frame into single agents and combos
  singles_df = cps_l2fc %>% dplyr::filter(pert_type == singles_type)
  combos_df = cps_l2fc %>% dplyr::filter(pert_type == combos_type) %>%
    dplyr::rename_with(~ paste0("combo_", .x), contains(l2fc_col))
  
  # Add single agents to pert1 of the combos
  join_cols = sig_cols[sig_cols != "pert_type" & !grepl('pert2_', sig_cols)]
  
  restructured_l2fc = combos_df %>% 
    dplyr::inner_join(singles_df %>% dplyr::rename_with(~ paste0("pert1_", .x), contains(l2fc_col)),
                      by = c(cell_line_cols, join_cols), suffix = c('', '.y')) %>%
    dplyr::relocate(all_of(c(cell_line_cols, sig_cols, 'LFC'))) %>% dplyr::select(!contains('.y'))
  
  # Add single agents to pert2 of the combos
  join_cols = sig_cols[!grepl('pert', sig_cols)]
  cross_join_cols = setNames(sig_cols[grepl('pert1_', sig_cols)], 
                             sig_cols[grepl('pert2_', sig_cols)])
  
  restructured_l2fc = restructured_l2fc %>% 
    dplyr::inner_join(singles_df %>% dplyr::rename_with(~ paste0("pert2_", .x), contains(l2fc_col)),
                      by = c(cell_line_cols, join_cols, cross_join_cols), suffix = c('', '.y')) %>%
    dplyr::relocate(all_of(c(cell_line_cols, sig_cols, 'LFC'))) %>% dplyr::select(!contains('.y'))
  
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
                           l2fc_root = 'median_l2fc',
                           cap_for_viability = 1) {
  # From l2fc calculate FCs
  synergy_scores = restructured_l2fc %>% 
    dplyr::mutate(pert1_capped_viab = pmin(2^get(paste0('pert1_', l2fc_root)), cap_for_viability),
                  pert2_capped_viab = pmin(2^get(paste0('pert2_', l2fc_root)), cap_for_viability),
                  combo_capped_viab = pmin(2^get(paste0('combo_', l2fc_root)), cap_for_viability),
                  bliss = pert1_capped_viab * pert2_capped_viab) %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(hsa = min(pert1_capped_viab, pert2_capped_viab)) %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(synergy = 
                    dplyr::case_when(combo_capped_viab > hsa ~ combo_capped_viab - hsa,
                                     combo_capped_viab < bliss ~ bliss - combo_capped_viab,
                                     .default = 0))
  
  return(synergy_scores)
}

# Sample DMSO l2fc values
sample_dmso_l2fc = function(dmso_l2fc, l2fc_col = 'l2fc',
                            grouping_cols,
                            n_samples) {
  LFC_sample <- list()
  for (i in 1:n_samples){
    LFC_sample[[i]] <-  dmso_l2fc %>%
      group_by(pick(all_of(grouping_cols))) %>%
      slice_sample(n = 3) %>%
      summarise(pert1_LFC = median(.[[l2fc_col]]), .groups = 'drop')
    if(i %% 50 == 0) { print(i) }
  }
  LFC_sample <- LFC_sample %>%
    dplyr::bind_rows()
  
  LFC_sampled = LFC_sample %>% 
    dplyr::group_by(pick(all_of(grouping_cols))) %>%
    dplyr::mutate(pert2_LFC = sample(pert1_LFC, n_samples, replace = TRUE),
                  combo_LFC = sample(pert1_LFC, n_samples, replace = TRUE)) %>%
    dplyr::ungroup()
  return(LFC_sampled)
}

# 
get_cdf_value = function(group_name, synergy, h5_file) {
  ecdf_obj = stats::ecdf(h5read(h5_file, group_name))
  out_value = ecdf_obj(synergy)
  return(out_value)
}
