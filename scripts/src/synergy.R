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

# bliss synergy function
# uses the first restructuring function
# l2fc_root is set to l2fc_col
calculate_bliss = function(restructured_l2fc,
                           sig_cols,
                           cell_line_cols,
                           l2fc_root = 'median_l2fc') {
  # function testing - to be removed
  restructured_l2fc = restructured_l2fc
  sig_cols = c('x_project_id', 'pert_type', 'pert1_iname', 'pert1_dose', 'pert2_iname', 'pert2_dose')
  cell_line_cols = c('culture', 'pool_id', 'ccle_name')
  l2fc_root= 'LFC_cb'
  # -
  
  # Calculate capped viability
  capped_viab = restructured_l2fc %>% 
    dplyr::mutate(pert1_capped_viab = pmin(2^get(paste0('pert1_', l2fc_root)), 1),
                  pert2_capped_viab = pmin(2^get(paste0('pert2_', l2fc_root)), 1),
                  Bliss.NullFC = pert1_capped_viab * pert2_capped_viab)
  
  # Bliss calcs axis 1 - on pert1, pert1_dose, pert2
  temp1 = capped_viab %>% 
    dplyr::group_by(pick(all_of(c(cell_line_cols, sig_cols[!grepl('pert2_[a-z]?dose', sig_cols)])))) %>%
    dplyr::filter(length(unique(pert2_dose)) > 3) %>% 
    dplyr::summarise(riem_AUC_med_combination = mean(pmin(2^get(paste0('combo_', l2fc_root)), 1)),
                     riem_AUC_med_bliss = mean(Bliss.NullFC),
                     n_pts = dplyr::n()) %>% dplyr::ungroup()
  
  # Bliss calcs axis 2 - on pert1, pert2, pert2_dose
  temp2 = capped_viab %>% 
    dplyr::group_by(pick(all_of(c(cell_line_cols, sig_cols[!grepl('pert1_[a-z]?dose', sig_cols)])))) %>%
    dplyr::filter(length(unique(pert1_dose)) > 3) %>% 
    dplyr::summarise(riem_AUC_med_combination = mean(pmin(2^get(paste0('combo_', l2fc_root)), 1)),
                     riem_AUC_med_bliss = mean(Bliss.NullFC),
                     n_pts = dplyr::n()) %>% dplyr::ungroup()
  
  AUC_bliss.df = dplyr::bind_rows(temp1,temp2)
  
  return(AUC_bliss.df)
}
