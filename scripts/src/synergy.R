# Collection of functions used by the synergy module

#' restructure_l2fc
#'
#' In the context of a CPS, this function restructures a l2fc file by creating additional
#' columns containing the l2fc values for the single agents.
#'
#' This function splits the l2fc dataframe into combination entries and single agent entries.
#' The single agent entries are merged onto the combonation entries - once for each member of the
#' combination.
#'
#' @import tidyverse
#' @param cps_l2fc
#' @param sig_cols A vector of column names defining unique profiles.
#' @param cell_line_cols A vector of column names defining unique cell lines.
#' @param singles_type
#' @param combos_type
#' @param l2fc_col A string name of the column containing the l2fc values.
#' @param names_prefix A vector of three prefixes used to identify perturbation 1, perturbation 2, and the combination.
#' @param names_sep A string used to create new column names with `names_prefix` and `l2fc_col`.
#' @return A dataframe
restructure_l2fc = function(cps_l2fc, sig_cols, cell_line_cols,
                            singles_type = "trt_cp",
                            combos_type = "trt_combo",
                            l2fc_col = "median_l2fc",
                            names_prefix = c("pert1", "pert2", "combo"),
                            names_sep = "_") {
  # Create new column names using prefix, l2fc_col, and names_sep
  new_names = paste(names_prefix, l2fc_col, sep = names_sep)

  # Slice data frame into single agents and combos
  singles_df = cps_l2fc %>% dplyr::filter(pert_type == singles_type)
  combos_df = cps_l2fc %>% dplyr::filter(pert_type == combos_type) %>% dplyr::rename(setNames(l2fc_col, new_names[3]))
  # Change the old l2fc_col name into combo + old l2fc_col

  # Add single agents to pert1 of the combos dataframe
  join_cols = sig_cols[sig_cols != "pert_type" & !grepl(paste0(names_prefix[2], names_sep), sig_cols)]

  restructured_l2fc = combos_df %>%
    dplyr::left_join(singles_df %>% dplyr::rename(setNames(l2fc_col, new_names[1])),
                     by = c(cell_line_cols, join_cols),
                     suffix = c("", ".y")) %>%
    dplyr::relocate(tidyselect::any_of(c(cell_line_cols, sig_cols, new_names))) %>%
    dplyr::select(!tidyselect::contains(".y"))

  # Add single agents to pert2 of the combos dataframe
  join_cols = sig_cols[!grepl("pert", sig_cols)]
  cross_join_cols = setNames(sig_cols[grepl(paste0(names_prefix[1], names_sep), sig_cols)],
                             sig_cols[grepl(paste0(names_prefix[2], names_sep), sig_cols)])

  restructured_l2fc = restructured_l2fc %>%
    dplyr::left_join(singles_df %>% dplyr::rename(setNames(l2fc_col, new_names[2])),
                     by = c(cell_line_cols, join_cols, cross_join_cols), suffix = c("", ".y")) %>%
    dplyr::relocate(tidyselect::any_of(c(cell_line_cols, sig_cols, new_names))) %>%
    dplyr::select(!tidyselect::contains(".y"))

  return(restructured_l2fc)
}

# Restructure 2 - grouping single agents with combinations
restructure2_l2fc = function(cps_l2fc, sig_cols, cell_line_cols,
                             singles_type = "trt_cp",
                             combos_type = "trt_combo",
                             l2fc_col = "median_l2fc") {
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

#' Calculate synergy scores
#'
#' In the context of a CPS, this function calculates synergy scores from l2fc values.
#'
#' From the restructure l2fc dataframe, l2fc columns are converted to viabilities which are then used
#' to calculate HSA, Bliss, and synergy scores. This function relies on the `data.table` package for speed.
#'
#' @import data.table
#' @param restructure_l2fc
#' @param l2fc_root
#' @param viab_cap Upper bound for the viability. Defaults to 1.
#' @param names_prefix A vector of three prefixes used to identify perturbation 1, perturbation 2, and the combination.
#' @param names_sep A string used to create new column names with `names_prefix` and `l2fc_root`.
#' @retun A dataframe
calculate_synergy = function(restructured_l2fc,
                             l2fc_root = "median_l2fc",
                             viab_root = "viab",
                             viab_cap = 1,
                             names_prefix = c("pert1", "pert2", "combo"),
                             names_sep = "_") {

  # Convert data frame to data.table
  if (!data.table::is.data.table(restructured_l2fc)) {
    restructured_l2fc = data.table::as.data.table(restructured_l2fc)
  }

  # Create vector of names used to create new columns
  l2fc_names = paste(names_prefix, l2fc_root, sep = names_sep)
  viab_names = paste(names_prefix, viab_root, sep = names_sep)

  # Convert l2fcs to viabilities
  restructured_l2fc[, (l2fc_names) := .(pmin(2^get(l2fc_names[1]), viab_cap),
                                        pmin(2^get(l2fc_names[2]), viab_cap),
                                        pmin(2^get(l2fc_names[3]), viab_cap))]
  data.table::setnames(restructured_l2fc, l2fc_names, viab_names)

  # Calculate synergies
  restructured_l2fc[, hsa := pmin(get(viab_names[1]), get(viab_names[2]))]
  restructured_l2fc[, bliss := get(viab_names[1]) * get(viab_names[2])]
  restructured_l2fc[, synergy := data.table::fifelse(get(viab_names[3]) > hsa, hsa - get(viab_names[3]),
                                 data.table::fifelse(get(viab_names[3]) < bliss, bliss - get(viab_names[3]), 0))]

  # Testing non viab columns
  #restructured_l2fc[, (viab_names[3]) := .(pmin(2^get(l2fc_names[3]), viab_cap))]
  #restructured_l2fc[, hsa := pmin(pmin(2^get(l2fc_names[1]), viab_cap),
  #                                pmin(2^get(l2fc_names[2]), viab_cap))]
  #restructured_l2fc[, bliss := pmin(2^get(l2fc_names[1]), viab_cap) * pmin(2^get(l2fc_names[2]), viab_cap)]
  #restructured_l2fc[, synergy := data.table::fcase(get(viab_names[3]) > hsa, hsa - get(viab_names[3]),
  #                                                 get(viab_names[3]) < bliss, bliss - get(viab_names[3]),
  #                                                 default = 0)]
  #restructured_l2fc[, (viab_names[3]) := NULL]

  return(restructured_l2fc)
}

# Sample DMSO l2fc values
median_sample = function(x, n_samples,
                         size = 3, seed = 2,
                         replace = FALSE, prob = NULL) {
  
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
