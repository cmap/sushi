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
#' @import data.table
#' @param cps_l2fc
#' @param cell_line_cols A vector of column names defining unique cell lines.
#' @param sig_cols A vector of column names defining unique profiles.
#' @param l2fc_col A string name of the column containing the l2fc values.
#' @param singles_type
#' @param combos_type
#' @param names_prefix A vector of three prefixes used to identify perturbation 1, perturbation 2, and the combination.
#' @param names_sep A string used to create new column names with `names_prefix` and `l2fc_col`.
#' @param ignore_cols A vector of column names to ignore. These columns might be in `sig_cols` and are ignored to prevent merging issues.
#' @return A dataframe
restructure_l2fc = function(cps_l2fc,
                            cell_line_cols,
                            sig_cols,
                            l2fc_col = "median_l2fc",
                            singles_type = "trt_cp",
                            combos_type = "trt_combo",
                            names_prefix = c("pert", "pert2", "combo"),
                            names_sep = "_",
                            ignore_cols = c("pert_type", "pert_plate")) {
  # Create new column names using prefix, l2fc_col, and names_sep
  new_names = paste(names_prefix, l2fc_col, sep = names_sep)
  
  pert1_cols = sig_cols[grepl(paste0(names_prefix[1], names_sep), sig_cols) & !sig_cols %in% ignore_cols]
  pert2_cols = sig_cols[grepl(paste0(names_prefix[2], names_sep), sig_cols) & !sig_cols %in% ignore_cols]

  # Check that pert1_cols and pert2_cols are of the same size
  if (length(pert1_cols) != length(pert2_cols)) {
    print(pert1_cols)
    print(pert2_cols)
    stop("Number of columns describing pert1 does not match number of columns describing pert2.")
  }
  
  # Slice dataframe into single agents and combos
  singles_df = cps_l2fc[pert_type == singles_type, ]
  combos_df = cps_l2fc[pert_type == combos_type, data.table::setnames(.SD, l2fc_col, new_names[3])]
  # Change the old l2fc_col name into combo + old l2fc_col
  
  # Add single agents to pert1 of the combos dataframe
  combos_df[singles_df, new_names[1] := get(l2fc_col), on = c(cell_line_cols, pert1_cols)]
  
  # Add single agents to pert2 of the combos dataframe
  cross_join_cols = setNames(pert1_cols, pert2_cols)
  combos_df[singles_df, new_names[2] := get(l2fc_col), on = c(cell_line_cols, cross_join_cols)]
  
  return(combos_df)
}

#' Calculate synergy scores
#'
#' In the context of a CPS, this function calculates synergy scores from l2fc values.
#'
#' From the restructure l2fc dataframe, l2fc columns are converted to viabilities which are then used
#' to calculate HSA, Bliss, and synergy scores. This function relies on the `data.table` package for speed
#' and thus will modify the restructured_l2fcl input.
#'
#' @import data.table
#' @param restructure_l2fc
#' @param l2fc_cols Vector of column names containing l2fc values.
#' @param viab_cap Upper bound for the viability. Defaults to 1.
#' @return A dataframe
calculate_synergy = function(restructured_l2fc, l2fc_cols, viab_cap = 1) {
  # Convert data frame to data.table
  if (!data.table::is.data.table(restructured_l2fc)) {
    restructured_l2fc = data.table::as.data.table(restructured_l2fc)
  }

  # Testing non viab columns
  restructured_l2fc[, hsa := pmin(2^pmin(get(l2fc_cols[1]), get(l2fc_cols[2])), viab_cap)]
  restructured_l2fc[, bliss := pmin(2^get(l2fc_cols[1]), viab_cap) * pmin(2^get(l2fc_cols[2]), viab_cap)]
  restructured_l2fc[, temp_viab := pmin(2^get(l2fc_cols[3]), viab_cap)]
  restructured_l2fc[, synergy := data.table::fcase(temp_viab < bliss, bliss - temp_viab,
                                                   temp_viab >= bliss & temp_viab <= hsa, 0,
                                                   temp_viab > hsa, hsa - temp_viab,
                                                   default = NA)]
  restructured_l2fc[, temp_viab := NULL] # drop temp viab column

  return(restructured_l2fc)
}

# Sample DMSO l2fc values
median_sample = function(x, n_samples, size = 3, seed = 2, replace = FALSE, prob = NULL) {
  
  # Check that there are enough entries to sample to n
  num_pick_combinations = base::choose(length(x), size)
  print(paste0(length(x), " choose ", size, " = ", num_pick_combinations))
  
  # Stop if the number to sample to is too low
  if(num_pick_combinations < n_samples) {
    stop("ERROR: Cannot sample up to n.")
  }
  
  # Set seed
  base::set.seed(seed)
  # Resample values
  resampled_values = base::replicate(n_samples, median(base::sample(x, size, replace = replace, prob = prob)))
  
  return(resampled_values)
}

#' Get cdf to change into pvalue
#' 
#' @import rhdf5
#' @param group_name
#' @param synergy_value
#' @param h5_file
#' @param n_samples
get_pvalue = function(group_name, synergy_value, h5_file, n_samples = 10000) {
  ecdf_obj = stats::ecdf(rhdf5::H5Dread(rhdf5::H5Dopen(h5_file, name = group_name)))
  cdf_value = smooth_pvalue(ecdf_obj(synergy_value), n_samples)
  pvalue = 2 * pmin(1 - cdf_value, cdf_value)
  
  return(pvalue)
}

#' Smooth pvalue
#'
#' Use Laplace smoothing on pvalues to prevent zero values.
#'
#' @param pval_naive description
#' @param n_samples Number of samples
#' @return Smoothed pvalue
smooth_pvalue = function(pval_naive, n_samples) {
  return((pval_naive * n_samples + 1) / (n_samples + 2))
}
