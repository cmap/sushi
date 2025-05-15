# Collection of functions used by the synergy module

#' restructure_l2fc
#'
#' In the context of a CPS, this function restructures a l2fc file by creating additional
#' columns containing the l2fc values for the single agents.
#'
#' This function splits the l2fc data.table into rows with single agents and rows with combinaitons using
#' `single_type` and `combo_type` from the column `pert_type`. The single agent rows are joined onto the
#' combination rows - once for each member of the combination. The new l2fc columns created by the joins
#' are named following `new_col_names`.
#'
#' @import data.table
#' @param cps_l2fc
#' @param join_cols A vector of column names used for joining.
#' @param pert_cols_list A list of vectors where each item is a vector of column names describing a perturbation.
#' @param l2fc_col A string name of the column containing the l2fc values.
#' @param single_type String value in the pert_type column that denotes single agent conditions.
#' @param combo_type String value in the pert_type column that denotes combination conditions.
#' @param new_col_names A vector of l2fc column names to be created by this function.
#' @return A data.table
restructure_l2fc = function(cps_l2fc, join_cols, pert_cols_list,
                            l2fc_col = "median_l2fc",
                            single_type = "trt_cp",
                            combo_type = "trt_combo",
                            new_col_names = c("pert1_l2fc", "pert2_l2fc", "combo_l2fc")) {

  # Check that pert_type is a column in cps_l2fc
  if (!"pert_type" %in% colnames(cps_l2fc)) {
    stop("The column pert_type is not in the provided data.table.")
  }

  # Check that the specified types are in the pert_type column
  if (!all(c(single_type, combo_type) %in% unique(cps_l2fc$pert_type))) {
    stop(paste0(single_type, " and ", combo_type, " are not both present in the pert_type column."))
  }

  # Check that vectors in pert_cols_list are of the same size
  if (length(unique(lengths(pert_cols_list))) !=  1) {
    print(pert_cols_list)
    stop("Number of columns describing each perturbation does not match.")
  }

  # Slice dataframe into single agents and combos
  singles_df = cps_l2fc[pert_type == single_type, ]
  combos_df = cps_l2fc[pert_type == combo_type,
                       data.table::setnames(.SD, l2fc_col, new_col_names[length(new_col_names)])]
  # In combos_df rename l2fc_col to last item of new_col_names

  # For each single agent, add its l2fc column onto combos_df
  # Loop over new_col_names ignoring last item
  for (idx in seq_along(new_col_names[-length(new_col_names)])) {
    # The first iterations will not be a cross join, but subsequent joins will be
    cross_join_cols = setNames(pert_cols_list[[1]], pert_cols_list[[idx]])
    combos_df[singles_df, new_col_names[idx] := get(l2fc_col), on = c(join_cols, cross_join_cols)]
  }

  # Rorder columns - move combo l2fc to after single l2fcs
  data.table::setcolorder(combos_df, neworder = new_col_names[length(new_col_names)],
                          after = new_col_names[length(new_col_names) - 1], skip_absent = TRUE)

  return(combos_df)
}

#' Calculate synergy scores
#'
#' In the context of a CPS, this function calculates synergy scores from l2fc values.
#'
#' From the restructure l2fc data.table, l2fc columns are converted to viabilities which are then used
#' to calculate HSA, Bliss, and synergy scores. This function relies on the `data.table` package for speed
#' and thus will modify the restructured_l2fc input.
#'
#' @import data.table
#' @param restructure_l2fc
#' @param l2fc_cols Vector of column names containing l2fc values.
#' @param viab_cap Upper bound for the viability. Defaults to 1.
#' @return A data.table
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
  restructured_l2fc[, temp_viab := NULL] # drop temp_viab column

  return(restructured_l2fc)
}

# Sample DMSO l2fc values
median_sample = function(x, n_samples, size = 3, seed = 2, replace = FALSE, prob = NULL) {
  # Check that there are enough entries to sample to n
  num_pick_combinations = base::choose(length(x), size)
  print(paste0(length(x), " choose ", size, " = ", num_pick_combinations))

  # Stop if the number to sample up to is too high
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
