
# Function expects annotation column to be "growth_pattern"
apply_growth_correction = function(df, raw_l2fc_col = "l2fc", growth_pattern_col = "growth_pattern") {
  centered_df = df |>
    dplyr::filter(is.finite(.data[[raw_l2fc_col]]), !is.na(.data[[growth_pattern_col]]),
                  !is.na(cell_set)) |>
    dplyr::mutate(growth_pattern = as.factor(.data[[growth_pattern_col]]),
                  centered_l2fc = .data[[raw_l2fc_col]] - mean(.data[[raw_l2fc_col]]),
                  cell_set = as.factor(cell_set))

  # Identify number of growth patterns
  num_growth_annots = length(unique(centered_df[[growth_pattern_col]]))

  # Create right model depending on combinations
  if (num_growth_annots > 1) {
    if (length(unique(centered_df$cell_set)) > 1) {
      # Use all three columns in model formula
      model_formula = sprintf("centered_l2fc ~ %s + median_log_normalized_ctl_vehicle * cell_set",
                              growth_pattern_col)
    } else {
      # Use just two
      model_formula = sprintf("centered_l2fc ~ %s + median_log_normalized_ctl_vehicle",
                              growth_pattern_col)
    }
  } else {
    if (length(unique(centered_df$cell_set)) > 1) {
      # dont fit with growth annots
      model_formula = sprintf("centered_l2fc ~ median_log_normalized_ctl_vehicle * cell_set")
    } else {
      # just use negcon meds
      model_formula = sprintf("centered_l2fc ~ median_log_normalized_ctl_vehicle")
    }
  }

  fit = lm(formula = as.formula(model_formula), data = centered_df)
  centered_df$l2fc_uncorrected = centered_df$l2fc
  centered_df$l2fc = fit$residuals + mean(centered_df$l2fc)

  # Drop columns created by this function and return output
  return(centered_df |> dplyr::select(-centered_l2fc))
}