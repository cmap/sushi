
# Function expects annotation column to be "growth_pattern"
apply_growth_correction = function(df, raw_l2fc_col = "l2fc", growth_pattern_col = "growth_pattern") {
  centered_df = df |>
    dplyr::filter(is.finite(.data[[raw_l2fc_col]]), !is.na(.data[[growth_pattern_col]])) |>
    dplyr::mutate(growth_pattern = as.factor(.data[[growth_pattern_col]]),
                  centered_l2fc = .data[[raw_l2fc_col]] - mean(.data[[raw_l2fc_col]]))

  # Identify number of growth patterns
  num_growth_annots = length(unique(centered_df[[growth_pattern_col]]))

  # Run correction if there are more than 1 growth pattern
  if (num_growth_annots > 1) {
    model_formula = sprintf("centered_l2fc ~ %s", growth_pattern_col)
    fit = lm(formula = as.formula(model_formula), data = centered_df)
    centered_df$l2fc_uncorrected = centered_df$l2fc
    centered_df$l2fc = fit$residuals + mean(centered_df$l2fc)
  } else {
    message("Only one growth pattern was detected. Growth pattern correction is not applied.")
    centered_df$l2fc_uncorrected = centered_df$l2fc
  }

  # Drop columns created by this function and return output
  return(centered_df |> dplyr::select(-centered_l2fc))
}