
# Initial function from MK
apply_growth_correction = function(df) {
  df = df |>
    dplyr::mutate(culture = as.factor(culture)) |>
    dplyr::filter(is.finite(LFC_corrected), !is.na(culture)) %>%
    dplyr::mutate(y = LFC_corrected - mean(LFC_corrected))

  fit <- lm(y ~ culture, df)
  df$LFC_regressed <- fit$residuals + mean(df$LFC_corrected)

  df
}