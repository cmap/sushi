#'  normalize
#'
#'  takes a filtered dataframe of raw read counts and normalizes
#'  counts using control barcodes
#'
#' @param X - dataframe of annotated readcounts that must include the following columns:
#'           log_n or n: raw readcounts or log10(n) of read counts. Computes log_n if not present
#'           log_dose: log10 of dose at which control barcode was spiked in, if applicable
#'           sample_ID: some identifier that distinguishes between each sample
#'           Name: contains the name of the control barcode that the read corresponds to, or NA
#' @param barcodes - a vector of control barcode Name identifiers
#' @return Table with counts normalized to control barcodes
#' @export
normalize <- function(X, barcodes) {
  if (!('log_n' %in% colnames(X)) & 
      ('n' %in% colnames(X))) {
    X <- X %>% dplyr::mutate(log_n = log10(n))
  }
  
  normalized <- X %>%
    dplyr::group_by(profile_id) %>%
    dplyr::mutate(log_normalized_n = glm(y ~ x,
                                         data = dplyr::tibble(
                                           y = log_dose[Name %in% barcodes],
                                           x = log_n[Name %in% barcodes])) %>%
                    predict(newdata = dplyr::tibble(x = log_n))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(normalized_n = 10^log_normalized_n)

  return(normalized)
}