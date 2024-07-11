#' generate_biomarkers
#' 
#' collapses filtered normalized counts and computes MAD/sqrt(n) metrics.
#' cell lines with MAD/sqrt(n) > 0.5/log10(2) are filtered out, and file with filtered out cell lines is written. 
#' log10(median counts) vs. MAD/sqrt(n) graph is saved, and collapsed filtered count table is returned
#'
#' 
#' @param collapsed_values: collapsed_l2fc values
#' @return Biomarker object containing lin_out, rf_out and disc_out tables
#' @export
generate_biomarkers = function(collapsed_values) {
  bio_in = collapsed_values %>%
    dplyr::filter(trt_pass_QC) %>%
    reshape2::dcast(DepMap_ID~sig_id, value.var="median_l2fc") %>%
    tibble::column_to_rownames("DepMap_ID")

  bio_out = cdsrbiomarker::get_biomarkers(bio_in)

  return(bio_out)
}
