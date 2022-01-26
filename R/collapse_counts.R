#' collapse_counts
#' 
#' collapses filtered normalized counts and computes MAD/sqrt(n) metrics.
#' cell lines with MAD/sqrt(n) > 0.5 are filtered out, and file with filtered out cell lines is written. 
#' log10(median counts) vs. MAD/sqrt(n) graph is saved, and collapsed filtered count table is returned
#'
#'  @param l2fc - normalized counts with bad replicates already filtered out
#'  @returns collapsed_counts 
#'  @export 

collapse_counts = function(l2fc) {
  collapsed_counts = l2fc %>% 
    dplyr::filter(control_pass_QC) %>% 
    dplyr::group_by_at(setdiff(names(.), c("bio_rep", "sum_normalized_n", "control_mad_sqrtN", "l2fc", "control_pass_QC", "control_median_normalized_n"))) %>% 
    dplyr::summarise(trt_median_normalized_n = median(sum_normalized_n),
                     median_l2fc = median(l2fc),
                     trt_mad_sqrtN = mad(log10(sum_normalized_n))/sqrt(dplyr::n())) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(trt_pass_QC = ifelse(trt_mad_sqrtN > 0.5, F, T)) %>% 
    dplyr::relocate(trt_median_normalized_n, trt_mad_sqrtN, trt_pass_QC, median_l2fc, .after=last_col())
  
  return(collapsed_counts)
}