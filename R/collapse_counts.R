#' collapse_counts
#' 
#' collapses l2fc values and computes MAD/sqrt(n) metrics for treatment conditions
#'
#'  @param l2fc - l2fc table with MAD/sqrt(n) metric for control condition
#'  @return - collapsed_counts 
#'  @export 
collapse_counts = function(l2fc) {
  collapsed_counts = l2fc %>% 
    dplyr::filter(control_pass_QC) %>% 
    dplyr::group_by_at(setdiff(names(.), c('bio_rep', 'mean_normalized_n', 'num_tech_reps', 
                                           'control_median_normalized_n', 'control_mad_sqrtN', 'num_ctrl_bio_reps', 'control_pass_QC',
                                           'l2fc'))) %>% 
    dplyr::summarise(trt_median_normalized_n= median(mean_normalized_n),
                     median_l2fc= median(l2fc), num_bio_reps= dplyr::n(),
                     trt_mad_sqrtN= mad(log2(mean_normalized_n)) / sqrt(dplyr::n())) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(trt_pass_QC= ifelse(trt_mad_sqrtN > 0.5/log10(2), F, T)) %>% # New: adjusted cut off to log2
    dplyr::relocate(trt_median_normalized_n, trt_mad_sqrtN, num_bio_reps, median_l2fc, trt_pass_QC, .after=last_col())
  
  return(collapsed_counts)
}
