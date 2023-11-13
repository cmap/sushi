#' collapse_counts
#' 
#' collapses l2fc values and computes MAD/sqrt(n) metrics for treatment conditions
#'
#'  @param l2fc - l2fc table with MAD/sqrt(n) metric for control condition
#'  @return - collapsed_counts 
#'  @export 
collapse_counts = function(l2fc) {
  collapsed_counts = l2fc %>% 
    dplyr::filter(is.na(counts_flag)) %>% 
    dplyr::group_by_at(setdiff(names(.), c('bio_rep', 'mean_n','mean_normalized_n', 'num_tech_reps', 'control_median_n',
                                           'control_median_normalized_n', 'control_mad_sqrtN', 'num_ctrl_bio_reps', 
                                           'control_MAD_QC','l2fc', 'counts_flag'))) %>% 
    dplyr::summarise(trt_median_n= median(mean_n), trt_median_normalized_n= median(mean_normalized_n),
                     trt_mad_sqrtN= mad(log2(mean_normalized_n)) / sqrt(dplyr::n()),
                     median_l2fc= median(l2fc), num_bio_reps= dplyr::n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(trt_MAD_QC= ifelse(trt_mad_sqrtN > 0.5/log10(2), F, T)) %>% # New: adjusted cut off to log2
    dplyr::relocate(trt_median_n, trt_median_normalized_n, trt_mad_sqrtN, 
                    num_bio_reps, median_l2fc, trt_MAD_QC, .after=last_col())
  
  return(collapsed_counts)
}
