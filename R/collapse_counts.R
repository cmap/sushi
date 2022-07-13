#' collapse_counts
#' 
#' collapses l2fc values and computes MAD/sqrt(n) metrics for treatment conditions
#'
#'  @param l2fc - l2fc table with MAD/sqrt(n) metric for control condition
#'  @return - collapsed_counts 
#'  @export 
collapse_counts = function(l2fc, sig_cols=c('cell_set','treatment','dose','dose_unit','day'), cell_id_cols=c('CCLE_name', 'DepMap_ID', 'prism_cell_set') ) {
  l2fc$sig_id = do.call(paste,c(l2fc[sig_cols], sep=':'))
  
  collapsed_counts = l2fc %>% 
    dplyr::filter(control_pass_QC) %>% 
    dplyr::group_by_at(c('sig_id', cell_id_cols)) %>% 
    dplyr::summarise(trt_median_normalized_n = median(mean_normalized_n),
                     median_l2fc = median(l2fc),
                     trt_mad_sqrtN = mad(log10(mean_normalized_n))/sqrt(dplyr::n())) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(trt_pass_QC = ifelse(trt_mad_sqrtN > 0.5, F, T)) %>% 
    dplyr::relocate(trt_median_normalized_n, trt_mad_sqrtN, trt_pass_QC, median_l2fc, .after=last_col())
  
  return(collapsed_counts)
}
