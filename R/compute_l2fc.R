#' compute_l2fc
#' 
#' takes normalized counts and computes log-fold change values as compared to the designated control condition
#'
#' @param normalized_counts - table with normalized_n column and trt_type column that designates the 
#'          the control sample
#' @param control_type - string that denotes which samples to compute log fold change against. Matches trt_type field. negcon by default.
#' @param sig_cols sig_cols - a vector of column names denoting which values specify each individual signature
#'                    cell_set,treatment,dose,dose_unit,day by default.
#' @return - l2fc data.frame with l2fc column
#' @export
compute_l2fc = function(normalized_counts, 
                        control_type = "negcon",
                        sig_cols=c('cell_set','treatment','dose','dose_unit','day')) {
  normalized_counts$sig_id = do.call(paste,c(normalized_counts[sig_cols], sep=':'))
  
  treatments = normalized_counts %>% 
    dplyr::filter(trt_type!=control_type, trt_type!="day_0",
           is.na(Name)) %>% 
    dplyr::select(-Name, -log_dose, -n, -log_n, -log_normalized_n) %>% 
    dplyr::group_by_at(setdiff(names(.), c("normalized_n", "tech_rep", "profile_id"))) %>% 
    dplyr::summarise(sum_normalized_n = sum(normalized_n)) %>% 
    dplyr::ungroup()
  
  controls = normalized_counts %>% 
    dplyr::filter(trt_type==control_type,
           is.na(Name)) %>% 
    dplyr::select(-Name, -log_dose, -n, -log_n, -log_normalized_n) %>% 
    dplyr::group_by_at(setdiff(names(.), c("normalized_n", "tech_rep", "profile_id"))) %>% 
    dplyr::summarise(sum_normalized_n = sum(normalized_n)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(CCLE_name, DepMap_ID, prism_cell_set) %>% 
    dplyr::summarise(control_median_normalized_n = median(sum_normalized_n),
                     control_mad_sqrtN = mad(log10(sum_normalized_n))/sqrt(dplyr::n())) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(control_pass_QC = ifelse(control_mad_sqrtN > 0.5, F, T)) %>% 
    dplyr::select(CCLE_name, DepMap_ID, prism_cell_set, control_median_normalized_n, control_mad_sqrtN, control_pass_QC)
  
  if(nrow(controls)==0) {
    print("No samples found for indicated control type.")
    stop()
  }
  
  l2fc = treatments %>% 
    merge(controls, by=c("CCLE_name", "DepMap_ID", "prism_cell_set"), all.x=T, all.y=T) %>% 
    dplyr::mutate(l2fc=log2(sum_normalized_n/control_median_normalized_n)) %>% 
    dplyr::relocate(project_code, CCLE_name, DepMap_ID, prism_cell_set, trt_type, control_barcodes, sig_id,
                    bio_rep) 
  
  return(l2fc)
}
