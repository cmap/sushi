#' compute_l2fc
#' 
#' takes normalized counts and computes log-fold change values as compared to the annotated control columns
#'
#' @param normalized_counts - table with normalized_n column and control_sample column that designates the 
#'          name of the control sample for each treatment sample
#' @param control_type - string that denotes which compute samples to compute log fold change. Matches trt_type field
#' @return l2fc data.frame with l2fc column
#' 
compute_l2fc = function(normalized_counts, control_type) {
  treatments = normalized_counts %>% 
    dplyr::filter(trt_type!=control_type, trt_type!="day_0",
           is.na(Name)) %>% 
    dplyr::select(-Name, -log_dose, -n, -log_n, -log_normalized_n) %>% 
    dplyr::group_by_at(setdiff(names(.), c("normalized_n", "tech_rep"))) %>% 
    dplyr::summarise(sum_normalized_n = sum(normalized_n)) %>% 
    dplyr::ungroup()
  
  controls = normalized_counts %>% 
    dplyr::filter(trt_type==control_type,
           is.na(Name)) %>% 
    dplyr::select(-Name, -log_dose, -n, -log_n, -log_normalized_n) %>% 
    dplyr::group_by_at(setdiff(names(.), c("normalized_n", "tech_rep"))) %>% 
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
    dplyr::relocate(project_code, CCLE_name, DepMap_ID, prism_cell_set, profile_id, trt_type, control_barcodes,
                    bio_rep) 
  
  return(l2fc)
}