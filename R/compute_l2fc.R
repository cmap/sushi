#' compute_l2fc
#' 
#' takes normalized counts and computes log-fold change values as compared to the designated control condition
#'
#' @param normalized_counts - table with normalized_n column and trt_type column that designates the 
#'          the control sample
#' @param control_type - string that denotes which samples to compute log fold change against. Matches trt_type field. negcon by default.
#' @param sig_cols - a vector of column names denoting which values specify each individual signature
#'                    cell_set,treatment,dose,dose_unit,day by default.
#' @param control_sigs - a vector of column names denoting which values specify each individual control condition
#'                    cell_set,day by default.
#' @param count_col_name - a string containing the name of the column to use as counts to calculate l2fc values. Generally normalized_n if
#'                           running on normalied_counts or n if running on filtered_counts
#' @return - l2fc data.frame with l2fc column
#' @export
compute_l2fc = function(normalized_counts, 
                        control_type = "negcon",
                        sig_cols=c('cell_set','treatment','dose','dose_unit','day'),
                        control_sigs= c('cell_set', 'day'), # will probably be a subset of sig_cols
                        count_col_name="normalized_n") {
  normalized_counts$sig_id = do.call(paste,c(normalized_counts[sig_cols], sep=':'))
  
  normalized_counts = normalized_counts %>% 
    dplyr::filter(!(trt_type %in% c("empty", "", "CB_only")) & !is.na(trt_type))
  
  # collapse tech reps for treatments
  treatments = normalized_counts %>% 
    dplyr::filter(trt_type!=control_type, trt_type!="day_0",
           !is.na(CCLE_name)) %>% 
    dplyr::select(-any_of(c("Name", "log_dose", "n", "log_n", "log_normalized_n", "normalized_n")
                          [!(c("Name", "log_dose", "n", "log_n", "log_normalized_n", "normalized_n") %in% count_col_name)])) %>% 
    dplyr::group_by_at(setdiff(names(.), c(count_col_name, "tech_rep", "profile_id"))) %>% 
    dplyr::summarise(mean_normalized_n = mean(!! rlang::sym(count_col_name))) %>% 
    dplyr::ungroup()
  
  # collapse controls 
  controls = normalized_counts %>% 
    dplyr::filter(trt_type==control_type,
           !is.na(CCLE_name)) %>% 
    dplyr::select(-any_of(c("Name", "log_dose", "n", "log_n", "log_normalized_n", "normalized_n")
                          [!(c("Name", "log_dose", "n", "log_n", "log_normalized_n", "normalized_n") %in% count_col_name)])) %>% 
    dplyr::group_by_at(setdiff(names(.), c(count_col_name, "tech_rep", "profile_id"))) %>% 
    dplyr::summarise(mean_normalized_n = mean(!! rlang::sym(count_col_name))) %>% 
    dplyr::ungroup() %>% # tech rep collapse
    dplyr::group_by_at(c(control_sigs, "CCLE_name", "DepMap_ID", "prism_cell_set")) %>% 
    dplyr::summarise(control_median_normalized_n = median(mean_normalized_n),
                     control_mad_sqrtN = mad(log10(mean_normalized_n))/sqrt(dplyr::n())) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(control_pass_QC = ifelse(control_mad_sqrtN > 0.5, F, T)) %>% 
    dplyr::select(CCLE_name, DepMap_ID, prism_cell_set, any_of(control_sigs), control_median_normalized_n, control_mad_sqrtN, control_pass_QC)
  
  if(nrow(controls)==0) {
    print("No samples found for indicated control type.")
    stop()
  }
  
  l2fc = treatments %>% 
    merge(controls, by= c("CCLE_name", "DepMap_ID", "prism_cell_set", control_sigs), all.x=T, all.y=T) %>%
    dplyr::mutate(l2fc=log2(mean_normalized_n/control_median_normalized_n)) %>%
    dplyr::relocate(project_code, CCLE_name, DepMap_ID, prism_cell_set, trt_type, control_barcodes, sig_id, bio_rep)
  
  return(l2fc)
}
