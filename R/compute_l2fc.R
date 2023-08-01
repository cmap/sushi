#' compute_l2fc
#' 
#' takes normalized counts and computes log-fold change values as compared to the designated control condition
#'
#' @param normalized_counts - table with normalized_n column and trt_type column that designates the 
#'          the control sample
#' @param control_type - string that denotes which samples to compute log fold change against. Matches trt_type field. 
#'          negcon by default.
#' @param sig_cols - a vector of column names denoting which values specify each individual signature
#'                    cell_set,treatment,dose,dose_unit,day by default.
#' @param ctrl_cols - a vector of column names denoting which values specify each individual control condition
#'                    cell_set,day by default.
#' @param count_col_name - a string containing the name of the column to use as counts to calculate l2fc values. 
#'          Generally normalized_n if running on normalied_counts or n if running on filtered_counts
#' @return - l2fc data.frame with l2fc column
#' @export
compute_l2fc = function(normalized_counts, 
                        control_type = "negcon",
                        sig_cols=c('cell_set','treatment','dose','dose_unit','day'),
                        ctrl_cols= c('cell_set', 'day'), # will probably be a subset of sig_cols
                        count_col_name="normalized_n") {
  
  if(!all(ctrl_cols %in% sig_cols)) {
    print("Control columns are not a subset of sig columns.") # new
    stop()
  }
  
  normalized_counts$sig_id= do.call(paste,c(normalized_counts[sig_cols], sep=':'))
  
  # collapse tech reps
  collapsed_tech_rep= normalized_counts %>%
    dplyr::filter(!(trt_type %in% c("empty", "", "CB_only")) & !is.na(trt_type), !is.na(CCLE_name)) %>%
    dplyr::group_by_at(setdiff(names(.), c('pcr_plate','pcr_well', 'Name', 'log2_dose', 'cb_intercept',
                                           'profile_id', 'tech_rep', 'n', 'log2_n', 'normalized_n', 'log2_normalized_n',
                                           'low_counts', count_col_name))) %>% 
    dplyr::summarise(mean_normalized_n = mean(!! rlang::sym(count_col_name)), 
                     num_tech_reps= n()) %>% 
    dplyr::ungroup()
  
  # collapse controls
  controls= collapsed_tech_rep %>% 
    dplyr::filter(trt_type==control_type) %>% 
    dplyr::group_by_at(c('project_code', 'CCLE_name', 'DepMap_ID', 'prism_cell_set', ctrl_cols)) %>% 
    dplyr::summarise(control_median_normalized_n = median(mean_normalized_n),
                     control_mad_sqrtN = mad(log2(mean_normalized_n))/sqrt(dplyr::n()),
                     num_ctrl_bio_reps = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(control_pass_QC = ifelse(control_mad_sqrtN > 0.5/log10(2), F, T)) #%>% # New: adjusted cut off to log2
  
  if(nrow(controls)==0) {
    print("No samples found for indicated control type.")
    stop()
  }
  
  l2fc= collapsed_tech_rep %>% dplyr::filter(!trt_type %in% c(control_type, 'day_0')) %>% 
    merge(controls, by= c('project_code',"CCLE_name", "DepMap_ID", "prism_cell_set", ctrl_cols), all.x=T, all.y=T) %>%
    dplyr::mutate(l2fc= log2(mean_normalized_n/control_median_normalized_n)) %>%
    dplyr::relocate(project_code, CCLE_name, DepMap_ID, prism_cell_set, trt_type, control_barcodes, sig_id, bio_rep)
  
  return(l2fc)
}
