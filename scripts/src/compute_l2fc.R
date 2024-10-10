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
#' @param count_threshold - threshold for determining low counts, defaults to 40.
#' @param cell_line_cols - Vector of columns that define a cell line. Defaults to project_code, DepMap_ID, and CCLE_name
#' @return - l2fc data.frame with l2fc column
#' @export
compute_l2fc= function(normalized_counts,
                       control_type = "negcon",
                       sig_cols=c('cell_set','treatment','dose','dose_unit','day'),
                       ctrl_cols= c('cell_set', 'day'), # will probably be a subset of sig_cols
                       count_col_name="normalized_n", count_threshold= 40,
                       cell_line_cols= c('project_code', 'DepMap_ID', 'CCLE_name')) {
  
  # Validation: Check that sig_cols are in normalized_counts ----
  if(validate_columns_exist(sig_cols, normalized_counts) == FALSE) {
    print(sig_cols)
    stop('Not all sig_cols (printed above) are present in normalized_counts.')
  }
  
  # Validation: Check that cell_line_cols are in normalized_counts ----
  if(validate_columns_exist(cell_line_cols, normalized_counts) == FALSE) {
    print(cell_line_cols)
    stop('Not all cell_line_cols (printed above) are present in normalized_counts.')
  }
  
  # Collapsing technical replicates ----
  # Detect bio_rep column to be used to collapse technical replicates
  if('bio_rep' %in% colnames(normalized_counts)) {
    bio_rep_id_cols= c(sig_cols, 'bio_rep')
  } else {
    bio_rep_id_cols= sig_cols
    print('WARNING: bio_rep column not detected. Assuming that there are NO biological replicates.') 
    print('Technical replicate collapse will be performed across the sig_cols.')
  }
  
  # collapse tech reps
  print('Collapsing technical replicates on the following columns: ')
  print(unique(c(cell_line_cols, 'trt_type', bio_rep_id_cols, ctrl_cols)))
  collapsed_tech_rep= normalized_counts %>%
    dplyr::filter(!(trt_type %in% c("empty", "", "CB_only")) & !is.na(trt_type), !is.na(CCLE_name)) %>%
    dplyr::group_by(pick(all_of(c(cell_line_cols, 'trt_type', bio_rep_id_cols, ctrl_cols)))) %>%
    dplyr::summarise(mean_n= mean(n),
                     mean_normalized_n = mean(!!rlang::sym(count_col_name))) %>% dplyr::ungroup()
  
  # Print out the occurrence of each count of tech_reps
  # print('Number of technical replicate collapsed across all cell lines and biological replicates:')
  # print(collapsed_tech_rep %>% dplyr::group_by(num_tech_reps) %>% 
  #         dplyr::summarise(count= dplyr::n()) %>% dplyr::ungroup())
    
  # Pull out negative controls and collapse any biological replicates ----
  print('Collapsing control conditions on the following columns: ')
  print(unique(c(cell_line_cols, ctrl_cols)))
  controls= collapsed_tech_rep %>% dplyr::filter(trt_type== control_type) %>% 
    dplyr::group_by(pick(all_of(c(cell_line_cols, ctrl_cols)))) %>%
    dplyr::summarise(control_median_n= median(mean_n),
                     control_median_normalized_n = median(mean_normalized_n),
                     control_mad_sqrtN = mad(log2(mean_normalized_n))/sqrt(dplyr::n()),
                     num_ctrl_bio_reps = dplyr::n()) %>% dplyr::ungroup()
  
  # Validation: Check that negative controls were extracted ----
  if(nrow(controls)==0) {
    stop("No samples found for the indicated control_type.")
  }
  
  # Join neg_cons and compute l2fc ----
  l2fc= collapsed_tech_rep %>% dplyr::filter(!trt_type %in% c(control_type, 'day_0')) %>% 
    dplyr::inner_join(controls, by= c(cell_line_cols, ctrl_cols), relationship='many-to-one') %>%
    dplyr::mutate(l2fc= log2(mean_normalized_n/control_median_normalized_n),
                  counts_flag= ifelse(control_median_n < count_threshold, paste0('negcon<', count_threshold), NA))
  
  return(l2fc)
}
