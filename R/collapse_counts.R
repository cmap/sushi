#' validate_columns_exist
#' 
#' This function checks that a list of columns are present in a dataframe.
#' 
#' @param selected_columns A vector of strings each representing a column name
#' @param df A dataframe to check against
#' @return Boolean
validate_columns_exist= function(selected_columns, df) {
  # Check that all of selected_columns are in df
  if(any(!selected_columns %in% colnames(df))) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' validate_num_bio_reps
#' 
#' This function checks that all the expected flowcells are present in a table of detected flowcells.
#' There can be more detected flowcells than there are expected flowcells.
#' 
#' @param detected_flowcells A dataframe with the columns "flowcell_name" and "flowcell_lane".
#' @param expected_flowcells A dataframe with the columns "flowcell_name" and "flowcell_lane".
validate_num_bio_reps= function(max_bio_rep_id, max_bio_rep_count) {
  if(max_bio_rep_id > 1 & max_bio_rep_count == 1) {
    stop('Unable to collapse bio_reps over the specified sig_cols.')
  } else if(max_bio_rep_id < max_bio_rep_count) {
    stop('Bio_reps were incorrectly collapsed resulting in more replicates than specified.')
  } else if(max_bio_rep_id > 1 & max_bio_rep_id > max_bio_rep_count) {
    print('Warning - Number of replicates that were collapses is smaller than the expected number of replicates.')
    print('This could be due to a problem in processing or from poor data/sequencing quality.')
  } else {}
}

#' collapse_counts
#' 
#' Collapses the l2fc values of biological replicates for treatment condtions and
#' computes the MAD/sqrt(n).
#'
#'  @param l2fc Dataframe of l2fc values The following columns are required -
#'              DepMap_ID, CCLE_name, counts_flag, mean_n, mean_normalized_n, and l2fc.
#'  @param sig_cols List of columns that define an individual condition. This should not include any replicates.
#'                  The columns in this list should be present in the l2fc dataframe.
#'  @return - collapsed_counts
#'  @export 
collapse_counts = function(l2fc, sample_meta, sig_cols) {
  # Validation: Check that sig_cols are present in l2fc ----
  if(validate_columns_exist(sig_cols, l2fc) == FALSE) {
    print(sig_cols)
    stop('Not all sig_cols (printed above) are present in the l2fc file.')
  }
  
  # Validation: Check that sig_cols are present in sample meta ----
  if(validate_columns_exist(sig_cols, sample_meta) == FALSE) {
    print(sig_cols)
    stop('Not all sig_cols (printed above) are present in the sample meta.')
  }
  
  # Create sample meta annotations for each unique sig_id ----
  # Select ONLY columns in the sample meta that describe the sig_cols and are not sequencing related columns
  sig_id_metadata= sample_meta %>% tidyr::unite(col='sig_id', all_of(sig_cols), sep=':', remove=F) %>%
    dplyr::group_by(sig_id) %>% 
    dplyr::summarise(across(everything(), function(x) list(unique(x)))) %>% dplyr::ungroup() %>%
    dplyr::select(where(function(x) max(lengths(x)) == 1) & 
                    !any_of(c('flowcell_names', 'flowcell_lanes', 'index_1', 'index_2'))) %>% 
    tidyr::unnest(cols= colnames(.))
  
  # Median collapsing bio replicates ----
  # Potentially parameter? 
  cell_line_cols= c('DepMap_ID', 'CCLE_name')
  
  collapsed_counts= l2fc %>% dplyr::filter(is.na(counts_flag)) %>% 
    tidyr::unite(col='sig_id', all_of(sig_cols), sep=':', remove=T) %>%
    dplyr::group_by(pick(all_of(c('project_code', cell_line_cols, 'sig_id')))) %>%
    dplyr::summarise(trt_median_n= median(mean_n), trt_median_normalized_n= median(mean_normalized_n),
                     trt_mad_sqrtN= mad(log2(mean_normalized_n)) / sqrt(dplyr::n()),
                     median_l2fc= median(l2fc), num_bio_reps= dplyr::n()) %>% dplyr::ungroup() %>% 
    dplyr::mutate(trt_MAD_QC= ifelse(trt_mad_sqrtN > 0.5/log10(2), F, T)) %>% # Adjusted cut off from log10 to log2
    dplyr::left_join(sig_id_metadata, by= intersect(colnames(.), colnames(sig_id_metadata)), relationship='many-to-one') %>%
    dplyr::relocate(trt_median_n, trt_median_normalized_n, trt_mad_sqrtN, 
                    num_bio_reps, median_l2fc, trt_MAD_QC, .after=last_col())
  # Metadata joined on potentially project_code and sig_id so that equality is assessed on strings.
  # Joining with sig_cols could result in rounding/floating point equality issues like 0.1 + 0.2 != 0.3
  
  # Validation: Check that replicates were collapsed ----
  if('bio_rep' %in% colnames(l2fc)) {
    max_bio_rep_id= max(unique(l2fc$bio_rep))
    max_bio_rep_count= max(unique(collapsed_counts$num_bio_reps))
    validate_num_bio_reps(max_bio_rep_id, max_bio_rep_count)
  }
  
  return(collapsed_counts)
}
