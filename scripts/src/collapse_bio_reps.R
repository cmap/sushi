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

#' collapse_bio_reps
#' 
#' Collapses the l2fc values of biological replicates for treatment conditions and
#' computes the MAD/sqrt(n).
#'
#' @param l2fc Dataframe of l2fc values The following columns are required -
#'              DepMap_ID, CCLE_name, counts_flag, mean_n, mean_normalized_n, and l2fc.
#' @param sig_cols List of columns that define an individual condition. This should not include any replicates.
#'                  The columns in this list should be present in the l2fc dataframe.
#' @param cell_line_cols List of columns that define a cell line. Defaults to project_code, DepMap_ID, and CCLE_name
#' @returns - collapsed_counts
collapse_bio_reps = function(l2fc) {
  collapsed_counts = l2fc %>%
    dplyr::filter(is.na(counts_flag)) %>%
    dplyr::group_by_at(setdiff(names(.), c('bio_rep', 'mean_n','mean_normalized_n', 'num_tech_reps', 'control_median_n',
                                           'control_median_normalized_n', 'control_mad_sqrtN', 'num_ctrl_bio_reps',
                                           'control_MAD_QC','l2fc', 'counts_flag', 'lysate_plate', 'lysate_well', 'pert_plate','sig_id'))) %>%
    dplyr::summarise(trt_median_n= median(mean_n), trt_median_normalized_n= median(mean_normalized_n),
                     trt_mad_sqrtN= mad(log2(mean_normalized_n)) / sqrt(dplyr::n()),
                     median_l2fc= median(l2fc), num_bio_reps= dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(trt_MAD_QC= ifelse(trt_mad_sqrtN > 0.5/log10(2), F, T)) %>% # New: adjusted cut off to log2
    dplyr::relocate(trt_median_n, trt_median_normalized_n, trt_mad_sqrtN,
                    num_bio_reps, median_l2fc, trt_MAD_QC, .after=last_col())

}
  return(collapsed_counts)
