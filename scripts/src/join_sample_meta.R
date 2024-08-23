#' validate_columns_exist
#' 
#' This function checks that a list of columns are present in a dataframe.
#' Columns that were not found in the dataframe are printed out.
#' 
#' @param selected_columns A vector of strings each representing a column name
#' @param df A dataframe to check against
#' @return Boolean
validate_columns_exist= function(selected_cols, df) {
  # Check that all of selected_columns are in df - base::setdiff(A, B) = A[!A %in% B].
  unmatched_cols= base::setdiff(selected_cols, colnames(df))
  
  if(length(unmatched_cols) > 0) {
    print('The following columns are missing: ')
    print(unmatched_cols)
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' join_sample_meta
#' 
#' Joins a given data frame with the sample meta.
#' 
#' @param df Input dataframe that should contain the columns specified in the "key_cols" parameter and "cell_set".
#' @param sample_meta Dataframe of the sample meta used in the run.
#' @param cell_set_meta Datafrane of the cell set metadata used in the run. This should contain the "cell_set" column.
#' @param key_cols Vector of column names used as identifiers in the sample meta.
#' @returns Data frame with additional columns from the sample meta.
join_sample_meta= function(df, sample_meta, cell_set_meta, key_cols) {
  # Validation: Check that key_cols are present in df ----
  if(validate_columns_exist(key_cols, df) == FALSE) {
    stop('Not all key_cols (printed above) are present in the provided data frame.')
  }
  
  # Validation: Check that key_cols are present in the sample meta ----
  if(validate_columns_exist(key_cols, sample_meta) == FALSE) {
    stop('Not all key_cols (printed above) are present in the sample meta.')
  }
  
  # Validation: Check that cell_set exists in df and cell_set meta ----
  if(validate_columns_exist(c('cell_set'), sample_meta) == FALSE) {
    stop('The cell_set column is NOT present in the sample meta.')
  }
  
  if(validate_columns_exist(c('cell_set'), cell_set_meta) == FALSE) {
    stop('The cell_set column is NOT present in the cell set meta.')
  }
  
  # Collapse the sample meta using key_cols and join onto the input df ----
  # Collapse unique values into a single row and then filter out columns with the separator.
  # Columns with only one unique value in a group are selected.
  collapsed_metadata= sample_meta %>% dplyr::left_join(cell_set_meta, by= 'cell_set') %>%
    dplyr::group_by(pick(all_of(key_cols))) %>% 
    dplyr::summarise(across(everything(), function(x) paste(sort(unique(x)), collapse= ':::'))) %>% dplyr::ungroup() %>%
    dplyr::select(all_of(key_cols), where(function(x) base::any(!grepl(':::', x))))

  expanded_df= dplyr::left_join(df, collapsed_metadata, 
                                by= base::intersect(colnames(df), colnames(collapsed_metadata)), 
                                relationship='many-to-one')
  
  # Print out the sample meta columns that were added to the dataframe ----
  added_cols= base::setdiff(colnames(expanded_df), colnames(df))
  if(length(added_cols > 0)) {
    print('The following columns from the sample meta were added:')
    print(added_cols)
  } else {
    print('No additional columns from the sample meta were added.')
  }
  
  return(expanded_df)
}